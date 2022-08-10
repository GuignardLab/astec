import os
import sys
import copy
import statistics
import numpy as np
import math
import multiprocessing
import itertools
import pickle as pkl
from sklearn.neighbors import NearestNeighbors
from operator import itemgetter


import astec.utils.common as common
import astec.utils.atlas_division as uatlasd
import astec.utils.atlas_cell as uatlasc
import astec.utils.atlas_embryo as uatlase
import astec.utils.icp as icp
import astec.utils.ioproperties as ioproperties
import astec.utils.ascidian_name as uname
import astec.utils.neighborhood_distance as uneighborhood


#
#
#
#
#

monitoring = common.Monitoring()

_instrumented_ = False


########################################################################################
#
# transformations
#
########################################################################################

def _rotation_matrix_from_rotation_vector(rot):
    # Euler-Rodrigues formula
    # R = I + f(theta) X(r) + g(theta) [X(r) * X(r)]
    # theta: rotation angle (modulus of rotation vector),
    # g(theta) = (1 - cos(theta)) / (theta * theta)
    # f(theta) = sin(theta) / theta
    # X(r): matrix of the cross product by r
    mat = np.zeros((3, 3))

    t2 = rot[0] * rot[0] + rot[1] * rot[1] + rot[2] * rot[2]
    theta = math.sqrt(t2)

    if theta > 1e-8:
        f = math.sin(theta) / theta
        g = (1.0 - math.cos(theta)) / t2

        mat[0, 0] = 1.0 - g * (rot[1] * rot[1] + rot[2] * rot[2])
        mat[1, 1] = 1.0 - g * (rot[2] * rot[2] + rot[0] * rot[0])
        mat[2, 2] = 1.0 - g * (rot[0] * rot[0] + rot[1] * rot[1])

        mat[0, 1] = g * rot[0] * rot[1]
        mat[0, 2] = g * rot[0] * rot[2]
        mat[1, 2] = g * rot[2] * rot[1]

        mat[1, 0] = mat[0, 1]
        mat[2, 0] = mat[0, 2]
        mat[2, 1] = mat[1, 2]

        mat[0, 1] -= f * rot[2]
        mat[0, 2] += f * rot[1]
        mat[1, 2] -= f * rot[0]

        mat[1, 0] += f * rot[2]
        mat[2, 0] -= f * rot[1]
        mat[2, 1] += f * rot[0]

    else:
        mat[0, 0] = 1.0
        mat[1, 1] = 1.0
        mat[2, 2] = 1.0

    return mat


########################################################################################
#
# classes
# - computation parameters
#
########################################################################################

class EmbryoRegistrationParameters(common.PrefixedParameter):
    def __init__(self, prefix=None):

        common.PrefixedParameter.__init__(self, prefix=prefix)

        if "doc" not in self.__dict__:
            self.doc = {}

        #
        # initialization of the 3D rigid transformation
        #
        doc = "\t To coregister two embryos, a first 3D rotation is done that align vectors issued from\n"
        doc += "\t the floating embryo (the embryo to be named) onto vectors issued from the reference\n"
        doc += "\t embryo (an already named embryo)\n"
        doc += "\t - 'sphere_wrt_z': an uniform sampling of 3D directions is done for the floating\n"
        doc += "\t    embryo (parameter 'direction_sphere_radius') while the 'z' direction is used\n"
        doc += "\t    for the reference embryo\n"
        doc += "\t - 'sphere_wrt_symaxis': an uniform sampling of 3D directions is done for the floating\n"
        doc += "\t    embryo (parameter 'direction_sphere_radius') while one (out of two) vector\n"
        doc += "\t    defining the symmetry axis direction is used for the reference embryo\n"
        doc += "\t    (to be used for test purposes)\n"
        doc += "\t - 'symaxis_wrt_symaxis': symmetry direction candidates are used for the floating \n"
        doc += "\t    embryo (parameters 'distribution_sphere_radius' and 'distribution_sigma') while\n"
        doc += "\t    one vector defining the symmetry axis (estimated from the cell name) is used for\n"
        doc += "\t    the reference embryo. Due to computational performances, this option should be\n"
        doc += "\t    considered first. \n"
        self.doc['rotation_initialization'] = doc
        self.rotation_initialization = 'symaxis_wrt_symaxis'

        doc = "\t To get an uniform sampling of the 3d directions in space (options 'sphere_wrt_z'\n"
        doc += "\t and 'sphere_wrt_symaxis' of 'rotation_initialization'), a discrete sphere is build\n"
        doc += "\t and each point of the outer surface gives a sample. The larger the radius, the more\n"
        doc += "\t vectors and the higher computational time\n"
        doc += "\t radius = 2.0: 26 vectors, angle between neighboring vectors in [36.26, 60.0] degrees\n"
        doc += "\t radius = 2.3: 38 vectors, angle between neighboring vectors in [26.57, 54.74] degrees\n"
        doc += "\t radius = 2.5: 54 vectors, angle between neighboring vectors in [24.09, 43.09] degrees\n"
        doc += "\t radius = 2.9: 66 vectors, angle between neighboring vectors in [18.43, 43.09] degrees\n"
        doc += "\t radius = 3.0: 90 vectors, angle between neighboring vectors in [17.72, 43.09] degrees\n"
        doc += "\t radius = 3.5: 98 vectors, angle between neighboring vectors in [15.79, 32.51] degrees\n"
        doc += "\t radius = 3.7: 110 vectors, angle between neighboring vectors in [15.26, 32.51] degrees\n"
        doc += "\t radius = 3.8: 134 vectors, angle between neighboring vectors in [14.76, 29.50] degrees\n"
        doc += "\t radius = 4.0: 222 vectors, angle between neighboring vectors in [10.31, 22.57] degrees\n"
        doc += "\t radius = 5.0: 222 vectors, angle between neighboring vectors in [10.31, 22.57] degrees\n"
        self.doc['direction_sphere_radius'] = doc
        self.direction_sphere_radius = 2.5

        #
        # 2D rotation
        #
        doc = "\t Increment (in degrees) between two successive angles when enumerating\n"
        doc += "\t rotations along the z or the symmetry axis\n"
        self.doc['z_rotation_angle_increment'] = doc
        self.z_rotation_angle_increment = 15

        #
        # save transformations into a file
        #
        doc = "\t File name to save or to read the transformations between the embryos.\n"
        doc += "\t Useful when parsing name choice parameters (for test purposes then).\n"
        self.doc['transformation_filename'] = doc
        self.transformation_filename = None

        #
        #
        #
        doc = "\t The optimal transformation between two embryos is the one giving the smallest\n"
        doc += "\t residual value. The residual value of a transformation is the sum of the smallest\n"
        doc += "\t distances between paired points. This parameter gives the ratio of\n"
        doc += "\t points to be retained\n"
        self.doc['pair_ratio_for_residual'] = doc
        self.pair_ratio_for_residual = 0.8

        #
        #
        #
        doc = "\t Number of processors for parallelization\n"
        self.doc['processors'] = doc
        self.processors = 10

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('#')
        print('# EmbryoRegistrationParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        self.varprint('rotation_initialization', self.rotation_initialization,
                      self.doc.get('rotation_initialization', None))

        self.varprint('direction_sphere_radius', self.direction_sphere_radius,
                      self.doc.get('direction_sphere_radius', None))

        self.varprint('z_rotation_angle_increment', self.z_rotation_angle_increment,
                      self.doc.get('z_rotation_angle_increment', None))

        self.varprint('transformation_filename', self.transformation_filename,
                      self.doc.get('transformation_filename', None))

        self.varprint('pair_ratio_for_residual', self.pair_ratio_for_residual,
                      self.doc.get('pair_ratio_for_residual', None))

        self.varprint('processors', self.processors, self.doc['processors'])
        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# EmbryoRegistrationParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'rotation_initialization', self.rotation_initialization,
                      self.doc.get('rotation_initialization', None))

        self.varwrite(logfile, 'direction_sphere_radius', self.direction_sphere_radius,
                      self.doc.get('direction_sphere_radius', None))

        self.varwrite(logfile, 'z_rotation_angle_increment', self.z_rotation_angle_increment,
                      self.doc.get('z_rotation_angle_increment', None))

        self.varwrite(logfile, 'transformation_filename', self.transformation_filename,
                      self.doc.get('transformation_filename', None))

        self.varwrite(logfile, 'pair_ratio_for_residual', self.pair_ratio_for_residual,
                      self.doc.get('pair_ratio_for_residual', None))

        self.varwrite(logfile, 'processors', self.processors, self.doc['processors'])
        logfile.write("\n")

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            self.write_parameters_in_file(logfile)
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameters):

        self.rotation_initialization = self.read_parameter(parameters, 'rotation_initialization',
                                                           self.rotation_initialization)

        self.direction_sphere_radius = self.read_parameter(parameters, 'direction_sphere_radius',
                                                           self.direction_sphere_radius)

        self.z_rotation_angle_increment = self.read_parameter(parameters, 'z_rotation_angle_increment',
                                                              self.z_rotation_angle_increment)

        self.transformation_filename = self.read_parameter(parameters, 'transformation_filename',
                                                           self.transformation_filename)

        self.pair_ratio_for_residual = self.read_parameter(parameters, 'pair_ratio_for_residual',
                                                           self.pair_ratio_for_residual)

        self.processors = self.read_parameter(parameters, 'processors', self.processors)

    def update_from_parameter_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            monitoring.to_log_and_console("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = common.load_source(parameter_file)
        self.update_from_parameters(parameters)


class InitNamingParameters(EmbryoRegistrationParameters, uatlase.AtlasParameters):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix=None):

        if "doc" not in self.__dict__:
            self.doc = {}

        EmbryoRegistrationParameters.__init__(self, prefix=prefix)
        uatlase.AtlasParameters.__init__(self, prefix=prefix)

        doc = "\t if True, generate python files (prefixed by 'figures_') that generate figures.\n"
        doc += "\t Those files will be saved into the 'outputDir' directory.\n"
        doc += "\t 'generate_figure' can be\n"
        doc += "\t - a boolean value: if True, all figure files are generated; if False, none of them\n"
        doc += "\t - a string: if 'all', all figure files are generated; else, only the specified\n"
        doc += "\t   figure file is generated (see below for the list)\n"
        doc += "\t - a list of strings: if 'all' is in the list, all figure files are generated;\n"
        doc += "\t   else, only the specified figure files are generated (see below for the list)\n"
        doc += "\t List of figures:\n"
        doc += "\t 'success-wrt-atlasnumber':\n"
        doc += "\t    plot leave-one-out results when naming a time points with [1, n-1] atlases\n"
        self.doc['generate_figure'] = doc
        self.generate_figure = False

        doc = "\t Input property file to be named. Must contain lineage, volumes and contact surfaces\n"
        self.doc['inputFile'] = doc
        self.inputFile = None
        doc = "\t Output property file."
        self.doc['outputFile'] = doc
        self.outputFile = None

        #
        doc = "\t Cell number of the developmental stade to be named\n"
        self.doc['cell_number'] = doc
        self.cell_number = 64

        doc = "\t Maximal number of iterations to be done when naming with unanimity rule\n"
        doc += "\t None means until stability is reached.\n"
        self.doc['unanimity_iterations'] = doc
        self.unanimity_iterations = None

        doc = "\t After naming, check whether some names are duplicated, and, if yes, remove them\n"
        self.doc['check_duplicate'] = doc
        self.check_duplicate = True

        doc = "\t After naming, check whether the named cell volume is coherent with its mother/daughters volumes\n"
        self.doc['check_volume'] = doc
        self.check_volume = True

        #
        # for test:
        # names will be deleted, and tried to be rebuilt
        doc = "\t Input property file to be tested (must include cell names).\n"
        doc += "\t If given, 'inputFile' is ignored.\n"
        self.doc['testFile'] = doc
        self.testFile = None

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('#')
        print('# InitNamingParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        EmbryoRegistrationParameters.print_parameters(self)
        uatlase.AtlasParameters.print_parameters(self)

        self.varprint('inputFile', self.inputFile, self.doc.get('inputFile', None))
        self.varprint('outputFile', self.outputFile, self.doc.get('outputFile', None))

        self.varprint('cell_number', self.cell_number, self.doc.get('cell_number', None))
        self.varprint('unanimity_iterations', self.unanimity_iterations, self.doc.get('unanimity_iterations', None))

        self.varprint('check_duplicate', self.check_duplicate, self.doc.get('cell_ncheck_duplicateumber', None))
        self.varprint('check_volume', self.check_volume, self.doc.get('check_volume', None))

        self.varprint('testFile', self.testFile, self.doc.get('testFile', None))
        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# InitNamingParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        EmbryoRegistrationParameters.write_parameters_in_file(self, logfile)
        uatlase.AtlasParameters.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'inputFile', self.inputFile, self.doc.get('inputFile', None))
        self.varwrite(logfile, 'outputFile', self.outputFile, self.doc.get('outputFile', None))

        self.varwrite(logfile, 'cell_number', self.cell_number, self.doc.get('cell_number', None))
        self.varwrite(logfile, 'unanimity_iterations', self.unanimity_iterations, self.doc.get('unanimity_iterations',
                                                                                               None))

        self.varwrite(logfile, 'check_duplicate', self.check_duplicate, self.doc.get('check_duplicate', None))
        self.varwrite(logfile, 'check_volume', self.check_volume, self.doc.get('check_volume', None))

        self.varwrite(logfile, 'testFile', self.testFile, self.doc.get('testFile', None))

        logfile.write("\n")

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            self.write_parameters_in_file(logfile)
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameters):

        EmbryoRegistrationParameters.update_from_parameters(self, parameters)
        uatlase.AtlasParameters.update_from_parameters(self, parameters)

        self.inputFile = self.read_parameter(parameters, 'inputFile', self.inputFile)
        self.outputFile = self.read_parameter(parameters, 'outputFile', self.outputFile)

        self.cell_number = self.read_parameter(parameters, 'cell_number', self.cell_number)
        self.unanimity_iterations = self.read_parameter(parameters, 'unanimity_iterations', self.unanimity_iterations)

        self.check_duplicate = self.read_parameter(parameters, 'check_duplicate', self.check_duplicate)
        self.check_volume = self.read_parameter(parameters, 'check_volume', self.check_volume)

        self.testFile = self.read_parameter(parameters, 'testFile', self.testFile)

    def update_from_parameter_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            monitoring.to_log_and_console("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = common.load_source(parameter_file)
        self.update_from_parameters(parameters)


########################################################################################
#
# sphere discretization
#
########################################################################################

def _build_voxel_sphere(r=3, ir=3):
    # r is radius
    # matrix size will be 2*r + 1 (center) + 2 (margins)
    s = int(2 * ir + 1 + 2)
    m = np.zeros((s, s, s), dtype=np.int8)
    # center is r+1
    c = ir + 1
    # fill the sphere
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(m.shape[2]):
                di = float(i) - float(c)
                dj = float(j) - float(c)
                dk = float(k) - float(c)
                if math.sqrt(di * di + dj * dj + dk * dk) <= r:
                    m[i][j][k] = 1
    # set the border voxel to 2
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(m.shape[2]):
                if m[i][j][k] < 1:
                    continue
                if i > 0 and m[i - 1][j][k] < 1:
                    m[i][j][k] = 2
                    continue
                if i < m.shape[0] - 1 and m[i + 1][j][k] < 1:
                    m[i][j][k] = 2
                    continue
                if j > 0 and m[i][j - 1][k] < 1:
                    m[i][j][k] = 2
                    continue
                if j < m.shape[1] - 1 and m[i][j + 1][k] < 1:
                    m[i][j][k] = 2
                    continue
                if k > 0 and m[i][j][k - 1] < 1:
                    m[i][j][k] = 2
                    continue
                if k < m.shape[2] - 1 and m[i][j][k + 1] < 1:
                    m[i][j][k] = 2
                    continue
    # erase the inner sphere
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(m.shape[2]):
                if m[i][j][k] == 1:
                    m[i][j][k] = 0

    return m


def _build_vector_sphere(r=3):
    # m will have non-zero value at sphere border
    ir = math.ceil(r)
    m = _build_voxel_sphere(r=r, ir=ir)
    # center is ir+1
    c = ir + 1
    xv = []
    yv = []
    zv = []
    angles = []
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(m.shape[2]):
                if m[i][j][k] == 0:
                    continue
                di = float(i) - float(c)
                dj = float(j) - float(c)
                dk = float(k) - float(c)
                norm = math.sqrt(di * di + dj * dj + dk * dk)
                u = np.array([di / norm, dj / norm, dk / norm])
                xv += [u[0]]
                yv += [u[1]]
                zv += [u[2]]
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        for dk in range(-1, 2):
                            if di == 0 and dj == 0 and dk == 0:
                                continue
                            if m[i+di][j+dj][k+dk] == 0:
                                continue
                            v = np.array([(i + di - c), (j + dj - c), (k + dk - c)])
                            angles += [math.acos(np.dot(u, v / np.linalg.norm(v)))]
    v = np.zeros((3, len(xv)))
    v[0, :] = xv
    v[1, :] = yv
    v[2, :] = zv

    monitoring.to_log_and_console("      ... direction distribution build with r = " + str(r), 2)
    monitoring.to_log_and_console("          vectors = " + str(len(xv)), 2)
    min_angles = min(angles)*180.0/np.pi
    max_angles = max(angles) * 180.0 / np.pi
    msg = "          angles between adjacent vectors in [{:.2f}, {:.2f}]".format(min_angles, max_angles)
    monitoring.to_log_and_console(msg, 2)

    return v


########################################################################################
#
# procedures for embryos
#
########################################################################################

def _get_times(prop, n_cells=64, time_digits_for_cell_id=4, mode='floor'):
    #
    # get the time range with the specified number of cells
    # if it does not exist, get the time range with the immediate upper ('ceil') or lower ('floor') number of cells
    #
    proc = "_get_times"
    #
    # count the cells per time
    #
    cells = list(prop.keys())
    cells = sorted(cells)
    cells_per_time = {}
    div = 10 ** time_digits_for_cell_id
    for c in cells:
        t = int(c) // div
        cells_per_time[t] = cells_per_time.get(t, 0) + 1
    #
    # sort the cells_per_time dictionary by time
    # assume that cell number will also be in increasing order
    #
    cells_per_time = dict(sorted(cells_per_time.items(), key=lambda item: item[0]))

    #
    # if there are times with the specified number of cells, return them
    #
    times = [t for t in cells_per_time if cells_per_time[t] == n_cells]
    if len(times) > 0:
        return times, n_cells

    #
    # if not, find times with immediate lower and upper number of cells
    #

    mincells = [cells_per_time[t] for t in cells_per_time if cells_per_time[t] < n_cells]
    maxcells = [cells_per_time[t] for t in cells_per_time if cells_per_time[t] > n_cells]

    if mode == 'floor':
        if len(mincells) > 0:
            n_cells = max(mincells)
            times = [t for t in cells_per_time if cells_per_time[t] == n_cells]
            return times, n_cells
        if len(maxcells) > 0:
            n_cells = min(maxcells)
            times = [t for t in cells_per_time if cells_per_time[t] == n_cells]
            return times, n_cells
    elif mode == 'ceil':
        if len(maxcells) > 0:
            n_cells = min(maxcells)
            times = [t for t in cells_per_time if cells_per_time[t] == n_cells]
            return times, n_cells
        if len(mincells) > 0:
            n_cells = max(mincells)
            times = [t for t in cells_per_time if cells_per_time[t] == n_cells]
            return times, n_cells

    msg = "unknown estimation mode '" + str(mode) + "'"
    monitoring.to_log_and_console(proc + ": " + msg)
    return None, None


def _get_barycenters(embryo, t, transformation=None):
    barycenters = embryo.cell_barycenter
    div = 10 ** embryo.time_digits_for_cell_id
    lc = [c for c in barycenters if (int(c) // div == t) and int(c) % div != 1]
    b = np.zeros((3, len(lc)))
    if transformation is None:
        for i, c in enumerate(lc):
            b[:, i] = barycenters[c]
    else:
        v = np.ones(4)
        for i, c in enumerate(lc):
            v[:3] = barycenters[c]
            b[:, i] = (np.matmul(transformation, v))[:3]
    return b


########################################################################################
#
#
#
########################################################################################

def _coregister_embryos_one_direction(parameters):
    n, vz, i, ni, scale_mat, flo_embryo_bar, ref_embryo_bar, flo_bars, ref_bars, registration_par, verbose = parameters

    if verbose:
        msg = "processing registration[{:3d}/{:3d}] = ".format(i, ni)
        msg += "({:5.2f} {:5.2f} {:5.2f})".format(n[0], n[1], n[2]) + " vs "
        msg += "({:5.2f} {:5.2f} {:5.2f})".format(vz[0], vz[1], vz[2])
        monitoring.to_log_and_console("      ... " + msg)

    # vz normalisation
    vz /= np.linalg.norm(vz)

    #
    rotz = np.identity(4)
    rotxy = np.identity(4)
    sp = np.dot(n, vz)

    #
    # rotation that put n along vz
    #
    if sp > 0.999:
        # n is aligned with vz
        # do nothing
        pass
    elif sp < -0.999:
        # n is aligned with -vz
        # rotation of pi wrt x x vz (if x has the smallest coordinate)
        if abs(vz[0]) <= abs(vz[1]) and abs(vz[0]) <= abs(vz[2]):
            vrot = np.cross(vz, np.array([1, 0, 0]))
        elif abs(vz[1]) <= abs(vz[0]) and abs(vz[1]) <= abs(vz[2]):
            vrot = np.cross(vz, np.array([0, 1, 0]))
        else:
            vrot = np.cross(vz, np.array([0, 0, 1]))
        vrot *= np.pi / np.linalg.norm(vrot)
        rotz[:3, :3] = _rotation_matrix_from_rotation_vector(vrot)
    else:
        # rotation of math.acos(n.z) wrt vector n x z
        # that puts n along z
        vrot = np.cross(n, vz)
        vrot *= math.acos(sp) / np.linalg.norm(vrot)
        rotz[:3, :3] = _rotation_matrix_from_rotation_vector(vrot)

    # for 64 cells, the smallest angle (from the barycenter) between 2 adjacent cells
    # is 7 degrees
    rigid_matrice = []
    affine_matrice = []
    average_residual = []
    angles = np.linspace(0, 2 * np.pi, int(360 / registration_par.z_rotation_angle_increment), endpoint=False)
    for a in angles:
        # compute initial transformation
        rotxy[:3, :3] = _rotation_matrix_from_rotation_vector(a * vz)
        init_mat = np.matmul(rotxy, np.matmul(rotz, scale_mat))
        trs = ref_embryo_bar - np.matmul(init_mat, flo_embryo_bar)
        init_mat[:3, 3] = trs[:3]

        # apply initial transformation
        # compute a rigid transformation
        flo = np.zeros((3, flo_bars.shape[1]))
        v = np.ones(4)
        for j in range(flo_bars.shape[1]):
            v[:3] = flo_bars[:, j]
            flo[:, j] = (np.matmul(init_mat, v))[:3]
        rigid_mat = icp.icp(ref=ref_bars, flo=flo, transformation_type="rigid")
        res_rigid_mat = np.matmul(rigid_mat, init_mat)
        # compute an affine transformation
        v = np.ones(4)
        for j in range(flo_bars.shape[1]):
            v[:3] = flo_bars[:, j]
            flo[:, j] = (np.matmul(res_rigid_mat, v))[:3]
        affine_mat = icp.icp(ref=ref_bars, flo=flo, transformation_type="affine")
        res_affine_mat = np.matmul(affine_mat, res_rigid_mat)

        # compute residuals
        v = np.ones(4)
        for j in range(flo_bars.shape[1]):
            v[:3] = flo_bars[:, j]
            flo[:, j] = (np.matmul(res_affine_mat, v))[:3]
        nbrs = NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(ref_bars.T)
        distances, indices = nbrs.kneighbors(flo.T)
        distances = distances.flatten()
        sort_index = np.argsort(distances)
        # 0.8 is the value used by Gael
        retained_points = int(registration_par.pair_ratio_for_residual * len(distances))
        lr = 0.0
        for j in range(retained_points):
            lr += distances[sort_index[j]]
        lr /= retained_points

        rigid_matrice += [res_rigid_mat]
        affine_matrice += [res_affine_mat]
        average_residual += [lr]

    return rigid_matrice, affine_matrice, average_residual


def _coregister_embryos_one_timepoint(embryo, info_embryo, atlas, info_atlas, parameters, verbose=True):

    # reference embryo
    # ref_barycenter is the barycenter in homogeneous coordinates
    # NB: it is in voxel coordinates
    # ref_barycenters is an np.array of shape (3, n) of the n cell barycenters at time reference_time
    # ref_barycenters[0][i] is the the x coordinate of the ith barycenter

    flo_embryo_barycenter = np.ones(4)
    flo_embryo_barycenter[:3] = info_embryo['barycenter']
    flo_barycenters = _get_barycenters(embryo, info_embryo['first_time'])

    ref_embryo_barycenter = np.ones(4)
    ref_embryo_barycenter[:3] = info_atlas['barycenter']
    ref_barycenters = _get_barycenters(atlas, info_atlas['first_time'])

    # scale to correct for volume differences
    # a voxel size correction of scale = np.cbrt(volref / volflo), applied to the floating points
    # makes the floating embryo volume equal to the reference one
    #
    # since we also have the global correction embryo.get_voxelsize_correction() with
    # volref * (ref.get_voxelsize_correction())^3 = volflo * (flo.get_voxelsize_correction())^3
    # it somes that scale = flo.get_voxelsize_correction() / ref.get_voxelsize_correction()

    scale = embryo.get_voxelsize_correction(info_atlas['first_time']) / \
            atlas.get_voxelsize_correction(info_embryo['first_time'])
    scale_mat = scale * np.identity(4)
    scale_mat[3, 3] = 1.0

    #
    # we could improve the computational time by parallelising the calculation
    # with multiprocessing
    #
    pool = multiprocessing.Pool(processes=parameters.processors)
    mapping = []

    #
    # 2D registrations that put info_atlas['vectors'][:, j] on  info_embryo['vectors'][:, i]
    #
    if True:
        #
        # parallel
        #
        for j in range(info_atlas['vectors'].shape[1]):
            for i in range(info_embryo['vectors'].shape[1]):
                mapping.append((info_embryo['vectors'][:, i], info_atlas['vectors'][:, j],
                                j*info_embryo['vectors'].shape[1]+i+1,
                                info_atlas['vectors'].shape[1]*info_embryo['vectors'].shape[1], scale_mat,
                                flo_embryo_barycenter, ref_embryo_barycenter, flo_barycenters, ref_barycenters,
                                parameters, verbose))

        outputs = pool.map(_coregister_embryos_one_direction, mapping)
        pool.close()
        pool.terminate()
    else:
        #
        # sequential
        #
        outputs = []
        for j in range(info_atlas['vectors'].shape[1]):
            for i in range(info_embryo['vectors'].shape[1]):
                par = (info_embryo['vectors'][:, i], info_atlas['vectors'][:, j],
                       j*info_embryo['vectors'].shape[1]+i+1,
                       info_atlas['vectors'].shape[1]*info_embryo['vectors'].shape[1], scale_mat,
                       flo_embryo_barycenter, ref_embryo_barycenter, flo_barycenters, ref_barycenters,
                       parameters, verbose)
                outputs += [_coregister_embryos_one_direction(par)]

    rigid_matrice = []
    affine_matrice = []
    average_residual = []
    #
    # outputs is a list of tuples
    #
    for rig_mat, aff_mat, lr in outputs:
        rigid_matrice += rig_mat
        affine_matrice += aff_mat
        average_residual += lr

    # residuals study
    index_min = min(range(len(average_residual)), key=average_residual.__getitem__)
    if verbose:
        min_residual = average_residual[index_min]
        ave_residual = statistics.mean(average_residual)
        stddev_residual = statistics.stdev(average_residual)
        msg = "residuals:  min = {:5.2f}".format(min_residual)
        msg += " - average = {:5.2f}".format(ave_residual)
        msg += " (+/- stddev = {:5.2f}".format(stddev_residual) + ")"
        monitoring.to_log_and_console("   ... " + msg)
    return rigid_matrice[index_min], affine_matrice[index_min]


def _coregister_embryos(transformations, embryo, info_embryo, atlas, info_atlas, atlas_name, parameters, verbose=True):
    proc = "_coregister_embryos"

    first_flotime = info_embryo['first_time']
    first_reftime = info_atlas['first_time']

    if verbose:
        msg = "co-registration of time #" + str(first_flotime) + " of embryo to be named"
        msg += " with time #" + str(first_reftime) + " of atlas '" + str(atlas_name) + "'"
        monitoring.to_log_and_console("   ... " + msg)

    #
    #
    #
    if first_flotime not in transformations:
        transformations[first_flotime] = {}
    if atlas_name not in transformations[first_flotime]:
        transformations[first_flotime][atlas_name] = {}

    if first_reftime in transformations[first_flotime][atlas_name]:
        if verbose:
            msg = "transformation time #" + str(first_flotime) + " of embryo to be named"
            msg += " with time #" + str(first_reftime) + " of atlas '" + str(atlas_name) + "'"
            msg += " already exists"
            monitoring.to_log_and_console(proc + ": " + msg)
    else:
        #
        # the rigid transformation can be used to initialize rigid transformations when
        # registering other times of embryo to be named to other times of atlas
        #
        rigid_mat, affine_mat = _coregister_embryos_one_timepoint(embryo, info_embryo, atlas, info_atlas, parameters,
                                                                  verbose=verbose)
        transformations[first_flotime][atlas_name][first_reftime] = affine_mat

    return transformations


#
# register one embryo with respect to several named embryos
#
def _register_embryos(transformations, embryo, atlases, parameters, time_digits_for_cell_id=4, verbose=True):
    proc = "_register_embryos"

    #
    # get the time range of the floating embryo with the specified number of cells
    #
    lowerflotimerange, lowerncells = _get_times(embryo.cell_contact_surface, parameters.cell_number,
                                                time_digits_for_cell_id=time_digits_for_cell_id, mode='floor')
    upperflotimerange, upperncells = _get_times(embryo.cell_contact_surface, parameters.cell_number,
                                                time_digits_for_cell_id=time_digits_for_cell_id, mode='ceil')
    if lowerncells <= parameters.cell_number <= upperncells:
        cell_number = parameters.cell_number
        if cell_number - lowerncells <= upperncells - cell_number:
            flotimerange = lowerflotimerange
        else:
            flotimerange = upperflotimerange
    elif lowerncells == upperncells:
        cell_number = lowerncells
        flotimerange = lowerflotimerange
        msg = "will register embryos at " + str(cell_number) + " instead of " + str(parameters.cell_number) + " cells"
        monitoring.to_log_and_console(proc + ": " + msg)
    else:
        msg = "weird, requested cell number was " + str(parameters.cell_number) + ". "
        msg += "Found cell range was [" + str(lowerncells) + ", " + str(upperncells) + "]"
        monitoring.to_log_and_console(proc + ": " + msg)
        return transformations, [], False

    flo_embryo = {'time_range': flotimerange, 'cell_number': cell_number}
    flo_embryo['first_time'] = statistics.median_low(flo_embryo['time_range'])
    flo_embryo['barycenter'] = embryo.get_embryo_barycenter(flo_embryo['first_time'])

    msg = "   ... embryo to be named: "
    msg += "chosen time is " + str(flo_embryo['first_time']) + " out of "
    msg += str(flo_embryo['time_range']) + "\n"
    msg += "                    number of cells is " + str(flo_embryo['cell_number'])
    monitoring.to_log_and_console(msg)

    #
    # get the time range of the reference embryos with the specified number of cells
    #
    ref_atlases = {}
    for a in atlases.keys():
        lowerreftimerange, lowerncells = _get_times(atlases[a].cell_contact_surface, cell_number,
                                                    time_digits_for_cell_id=time_digits_for_cell_id, mode='floor')
        upperreftimerange, upperncells = _get_times(atlases[a].cell_contact_surface, cell_number,
                                                    time_digits_for_cell_id=time_digits_for_cell_id, mode='ceil')
        if lowerncells <= cell_number <= upperncells:
            if cell_number - lowerncells <= upperncells - cell_number:
                ref_atlases[a] = {'time_range': lowerreftimerange, 'cell_number': lowerncells}
            else:
                ref_atlases[a] = {'time_range': upperreftimerange, 'cell_number': upperncells}
            ref_atlases[a]['first_time'] = statistics.median_low(ref_atlases[a]['time_range'])
            ref_atlases[a]['barycenter'] = atlases[a].get_embryo_barycenter(ref_atlases[a]['first_time'])
            msg = "   ... reference embryo '" + str(a) + "': "
            msg += "chosen time is " + str(ref_atlases[a]['first_time']) + " out of "
            msg += str(ref_atlases[a]['time_range']) + "\n"
            msg += "                    number of cells is " + str(ref_atlases[a]['cell_number'])
            monitoring.to_log_and_console(msg)
            continue
        msg = "requested cell number of " + str(parameters.cell_number) + " is not in the found "
        msg += "cell range of atlas '" + str(a) + "' which is [" + str(lowerncells) + ", " + str(upperncells) + "].\n"
        msg += "\t skip atlas '" + str(a) + "'"
        monitoring.to_log_and_console(proc + ": " + msg)

    #
    # is there something to do ?
    #
    are_transformations_already_computed = True
    if flo_embryo['first_time'] not in transformations:
        are_transformations_already_computed = False
    else:
        for a in ref_atlases:
            if a not in transformations[flo_embryo['first_time']]:
                are_transformations_already_computed = False
                break
            if ref_atlases[a]['first_time'] not in transformations[flo_embryo['first_time']][a]:
                are_transformations_already_computed = False
                break
    if are_transformations_already_computed:
        msg = "   ... co-registrations of embryo to be named already done"
        monitoring.to_log_and_console(msg)
        return transformations, flo_embryo['time_range'], False

    #
    #
    #
    if parameters.rotation_initialization == 'sphere_wrt_z' or \
            parameters.rotation_initialization == 'sphere_wrt_symaxis':
        flo_embryo['vectors'] = _build_vector_sphere(parameters.direction_sphere_radius)
    elif parameters.rotation_initialization == 'symaxis_wrt_symaxis':
        candidates = embryo.get_direction_distribution_candidates(flo_embryo['first_time'], parameters)
        nvectors = len(candidates)
        if isinstance(parameters.maxima_number, int) and 0 < parameters.maxima_number < nvectors:
            nvectors = parameters.maxima_number
        flo_embryo['vectors'] = np.zeros((3, nvectors))
        for i, p in enumerate(candidates):
            if i >= nvectors:
                break
            flo_embryo['vectors'][:, i] = p['vector']
    else:
        msg = "unknown rotation initialization '" + str(parameters.rotation_initialization) + "'"
        monitoring.to_log_and_console(proc + ": " + msg)
        return transformations, flo_embryo['time_range'], False

    if parameters.rotation_initialization == 'sphere_wrt_z':
        for a in ref_atlases:
            ref_atlases[a]['vectors'] = np.zeros((3, 1))
            ref_atlases[a]['vectors'][2, 0] = 1
    elif parameters.rotation_initialization == 'sphere_wrt_symaxis' or \
            parameters.rotation_initialization == 'symaxis_wrt_symaxis':
        for a in ref_atlases:
            ref_atlases[a]['vectors'] = np.zeros((3, 1))
            ref_atlases[a]['vectors'][:, 0] = atlases[a].get_symmetry_axis_from_names(ref_atlases[a]['first_time'])
    else:
        msg = "unknown rotation initialization '" + str(parameters.rotation_initialization) + "'"
        monitoring.to_log_and_console(proc + ": " + msg)
        return transformations, flo_embryo['time_range'], False

    for a in ref_atlases:
        transformations = _coregister_embryos(transformations, embryo, flo_embryo, atlases[a], ref_atlases[a], a,
                                              parameters, verbose=verbose)

    return transformations, flo_embryo['time_range'], True


#########################################################################################
#
# assessment wrt ground-truth
#
########################################################################################

def _get_named_timepoints(prop, time_digits_for_cell_id=4):
    cells = set(list(prop['cell_name'].keys()))
    div = 10 ** time_digits_for_cell_id
    times = list(set([int(c) // div for c in cells]))
    return times


def _test_naming(prop, reference_prop, time_digits_for_cell_id=4, verbose=True):
    #
    #
    #
    times = _get_named_timepoints(prop, time_digits_for_cell_id=time_digits_for_cell_id)
    if len(times) == 0:
        if verbose:
            msg = "no names were set"
            monitoring.to_log_and_console("... " + msg)
    elif len(times) > 1:
        if verbose:
            msg = "weird, names were set on several time points " + str(times)
            monitoring.to_log_and_console("... " + msg)

    cells = set(list(prop['cell_name'].keys()))
    div = 10 ** time_digits_for_cell_id
    times = list(set([int(c) // div for c in cells]))

    refcells = set([c for c in reference_prop['cell_name'] if int(c) // div in times])

    unnamed_cells = refcells.difference(cells)
    unnamed_names = sorted([reference_prop['cell_name'][c] for c in unnamed_cells])
    not_in_reference_cells = cells.difference(refcells)

    common_cells = refcells.intersection(cells)
    right_cells = []
    wrong_cells = []
    for c in common_cells:
        if prop['cell_name'][c] != reference_prop['cell_name'][c]:
            wrong_cells += [c]
            continue
        right_cells += [c]
    wrong_names = [(reference_prop['cell_name'][c], prop['cell_name'][c]) for c in wrong_cells]
    wrong_names = sorted(wrong_names, key=itemgetter(1))

    if verbose:
        msg = "ground-truth cells = " + str(len(refcells)) + " for times " + str(times) + "\n"
        msg += "\t tested cells = " + str(len(cells)) + "\n"
        msg += "\t right retrieved names = {:d}/{:d}\n".format(len(right_cells), len(refcells))
        if len(wrong_cells) > 0:
            msg += "\t ... error in naming = {:d}/{:d}\n".format(len(wrong_cells), len(refcells))
            msg += "\t     reference cells wrongly named (right name, given name): " + str(wrong_names) + "\n"
        if len(unnamed_cells) > 0:
            msg += "\t ... unnamed cells = {:d}/{:d}\n".format(len(unnamed_cells), len(refcells))
            msg += "\t     reference cells unnamed: " + str(unnamed_names) + "\n"
        if len(not_in_reference_cells) > 0:
            msg += "\t ... names not in reference = {:d}/{:d}\n".format(len(not_in_reference_cells), len(refcells))
        monitoring.to_log_and_console("summary" + ": " + msg)

    return len(right_cells), len(wrong_cells), len(unnamed_cells), len(not_in_reference_cells)


def _test_distance(prop, reference_prop, time_digits_for_cell_id=4, verbose=True):
    #
    #
    #
    times = _get_named_timepoints(prop, time_digits_for_cell_id=time_digits_for_cell_id)
    if len(times) == 0:
        if verbose:
            msg = "no names were set"
            monitoring.to_log_and_console("... " + msg)
    elif len(times) > 1:
        if verbose:
            msg = "weird, names were set on several time points " + str(times)
            monitoring.to_log_and_console("... " + msg)

    cells = set(list(prop['cell_name'].keys()))
    div = 10 ** time_digits_for_cell_id
    times = list(set([int(c) // div for c in cells]))

    refcells = set([c for c in reference_prop['cell_name'] if int(c) // div in times])

    common_cells = refcells.intersection(cells)
    right_cells = []
    wrong_cells = []
    for c in common_cells:
        if prop['cell_name'][c] != reference_prop['cell_name'][c]:
            wrong_cells += [c]
            continue
        right_cells += [c]

    keydistance = 'morphonet_float_init_naming_neighborhood_assessment'
    right_distances = []
    for c in right_cells:
        right_distances += [prop[keydistance][c]]
    wrong_distances = []
    for c in wrong_cells:
        wrong_distances += [prop[keydistance][c]]

    return right_distances, wrong_distances

########################################################################################
#
#
#
########################################################################################

def _get_correspondence(flo_time_range, embryo):
    #
    # for the cells of each time where the number of cells is equal to parameters.cell_number
    # gives the cell of the first time where the number of cells is equal to parameters.cell_number
    # it will allows to name only this first time while names can be issued from any time
    #
    correspondence = {}
    flotimesmin = min(flo_time_range)
    lineage = embryo.cell_lineage
    reverse_lineage = {v: k for k, values in lineage.items() for v in values}
    cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
    div = 10 ** embryo.time_digits_for_cell_id
    for c in cells:
        if (int(c) // div) not in flo_time_range or int(c) % div == 1:
            continue
        if int(c) // div == flotimesmin:
            correspondence[c] = c
            continue
        d = c
        while True:
            correspondence[c] = reverse_lineage[d]
            d = reverse_lineage[d]
            if int(d) // div == flotimesmin:
                break
    return correspondence


def _set_candidate_name(name_from_both, name_from_flo, name_from_ref, embryo, reference, tflo, tref,
                        transformation, time_digits_for_cell_id=4):

    # get cell ids of both the reference and floating cells
    div = 10 ** time_digits_for_cell_id
    ref_barycenters = reference.cell_barycenter
    cell_tref = [c for c in ref_barycenters if (int(c) // div == tref) and int(c) % div != 1]
    flo_barycenters = embryo.cell_barycenter
    cell_tflo = [c for c in flo_barycenters if (int(c) // div == tflo) and int(c) % div != 1]

    # get reference barycenter coordinates
    bref = np.zeros((3, len(cell_tref)))
    for i, c in enumerate(cell_tref):
        bref[:, i] = ref_barycenters[c]

    # get floating barycenter coordinates
    bflo = np.zeros((3, len(cell_tflo)))
    v = np.ones(4)
    for i, c in enumerate(cell_tflo):
        v[:3] = flo_barycenters[c]
        bflo[:, i] = (np.matmul(transformation, v))[:3]

    # find reference barycenters that are closest to floating ones
    # i_flo_to_ref[i] is the index (in the list) of the reference cell closest
    # to the ith floating cell in cell_tflo list
    nghbref = NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(bref.T)
    d_flo_to_ref, i_flo_to_ref = nghbref.kneighbors(bflo.T)
    i_flo_to_ref = i_flo_to_ref.flatten()

    # find floating barycenters that are closest to reference ones
    nghbflo = NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(bflo.T)
    d_ref_to_flo, i_ref_to_flo = nghbflo.kneighbors(bref.T)
    i_ref_to_flo = i_ref_to_flo.flatten()

    # find cell/barycenters which are closest to each other
    ref_name = reference.cell_name
    reciprocal_n = 0
    # i : index of cell c in cell_tflo
    for i, cf in enumerate(cell_tflo):
        j = i_flo_to_ref[i]
        cr = cell_tref[j]
        if cr not in ref_name:
            continue
        if i_ref_to_flo[j] == i:
            name_from_both[cf] = name_from_both.get(cf, []) + [ref_name[cr]]
            reciprocal_n += 1
        else:
            name_from_flo[cf] = name_from_flo.get(cf, []) + [ref_name[cr]]

    for j, cr in enumerate(cell_tref):
        i = i_ref_to_flo[j]
        cf = cell_tflo[i]
        if cr not in ref_name:
            continue
        if i_flo_to_ref[i] == j:
            pass
        else:
            name_from_ref[cf] = name_from_ref.get(cf, []) + [ref_name[cr]]

    return name_from_both, name_from_flo, name_from_ref, reciprocal_n


def _get_floating_names(transformations, embryo, atlases, time_digits_for_cell_id=4, verbose=True):
    #
    # correspondences
    # build dictionaries indexed by cell ids (from embryo to be named)
    # each dictionary entry being a list of names issued  from atlases
    #

    name_both = {}
    name_flo = {}
    name_ref = {}

    ref_atlases = atlases.get_atlases()
    # parse time points of the floating embryo
    for tflo in transformations:
        # parse atlases of the reference atlases
        # there may be more references in the transformation dictionary
        for a in ref_atlases:
            if a not in transformations[tflo]:
                msg = "weird, '" + str(a) + "' is not in computed transformations at time #" + str(tflo)
                monitoring.to_log_and_console("      ... " + msg, 2)
                continue
            for tref in transformations[tflo][a]:
                output = _set_candidate_name(name_both, name_flo, name_ref, embryo, ref_atlases[a],
                                             tflo, tref, transformations[tflo][a][tref],
                                             time_digits_for_cell_id=time_digits_for_cell_id)
                name_both, name_flo, name_ref, n = output
                if verbose:
                    msg = "set " + str(n) + " names from time #" + str(tref) + " of '" + str(a) + "'"
                    msg += " on time #" + str(tflo)
                    monitoring.to_log_and_console("      ... " + msg, 2)

    return name_both, name_flo, name_ref


def _get_unanimity_names(prop, name):
    cells = sorted(list(name.keys()))
    n = 0
    for c in cells:
        if len(set(name[c])) == 1:
            prop['cell_name'][c] = name[c][0]
            del name[c]
            n += 1
    return prop, name, n


def _iterate_unanimity_names(prop, name_both, parameters, verbose=False):
    iteration = 0
    prop['cell_name'] = {}
    while True:
        prop, name_both, n = _get_unanimity_names(prop, name_both)
        if n == 0:
            break
        if verbose:
            monitoring.to_log_and_console("   ... found " + str(n) + " consensual names")
        used_names = list(prop['cell_name'].values())
        for c in name_both:
            l = [n for n in name_both[c] if n not in used_names]
            name_both[c] = l
        iteration += 1
        if isinstance(parameters.unanimity_iterations, int) and iteration >= parameters.unanimity_iterations:
            break

    return prop


########################################################################################
#
#
#
########################################################################################

def _remove_duplicate(prop, embryo, verbose=True):
    proc = "_remove_duplicate"

    if 'cell_name' not in prop:
        monitoring.to_log_and_console(str(proc) + ": 'cell_name' was not in dictionary")
        return prop

    times = _get_named_timepoints(prop, time_digits_for_cell_id=embryo.time_digits_for_cell_id)
    if len(times) == 0:
        msg = "no names were set"
        monitoring.to_log_and_console("... " + msg)
        return prop
    elif len(times) > 1:
        msg = "weird, names were set on several time points " + str(times)
        monitoring.to_log_and_console("... " + msg)
        return prop

    cell_per_name = {}
    for c in prop['cell_name']:
        cell_per_name[prop['cell_name'][c]] = cell_per_name.get(prop['cell_name'][c], []) + [c]

    duplicates = [n for n in cell_per_name if len(cell_per_name[n]) > 1]
    if len(duplicates) > 0:
        for n in duplicates:
            for c in cell_per_name[n]:
                del prop['cell_name'][c]
            if verbose:
                msg = "   ... name " + str(n) + " was given " + str(len(cell_per_name[n])) + "  times"
                monitoring.to_log_and_console(msg)

    return prop


########################################################################################
#
#
#
########################################################################################

def _get_transformation_filename(parameters):
    if parameters.transformation_filename is None:
        return None
    if isinstance(parameters.transformation_filename, str):
        if parameters.transformation_filename.endswith(".pkl") is True:
            return parameters.transformation_filename
        return parameters.transformation_filename + '.pkl'
    return None


########################################################################################
#
#
#
########################################################################################

def _get_atlases_branch_volume(name, atlases):
    ref_atlases = atlases.get_atlases()
    vols = []
    for a in ref_atlases:
        ref_name = ref_atlases[a].cell_name
        ref_volume = ref_atlases[a].cell_volume
        for c in ref_name:
            if ref_name[c] != name:
                continue
            if c not in ref_volume:
                msg = "weird, cell " + str(c) + " is not in volumes of atlas '" + str(a) + "'"
                monitoring.to_log_and_console("      ... " + msg)
                continue
            #
            # get volume correction for the curr]ent time
            # to get homogeneized volumes
            #
            vols += [ref_atlases[a].get_rectified_cell_volume(c)]
    return vols


def _get_embryo_branch_volume(embryo, cell, lineage, reverse_lineage):
    vols = []
    c = cell
    cell_volume = embryo.cell_volume
    while True:
        if c in cell_volume:
            vols += [embryo.get_rectified_cell_volume(c)]
        if c not in lineage:
            break
        if len(lineage[c]) != 1:
            break
        c = lineage[c][0]

    c = cell
    while True:
        if c not in reverse_lineage:
            break
        c = reverse_lineage[c]
        if len(lineage[c]) != 1:
            break
        if c in cell_volume:
            vols += [embryo.get_rectified_cell_volume(c)]
    return vols


def _check_initial_naming_volume(prop, embryo, atlases, time_digits_for_cell_id=4, verbose=True):

    times = _get_named_timepoints(prop, time_digits_for_cell_id=time_digits_for_cell_id)
    if len(times) == 0:
        if verbose:
            msg = "no names were set"
            monitoring.to_log_and_console("... " + msg)
    elif len(times) > 1:
        if verbose:
            msg = "weird, names were set on several time points " + str(times)
            monitoring.to_log_and_console("... " + msg)

    lineage = embryo.cell_lineage
    reverse_lineage = {v: k for k, values in lineage.items() for v in values}

    cells = sorted(list(prop['cell_name'].keys()))
    unnamed_cells = []
    for c in cells:
        branch_volume = _get_embryo_branch_volume(embryo, c, lineage, reverse_lineage)
        cell_volume = np.mean(branch_volume)
        name = prop['cell_name'][c]

        #
        # get volumes of cell with the same name in reference
        #
        volumes = _get_atlases_branch_volume(name, atlases)
        mean_volume = np.mean(volumes)
        std_volume = np.std(volumes)
        if cell_volume >= mean_volume:
            #
            # get volumes of the mother cell in reference
            #
            mname = uname.get_mother_name(name)
            m_volumes = _get_atlases_branch_volume(mname, atlases)
            if len(m_volumes) == 0:
                if verbose:
                    msg = "   ... can not check " + str(c) + " (given name was '" + str(name) + "')."
                    msg += " No mother cell volumes."
                    monitoring.to_log_and_console(msg)
                continue
            m_mean_volume = np.mean(m_volumes)
            if cell_volume >= m_mean_volume:
                unnamed_cells += [(c, name)]
                del prop['cell_name'][c]
                if verbose:
                    msg = "   ... unnamed cell " + str(c) + " (given name was '" + str(name) + "').\n"
                    msg += "       Cell volume larger than reference mother cell volumes."
                    monitoring.to_log_and_console(msg)
                continue
            #   emb_volume[c] < m_mean_volume:
            if len(m_volumes) == 1:
                diff_to_cell = (cell_volume - mean_volume)
                diff_to_mcell = (m_mean_volume - cell_volume)
            else:
                m_std_volume = np.std(m_volumes)
                if std_volume > 0 and m_std_volume > 0:
                    diff_to_cell = (cell_volume - mean_volume) / std_volume
                    diff_to_mcell = (m_mean_volume - cell_volume) / m_std_volume
                else:
                    diff_to_cell = (cell_volume - mean_volume)
                    diff_to_mcell = (m_mean_volume - cell_volume)
            if diff_to_cell >= diff_to_mcell:
                unnamed_cells += [(c, name)]
                del prop['cell_name'][c]
                if verbose:
                    msg = "   ... unnamed cell " + str(c) + " (given name was '" + str(name) + "').\n"
                    msg += "       Cell volume closer to reference mother cell volumes."
                    monitoring.to_log_and_console(msg)
                continue
        else:
            dnames = uname.get_daughter_names(name)
            d_volumes = [0, 0]
            d_mean_volume = [0, 0]
            d_volumes[0] = _get_atlases_branch_volume(dnames[0], atlases)
            d_volumes[1] = _get_atlases_branch_volume(dnames[1], atlases)
            if len(d_volumes[0]) == 0 and len(d_volumes[1]) == 0:
                if verbose:
                    msg = "   ... can not check " + str(c) + " (given name was '" + str(name) + "')."
                    msg += " No daughter cell volumes."
                    monitoring.to_log_and_console(msg)
                continue
            elif len(d_volumes[0]) == 0:
                d_mean_volume[1] = np.mean(d_volumes[1])
                idaughter = 1
            elif len(d_volumes[1]) == 0:
                d_mean_volume[0] = np.mean(d_volumes[0])
                idaughter = 0
            else:
                d_mean_volume = [0, 0]
                d_mean_volume[0] = np.mean(d_volumes[0])
                d_mean_volume[1] = np.mean(d_volumes[1])
                if d_mean_volume[1] > d_mean_volume[0]:
                    idaughter = 1
                else:
                    idaughter = 0
            if cell_volume <= d_mean_volume[idaughter]:
                unnamed_cells += [(c, name)]
                del prop['cell_name'][c]
                if verbose:
                    msg = "   ... unnamed cell " + str(c) + " (given name was '" + str(name) + "').\n"
                    msg += "       Cell volume smaller than the largest reference daughter ('"
                    msg += str(dnames[idaughter]) + "') cell volumes."
                    monitoring.to_log_and_console(msg)
                continue
            # emb_volume[c] > d_mean_volume[idaughter]
            if len(d_volumes[idaughter]) == 1:
                diff_to_cell = (mean_volume - cell_volume)
                diff_to_dcell = (cell_volume - d_mean_volume[idaughter])
            else:
                d_std_volume = np.std(d_volumes[idaughter])
                if std_volume > 0 and d_std_volume > 0:
                    diff_to_cell = (mean_volume - cell_volume) / std_volume
                    diff_to_dcell = (cell_volume - d_mean_volume[idaughter]) / d_std_volume
                else:
                    diff_to_cell = (mean_volume - cell_volume)
                    diff_to_dcell = (cell_volume - d_mean_volume[idaughter])
            if diff_to_cell >= diff_to_dcell:
                unnamed_cells += [(c, name)]
                del prop['cell_name'][c]
                if verbose:
                    msg = "   ... unnamed cell " + str(c) + " (given name was '" + str(name) + "').\n"
                    msg += "       Cell volume closer to the largest reference daughter ('"
                    msg += str(dnames[idaughter]) + "') cell volumes."
                    monitoring.to_log_and_console(msg)
                continue
    if verbose:
        msg = "   unname " + str(len(unnamed_cells)) + " cells"
        monitoring.to_log_and_console(msg)
    return prop


#######################################################################################
#
# leave-one-out tests
#
########################################################################################

def _get_success(atlases, parameters, time_digits_for_cell_id=4, verbose=True):
    monitoring.to_log_and_console(" - compute transformations ")
    #
    # compute transformations
    #
    all_atlases = atlases.get_atlases()
    transformations = {}
    flo_time_range = {}
    transformation_filename = _get_transformation_filename(parameters)
    if isinstance(transformation_filename, str) and os.path.isfile(transformation_filename):
        inputfile = open(transformation_filename, 'rb')
        transformations = pkl.load(inputfile)
        inputfile.close()
    for a in all_atlases:
        #
        # compute transformations with a as the floating embryo
        #
        tmp_atlases = copy.deepcopy(all_atlases)
        del tmp_atlases[a]
        if a not in transformations:
            transformations[a] = {}
        output = _register_embryos(transformations[a], all_atlases[a], tmp_atlases, parameters,
                                   time_digits_for_cell_id=time_digits_for_cell_id, verbose=verbose)
        transformations[a], flo_time_range[a], transformations_have_been_computed = output
        if isinstance(transformation_filename, str) and transformations_have_been_computed:
            outputfile = open(transformation_filename, 'wb')
            pkl.dump(transformations, outputfile)
            outputfile.close()

    #
    # set names
    #
    right_naming = {}
    wrong_naming = {}
    unnamed_naming = {}
    right_distance = {}
    wrong_distance = {}
    for a in all_atlases:
        # atlas to be tested
        ref_test_atlas = copy.deepcopy(all_atlases[a])
        test_atlas = copy.deepcopy(ref_test_atlas)
        test_atlas.del_property('cell_name')

        # reference atlas list
        ref_atlas_list = list(all_atlases.keys())
        ref_atlas_list.remove(a)

        #
        # pick combinations of reference atlas
        #
        monitoring.to_log_and_console("*** compare " + str(a) + " versus " + str(ref_atlas_list))
        for i in range(1, len(ref_atlas_list) + 1):
            subset_atlas = list(itertools.combinations(ref_atlas_list, i))
            msg = " ** compare " + str(a)
            msg += " versus " + str(len(subset_atlas)) + " combinations of " + str(i) + " atlases"
            monitoring.to_log_and_console(msg)
            for j, s in enumerate(subset_atlas):
                msg = "  * " + str(j+1) + "/" + str(len(subset_atlas)) + " compare " + str(a) + " versus " + str(s)
                monitoring.to_log_and_console(msg)
                #
                # build references atlases
                #
                tmp_atlases = uatlasd.DivisionAtlases(parameters=parameters)
                tmp = {}
                for b in s:
                    tmp[b] = all_atlases[b]
                tmp_atlases.set_atlases(tmp)
                #
                # collect names from each co-registration
                #
                output = _get_floating_names(transformations[a], test_atlas, tmp_atlases,
                                             time_digits_for_cell_id=time_digits_for_cell_id, verbose=verbose)
                name_both, name_flo, name_ref = output

                #
                # get consensual names
                #
                prop = copy.deepcopy(test_atlas.get_property())
                prop = _iterate_unanimity_names(prop, name_both, parameters, verbose=verbose)

                #
                # eventually remove some names
                #
                if parameters.check_duplicate:
                    if verbose:
                        monitoring.to_log_and_console("... removing duplicates")
                    prop = _remove_duplicate(prop, test_atlas, verbose=verbose)

                if parameters.check_volume:
                    if verbose:
                        monitoring.to_log_and_console("... checking names wrt volumes")
                    prop = _check_initial_naming_volume(prop, test_atlas, tmp_atlases,
                                                        time_digits_for_cell_id=time_digits_for_cell_id,
                                                        verbose=verbose)
                #
                #
                #
                prop = _evaluate_initial_naming(prop, transformations[a], test_atlas, tmp_atlases)

                #
                #
                #
                right, wrong, unnamed, not_in = _test_naming(prop, ref_test_atlas.get_property(),
                                                             time_digits_for_cell_id=time_digits_for_cell_id,
                                                             verbose=False)
                right_naming[i] = right_naming.get(i, []) + [right]
                wrong_naming[i] = wrong_naming.get(i, []) + [wrong]
                unnamed_naming[i] = unnamed_naming.get(i, []) + [unnamed]

                right, wrong = _test_distance(prop, ref_test_atlas.get_property(),
                                              time_digits_for_cell_id=time_digits_for_cell_id, verbose=False)
                right_distance[i] = right_distance.get(i, []) + right
                wrong_distance[i] = wrong_distance.get(i, []) + wrong


    # print("right_naming = " + str(right_naming))
    # print("wrong_naming = " + str(wrong_naming))
    return right_naming, wrong_naming, unnamed_naming, right_distance, wrong_distance


def _print_dict_of_list(f, d, name):
    f.write(str(name) + " = ")
    f.write("{")
    ks = list(d.keys())
    for j, k in enumerate(ks):
        f.write(str(k) + ": [")
        for i, e in enumerate(d[k]):
            f.write("{:.3f}".format(e))
            if i+1 < len(d[k]):
                f.write(", ")
        f.write("]")
        if j + 1 < len(ks):
            f.write(", ")
    f.write("}")
    f.write("\n")


def _figure_atlas_init_naming(atlases, parameters, time_digits_for_cell_id=4):
    proc = "_figure_atlas_init_naming"
    output = _get_success(atlases, parameters, time_digits_for_cell_id=time_digits_for_cell_id, verbose=False)
    right_naming, wrong_naming, unnamed_naming, right_distance, wrong_distance = output

    filename = 'figure_atlas_init_naming'
    file_suffix = None
    if parameters.figurefile_suffix is not None and isinstance(parameters.figurefile_suffix, str) and \
            len(parameters.figurefile_suffix) > 0:
        file_suffix = '_' + parameters.figurefile_suffix
    if file_suffix is not None:
        filename += file_suffix
    filename += '.py'

    if parameters.outputDir is not None and isinstance(parameters.outputDir, str):
        if not os.path.isdir(parameters.outputDir):
            if not os.path.exists(parameters.outputDir):
                os.makedirs(parameters.outputDir)
            else:
                monitoring.to_log_and_console(proc + ": '" + str(parameters.outputDir) + "' is not a directory ?!")
        if os.path.isdir(parameters.outputDir):
            filename = os.path.join(parameters.outputDir, filename)
    #
    #
    #
    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    f.write("right_naming = " + str(right_naming) + "\n")
    f.write("wrong_naming = " + str(wrong_naming) + "\n")
    f.write("unnamed_naming = " + str(unnamed_naming) + "\n")
    f.write("labels = sorted(list(right_naming.keys()))\n")

    f.write("\n")
    _print_dict_of_list(f, right_distance, "right_distance")
    _print_dict_of_list(f, wrong_distance, "wrong_distance")

    f.write("\n")
    f.write("right_list = []\n")
    f.write("wrong_list = []\n")
    f.write("unnamed_list = []\n")
    f.write("for i, l in enumerate(labels):\n")
    f.write("    right_list += [right_naming[l]]\n")
    f.write("    wrong_list += [wrong_naming[l]]\n")
    f.write("    unnamed_list += [unnamed_naming.get(l, [])]\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(16, 6.5))\n")
    f.write("rbox = ax.boxplot(right_list, patch_artist=False)\n")
    f.write("wbox = ax.boxplot(wrong_list, patch_artist=False)\n")

    f.write("\n")
    f.write("rmed = [b.get_ydata()[0] for b in rbox['medians']]\n")
    f.write("wmed = [b.get_ydata()[0] for b in wbox['medians']]\n")
    f.write("ax.plot(labels, rmed, color='blue', label='right names')\n")
    f.write("ax.plot(labels, wmed, color='cyan', label='wrong names')\n")

    f.write("\n")
    f.write("ax.yaxis.grid(True)\n")
    f.write("ax.set_xticks([y + 1 for y in range(len(right_naming))])\n")
    f.write("ax.set_xlabel('Number of atlases')\n")
    f.write("ax.set_ylabel('#cells')\n")
    f.write("ax.legend()\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('right_wrong_wrt_atlasnumber")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")
    f.write("\n")

    f.write("fig, ax = plt.subplots(figsize=(16, 6.5))\n")
    f.write("rbox = ax.boxplot(right_list, patch_artist=False)\n")
    f.write("ubox = ax.boxplot(unnamed_list, patch_artist=False)\n")

    f.write("\n")
    f.write("rmed = [b.get_ydata()[0] for b in rbox['medians']]\n")
    f.write("umed = [b.get_ydata()[0] for b in ubox['medians']]\n")
    f.write("ax.plot(labels, rmed, color='blue', label='right names')\n")
    f.write("ax.plot(labels, umed, color='cyan', label='unnamed cells')\n")

    f.write("\n")
    f.write("ax.yaxis.grid(True)\n")
    f.write("ax.set_xticks([y + 1 for y in range(len(right_naming))])\n")
    f.write("ax.set_xlabel('Number of atlases')\n")
    f.write("ax.set_ylabel('#cells')\n")
    f.write("ax.legend()\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('right_unamed_wrt_atlasnumber")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    f.write("fig, axs = plt.subplots(2, 3, sharex=True, figsize=(16, 6.5))\n")
    f.write("for i in right_distance:\n")
    f.write("    y = (i - 1) // 3\n")
    f.write("    x = (i - 1) % 3\n")
    f.write("    if i in wrong_distance:\n")
    f.write("        labels = ['right name', 'wrong name']\n")
    f.write("        axs[y, x].hist([right_distance[i], wrong_distance[i]], 25, histtype='bar', label=labels)\n")
    f.write("    axs[y, x].set_title('#atlases = ' + str(i))\n")
    f.write("    axs[y, x].set(xlim=(0, 1))\n")
    f.write("fig.suptitle('Assessment score (1 - neighborhood distance)', fontsize=16)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('right_wrong_distances")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.close()


def generate_figure(atlases, parameters, time_digits_for_cell_id=4):
    do_generate_figure = (isinstance(parameters.generate_figure, bool) and parameters.generate_figure) or \
                         (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'all') or \
                         (isinstance(parameters.generate_figure, list) and 'all' in parameters.generate_figure)

    #
    # plot cell number wrt time without and with temporal registration
    #
    if (isinstance(parameters.generate_figure, str) and
        parameters.generate_figure == 'atlas-init-naming-leave-one-out') \
            or (isinstance(parameters.generate_figure, list)
                and 'atlas-init-naming-leave-one-out' in parameters.generate_figure) \
            or do_generate_figure:
        monitoring.to_log_and_console("... generate atlas init naming leave one out figure", 1)
        _figure_atlas_init_naming(atlases, parameters, time_digits_for_cell_id=time_digits_for_cell_id)
        monitoring.to_log_and_console("... done", 1)


########################################################################################
#
#
#
########################################################################################

def _evaluate_initial_naming(prop, transformations, embryo, atlases):
    proc = "_evaluate_initial_naming"

    if 'cell_name' not in prop:
        monitoring.to_log_and_console(str(proc) + ": 'cell_name' was not in dictionary")
        return prop

    times = _get_named_timepoints(prop, time_digits_for_cell_id=embryo.time_digits_for_cell_id)
    if len(times) == 0:
        msg = "no names were set"
        monitoring.to_log_and_console("... " + msg)
        return prop
    elif len(times) > 1:
        msg = "weird, names were set on several time points " + str(times)
        monitoring.to_log_and_console("... " + msg)
    elif len(times) == 1 and times[0] != list(transformations.keys())[0]:
        msg = "named time point is " + str(times[0]) + " while transformations were computed at time "
        msg += str(list(transformations.keys()))
        monitoring.to_log_and_console("... " + msg)

    ref_atlases = atlases.get_atlases()
    keydistance = 'morphonet_float_init_naming_neighborhood_assessment'
    keyassessment = 'morphonet_float_init_naming_neighborhood_assessment_quality'
    prop[keydistance] = {}
    prop[keyassessment] = {}
    div = 10 ** embryo.time_digits_for_cell_id

    cell_contact_surface = embryo.cell_contact_surface
    cells = sorted(list(prop['cell_name'].keys()))
    for c in cells:
        if c not in cell_contact_surface:
            msg = "weird, cell '" + str(c) + "' is not in 'cell_contact_surface' dictionary "
            monitoring.to_log_and_console("... " + msg)
            continue
        #
        # get the cell neighborhood
        # normalize it
        # replace cell id with names
        #
        t_flo = c // div
        if t_flo not in transformations:
            msg = "weird, time point '" + str(t_flo) + "' is not in transformation dictionary "
            monitoring.to_log_and_console("... " + msg)
            continue
        neighborhood = {'foo': {}}
        neighborhood['foo']['embryo'] = copy.deepcopy(cell_contact_surface[c])
        neighborhood['foo']['embryo'] = uatlasc.neighborhood_normalization(neighborhood['foo']['embryo'], embryo,
                                                                           t_flo, cell_normalization='global')
        surface = 0.0
        named_surface = 0.0
        cellids = list(neighborhood['foo']['embryo'].keys())
        for d in cellids:
            surface += neighborhood['foo']['embryo'][d]
            if int(d) % div == 1 or int(d) % div == 0:
                named_surface += neighborhood['foo']['embryo'][d]
                neighborhood['foo']['embryo']['background'] = neighborhood['foo']['embryo'][d]
                del neighborhood['foo']['embryo'][d]
            elif d in prop['cell_name']:
                named_surface += neighborhood['foo']['embryo'][d]
                neighborhood['foo']['embryo'][prop['cell_name'][d]] = neighborhood['foo']['embryo'][d]
                del neighborhood['foo']['embryo'][d]
            else:
                #
                # we kept the cell id as key
                # it may be problematic if the same cell id correspond to an unnamed cell
                # in the neighborhood of the cell with the same name in an atlas
                #
                pass

        prop[keyassessment][c] = named_surface / surface

        #
        # get other neighborhoods
        # transformations[t_flo][atlas][t_ref] is the affine transformation that maps the floating frame
        # (at t_flo) into the reference one (at t_ref)
        #
        for a in ref_atlases:
            if a not in transformations[t_flo]:
                msg = "weird, atlas '" + str(a) + "' is not in transformation dictionary at time point " + str(t_flo)
                monitoring.to_log_and_console("... " + msg)
                continue
            t_ref = list(transformations[t_flo][a].keys())[0]
            #
            # find cells at t_ref in atlas 'a' with the same name than c
            #
            adiv = 10 ** ref_atlases[a].time_digits_for_cell_id
            acell = [ac for ac in ref_atlases[a].cell_name if ac // adiv == t_ref and ref_atlases[a].cell_name[ac] ==
                     prop['cell_name'][c]]
            if len(acell) == 0:
                # msg = "weird, atlas '" + str(a) + "' has no cell named " + str(prop['cell_name'][c])
                # msg += " at time point " + str(t_ref)
                # monitoring.to_log_and_console("... " + msg)
                continue
            elif len(acell) > 1:
                msg = "weird, atlas '" + str(a) + "' has several cells named " + str(prop['cell_name'][c])
                msg += " at time point " + str(t_ref)
                monitoring.to_log_and_console("... " + msg)
                continue
            neighborhood['foo'][a] = copy.deepcopy(ref_atlases[a].cell_contact_surface[acell[0]])
            neighborhood['foo'][a] = uatlasc.neighborhood_normalization(neighborhood['foo'][a], ref_atlases[a],
                                                                        t_ref, cell_normalization='global')
            cellids = list(neighborhood['foo'][a].keys())
            for d in cellids:
                if int(d) % adiv == 1 or int(d) % adiv == 0:
                    neighborhood['foo'][a]['background'] = neighborhood['foo'][a][d]
                    del neighborhood['foo'][a][d]
                elif d in ref_atlases[a].cell_name:
                    neighborhood['foo'][a][ref_atlases[a].cell_name[d]] = neighborhood['foo'][a][d]
                    del neighborhood['foo'][a][d]
                else:
                    #
                    # we kept the cell id as key
                    # it may be problematic if the same cell id correspond to an unnamed cell
                    # in the neighborhood of the cell from the embryo to be tested
                    #
                    pass
        #
        # for named cell c, we built neighborhood['foo'][embryo, a in atlases]
        # which is a dictionary/vector of contact surfaces whose key are cell names (when it is possible)
        # build same contact surfaces
        #
        neighborhood = uneighborhood.build_same_contact_surfaces(neighborhood, ['foo'])
        distances = []
        for a in neighborhood['foo']:
            if a == 'embryo':
                continue
            nm, n0, n1 = uneighborhood.cell_distance_elements(neighborhood['foo']['embryo'], neighborhood['foo'][a])
            distances += [1.0 - nm / (n0 + n1)]
        if len(distances) == 0:
            msg = "weird, no atlases are found to name cell " + str(c) + " as " + str(prop['cell_name'][c])
            monitoring.to_log_and_console("... " + msg)
        else:
            prop[keydistance][c] = np.mean(distances)
    return prop


########################################################################################
#
#
#
########################################################################################

def _initial_naming(prop, embryo, atlases, parameters, time_digits_for_cell_id=4, verbose=True):
    proc = "_initial_naming"

    properties = embryo.property_list()
    if 'cell_lineage' not in properties:
        monitoring.to_log_and_console(str(proc) + ": 'cell_lineage' was not in dictionary")
        return None

    if 'cell_contact_surface' not in properties:
        monitoring.to_log_and_console(str(proc) + ": 'cell_contact_surface' was not in dictionary")
        return None

    if 'cell_volume' not in properties:
        monitoring.to_log_and_console(str(proc) + ": 'cell_volume' was not in dictionary")
        return None

    if 'cell_barycenter' not in properties:
        monitoring.to_log_and_console(str(proc) + ": 'cell_barycenter' was not in dictionary")
        return None

    #
    # compute transformations between the (floating) embryo to be named and (reference) embryos
    # transformations[t_flo][atlas][t_ref] is the affine transformation that maps the floating frame
    # (at t_flo) into the reference one (at t_ref)
    #
    if verbose:
        monitoring.to_log_and_console("... co-registrations")

    ref_atlases = atlases.get_atlases()
    transformations = {}

    transformation_filename = _get_transformation_filename(parameters)
    if isinstance(transformation_filename, str) and os.path.isfile(transformation_filename):
        inputfile = open(transformation_filename, 'rb')
        transformations = pkl.load(inputfile)
        inputfile.close()

    output = _register_embryos(transformations, embryo, ref_atlases, parameters,
                               time_digits_for_cell_id=time_digits_for_cell_id, verbose=verbose)
    transformations, flo_time_range, transformations_have_been_computed = output

    if isinstance(transformation_filename, str) and transformations_have_been_computed:
        outputfile = open(transformation_filename, 'wb')
        pkl.dump(transformations, outputfile)
        outputfile.close()

    if verbose:
        monitoring.to_log_and_console("... setting names")

    #
    # collect names from each co-registration
    # name_both, name_flo, name_ref are dictionaries indexed by cell ids
    # name_both
    output = _get_floating_names(transformations, embryo, atlases,
                                 time_digits_for_cell_id=time_digits_for_cell_id, verbose=verbose)
    name_both, name_flo, name_ref = output

    #
    # get consensual names
    #
    prop = _iterate_unanimity_names(prop, name_both, parameters, verbose=verbose)

    #
    # eventually remove some names
    #
    if parameters.check_duplicate:
        monitoring.to_log_and_console("... removing duplicates")
        prop = _remove_duplicate(prop, embryo)

    if parameters.check_volume:
        monitoring.to_log_and_console("... checking names wrt volumes")
        prop = _check_initial_naming_volume(prop, embryo, atlases, time_digits_for_cell_id=time_digits_for_cell_id,
                                            verbose=True)
    #
    #
    #
    prop = _evaluate_initial_naming(prop, transformations, embryo, atlases)

    return prop


########################################################################################
#
#
#
########################################################################################

def init_naming_process(experiment, parameters):
    proc = "init_naming_process"

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, InitNamingParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()

    uatlasd.monitoring.copy(monitoring)
    uatlasc.monitoring.copy(monitoring)
    icp.monitoring.copy(monitoring)

    #
    # should we clean reference here?
    #
    atlases = uatlasd.DivisionAtlases(parameters=parameters)
    atlases.add_atlases(parameters.atlasFiles, parameters, time_digits_for_cell_id=time_digits_for_cell_id)
    generate_figure(atlases, parameters, time_digits_for_cell_id=time_digits_for_cell_id)

    #
    # read input properties to be named
    #
    inputFile = None
    prop = {}
    reference_prop = {}

    if parameters.testFile is not None:
        reference_prop = ioproperties.read_dictionary(parameters.testFile, inputpropertiesdict={})
        if 'cell_name' not in reference_prop:
            monitoring.to_log_and_console(str(proc) + ": 'cell_name' is not in '" + str(parameters.testFile) + "'")
            sys.exit(1)
        inputFile = parameters.testFile
        prop = copy.deepcopy(reference_prop)
        del prop['cell_name']

    elif parameters.inputFile is not None:
        inputFile = parameters.inputFile
        prop = ioproperties.read_dictionary(parameters.inputFile, inputpropertiesdict={})

    if inputFile is None:
        return prop

    name = inputFile.split(os.path.sep)[-1]
    if name.endswith(".xml") or name.endswith(".pkl"):
        name = name[:-4]

    if prop == {}:
        monitoring.to_log_and_console(str(proc) + ": no properties?!")
        sys.exit(1)

    #
    # build an atlas from the embryo to be named,
    # temporally align it with the reference embryo
    #
    embryo = uatlase.Atlas(prop, time_digits_for_cell_id=time_digits_for_cell_id)
    ref_atlases = atlases.get_atlases()
    ref_atlas = atlases.get_reference_atlas()
    embryo.temporally_align_with(ref_atlases[ref_atlas])
    msg = "   ... "
    msg += "linear time warping of '" + str(name) + "' wrt '" + str(ref_atlas) + "' is "
    msg += "({:.3f}, {:.3f})".format(embryo.temporal_alignment[0], embryo.temporal_alignment[1])
    monitoring.to_log_and_console(msg, 1)

    prop = _initial_naming(prop, embryo, atlases, parameters, time_digits_for_cell_id=time_digits_for_cell_id)
    #
    #
    #

    if parameters.testFile is not None:
        output = _test_naming(prop, reference_prop, time_digits_for_cell_id=time_digits_for_cell_id)
        right_naming, wrong_naming, unnamed_cells, not_in_reference = output

    if isinstance(parameters.outputFile, str):
        ioproperties.write_dictionary(parameters.outputFile, prop)

    return prop
