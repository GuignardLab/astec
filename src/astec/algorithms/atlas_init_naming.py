import os
import sys
import copy
import statistics
import numpy as np
import math
import multiprocessing
import itertools
from sklearn.neighbors import NearestNeighbors

import astec.utils.common as common
import astec.utils.atlas_division as uatlasd
import astec.utils.atlas_cell as uatlasc
import astec.utils.atlas_embryo as uatlase
import astec.utils.ioproperties as ioproperties

#
#
#
#
#

monitoring = common.Monitoring()

_instrumented_ = False

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

        doc = "\t To get an uniform sampling of the 3d directions in space, a discrete sphere is build\n"
        doc += "\t and each point of the outer surface gives a sample. The larger the radius, the more\n"
        doc += "\t vectors and the higher computational time\n"
        doc += "\t radius = 2: 26 vectors\n"
        doc += "\t radius = 2.5: 54 vectors\n"
        doc += "\t radius = 2.9: 66 vectors\n"
        doc += "\t radius = 3: 90 vectors\n"
        self.doc['direction_sphere_radius'] = doc
        self.direction_sphere_radius = 2.5

        doc = "\t Increment (in degrees) between two successive angles when enumerating\n"
        doc += "\t rotation along the z axe\n"
        self.doc['z_rotation_angle_increment'] = doc
        self.z_rotation_angle_increment = 15

        doc = "\t The residual value of a transformation is the sum of the smallest\n"
        doc += "\t distances between paired points. This parameter gives the ratio of\n"
        doc += "\t points to be retained\n"
        self.doc['pair_ratio_for_residual'] = doc
        self.pair_ratio_for_residual = 0.8

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

        self.varprint('direction_sphere_radius', self.direction_sphere_radius,
                      self.doc.get('direction_sphere_radius', None))
        self.varprint('z_rotation_angle_increment', self.z_rotation_angle_increment,
                      self.doc.get('z_rotation_angle_increment', None))
        self.varprint('pair_ratio_for_residual', self.pair_ratio_for_residual,
                      self.doc.get('pair_ratio_for_residual', None))

        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# InitNamingParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'direction_sphere_radius', self.direction_sphere_radius,
                      self.doc.get('direction_sphere_radius', None))
        self.varwrite(logfile, 'z_rotation_angle_increment', self.z_rotation_angle_increment,
                      self.doc.get('z_rotation_angle_increment', None))
        self.varwrite(logfile, 'pair_ratio_for_residual', self.pair_ratio_for_residual,
                      self.doc.get('pair_ratio_for_residual', None))

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

        self.direction_sphere_radius = self.read_parameter(parameters, 'direction_sphere_radius',
                                                           self.direction_sphere_radius)
        self.z_rotation_angle_increment = self.read_parameter(parameters, 'z_rotation_angle_increment',
                                                              self.z_rotation_angle_increment)
        self.pair_ratio_for_residual = self.read_parameter(parameters, 'pair_ratio_for_residual',
                                                           self.pair_ratio_for_residual)

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

    def __init__(self, prefix='init_naming_'):

        if "doc" not in self.__dict__:
            self.doc = {}

        EmbryoRegistrationParameters.__init__(self, prefix=prefix)
        uatlase.AtlasParameters.__init__(self, prefix=prefix)

        doc = "\t Input property file to be named. Must contain lineage, volumes and contact surfaces\n"
        doc += "\t as well as some input names (one time point should be entirely named).\n"
        self.doc['inputFile'] = doc
        self.inputFile = []
        doc = "\t Output property file."
        self.doc['outputFile'] = doc
        self.outputFile = None

        #
        doc = "\t Cell number of the developmental stade to be named\n"
        self.doc['cell_number'] = doc
        self.cell_number = 64

        doc = "\t Maximal number of iterations to be done when naming with unanimity rule\n"
        self.doc['unanimity_iterations'] = doc
        self.unanimity_iterations = None

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
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(m.shape[2]):
                if m[i][j][k] == 0:
                    continue
                di = float(i) - float(c)
                dj = float(j) - float(c)
                dk = float(k) - float(c)
                norm = math.sqrt(di * di + dj * dj + dk * dk)
                xv += [di / norm]
                yv += [dj / norm]
                zv += [dk / norm]
    v = np.zeros((3, len(xv)))
    v[0, :] = xv
    v[1, :] = yv
    v[2, :] = zv
    return v


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


def _rotation_matrix_from_quaternion(qr):
    rot = np.zeros((3, 3))
    n = qr[3] * qr[3] + qr[2] * qr[2] + qr[1] * qr[1] + qr[0] * qr[0]
    rot[0, 0] = qr[0] * qr[0] + qr[1] * qr[1] - qr[2] * qr[2] - qr[3] * qr[3]
    rot[1, 1] = qr[0] * qr[0] - qr[1] * qr[1] + qr[2] * qr[2] - qr[3] * qr[3]
    rot[2, 2] = qr[0] * qr[0] - qr[1] * qr[1] - qr[2] * qr[2] + qr[3] * qr[3]
    rot[0, 1] = (qr[1] * qr[2] - qr[0] * qr[3]) * 2.0
    rot[0, 2] = (qr[1] * qr[3] + qr[0] * qr[2]) * 2.0
    rot[1, 0] = (qr[1] * qr[2] + qr[0] * qr[3]) * 2.0
    rot[1, 2] = (qr[2] * qr[3] - qr[0] * qr[1]) * 2.0
    rot[2, 0] = (qr[1] * qr[3] - qr[0] * qr[2]) * 2.0
    rot[2, 1] = (qr[2] * qr[3] + qr[0] * qr[1]) * 2.0
    rot /= n
    return rot


########################################################################################
#
# least-square estimation
#
########################################################################################

def _ls_rigid_transformation(ref=None, flo=None):
    # compute the transformation that brings floating points onto the reference ones
    proc = "_ls_rigid_transformation"
    if ref.shape[0] != 4 or flo.shape[0] != 4:
        monitoring.to_log_and_console(str(proc) + ": input matrices should be of shape[0]==4 (array of 4D points)")
        return None
    if ref.shape[1] != flo.shape[1]:
        monitoring.to_log_and_console(str(proc) + ": input matrices should have equal shape[1]")
        return None

    # barycenters
    ref_barycenter = np.average(ref, axis=1)
    flo_barycenter = np.average(flo, axis=1)

    # centered point clouds
    cref = copy.deepcopy(ref)
    for i in range(ref.shape[1]):
        cref[:, i] = ref[:, i] - ref_barycenter
    cflo = copy.deepcopy(flo)
    for i in range(flo.shape[1]):
        cflo[:, i] = flo[:, i] - flo_barycenter

    #
    # looking for the best rotation (expressed as a quaternion)
    #
    sum_ata = np.zeros((4, 4))
    a = np.zeros((4, 4))
    for i in range(ref.shape[1]):
        a[0, 1:] = cflo[:3, i] - cref[:3, i]
        a[1, 2] = - (cflo[2, i] + cref[2, i])
        a[1, 3] = cflo[1, i] + cref[1, i]
        a[2, 3] = - (cflo[0, i] + cref[0, i])
        a[1:, 0] = - a[0, 1:]
        a[2, 1] = - a[1, 2]
        a[3, 1] = - a[1, 3]
        a[3, 2] = - a[2, 3]
        sum_ata += np.matmul(a.T, a)

    # eigen values and vectors
    # get the eigenvector corresponding to the smallest eigenvalue
    evalues, evectors = np.linalg.eig(sum_ata)
    i_min = np.argmin(evalues)

    mat = np.zeros((4, 4))
    mat[:3, :3] = _rotation_matrix_from_quaternion(evectors[:, i_min])
    mat[3, 3] = 1
    # translation part
    rog = np.matmul(mat, flo_barycenter)
    mat[:3, 3] = ref_barycenter[:3] - rog[:3]

    return mat


def _ls_similitude_transformation(ref=None, flo=None):
    # compute the transformation that brings floating points onto the reference ones
    proc = "_ls_rigid_transformation"
    if ref.shape[0] != 4 or flo.shape[0] != 4:
        monitoring.to_log_and_console(str(proc) + ": input matrices should be of shape[0]==4 (array of 4D points)")
        return None
    if ref.shape[1] != flo.shape[1]:
        monitoring.to_log_and_console(str(proc) + ": input matrices should have equal shape[1]")
        return None

    # barycenters
    ref_barycenter = np.average(ref, axis=1)
    flo_barycenter = np.average(flo, axis=1)

    # centered point clouds
    cref = copy.deepcopy(ref)
    for i in range(ref.shape[1]):
        cref[:, i] = ref[:, i] - ref_barycenter
    cflo = copy.deepcopy(flo)
    for i in range(flo.shape[1]):
        cflo[:, i] = flo[:, i] - flo_barycenter

    #
    # sum of squared modulus
    #
    sum_cref = np.sum(cref[:3, :] * cref[:3, :])
    sum_cflo = np.sum(cflo[:3, :] * cflo[:3, :])

    #
    # looking for the best rotation (expressed as a quaternion)
    #
    sum_ata = np.zeros((4, 4))
    a = np.zeros((4, 4))
    for i in range(ref.shape[1]):
        a[0, 1:] = cflo[:3, i] - cref[:3, i]
        a[1, 2] = - (cflo[2, i] + cref[2, i])
        a[1, 3] = cflo[1, i] + cref[1, i]
        a[2, 3] = - (cflo[0, i] + cref[0, i])
        a[1:, 0] = - a[0, 1:]
        a[2, 1] = - a[1, 2]
        a[3, 1] = - a[1, 3]
        a[3, 2] = - a[2, 3]
        sum_ata += np.matmul(a.T, a)

    # eigen values and vectors
    # get the eigenvector corresponding to the smallest eigenvalue
    evalues, evectors = np.linalg.eig(sum_ata)
    i_min = np.argmin(evalues)

    mat = np.zeros((4, 4))
    mat[:3, :3] = _rotation_matrix_from_quaternion(evectors[:, i_min])
    mat[3, 3] = 1

    #
    # scaling
    #
    diff = cref - np.matmul(mat, cflo)
    sum_diff = np.sum(diff[:3, :] * diff[:3, :])
    scale = (sum_cflo + sum_cref - sum_diff) / (2 * sum_cflo)

    # translation part
    rog = np.matmul(mat, flo_barycenter)
    mat[:3, :3] *= scale
    mat[:3, 3] = ref_barycenter[:3] - scale * rog[:3]

    return mat


def _ls_affine_transformation(ref=None, flo=None):
    # compute the transformation that brings floating points onto the reference ones
    proc = "_ls_rigid_transformation"
    if ref.shape[0] != 4 or flo.shape[0] != 4:
        monitoring.to_log_and_console(str(proc) + ": input matrices should be of shape[0]==4 (array of 4D points)")
        return None
    if ref.shape[1] != flo.shape[1]:
        monitoring.to_log_and_console(str(proc) + ": input matrices should have equal shape[1]")
        return None

    # barycenters
    ref_barycenter = np.average(ref, axis=1)
    flo_barycenter = np.average(flo, axis=1)

    # centered point clouds
    cref = copy.deepcopy(ref)
    for i in range(ref.shape[1]):
        cref[:, i] = ref[:, i] - ref_barycenter
    cflo = copy.deepcopy(flo)
    for i in range(flo.shape[1]):
        cflo[:, i] = flo[:, i] - flo_barycenter

    #
    # covariance matrices
    #
    cov = np.zeros((3, 3))
    cov_flo = np.zeros((3, 3))
    for i in range(ref.shape[1]):
        cov[0, :] += cref[0, i] * cflo[:3, i]
        cov[1, :] += cref[1, i] * cflo[:3, i]
        cov[2, :] += cref[2, i] * cflo[:3, i]
        cov_flo[0, :] += cflo[0, i] * cflo[:3, i]
        cov_flo[1, :] += cflo[1, i] * cflo[:3, i]
        cov_flo[2, :] += cflo[2, i] * cflo[:3, i]

    mat = np.zeros((4, 4))
    mat[:3, :3] = np.matmul(cov, np.linalg.inv(cov_flo))
    mat[3, 3] = 1

    # translation part
    rog = np.matmul(mat, flo_barycenter)
    mat[:3, 3] = ref_barycenter[:3] - rog[:3]

    return mat


def _ls_transformation(ref=None, flo=None, transformation_type="affine"):
    proc = "_ls_transformation"

    if ref.shape[0] != 4 or flo.shape[0] != 4:
        monitoring.to_log_and_console(str(proc) + ": input matrices should be of shape[0]==4 (array of 4D points)")
        return None
    if ref.shape[1] != flo.shape[1]:
        monitoring.to_log_and_console(str(proc) + ": input matrices should have equal shape[1]")
        return None

    if transformation_type == "rigid":
        return _ls_rigid_transformation(ref=ref, flo=flo)
    elif transformation_type == "similitude":
        return _ls_similitude_transformation(ref=ref, flo=flo)
    elif transformation_type == "affine":
        return _ls_affine_transformation(ref=ref, flo=flo)
    monitoring.to_log_and_console(str(proc) + ": unknown transformation type '" + str(transformation_type) + "'")
    return None


def _lts_transformation(ref=None, flo=None, transformation_type="affine", retained_fraction=0.75, verbose=False):

    proc = "_lts_transformation"
    tol_r = 1e-4
    tol_t = 1e-4
    max_iterations = 100
    iteration = 0

    if ref.shape[0] != 4 or flo.shape[0] != 4:
        monitoring.to_log_and_console(str(proc) + ": input matrices should be of shape[0]==4 (array of 4D points)")
        return None
    if ref.shape[1] != flo.shape[1]:
        monitoring.to_log_and_console(str(proc) + ": input matrices should have equal shape[1]")
        return None

    mat = _ls_transformation(ref=ref, flo=flo, transformation_type=transformation_type)

    while True:
        pmat = copy.deepcopy(mat)
        # residuals
        diff = ref - np.matmul(mat, flo)
        residuals = np.linalg.norm(diff, axis=0)
        # print("residuals = " + str(residuals))
        sort_index = np.argsort(residuals)
        # print("sort_index = " + str(sort_index))
        #
        retained_points = int(retained_fraction * ref.shape[1])

        # transformation on selected points
        nref = np.zeros((4, retained_points))
        nflo = np.zeros((4, retained_points))
        for i in range(retained_points):
            nref[:, i] = ref[:, sort_index[i]]
            nflo[:, i] = flo[:, sort_index[i]]
        mat = _ls_transformation(ref=nref, flo=nflo, transformation_type=transformation_type)
        del nref
        del nflo

        # end condition
        dmat = (pmat - mat) * (pmat - mat)
        eps_r = math.sqrt(np.sum(dmat[:3, :3]))
        eps_t = math.sqrt(np.sum(dmat[:3, 3]))
        if verbose:
            msg = "   lts #" + str(iteration) + ": " + "eps_r = " + str(eps_r) + " - " + "eps_t = " + str(eps_t)
            monitoring.to_log_and_console(msg, 2)
        if (eps_r < tol_r and eps_t < tol_t) or iteration >= max_iterations:
            return mat
        iteration += 1

    # should not occur
    return None


########################################################################################
#
# icp
#
########################################################################################

def icp(ref=None, flo=None, transformation_type="affine", estimation="lts", verbose=False):
    proc = "icp"
    tol_r = 1e-4
    tol_t = 1e-4
    max_iterations = 100
    iteration = 0

    if ref.shape[0] != 3 or flo.shape[0] != 3:
        print(proc + ": input matrices should be of shape[0]==3 (array of 3D points)")
        return None

    mat = np.identity(4)

    # transform to 4D points
    lref = np.ones((4, ref.shape[1]))
    lref[:3, :] = ref
    lflo = np.ones((4, flo.shape[1]))
    lflo[:3, :] = flo

    #
    nbrs = NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(lref.T)

    while True:
        pmat = copy.deepcopy(mat)
        # update floating points
        tflo = np.matmul(mat, lflo)
        # get closest ref point for each floating point
        # closest = _pairing(lref, tflo)
        distances, indices = nbrs.kneighbors(tflo.T)
        closest = indices.flatten()

        # in case for the future
        # anticipate that None can be a pairing
        # get indices of flo that are paired
        i_flo = [i for i, c in enumerate(closest) if c is not None]
        npairings = len(i_flo)

        # build pairs
        sref = np.ones((4, npairings))
        sflo = np.ones((4, npairings))
        for i, j in enumerate(i_flo):
            sflo[:, i] = tflo[:, j]
            sref[:, i] = lref[:, closest[j]]

        # compute incremental transformation and compose it with previous one
        if estimation == "ls":
            imat = _ls_transformation(sref, sflo, transformation_type=transformation_type)
        elif estimation == "lts":
            imat = _lts_transformation(sref, sflo, transformation_type=transformation_type, verbose=False)
        else:
            print(proc + ": unknown transformation type '" + str(transformation_type) + "'")
            return None
        mat = np.matmul(imat, pmat)

        # end condition
        dmat = (pmat - mat) * (pmat - mat)
        eps_r = math.sqrt(np.sum(dmat[:3, :3]))
        eps_t = math.sqrt(np.sum(dmat[:3, 3]))
        if verbose:
            print(" icp #" + str(iteration) + ": " + "eps_r = " + str(eps_r) + " - " + "eps_t = " + str(eps_t))
        if (eps_r < tol_r and eps_t < tol_t) or iteration >= max_iterations:
            return mat
        iteration += 1
    return None


########################################################################################
#
# procedures for embryos
#
########################################################################################

def _get_times(prop, n_cells=64, time_digits_for_cell_id=4):
    cells = list(prop.keys())
    cells = sorted(cells)
    cells_per_time = {}
    div = 10 ** time_digits_for_cell_id
    for c in cells:
        t = int(c) // div
        cells_per_time[t] = cells_per_time.get(t, 0) + 1
    times = [t for t in cells_per_time if cells_per_time[t] == n_cells]
    return times


def _embryo_barycenter(embryo, t, time_digits_for_cell_id=4):
    b = np.zeros(3)
    s = 0.0
    volumes = embryo.cell_volume
    barycenters = embryo.cell_barycenter
    div = 10 ** time_digits_for_cell_id
    for c in volumes:
        if int(c) // div != t:
            continue
        if int(c) % div == 1 or int(c) % div == 0:
            continue
        if c not in barycenters:
            continue
        b += volumes[c] * barycenters[c]
        s += volumes[c]
    b /= s
    return b


def _embryo_volume(embryo, t, time_digits_for_cell_id=4):
    s = 0.0
    volumes = embryo.cell_volume
    div = 10 ** time_digits_for_cell_id
    for c in volumes:
        if int(c) // div != t:
            continue
        if int(c) % div == 1 or int(c) % (10 ** 4) == 0:
            continue
        s += volumes[c]
    return s


def _get_barycenters(embryo, t, transformation=None, time_digits_for_cell_id=4):
    barycenters = embryo.cell_barycenter
    div = 10 ** time_digits_for_cell_id
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
    n, i, ni, scale_mat, flo_embryo_bar, ref_embryo_bar, flo_bars, ref_bars, registration_par, verbose = parameters

    if verbose:
        msg = "processing direction[" + str(i) + "/" + str(ni) + "] = " + str(n)
        monitoring.to_log_and_console("      ... " + msg)
    #
    rotz = np.identity(4)
    rotxy = np.identity(4)
    sp = np.dot(n, np.array([0, 0, 1]))

    if sp > 0.999:
        # n is aligned with z
        # do nothing
        pass
    elif sp < -0.999:
        # n is aligned with -z
        # rotation of pi wrt x
        rotz[:3, :3] = _rotation_matrix_from_rotation_vector(np.array([np.pi, 0, 0]))
    else:
        # rotation of math.acos(n.z) wrt vector n x z
        # that puts n along z
        vrot = np.cross(n, np.array([0, 0, 1]))
        vnor = np.linalg.norm(vrot * vrot)
        vrot *= math.acos(sp) / vnor
        rotz[:3, :3] = _rotation_matrix_from_rotation_vector(vrot)

    # for 64 cells, the smallest angle (from the barycenter) between 2 adjacent cells
    # is 7 degrees
    rigid_matrice = []
    affine_matrice = []
    average_residual = []
    angles = np.linspace(0, 2 * np.pi, int(360 / registration_par.z_rotation_angle_increment), endpoint=False)
    for a in angles:
        # compute initial transformation
        rotxy[:3, :3] = _rotation_matrix_from_rotation_vector(a * np.array([0, 0, 1]))
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
        rigid_mat = icp(ref=ref_bars, flo=flo, transformation_type="rigid")
        res_rigid_mat = np.matmul(rigid_mat, init_mat)
        # compute an affine transformation
        v = np.ones(4)
        for j in range(flo_bars.shape[1]):
            v[:3] = flo_bars[:, j]
            flo[:, j] = (np.matmul(res_rigid_mat, v))[:3]
        affine_mat = icp(ref=ref_bars, flo=flo, transformation_type="affine")
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


def _coregister_embryos_one_timepoint(embryo, atlas, parameters, floating_time, reference_time,
                                      time_digits_for_cell_id=4, verbose=True):
    proc = "_coregister_embryos_one_timepoint"

    # reference embryo
    # ref_barycenter is the barycenter in homogeneous coordinates
    # NB: it is in voxel coordinates
    # ref_barycenters is an np.array of shape (3, n) of the n cell barycenters at time reference_time
    # ref_barycenters[0][i] is the the x coordinate of the ith barycenter

    ref_embryo_barycenter = np.ones(4)
    ref_embryo_barycenter[:3] = _embryo_barycenter(atlas, reference_time,
                                                   time_digits_for_cell_id=time_digits_for_cell_id)
    ref_barycenters = _get_barycenters(atlas, reference_time, time_digits_for_cell_id=time_digits_for_cell_id)

    flo_embryo_barycenter = np.ones(4)
    flo_embryo_barycenter[:3] = _embryo_barycenter(embryo, floating_time,
                                                   time_digits_for_cell_id=time_digits_for_cell_id)
    flo_barycenters = _get_barycenters(embryo, floating_time, time_digits_for_cell_id=time_digits_for_cell_id)

    # scale to correct for volume differences
    # a voxel size correction of scale = np.cbrt(volref / volflo), applied to the floating points
    # makes the floating embryo volume equal to the reference one
    #
    # since we also have the global correction embryo.get_voxelsize_correction() with
    # volref * (ref.get_voxelsize_correction())^3 = volflo * (flo.get_voxelsize_correction())^3
    # it somes that scale = flo.get_voxelsize_correction() / ref.get_voxelsize_correction()

    scale = embryo.get_voxelsize_correction(floating_time) / atlas.get_voxelsize_correction(reference_time)
    scale_mat = scale * np.identity(4)
    scale_mat[3, 3] = 1.0
    if False:
        ref_volume = _embryo_volume(atlas, reference_time, time_digits_for_cell_id=time_digits_for_cell_id)
        flo_volume = _embryo_volume(embryo, floating_time, time_digits_for_cell_id=time_digits_for_cell_id)
        monitoring.to_log_and_console("ref volume = " + str(ref_volume) + " - flo volume = " + str(flo_volume))
        monitoring.to_log_and_console("\t corrected flo volume = " + str(flo_volume * scale * scale * scale))

    #
    # get a set of vectors
    #
    vectors = _build_vector_sphere(parameters.direction_sphere_radius)


    #
    # we could improve the computational time by parallelising the calculation
    # with multiprocessing
    #
    pool = multiprocessing.Pool(processes=7)
    mapping = []

    for i in range(vectors.shape[1]):
        n = vectors[:, i]
        mapping.append((vectors[:, i], i+1, vectors.shape[1], scale_mat, flo_embryo_barycenter, ref_embryo_barycenter,
                        flo_barycenters, ref_barycenters, parameters, verbose))
    outputs = pool.map(_coregister_embryos_one_direction, mapping)
    pool.close()
    pool.terminate()

    rigid_matrice = []
    affine_matrice = []
    average_residual = []
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


def _coregister_embryos(transformations, embryo, atlas, atlas_name, parameters, time_digits_for_cell_id=4,
                        verbose=True):
    proc = "_coregister_embryos"

    flotimes = _get_times(embryo.cell_contact_surface, parameters.cell_number,
                          time_digits_for_cell_id=time_digits_for_cell_id)
    if len(flotimes) == 0:
        msg = "there is no time points with " + str(parameters.cell_number) + " in the embryo to be named\n"
        msg += "\t skip registration with atlas '" + str(atlas_name) + "'"
        monitoring.to_log_and_console("   ... " + msg)
        return transformations
    first_flotime = statistics.median_low(flotimes)
    reftimes = _get_times(atlas.cell_contact_surface, parameters.cell_number,
                          time_digits_for_cell_id=time_digits_for_cell_id)
    if len(reftimes) == 0:
        msg = "there is no time points with " + str(parameters.cell_number) + " in the atlas\n"
        msg += "\t skip registration with atlas '" + str(atlas_name) + "'"
        monitoring.to_log_and_console("   ... " + msg)
        return transformations
    first_reftime = statistics.median_low(reftimes)

    if verbose:
        msg = "co-registration of time #" + str(first_flotime) + " of embryo to be named"
        msg += " with time #" + str(first_reftime) + " of atlas '" + str(atlas_name) + "'"
        monitoring.to_log_and_console("   ... " + msg)
    #
    # the rigid transformation can be used to initialize rigid transformations when
    # registering other times of embryo to be named to other times of atlas
    #
    rigid_mat, affine_mat = _coregister_embryos_one_timepoint(embryo, atlas, parameters, first_flotime, first_reftime,
                                                              time_digits_for_cell_id=time_digits_for_cell_id,
                                                              verbose=verbose)
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
        transformations[first_flotime][atlas_name][first_reftime] = affine_mat

    return transformations


########################################################################################
#
#
#
########################################################################################

def _get_success(atlases, parameters, time_digits_for_cell_id=4, verbose=True):
    monitoring.to_log_and_console(" - compute transformations ")
    #
    # compute transformations
    #
    all_atlases = atlases.get_atlases()
    transformations = {}
    for a in all_atlases:
        #
        # compute transformations with a as the floating embryo
        #
        transformations[a] = {}
        for b in all_atlases:
            if b == a:
                continue
            monitoring.to_log_and_console("   - compute transformation of '" + str(a) + "' versus '" + str(b) + "'")
            transformations[a] = _coregister_embryos(transformations[a], all_atlases[a], all_atlases[b], b, parameters,
                                                     time_digits_for_cell_id=time_digits_for_cell_id, verbose=verbose)

    #
    # set names
    #
    right_naming = {}
    wrong_naming = {}
    for a in all_atlases:
        # atlas to be tested
        ref_test_atlas = copy.deepcopy(all_atlases[a])
        test_atlas = copy.deepcopy(ref_test_atlas)
        test_atlas.del_cell_name()
        # reference atlas list
        ref_atlas_list = list(all_atlases.keys())
        ref_atlas_list.remove(a)

        #
        # pick combinations of reference atlas
        #
        monitoring.to_log_and_console(" - compare " + str(a) + " versus " + str(ref_atlas_list))
        for i in range(1, len(ref_atlas_list) + 1):
            subset_atlas = list(itertools.combinations(ref_atlas_list, i))
            monitoring.to_log_and_console("      " + str(len(subset_atlas)) + " combinations of " + str(i) + " atlases")
            for s in subset_atlas:
                monitoring.to_log_and_console("     - compare " + str(a) + " versus " + str(s))
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
                name_both, name_flo, name_ref = _get_floating_names(transformations[a], test_atlas, tmp_atlases,
                                                                    parameters,
                                                                    time_digits_for_cell_id=time_digits_for_cell_id,
                                                                    verbose=verbose)

                #
                # get consensual names
                #
                prop = _iterate_unanimity_names(name_both, parameters, verbose=verbose)

                #
                #
                #
                right, wrong, not_in = _test_naming(prop, ref_test_atlas.get_property(),
                                                    time_digits_for_cell_id=time_digits_for_cell_id, verbose=False)
                right_naming[i] = right_naming.get(i, []) + [right]
                wrong_naming[i] = wrong_naming.get(i, []) + [wrong]

    # print("right_naming = " + str(right_naming))
    # print("wrong_naming = " + str(wrong_naming))
    return right_naming, wrong_naming


def _figure_success_wrt_atlasnumber(atlases, parameters, time_digits_for_cell_id=4):
    right_naming, wrong_naming = _get_success(atlases, parameters, time_digits_for_cell_id=time_digits_for_cell_id,
                                              verbose=False)
    proc = "_figure_success_wrt_atlasnumber"

    filename = 'figure_success_wrt_atlasnumber'
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
    f.write("labels = sorted(list(right_naming.keys()))\n")

    f.write("\n")
    f.write("right_list = []\n")
    f.write("wrong_list = []\n")
    f.write("for i, l in enumerate(labels):\n")
    f.write("    right_list += [right_naming[l]]\n")
    f.write("    wrong_list += [wrong_naming[l]]\n")

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
    f.write("    plt.savefig('success_wrt_atlasnumber")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")

    f.close()


def _generate_figure(atlases, parameters, time_digits_for_cell_id=4):

    generate_figure = (isinstance(parameters.generate_figure, bool) and parameters.generate_figure) or \
                      (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'all') or \
                      (isinstance(parameters.generate_figure, list) and 'all' in parameters.generate_figure)

    if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'success-wrt-atlasnumber') \
            or (isinstance(parameters.generate_figure, list) and 'success-wrt-atlasnumber' in parameters.generate_figure) \
            or generate_figure:
        monitoring.to_log_and_console("... generate figure success versus atlas number", 1)
        _figure_success_wrt_atlasnumber(atlases, parameters, time_digits_for_cell_id=time_digits_for_cell_id)
        monitoring.to_log_and_console("... done", 1)


########################################################################################
#
#
#
########################################################################################

def _test_naming(prop, reference_prop, time_digits_for_cell_id=4, verbose=True):
    proc = "_test_naming"
    #
    #
    #
    cells = list(prop['cell_name'].keys())
    div = 10 ** time_digits_for_cell_id
    times = set([int(c) // div for c in cells])
    if len(times) == 0:
        if verbose:
            msg = "no names were set"
            monitoring.to_log_and_console("... " + msg)
    elif len(times) > 1:
        if verbose:
            msg = "weird, names were set on several time points " + str(times)
            monitoring.to_log_and_console("... " + msg)

    refcells = [int(c) // div for c in reference_prop['cell_name'] if int(c) // div in times]
    right_naming = 0
    wrong_naming = 0
    not_in_reference = 0
    for c in cells:
        if c not in reference_prop['cell_name']:
            not_in_reference += 1
            continue
        if prop['cell_name'][c] != reference_prop['cell_name'][c]:
            wrong_naming += 1
            continue
        right_naming += 1

    if verbose:
        msg = "ground-truth cells = " + str(len(refcells)) + " for times " + str(times) + "\n"
        msg += "\t tested cells = " + str(len(cells)) + "\n"
        msg += "\t retrieved names = {:d}/{:d}\n".format(right_naming, len(cells))
        if wrong_naming > 0:
            msg += "\t ... error in naming = {:d}/{:d}\n".format(wrong_naming, len(cells))
        if not_in_reference > 0:
            msg += "\t ... names not in reference = {:d}/{:d}\n".format(not_in_reference, len(cells))
        monitoring.to_log_and_console("summary" + ": " + msg)

    return right_naming, wrong_naming, not_in_reference


########################################################################################
#
#
#
########################################################################################

def _get_correspondence(embryo, parameters, time_digits_for_cell_id=4, verbose=True):
    #
    # for the cells of each time where the number of cells is equal to parameters.cell_number
    # gives the cell of the first time where the number of cells is equal to parameters.cell_number
    # it will allows to name only this first time while names can be issued from any time
    #
    correspondence = {}
    flotimes = _get_times(embryo.cell_contact_surface, parameters.cell_number,
                          time_digits_for_cell_id=time_digits_for_cell_id)
    if len(flotimes) == 0:
        msg = "   ... there is no time points with 64 in the embryo to be named"
        monitoring.to_log_and_console(msg)
        return correspondence
    flotimesmin = min(flotimes)
    lineage = embryo.lineage
    reverse_lineage = {v: k for k, values in lineage.items() for v in values}
    cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
    div = 10 ** time_digits_for_cell_id
    for c in cells:
        if (int(c) // div) not in flotimes or int(c) % div == 1:
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


def _set_candidate_name(name_from_both, name_from_flo, name_from_ref, correspondence, embryo, reference, tflo, tref,
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
        if i_ref_to_flo[j] == i:
            name_from_both[correspondence[cf]] = name_from_both.get(correspondence[cf], []) + [ref_name[cr]]
            reciprocal_n += 1
        else:
            name_from_flo[correspondence[cf]] = name_from_flo.get(correspondence[cf], []) + [ref_name[cr]]

    for j, cr in enumerate(cell_tref):
        i = i_ref_to_flo[j]
        cf = cell_tflo[i]
        if i_flo_to_ref[i] == j:
            pass
        else:
            name_from_ref[correspondence[cf]] = name_from_ref.get(correspondence[cf], []) + [ref_name[cr]]

    return name_from_both, name_from_flo, name_from_ref, reciprocal_n


def _get_floating_names(transformations, embryo, atlases, parameters, time_digits_for_cell_id=4, verbose=True):
    #
    # correspondences
    #
    correspondence = _get_correspondence(embryo, parameters, time_digits_for_cell_id=time_digits_for_cell_id,
                                         verbose=verbose)

    name_both = {}
    name_flo = {}
    name_ref = {}
    if len(correspondence) == 0:
        return name_both, name_flo, name_ref

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
                name_both, name_flo, name_ref, n = _set_candidate_name(name_both, name_flo, name_ref, correspondence,
                                                                       embryo, ref_atlases[a], tflo, tref,
                                                                       transformations[tflo][a][tref],
                                                                       time_digits_for_cell_id=time_digits_for_cell_id)
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
            n +=1
    return prop, name, n


def _iterate_unanimity_names(name_both, parameters, verbose=False):
    iteration = 0
    prop = {'cell_name':{}}
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

def _initial_naming(embryo, atlases, parameters, time_digits_for_cell_id=4, verbose=True):
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
    msg = "setting names"
    if verbose:
        monitoring.to_log_and_console("... co-registrations")
    transformations = {}
    ref_atlases = atlases.get_atlases()
    for a in ref_atlases:
        transformations = _coregister_embryos(transformations, embryo, ref_atlases[a], a, parameters,
                                              time_digits_for_cell_id=time_digits_for_cell_id, verbose=verbose)
    if verbose:
        monitoring.to_log_and_console("... setting names")

    #
    # collect names from each co-registration
    # name_both, name_flo, name_ref are dictionaries indexed by cell ids
    # name_both
    name_both, name_flo, name_ref = _get_floating_names(transformations, embryo, atlases, parameters,
                                                        time_digits_for_cell_id=time_digits_for_cell_id,
                                                        verbose=verbose)

    #
    # get consensual names
    #
    prop = _iterate_unanimity_names(name_both, parameters, verbose=verbose)

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

    #
    # should we clean reference here?
    #
    atlases = uatlasd.DivisionAtlases(parameters=parameters)
    atlases.add_atlases(parameters.atlasFiles, parameters, time_digits_for_cell_id=time_digits_for_cell_id)

    if parameters.generate_figure is not None and parameters.generate_figure is not False:
        _generate_figure(atlases, parameters, time_digits_for_cell_id=time_digits_for_cell_id)
        sys.exit()
    #
    # read input properties to be named
    #
    prop = {}
    reference_prop = {}

    if parameters.testFile is not None:
        reference_prop = ioproperties.read_dictionary(parameters.testFile, inputpropertiesdict={})
        if 'cell_name' not in reference_prop:
            monitoring.to_log_and_console(str(proc) + ": 'cell_name' is not in '" + str(parameters.testFile) + "'")
            sys.exit(1)
        else:
            prop = copy.deepcopy(reference_prop)
            del prop['cell_name']
    elif parameters.inputFile is not None:
        prop = ioproperties.read_dictionary(parameters.inputFile, inputpropertiesdict={})

    if prop == {}:
        monitoring.to_log_and_console(str(proc) + ": no properties?!")
        sys.exit(1)

    #
    # build an atlas from the embryo to be named,
    # temporally align it with the reference embryo
    #
    embryo = uatlase.Atlas(prop)
    ref_atlases = atlases.get_atlases()
    ref_atlas = atlases.get_reference_atlas()
    embryo.temporally_align_with(ref_atlases[ref_atlas], time_digits_for_cell_id=time_digits_for_cell_id)

    prop = _initial_naming(embryo, atlases, parameters, time_digits_for_cell_id=time_digits_for_cell_id)
    #
    #
    #

    if parameters.testFile is not None:
        right_naming, wrong_naming, not_in_reference = _test_naming(prop, reference_prop,
                                                                    time_digits_for_cell_id=time_digits_for_cell_id)

    if isinstance(parameters.outputFile, str):
        ioproperties.write_dictionary(parameters.outputFile, prop)

    return prop
