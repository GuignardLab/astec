
import copy
import math
import os
import sys
from collections import Counter

import numpy as np
import sklearn.linear_model as sklm
import sklearn.neighbors as skn
import matplotlib.pyplot as plt

import astec.utils.common as common
import astec.utils.icp as icp
import astec.utils.properties as properties
import astec.utils.diagnosis as udiagnosis
import astec.utils.ioproperties as ioproperties

monitoring = common.Monitoring()


###########################################################
#
#
#
############################################################

def _build_vector_sphere(r=3):
    """

    Parameters
    ----------
    r : float
        sphere radius
        radius =  2.0:    26 vectors, angle between neighboring vectors in [36.26, 60.0] degrees
        radius =  2.3:    38 vectors, angle between neighboring vectors in [26.57, 54.74] degrees
        radius =  2.5:    54 vectors, angle between neighboring vectors in [24.09, 43.09] degrees
        radius =  2.9:    66 vectors, angle between neighboring vectors in [18.43, 43.09] degrees
        radius =  3.0:    90 vectors, angle between neighboring vectors in [17.72, 43.09] degrees
        radius =  3.5:    98 vectors, angle between neighboring vectors in [15.79, 32.51] degrees
        radius =  3.7:   110 vectors, angle between neighboring vectors in [15.26, 32.51] degrees
        radius =  3.8:   134 vectors, angle between neighboring vectors in [14.76, 29.50] degrees
        radius =  4.0:   222 vectors, angle between neighboring vectors in [10.31, 22.57] degrees
        radius =  5.0:   222 vectors, angle between neighboring vectors in [10.31, 22.57] degrees
        radius = 10.0:   978 vectors, angle between neighboring vectors in [4.40, 10.58] degrees
        radius = 15.0:  2262 vectors, angle between neighboring vectors in [2.98, 6.93] degrees
        radius = 20.0:  4026 vectors, angle between neighboring vectors in [2.25, 5.16] degrees
        radius = 25.0:  6366 vectors, angle between neighboring vectors in [1.73, 4.01] degrees
        radius = 30.0:  9194 vectors, angle between neighboring vectors in [1.46, 3.40] degrees
        radius = 35.0: 12542 vectors, angle between neighboring vectors in [1.26, 2.90] degrees
        radius = 40.0: 16418 vectors, angle between neighboring vectors in [1.08, 2.53] degrees
    Returns
    -------
    a 3D numpy ndarray where outer sphere points have non-zero value, and the center coordinates
    """
    #
    # integer radius: build a 3D matrix of size 2*(ir+1)+1
    #
    ir = math.ceil(r)
    s = int(2 * ir + 1 + 2)
    m = np.zeros((s, s, s), dtype=np.int8)
    #
    # center coordinate is r+1
    #
    c = ir + 1
    #
    # fill the sphere, mark points with distance <= r
    #
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(m.shape[2]):
                di = float(i) - float(c)
                dj = float(j) - float(c)
                dk = float(k) - float(c)
                if math.sqrt(di * di + dj * dj + dk * dk) <= r:
                    m[i][j][k] = 1
    #
    # set the border voxels to 2
    #
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
    #
    # erase the inner sphere
    #
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(m.shape[2]):
                if m[i][j][k] == 1:
                    m[i][j][k] = 0
    #
    # return the matrix and the center coordinates
    #
    return m, (c, c, c)


###########################################################
#
#
#
############################################################

class EmbryoSymmetryParameters(common.PrefixedParameter):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix=""):

        common.PrefixedParameter.__init__(self, prefix=prefix)

        if "doc" not in self.__dict__:
            self.doc = {}

        #
        #
        #

        doc = "\t Sphere radius to build the distribution support, for looking for symmetry direction\n"
        doc += "\t radius = 10:   978 vectors, angle between neighboring vectors in [4.40, 10.58] degrees\n"
        doc += "\t radius = 15:  2262 vectors, angle between neighboring vectors in [2.98, 6.93] degrees\n"
        doc += "\t radius = 20:  4026 vectors, angle between neighboring vectors in [2.25, 5.16] degrees\n"
        doc += "\t radius = 25:  6366 vectors, angle between neighboring vectors in [1.73, 4.01] degrees\n"
        doc += "\t radius = 30:  9194 vectors, angle between neighboring vectors in [1.46, 3.40] degrees\n"
        doc += "\t radius = 35: 12542 vectors, angle between neighboring vectors in [1.26, 2.90] degrees\n"
        doc += "\t radius = 40: 16418 vectors, angle between neighboring vectors in [1.08, 2.53] degrees\n"
        self.doc['sphere_radius'] = doc
        self.sphere_radius = 20

        doc = "\t Sigma (standard deviation) to build the direction distribution (in radian)\n"
        self.doc['sigma'] = doc
        self.sigma = 0.5

        doc = "\t Threshold on the distribution value. Only maxima above this threshold are\n"
        doc += "\t keep. Recall that the distribution values are normalize so that the maximum\n"
        doc += "\t is 1.\n"
        self.doc['maxima_threshold'] = doc
        self.maxima_threshold = 0.5

        doc = "\t Evaluation method to sort (and then select) distribution maxima\n"
        doc += "\t 'value':\n"
        doc += "\t 'single_pairing'\n"
        doc += "\t 'multiple_pairing'\n"
        doc += "\t 'symmetric_registration'\n"
        self.doc['maxima_evaluation'] = doc
        self.maxima_evaluation = 'value'

        doc = "\t Number of distribution maxima to be retained as candidates\n"
        doc += "\t 'None' or negative number means all of them\n"
        self.doc['maxima_number'] = doc
        self.maxima_number = 12

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('#')
        print('# EmbryoSymmetryParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        self.varprint('sphere_radius', self.sphere_radius, self.doc.get('sphere_radius', None))
        self.varprint('sigma', self.sigma, self.doc.get('sigma', None))
        self.varprint('maxima_threshold', self.maxima_threshold, self.doc.get('maxima_threshold', None))
        self.varprint('maxima_evaluation', self.maxima_evaluation, self.doc.get('maxima_evaluation', None))
        self.varprint('maxima_number', self.maxima_number, self.doc.get('maxima_number', None))
        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# EmbryoSymmetryParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'sphere_radius', self.sphere_radius, self.doc.get('sphere_radius', None))
        self.varwrite(logfile, 'sigma', self.sigma, self.doc.get('sigma', None))
        self.varwrite(logfile, 'maxima_threshold', self.maxima_threshold, self.doc.get('maxima_threshold', None))
        self.varwrite(logfile, 'maxima_evaluation', self.maxima_evaluation, self.doc.get('maxima_evaluation', None))
        self.varwrite(logfile, 'maxima_number', self.maxima_number, self.doc.get('maxima_number', None))

        logfile.write("\n")
        return

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

        self.sphere_radius = self.read_parameter(parameters, 'sphere_radius', self.sphere_radius)
        self.sigma = self.read_parameter(parameters, 'sigma', self.sigma)
        self.maxima_threshold = self.read_parameter(parameters, 'maxima_threshold', self.maxima_threshold)
        self.maxima_evaluation = self.read_parameter(parameters, 'maxima_evaluation', self.maxima_evaluation)
        self.sigma = self.read_parameter(parameters, 'maxima_number', self.maxima_number)

    def update_from_parameter_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = common.load_source(parameter_file)
        self.update_from_parameters(parameters)


###########################################################
#
#
#
############################################################

def _distribution_build_support(parameters):
    proc = "_distribution_build_support"

    if not isinstance(parameters, EmbryoSymmetryParameters):
        msg = "parameters expected type is 'EmbryoSymmetryParameters' but was '"
        msg += type(parameters) + "' instead"
        monitoring.to_log_and_console(proc + ": " + msg)
        return None

    m, c = _build_vector_sphere(parameters.sphere_radius)

    #
    # points of the outer sphere have a value of 2, all other points are at 0
    #

    direction_distribution = []

    angles = []
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(m.shape[2]):
                if m[i][j][k] == 0:
                    continue
                p = {'voxel': (i, j, k)}
                v = np.array([(i - c[0]), (j - c[1]), (k - c[2])])
                p['vector'] = v / np.linalg.norm(v)
                p['value'] = 0.0
                direction_distribution += [p]
                #
                # error evaluation: angles between neighboring vectors
                #
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        for dk in range(-1, 2):
                            if di == 0 and dj == 0 and dk == 0:
                                continue
                            if m[i+di][j+dj][k+dk] == 0:
                                continue
                            v = np.array([(i + di - c[0]), (j + dj - c[1]), (k + dk - c[2])])
                            angles += [math.acos(np.dot(p['vector'], v / np.linalg.norm(v)))]

    monitoring.to_log_and_console("      ... direction distribution build with r = " + str(parameters.sphere_radius), 2)
    monitoring.to_log_and_console("          vectors = " + str(len(direction_distribution)), 2)
    min_angles = min(angles)*180.0/np.pi
    max_angles = max(angles) * 180.0 / np.pi
    msg = "          angles between adjacent vectors in [{:.2f}, {:.2f}]".format(min_angles, max_angles)
    monitoring.to_log_and_console(msg, 2)

    return direction_distribution


def _distribution_direction_build(embryo, timepoint, parameters):
    proc = "_distribution_direction_build"
    if not isinstance(parameters, EmbryoSymmetryParameters):
        msg = "parameters expected type is 'EmbryoSymmetryParameters' but was '"
        msg += type(parameters) + "' instead"
        monitoring.to_log_and_console(proc + ": " + msg)
        return

    #
    # embryo properties
    #
    cell_contact_surface = embryo.cell_contact_surface
    barycenters = embryo.cell_barycenter
    time_digits_for_cell_id = embryo.time_digits_for_cell_id
    distribution = embryo.direction_distribution[timepoint]

    #
    # cells to be processed
    #
    div_cell = 10 ** time_digits_for_cell_id
    cells = [c for c in cell_contact_surface if (int(c) // div_cell == timepoint) and int(c) % div_cell != 1]

    #
    # initialisation
    #
    div = 2 * parameters.sigma * parameters.sigma
    # do not compute contribution for angles > 3 * sigma
    angle_threshold = 3.0 * parameters.sigma

    contribution_list = np.zeros(len(distribution))
    v = np.zeros(3)

    for c in cells:
        for d in cell_contact_surface[c]:
            if d <= c:
                continue

            #
            # approximation of the normal vector of the contact surface between c and d
            #
            v[0] = barycenters[c][0] - barycenters[d][0]
            v[1] = barycenters[c][1] - barycenters[d][1]
            v[2] = barycenters[c][2] - barycenters[d][2]
            v /= np.linalg.norm(v)

            #
            # get contributions for vector
            #

            plist = []
            sum_positive_contributions = 0.0
            nlist = []
            sum_negative_contributions = 0.0

            for i, p in enumerate(distribution):
                sp = p['vector'][0] * v[0] + p['vector'][1] * v[1] + p['vector'][2] * v[2]
                sp = max(min(sp, 1.0), -1.0)
                #
                # positive contribution
                # math.acos(x) return the arc cosine of x, in radians. The result is between 0 and pi.
                #
                angle = math.acos(sp)
                if angle < angle_threshold:
                    contrib = math.exp(- angle / div)
                    sum_positive_contributions += contrib
                    plist += [(i, contrib)]
                #
                # negative contribution
                # cos (pi-a) = - cos a
                # so acos( -sp ) will be pi -  acos( sp )
                angle = np.pi - angle
                if angle < angle_threshold:
                    contrib = math.exp(- angle / div)
                    sum_negative_contributions += contrib
                    nlist += [(i, contrib)]

            for (i, contrib) in plist:
                contribution_list[i] += cell_contact_surface[c][d] * contrib / sum_positive_contributions
            for (i, contrib) in nlist:
                contribution_list[i] += cell_contact_surface[c][d] * contrib / sum_negative_contributions

    #
    # set distribution values on distribution support
    #
    max_contribution = max(contribution_list)
    for i, v in enumerate(contribution_list):
        distribution[i]['value'] = v / max_contribution

    embryo.direction_distribution[timepoint] = distribution
    return


def _distribution_direction_maxima(embryo, timepoint, parameters):

    #
    # build a voxel array and set the distribution values
    # see _build_vector_sphere()
    #
    ir = math.ceil(parameters.sphere_radius)
    s = int(2 * ir + 1 + 2)
    values = np.zeros((s, s, s))
    distribution = embryo.direction_distribution[timepoint]

    for p in distribution:
        values[p['voxel'][0]][p['voxel'][1]][p['voxel'][2]] = p['value']

    maxima = []
    for p in distribution:
        #
        # hard threshold on maxima
        # recall the maximum of the distribution has to set to 1.0
        #
        if p['value'] <= parameters.maxima_threshold:
            continue
        i = p['voxel'][0]
        j = p['voxel'][1]
        k = p['voxel'][2]
        if values[i][j][k] == 0.0:
            continue
        if values[i][j][k] < values[i-1][j-1][k-1]:
            continue
        if values[i][j][k] < values[i][j-1][k-1]:
            continue
        if values[i][j][k] < values[i+1][j-1][k-1]:
            continue
        if values[i][j][k] < values[i-1][j][k-1]:
            continue
        if values[i][j][k] < values[i][j][k-1]:
            continue
        if values[i][j][k] < values[i+1][j][k-1]:
            continue
        if values[i][j][k] < values[i-1][j+1][k-1]:
            continue
        if values[i][j][k] < values[i][j+1][k-1]:
            continue
        if values[i][j][k] < values[i+1][j+1][k-1]:
            continue
        if values[i][j][k] < values[i-1][j-1][k]:
            continue
        if values[i][j][k] < values[i][j-1][k]:
            continue
        if values[i][j][k] < values[i+1][j-1][k]:
            continue
        if values[i][j][k] < values[i-1][j][k]:
            continue
        if values[i][j][k] < values[i+1][j][k]:
            continue
        if values[i][j][k] < values[i-1][j+1][k]:
            continue
        if values[i][j][k] < values[i][j+1][k]:
            continue
        if values[i][j][k] < values[i+1][j+1][k]:
            continue
        if values[i][j][k] < values[i-1][j-1][k+1]:
            continue
        if values[i][j][k] < values[i][j-1][k+1]:
            continue
        if values[i][j][k] < values[i+1][j-1][k+1]:
            continue
        if values[i][j][k] < values[i-1][j][k+1]:
            continue
        if values[i][j][k] < values[i][j][k+1]:
            continue
        if values[i][j][k] < values[i+1][j][k+1]:
            continue
        if values[i][j][k] < values[i-1][j+1][k+1]:
            continue
        if values[i][j][k] < values[i][j+1][k+1]:
            continue
        if values[i][j][k] < values[i+1][j+1][k+1]:
            continue
        maxima += [p]

    msg = "      ... found {:d} distribution maxima".format(len(maxima))
    monitoring.to_log_and_console(msg, 2)

    embryo.direction_distribution_maxima[timepoint] = maxima
    return


def _distribution_direction_evaluate_symicp(embryo, timepoint):

    icp.monitoring.copy(monitoring)

    embryo_barycenter = embryo.get_embryo_barycenter(timepoint)
    candidates = embryo.direction_distribution_maxima[timepoint]
    #
    # barycenters is a np.array of (3,n)
    # barycenters[:,i] is the ith barycenter
    #
    barycenters = embryo.get_cell_barycenter(timepoint)
    sym_barycenters = np.zeros_like(barycenters)

    for p in candidates:
        p['score_from_symmetry_icp'] = 0
        p['average_from_symmetry_icp'] = np.zeros(3)
        for i in range(barycenters.shape[1]):
            #
            # G: embryo barycenter, n: symmetry plane normal
            # search M' symmetric of M wrt (G, n)
            # GM =  (GM - (GM.n)n) + (GM.n)n
            # GM' = (GM - (GM.n)n) - (GM.n)n = GM - 2 (GM.n)n
            # M' = M - 2 (GM.n)n
            #
            sp = np.dot(p['vector'], barycenters[:, i] - embryo_barycenter)
            sym_barycenters[:, i] = barycenters[:, i] - 2 * sp * p['vector']

        #
        # register cell barycenters with their symmetric counterpart
        #
        res_rigid_mat = icp.icp(ref=barycenters, flo=sym_barycenters, transformation_type="rigid")

        # compute transformed floating points
        trsfs_barycenters = np.zeros((3, sym_barycenters.shape[1]))
        v = np.ones(4)
        for j in range(sym_barycenters.shape[1]):
            v[:3] = sym_barycenters[:, j]
            trsfs_barycenters[:, j] = (np.matmul(res_rigid_mat, v))[:3]

        # find reference barycenters that are closest to floating ones
        # i_flo_to_ref[i] is the index (in the list) of the reference cell closest
        # to the ith floating cell in cell_tflo list
        nghbref = skn.NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(barycenters.T)
        d_flo_to_ref, i_flo_to_ref = nghbref.kneighbors(trsfs_barycenters.T)
        i_flo_to_ref = i_flo_to_ref.flatten()

        # find floating barycenters that are closest to reference ones
        nghbflo = skn.NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(trsfs_barycenters.T)
        d_ref_to_flo, i_ref_to_flo = nghbflo.kneighbors(barycenters.T)
        i_ref_to_flo = i_ref_to_flo.flatten()

        # find cell/barycenters which are closest to each other
        reciprocal_n = 0
        # residual = 0.0
        for j in range(trsfs_barycenters.shape[1]):
            k = i_flo_to_ref[j]
            if i_ref_to_flo[k] != j:
                continue
            # there is a pairing between barycenters[:, j] and barycenters[:, k]
            # they should not be the same point
            if j == k:
                continue
            reciprocal_n += 1
            # residual += d_flo_to_ref[j] + d_ref_to_flo[k]
            # print("   " + str(barycenters[:, j]) + " <-> " + str(barycenters[:, k]))
            v = barycenters[:, j] - barycenters[:, k]
            v /= np.linalg.norm(v)
            sp = np.dot(p['vector'], v)
            if sp > 0:
                p['score_from_symmetry_icp'] += sp
                p['average_from_symmetry_icp'] += v
            else:
                p['score_from_symmetry_icp'] -= sp
                p['average_from_symmetry_icp'] -= v

        if reciprocal_n > 0:
            p['score_from_symmetry_icp'] /= reciprocal_n
            p['average_from_symmetry_icp'] /= np.linalg.norm(p['average_from_symmetry_icp'])

    embryo.direction_distribution_maxima[timepoint] = sorted(candidates, reverse=True,
                                                             key=lambda x: x['score_from_symmetry_icp'])


def _distribution_direction_candidates(embryo, timepoint, parameters):
    proc = "_distribution_direction_candidates"

    #
    # icp between embryo and its symmetric version
    #
    if parameters.maxima_evaluation == 'symmetric_registration':
        _distribution_direction_evaluate_symicp(embryo, timepoint)
        return

    time_digits_for_cell_id = embryo.time_digits_for_cell_id
    candidates = embryo.direction_distribution_maxima[timepoint]

    #
    # evaluation by distribution value
    #
    if parameters.maxima_evaluation == 'value':
        embryo.direction_distribution_maxima[timepoint] = sorted(candidates, reverse=True, key=lambda x: x['value'])
        return

    #
    # evaluation from pairings build with cell from each side of the symmetric plane
    #
    embryo_barycenter = embryo.get_embryo_barycenter(timepoint)
    cell_barycenters = embryo.cell_barycenter
    cell_contact_surface = embryo.cell_contact_surface
    div = 10 ** time_digits_for_cell_id

    if parameters.maxima_evaluation != 'single_pairing' and parameters.maxima_evaluation != 'multiple_pairing':
        msg = "unknown distribution maxima evaluation method '" + str(parameters.maxima_evaluation) + "'"
        monitoring.to_log_and_console(proc + ": " + msg)
        embryo.direction_distribution_maxima[timepoint] = sorted(candidates, reverse=True, key=lambda x: x['value'])
        return

    #
    # get cells for the time point
    #
    cells = [c for c in cell_barycenters if (int(c) // div == timepoint) and int(c) % div != 1]
    for p in candidates:
        #
        # initialisation
        #
        if parameters.maxima_evaluation == 'single_pairing':
            p['score_from_single_pairing'] = 0
            p['vectors_from_single_pairing'] = []
            p['average_from_single_pairing'] = np.zeros(3)
        if parameters.maxima_evaluation == 'multiple_pairing':
            p['score_from_multiple_pairing'] = 0
            p['vectors_from_multiple_pairing'] = []
            p['average_from_multiple_pairing'] = np.zeros(3)

        #
        # part cells into left and right parts wrt hyperplane defined by (embryo barycenter, candidate vector)
        #
        pcells = []
        ncells = []
        for c in cells:
            sp = np.dot(p['vector'], cell_barycenters[c] - embryo_barycenter)
            if sp >= 0.0:
                pcells += [c]
            else:
                ncells += [c]
        #
        # find symmetric cell for positive cells
        #
        p_symcells = {}
        for c in pcells:
            n_neighbors = [d for d in cell_contact_surface[c] if d in ncells]
            if len(n_neighbors) == 0:
                continue
            if parameters.maxima_evaluation == 'single_pairing':
                n_maxsurface = max([cell_contact_surface[c][d] for d in n_neighbors])
                p_symcells[c] = [k for k, v in cell_contact_surface[c].items() if v == n_maxsurface][0]
            if parameters.maxima_evaluation == 'multiple_pairing':
                v = np.zeros(3)
                s = 0.0
                for d in n_neighbors:
                    v += cell_contact_surface[c][d] * cell_barycenters[d]
                    s += cell_contact_surface[c][d]
                v /= s
                v -= cell_barycenters[c]
                v /= np.linalg.norm(v)
                p['vectors_from_multiple_pairing'] += [v]
        #
        # find symmetric cell for negative cells
        #
        n_symcells = {}
        for c in ncells:
            p_neighbors = [d for d in cell_contact_surface[c] if d in pcells]
            if len(p_neighbors) == 0:
                continue
            if parameters.maxima_evaluation == 'single_pairing':
                p_maxsurface = max([cell_contact_surface[c][d] for d in p_neighbors])
                n_symcells[c] = [k for k, v in cell_contact_surface[c].items() if v == p_maxsurface][0]
            if parameters.maxima_evaluation == 'multiple_pairing':
                v = np.zeros(3)
                s = 0.0
                for d in p_neighbors:
                    v += cell_contact_surface[c][d] * cell_barycenters[d]
                    s += cell_contact_surface[c][d]
                v /= s
                v -= cell_barycenters[c]
                v /= np.linalg.norm(v)
                v *= (-1)
                p['vectors_from_multiple_pairing'] += [v]

        #
        # evaluation = average of scalar products
        #
        if parameters.maxima_evaluation == 'single_pairing':
            n = 0
            for c in p_symcells:
                if p_symcells[c] not in n_symcells:
                    continue
                if n_symcells[p_symcells[c]] != c:
                    continue
                n += 1
                v = cell_barycenters[c] - cell_barycenters[p_symcells[c]]
                v /= np.linalg.norm(v)
                p['vectors_from_single_pairing'] += [v]
                sp = np.dot(p['vector'], v)
                if sp > 0:
                    p['score_from_single_pairing'] += sp
                    p['average_from_single_pairing'] += v
                else:
                    p['score_from_single_pairing'] -= sp
                    p['average_from_single_pairing'] -= v
            if n > 0:
                p['score_from_single_pairing'] /= n
                p['average_from_single_pairing'] /= np.linalg.norm(p['average_from_single_pairing'])
                del p['vectors_from_single_pairing']

        #
        #
        #
        if parameters.maxima_evaluation == 'multiple_pairing':
            for v in p['vectors_from_multiple_pairing']:
                sp = np.dot(p['vector'], v)
                if sp > 0:
                    p['score_from_multiple_pairing'] += sp
                    p['average_from_multiple_pairing'] += v
                else:
                    p['score_from_multiple_pairing'] -= sp
                    p['average_from_multiple_pairing'] -= v
            p['score_from_multiple_pairing'] /= len(p['vectors_from_multiple_pairing'])
            p['average_from_multiple_pairing'] /= np.linalg.norm(p['average_from_multiple_pairing'])
            del p['vectors_from_multiple_pairing']

    if parameters.maxima_evaluation == 'single_pairing':
        embryo.direction_distribution_maxima[timepoint] = sorted(candidates, reverse=True,
                                                                 key=lambda x: x['score_from_single_pairing'])
        return
    if parameters.maxima_evaluation == 'multiple_pairing':
        embryo.direction_distribution_maxima[timepoint] = sorted(candidates, reverse=True,
                                                                 key=lambda x: x['score_from_multiple_pairing'])
    embryo.direction_distribution_maxima[timepoint] = sorted(candidates, reverse=True, key=lambda x: x['value'])
    return


############################################################
#
#
#
############################################################

class AtlasParameters(udiagnosis.DiagnosisParameters, EmbryoSymmetryParameters):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix=None):

        if "doc" not in self.__dict__:
            self.doc = {}

        local_prefix = ["distribution_"]
        if prefix is not None:
            if isinstance(prefix, str):
                local_prefix = [prefix, "distribution_"]
            elif isinstance(prefix, list):
                local_prefix = prefix + ["distribution_"]
        udiagnosis.DiagnosisParameters.__init__(self, prefix=local_prefix)
        EmbryoSymmetryParameters.__init__(self, prefix=local_prefix)

        #
        # atlas general parameters
        #

        doc = "\t List of atlas files. An atlas file is a property file that contains lineage,\n"
        doc += "\t names, and contact surfaces for an embryo."
        self.doc['atlasFiles'] = doc
        self.atlasFiles = []

        doc = "\t Reference atlas. Use for time alignment. If not provide, the first atlas of\n"
        doc += "\t 'atlasFiles' is used as reference. Warning, the reference atlas has to be in\n"
        doc += "\t 'atlasFiles' list also."
        self.doc['referenceAtlas'] = doc
        self.referenceAtlas = None

        doc = "\t Output directory where to write atlas-individualized output files,\n"
        doc += "\t ie morphonet selection files or figure files."
        self.doc['outputDir'] = doc
        self.outputDir = "."

        doc = "\t True or False. Performs some diagnosis when reading an additional property file \n"
        doc += "\t into the atlases. Incrementing the verboseness ('-v' in the command line) may give\n"
        doc += "\t more details."
        self.doc['atlas_diagnosis'] = doc
        self.atlas_diagnosis = False

        doc = "\t if True, generate python files (prefixed by 'figures_') that generate figures.\n"
        doc += "\t Those files will be saved into the 'outputDir' directory.\n"
        doc += "\t 'generate_figure' can be\n"
        doc += "\t - a boolean value: if True, all figure files are generated; if False, none of them\n"
        doc += "\t - a string: if 'all', all figure files are generated; else, only the specified\n"
        doc += "\t   figure file is generated (see below for the list)\n"
        doc += "\t - a list of strings: if 'all' is in the list, all figure files are generated;\n"
        doc += "\t   else, only the specified figure files are generated (see below for the list)\n"
        doc += "\t List of figures:\n"
        doc += "\t 'cell-distance-along-branch': plot the cell-to-cell distance between successive\n"
        doc += "\t    along a branch (a cell without division) wrt the distance to the first cell.\n"
        doc += "\t    Cell neighborhoods are expressed with the neighbors of the first cell of the branch\n"
        doc += "\t    (thus it ignores divisions occurring in the cell neighborhood during the cell life).\n"
        doc += "\t 'cell-number-wrt-time': plot the number of cells wrt time point (ie image indices)\n"
        doc += "\t    without and with temporal registration (allows to assess the temporal registration)\n"
        doc += "\t 'neighbors-wrt-cell-number': plot the cell number in the cell neighborhood wrt\n"
        doc += "\t    the total cell number in the embryo\n"
        doc += "\t 'cell-distance-histograms': plot cell-to-cell distance histograms.\n"
        doc += "\t    warning: it may be long.\n"
        doc += "\t 'division-distance-histograms': plot division-to-division distance histograms.\n"
        doc += "\t 'distance-histograms': plot cell-to-cell distance histograms, \n"
        doc += "\t    as well as division-to-division distance histograms.\n"
        doc += "\t    warning: it may be long.\n"
        doc += "\t 'division-dendrograms': draw a dendrogram per division where atlases are grouped with\n"
        doc += "\t    distance between divisions\n"
        doc += "\t 'embryo-volume': plot the embryo volume (in voxel)\n"
        doc += "\t    without and with temporal registration (computed from cell number)\n"
        doc += "\t 'symmetry - axis':\n"
        self.doc['generate_figure'] = doc
        self.generate_figure = False

        doc = "\t suffix used to named the above python files as well as the generated figures."
        self.doc['figurefile_suffix'] = doc
        self.figurefile_suffix = ""

        #
        # parameters dedicated to extract neighborhoods
        #

        doc = "\t Delay from the division to extract the neighborhooods used for atlas building,\n"
        doc += "\t and thus for naming.\n"
        doc += "\t 0 means right after the division.\n"
        doc += "\t negative values means that the delay is counted backwards from the end of the branch.\n"
        self.doc['name_delay_from_division'] = doc
        self.name_delay_from_division = 0

        doc = "\t Delay from the division to extract the neighborhooods used for naming confidence.\n"
        doc += "\t 0 means right after the division.\n"
        doc += "\t negative values means that the delay is counted backwards from the end of the branch.\n"
        self.doc['confidence_delay_from_division'] = doc
        self.confidence_delay_from_division = None

        #
        # parameters dedicated to build neighborhoods
        #
        doc = "\t if 'True', add the symmetric neighborhood as additional exemplar.\n"
        self.doc['add_symmetric_neighborhood'] = doc
        self.add_symmetric_neighborhood = True

        doc = "\t if 'True', differentiate the cells of the symmetric half-embryo.\n"
        doc += "\t If 'False', consider all the cells of the symmetric half-embryo\n"
        doc += "\t as a single cell.\n"
        self.doc['differentiate_other_half'] = doc
        self.differentiate_other_half = True

        doc = "\t The same cell has different neighbors from an atlas to the other.\n"
        doc += "\t If 'True' build and keep an unique common neighborhood (set of\n"
        doc += "\t neighbors) for all atlases by keeping the closest ancestor for\n"
        doc += "\t neighboring cells. Eg, if a division has occurred in some embryos\n"
        doc += "\t and not in others, daughter cells will be fused so that all\n"
        doc += "\t neighborhoods only exhibit the parent cell."
        self.doc['use_common_neighborhood'] = doc
        self.use_common_neighborhood = True

        #
        #
        #
        doc = "\t - None\n"
        doc += "\t - 'local': normalization by the cell surface\n"
        doc += "\t   The normalization factor is then different from cell to cell within a embryo.\n"
        doc += "\t - 'global': normalization by embryo volume\n"
        doc += "\t   The normalization factor is the same from cell to cell within a embryo.\n"
        self.doc['cell_normalization'] = doc
        self.cell_normalization = 'global'

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('#')
        print('# AtlasParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        udiagnosis.DiagnosisParameters.print_parameters(self)
        EmbryoSymmetryParameters.print_parameters(self)

        self.varprint('atlasFiles', self.atlasFiles)
        self.varprint('referenceAtlas', self.referenceAtlas)
        self.varprint('outputDir', self.outputDir)

        self.varprint('atlas_diagnosis', self.atlas_diagnosis)

        self.varprint('generate_figure', self.generate_figure)
        self.varprint('figurefile_suffix', self.figurefile_suffix)

        self.varprint('name_delay_from_division', self.name_delay_from_division)
        self.varprint('confidence_delay_from_division', self.confidence_delay_from_division)

        self.varprint('add_symmetric_neighborhood', self.add_symmetric_neighborhood)
        self.varprint('differentiate_other_half', self.differentiate_other_half)
        self.varprint('use_common_neighborhood', self.use_common_neighborhood)

        self.varprint('cell_normalization', self.cell_normalization)
        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# AtlasParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        udiagnosis.DiagnosisParameters.write_parameters_in_file(self, logfile)
        EmbryoSymmetryParameters.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'atlasFiles', self.atlasFiles, self.doc.get('atlasFiles', None))
        self.varwrite(logfile, 'referenceAtlas', self.referenceAtlas, self.doc.get('referenceAtlas', None))
        self.varwrite(logfile, 'outputDir', self.outputDir, self.doc.get('outputDir', None))

        self.varwrite(logfile, 'atlas_diagnosis', self.atlas_diagnosis, self.doc.get('atlas_diagnosis', None))

        self.varwrite(logfile, 'generate_figure', self.generate_figure, self.doc.get('generate_figure', None))
        self.varwrite(logfile, 'figurefile_suffix', self.figurefile_suffix, self.doc.get('figurefile_suffix', None))

        self.varwrite(logfile, 'name_delay_from_division', self.name_delay_from_division,
                      self.doc.get('name_delay_from_division', None))
        self.varwrite(logfile, 'confidence_delay_from_division', self.confidence_delay_from_division,
                      self.doc.get('confidence_delay_from_division', None))

        self.varwrite(logfile, 'add_symmetric_neighborhood', self.add_symmetric_neighborhood,
                      self.doc.get('add_symmetric_neighborhood', None))
        self.varwrite(logfile, 'differentiate_other_half', self.differentiate_other_half,
                      self.doc.get('differentiate_other_half', None))
        self.varwrite(logfile, 'use_common_neighborhood', self.use_common_neighborhood,
                      self.doc.get('use_common_neighborhood', None))

        self.varwrite(logfile, 'cell_normalization', self.cell_normalization, self.doc.get('cell_normalization', None))
        logfile.write("\n")
        return

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

        udiagnosis.DiagnosisParameters.update_from_parameters(self, parameters)
        EmbryoSymmetryParameters.update_from_parameters(self, parameters)

        self.atlasFiles = self.read_parameter(parameters, 'atlasFiles', self.atlasFiles)
        self.atlasFiles = self.read_parameter(parameters, 'referenceFiles', self.atlasFiles)
        self.referenceAtlas = self.read_parameter(parameters, 'referenceAtlas', self.referenceAtlas)

        self.outputDir = self.read_parameter(parameters, 'outputDir', self.outputDir)

        self.atlas_diagnosis = self.read_parameter(parameters, 'atlas_diagnosis', self.atlas_diagnosis)
        self.atlas_diagnosis = self.read_parameter(parameters, 'diagnosis_properties', self.atlas_diagnosis)
        self.atlas_diagnosis = self.read_parameter(parameters, 'naming_diagnosis', self.atlas_diagnosis)
        self.atlas_diagnosis = self.read_parameter(parameters, 'diagnosis_naming', self.atlas_diagnosis)

        self.generate_figure = self.read_parameter(parameters, 'generate_figure', self.generate_figure)
        self.figurefile_suffix = self.read_parameter(parameters, 'figurefile_suffix', self.figurefile_suffix)

        self.name_delay_from_division = self.read_parameter(parameters, 'name_delay_from_division',
                                                            self.name_delay_from_division)
        self.name_delay_from_division = self.read_parameter(parameters, 'delay_from_division',
                                                            self.name_delay_from_division)
        self.confidence_delay_from_division = self.read_parameter(parameters, 'confidence_delay_from_division',
                                                                  self.confidence_delay_from_division)
        self.confidence_delay_from_division = self.read_parameter(parameters, 'delay_from_division',
                                                                  self.confidence_delay_from_division)

        self.add_symmetric_neighborhood = self.read_parameter(parameters, 'add_symmetric_neighborhood',
                                                              self.add_symmetric_neighborhood)
        self.differentiate_other_half = self.read_parameter(parameters, 'differentiate_other_half',
                                                            self.differentiate_other_half)
        self.use_common_neighborhood = self.read_parameter(parameters, 'use_common_neighborhood',
                                                           self.use_common_neighborhood)

        self.cell_normalization = self.read_parameter(parameters, 'cell_normalization', self.cell_normalization)

    def update_from_parameter_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = common.load_source(parameter_file)
        self.update_from_parameters(parameters)


############################################################
#
# Atlas = one property file
#
############################################################

class Atlas(object):
    def __init__(self, atlas_properties=None, time_digits_for_cell_id=4, verbose=False):
        proc = "Atlas.init"

        self.time_digits_for_cell_id = time_digits_for_cell_id
        self._properties = {'temporal_alignment': (1.0, 0.0), 'volume_local_estimation': (0.0, 1.0),
                            'target_volume': 6000000}

        if isinstance(atlas_properties, dict):
            if 'cell_lineage' in atlas_properties:
                self.cell_lineage = atlas_properties['cell_lineage']
            elif verbose:
                monitoring.to_log_and_console(str(proc) + ": 'cell_lineage' was not in dictionary")
            if 'cell_name' in atlas_properties:
                self.cell_name = atlas_properties['cell_name']
            elif verbose:
                monitoring.to_log_and_console(str(proc) + ": 'cell_name' was not in dictionary")
            if 'cell_contact_surface' in atlas_properties:
                self.cell_contact_surface = atlas_properties['cell_contact_surface']
            elif verbose:
                monitoring.to_log_and_console(str(proc) + ": 'cell_contact_surface' was not in dictionary")
            if 'cell_barycenter' in atlas_properties:
                self.cell_barycenter = atlas_properties['cell_barycenter']
            elif verbose:
                monitoring.to_log_and_console(str(proc) + ": 'cell_barycenter' was not in dictionary")
            if 'cell_volume' in atlas_properties:
                self.cell_volume = atlas_properties['cell_volume']
                self._volume_local_fitting()
            elif verbose:
                monitoring.to_log_and_console(str(proc) + ": 'cell_volume' was not in dictionary")

    ############################################################
    #
    # Property management
    #
    ############################################################

    def property_list(self):
        return self._properties.keys()

    def get_property(self):
        return self._properties

    def _del_one_property(self, property_name):
        if not isinstance(property_name, str):
            return
        if property_name in self._properties:
            del self._properties[property_name]
        if property_name == 'direction_distribution':
            self._del_one_property('direction_distribution_maxima')
            self._del_one_property('direction_distribution_candidates')
            if property_name == 'direction_distribution_maxima':
                self._del_one_property('direction_distribution_candidates')
        return

    def del_property(self, property_name):
        if isinstance(property_name, str):
            self._del_one_property(property_name)
        elif isinstance(property_name, list):
            for n in property_name:
                self._del_one_property(n)
        return

    ############################################################
    #
    # @property
    # From property file
    # - cell_lineage
    # - cell_name
    # - cell_volume
    # - cell_contact_surface
    # - cell_barycenter
    # Computed ones
    # - temporal_alignment: tuple (a, b). Time (at+b) of the embryo corresponds to
    #   the time t of a reference embryo
    # - volume_local_estimation: tuple (a, b). Embryo volume along time is fitted by (at+b)
    #   allows to compute a time-varying scaling factor to get a constant embryo volume
    #   along time (see self.get_voxelsize_correction())
    # - target_volume: this is the targeted contant volume to get a constant embryo volume
    #   along time (see self.get_voxelsize_correction())
    # - direction_distribution
    # - direction_distribution_maxima
    # - direction_distribution_candidates
    #
    ############################################################

    @property
    def cell_lineage(self):
        """
        The cell lineage, as in the property file
        Returns
        -------

        """
        if 'cell_lineage' in self._properties:
            return self._properties['cell_lineage']
        return None

    @cell_lineage.setter
    def cell_lineage(self, atlas_properties):
        self._properties['cell_lineage'] = copy.deepcopy(atlas_properties)
        return

    @property
    def cell_name(self):
        """
        The cell names, as in the property file
        Returns
        -------

        """
        if 'cell_name' in self._properties:
            return self._properties['cell_name']
        return None

    @cell_name.setter
    def cell_name(self, atlas_properties):
        self._properties['cell_name'] = copy.deepcopy(atlas_properties)
        return

    @property
    def cell_volume(self):
        """
        The cell volumes, as in the property file
        Returns
        -------

        """
        if 'cell_volume' in self._properties:
            return self._properties['cell_volume']
        return None

    @cell_volume.setter
    def cell_volume(self, atlas_properties):
        self._properties['cell_volume'] = copy.deepcopy(atlas_properties)
        return

    @property
    def cell_contact_surface(self):
        """
        The cell contact surfaces, as in the property file
        Returns
        -------

        """
        if 'cell_contact_surface' in self._properties:
            return self._properties['cell_contact_surface']
        return None

    @cell_contact_surface.setter
    def cell_contact_surface(self, atlas_properties):
        self._properties['cell_contact_surface'] = copy.deepcopy(atlas_properties)
        return

    @property
    def cell_barycenter(self):
        """
        The cell contact surfaces, as in the property file
        Returns
        -------

        """
        if 'cell_barycenter' in self._properties:
            return self._properties['cell_barycenter']
        return None

    @cell_barycenter.setter
    def cell_barycenter(self, atlas_properties):
        self._properties['cell_barycenter'] = copy.deepcopy(atlas_properties)
        return

    #
    # other properties
    #
    @property
    def temporal_alignment(self):
        if 'temporal_alignment' in self._properties:
            return self._properties['temporal_alignment']
        return None

    @property
    def volume_local_estimation(self):
        if 'volume_local_estimation' in self._properties:
            return self._properties['volume_local_estimation']
        return None

    @property
    def target_volume(self):
        if 'target_volume' in self._properties:
            return self._properties['target_volume']
        return None

    @target_volume.setter
    def target_volume(self, volume):
        self._properties['target_volume'] = volume
        return

    @property
    def direction_distribution(self):
        if 'direction_distribution' not in self._properties:
            self._properties['direction_distribution'] = {}
        return self._properties['direction_distribution']

    @property
    def direction_distribution_maxima(self):
        if 'direction_distribution_maxima' not in self._properties:
            self._properties['direction_distribution_maxima'] = {}
        return self._properties['direction_distribution_maxima']

    ############################################################
    #
    # Properties (end)
    #
    ############################################################

    #
    # temporal alignment of an atlas/embryo with the reference atlas/embryo
    # it consists at finding the time lineage warping based on the cell number
    # so that the cell number at (a * t + b) of the atlas is equal at the one
    # of the reference number at (t)
    #
    def temporally_align_with(self, reference):
        """

        Parameters
        ----------
        reference : Atlas
            reference atlas to be temporally aligned with.
        time_digits_for_cell_id : int

        Returns
        -------
            Tuple (a, b). Time point t of the atlas at hand is equivalent to the time
            point (at+b) of the reference atlas.

        """
        proc = "Atlas.temporally_align_with"
        if not isinstance(reference, Atlas):
            monitoring.to_log_and_console(str(proc) + ": 'reference' should be of 'Atlas' class")
            return
        a, b = properties.temporal_alignment(reference.cell_lineage, self.cell_lineage,
                                             reference.time_digits_for_cell_id, self.time_digits_for_cell_id)
        self._properties['temporal_alignment'] = (a, b)
        return

    #
    #
    #

    def _volume_local_fitting(self):
        div = 10 ** self.time_digits_for_cell_id
        volume = self.cell_volume
        volume_along_time = {}
        #
        # compute embryo volume for each timepoint 't'
        #
        for c in volume:
            t = int(c) // div
            volume_along_time[t] = volume_along_time.get(t, 0) + volume[c]

        #
        # get array of time point (x) and array of volunes (y)
        #
        x = list(volume_along_time.keys())
        x = sorted(x)
        y = [volume_along_time[i] for i in x]

        #
        # robust regression via ransac
        #
        xnp = np.array(x)[:, np.newaxis]
        ynp = np.array(y)[:, np.newaxis]
        ransac = sklm.RANSACRegressor()
        ransac.fit(xnp, ynp)
        self._properties['volume_local_estimation'] = (ransac.estimator_.coef_[0][0], ransac.estimator_.intercept_[0])

    #
    #
    #

    def get_voxelsize_correction(self, timepoint, target_volume=60000000):
        if 'volume_local_estimation' not in self._properties:
            self._volume_local_fitting()
        if target_volume != self.target_volume:
            self.target_volume = target_volume
        v_coefficients = self.volume_local_estimation
        t_volume = v_coefficients[0] * timepoint + v_coefficients[1]
        return np.cbrt(self.target_volume / t_volume)

    #
    #
    #
    def get_embryo_volume(self, timepoint):
        s = 0.0
        volumes = self.cell_volume
        div = 10 ** self.time_digits_for_cell_id
        for c in volumes:
            if int(c) // div != timepoint:
                continue
            if int(c) % div == 1 or int(c) % div == 0:
                continue
            s += volumes[c]
        return s

    #
    #
    #
    def get_symmetry_axis_from_names(self, timepoint):

        if 'symmetry_axis_from_names' in self._properties:
            if timepoint in self._properties['symmetry_axis_from_names']:
                return self._properties['symmetry_axis_from_names'][timepoint]
        else:
            self._properties['symmetry_axis_from_names'] = {}

        volumes = self.cell_volume
        barycenters = self.cell_barycenter
        names = self.cell_name
        if names is None:
            return np.zeros(3)
        #
        leftb = np.zeros(3)
        lefts = 0.0
        rightb = np.zeros(3)
        rights = 0.0
        div = 10 ** self.time_digits_for_cell_id
        #
        for c in volumes:
            if int(c) // div != timepoint:
                continue
            if int(c) % div == 1 or int(c) % div == 0:
                continue
            if c not in barycenters or c not in names:
                continue
            if names[c].split('.')[1][4] == '_':
                leftb += volumes[c] * barycenters[c]
                lefts += volumes[c]
            elif names[c].split('.')[1][4] == '*':
                rightb += volumes[c] * barycenters[c]
                rights += volumes[c]
        if lefts == 0.0 or rights == 0.0:
            return np.zeros(3)
        symdir = leftb/lefts - rightb/rights
        self._properties['symmetry_axis_from_names'][timepoint] = symdir / np.linalg.norm(symdir)
        return self._properties['symmetry_axis_from_names'][timepoint]

    #
    #
    #
    def get_embryo_barycenter(self, timepoint):
        if 'embryo_barycenter' in self._properties:
            if timepoint in self._properties['embryo_barycenter']:
                return self._properties['embryo_barycenter'][timepoint]
        else:
            self._properties['embryo_barycenter'] = {}

        b = np.zeros(3)
        s = 0.0
        volumes = self.cell_volume
        barycenters = self.cell_barycenter
        div = 10 ** self.time_digits_for_cell_id
        for c in volumes:
            if int(c) // div != timepoint:
                continue
            if int(c) % div == 1 or int(c) % div == 0:
                continue
            if c not in barycenters:
                continue
            b += volumes[c] * barycenters[c]
            s += volumes[c]
        self._properties['embryo_barycenter'][timepoint] = b / s
        return self._properties['embryo_barycenter'][timepoint]

    #
    #
    #
    def get_cell_barycenter(self, timepoint):
        barycenters = self.cell_barycenter
        div = 10 ** self.time_digits_for_cell_id
        lc = [c for c in barycenters if (int(c) // div == timepoint) and int(c) % div != 1]
        b = np.zeros((3, len(lc)))
        for i, c in enumerate(lc):
            b[:, i] = barycenters[c]
        return b

    #
    #
    #
    def get_direction_distribution_candidates(self, timepoint, parameters=None):
        proc = "get_direction_distribution_candidates"
        if parameters is not None:
            if not isinstance(parameters, EmbryoSymmetryParameters):
                msg = "parameters expected type is 'EmbryoSymmetryParameters' but was '"
                msg += type(parameters) + "' instead"
                monitoring.to_log_and_console(proc + ": " + msg)
                return None
            local_parameters = parameters
        else:
            local_parameters = EmbryoSymmetryParameters()

        #
        # build direction distribution support
        #
        # distribution is an array of dictionaries
        # 'voxel': (i, j, k) voxel coordinates in a sxsxs matrix with s computed by
        #           ir = math.ceil(parameters.distribution_sphere_radius)
        #           s = int(2 * ir + 1 + 2)
        # 'vector': a normalized vector v / np.linalg.norm(v)  with v = 'voxel' - center
        #           center coordinates are (c,c,c) with c = np.double(ir + 1)
        # 'value': the distribution values
        #          the maximal value is 1
        #
        if 'direction_distribution' not in self._properties:
            self._properties['direction_distribution'] = {}
        if timepoint not in self.direction_distribution:
            self.direction_distribution[timepoint] = _distribution_build_support(parameters=local_parameters)
            #
            # fill direction distribution from embryo properties
            #
            _distribution_direction_build(self, timepoint, parameters)

        #
        # array of dictionaries, with only the local maxima (with value > 0.5)
        #
        if 'direction_distribution_maxima' not in self._properties:
            self._properties['direction_distribution_maxima'] = {}
        if timepoint not in self.direction_distribution_maxima:
            _distribution_direction_maxima(self, timepoint, parameters)

        #
        # sort maxima
        #
        _distribution_direction_candidates(self, timepoint, parameters)

        return self.direction_distribution_maxima[timepoint]


############################################################
#
# Atlases: set of many atlas
#
############################################################

class Atlases(object):
    def __init__(self, parameters=None):
        # atlases
        self._atlases = {}
        # reference atlas for time alignment
        self._ref_atlas = None

        if parameters is not None:
            self.update_from_parameters(parameters)

    ############################################################
    #
    # getters / setters
    #
    ############################################################

    def get_atlases(self):
        return self._atlases

    def set_atlases(self, atlases):
        self._atlases = atlases

    def get_reference_atlas(self):
        proc = "get_reference_atlas"
        if self._ref_atlas is None:
            atlases = self.get_atlases()
            if len(atlases) == 0:
                monitoring.to_log_and_console(proc + ": no reference atlas nor registered atlases")
                return None
            names = sorted(list(atlases.keys()))
            monitoring.to_log_and_console("   ... set reference to '" + str(names[0]) + "'")
            self.set_reference_atlas(names[0])
        return self._ref_atlas

    def set_reference_atlas(self, reference_atlas):
        self._ref_atlas = reference_atlas

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameters):
        if not isinstance(parameters, AtlasParameters):
            return
        if parameters.referenceAtlas is not None:
            name = parameters.referenceAtlas.split(os.path.sep)[-1]
            if name.endswith(".xml") or name.endswith(".pkl"):
                name = name[:-4]
            self.set_reference_atlas(name)

    ############################################################
    #
    #
    #
    ############################################################

    def temporal_alignment(self):
        proc = "temporal_alignment"
        ref_atlas = self.get_reference_atlas()
        if ref_atlas is None:
            monitoring.to_log_and_console(proc + ": no reference atlas, can not perform temporal alignment")
            return
        atlases = self.get_atlases()
        if ref_atlas not in atlases:
            msg = "'" + str(ref_atlas) + "' is not in registered atlases, can not perform temporal alignment"
            monitoring.to_log_and_console(proc + ": " + msg)
            return
        for a in atlases:
            if a == ref_atlas:
                continue
            atlases[a].temporally_align_with(atlases[ref_atlas])
        return

    ############################################################
    #
    #
    #
    ############################################################

    def _add_atlas(self, prop, name, time_digits_for_cell_id=4):
        atlases = self.get_atlases()
        atlases[name] = Atlas(prop, time_digits_for_cell_id=time_digits_for_cell_id)
        self.set_atlases(atlases)

    def add_atlases(self, atlasfiles, parameters, time_digits_for_cell_id=4):
        """

        Parameters
        ----------
        atlasfiles: a file name or a liste of file names. Pkl or xml files containing embryo properties
        parameters:
        time_digits_for_cell_id: number of digits to encode cell label value (left padded with zero

        Returns
        -------
        Keeps a copy of 'cell_lineage', 'cell_name', 'cell_contact_surface' and 'cell_contact_surface'
        of each atlas/embryo of the list

        """
        proc = "add_atlases"

        if not isinstance(parameters, AtlasParameters):
            monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                          + str(type(parameters)))
            sys.exit(1)

        if isinstance(atlasfiles, str):
            prop = ioproperties.read_dictionary(atlasfiles, inputpropertiesdict={})
            name = atlasfiles.split(os.path.sep)[-1]
            if name.endswith(".xml") or name.endswith(".pkl"):
                name = name[:-4]
            if parameters.atlas_diagnosis:
                udiagnosis.diagnosis(prop, ['name', 'contact'], parameters,
                                     time_digits_for_cell_id=time_digits_for_cell_id)
            self._add_atlas(prop, name, time_digits_for_cell_id=time_digits_for_cell_id)
            del prop
        elif isinstance(atlasfiles, list):
            if len(atlasfiles) == 0:
                monitoring.to_log_and_console(str(proc) + ": empty atlas file list ?!")
                sys.exit(1)
            for f in atlasfiles:
                prop = ioproperties.read_dictionary(f, inputpropertiesdict={})
                name = f.split(os.path.sep)[-1]
                if name.endswith(".xml") or name.endswith(".pkl"):
                    name = name[:-4]
                if parameters.atlas_diagnosis:
                    udiagnosis.diagnosis(prop, ['name', 'contact'], parameters,
                                         time_digits_for_cell_id=time_digits_for_cell_id)
                self._add_atlas(prop, name, time_digits_for_cell_id=time_digits_for_cell_id)
                del prop

        #
        # temporal alignment (done from the cell number)
        #
        monitoring.to_log_and_console("... temporal alignment of lineages", 1)
        self.temporal_alignment()
        atlases = self.get_atlases()
        for n in atlases:
            msg = "   ... "
            msg += "linear time warping of '" + str(n) + "' wrt '" + str(self._ref_atlas) + "' is "
            msg += "({:.3f}, {:.3f})".format(atlases[n].temporal_alignment[0], atlases[n].temporal_alignment[1])
            monitoring.to_log_and_console(msg, 1)

    ############################################################
    #
    #
    #
    ############################################################

    def generate_figure(self, parameters):
        generate_figure = (isinstance(parameters.generate_figure, bool) and parameters.generate_figure) or \
                          (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'all') or \
                          (isinstance(parameters.generate_figure, list) and 'all' in parameters.generate_figure)

        #
        # plot cell number wrt time without and with temporal registration
        #
        if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'cell-number-wrt-time') \
                or (isinstance(parameters.generate_figure, list)
                    and 'cell-number-wrt-time' in parameters.generate_figure) \
                or generate_figure:
            monitoring.to_log_and_console("... generate cell number wrt time file", 1)
            _figures_temporal_registration(self, parameters)
            monitoring.to_log_and_console("... done", 1)

        #
        # plot embryo volume wrt time without and with temporal registration
        #
        if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'embryo-volume') \
                or (isinstance(parameters.generate_figure, list) and 'embryo-volume' in parameters.generate_figure) \
                or generate_figure:
            monitoring.to_log_and_console("... generate embryo volume figure file", 1)
            _figures_embryo_volume(self, parameters)
            monitoring.to_log_and_console("... done", 1)

        #
        # cell neighbors number wrt total number of cells in the embryo
        #
        if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'neighbors-wrt-cell-number') \
                or (isinstance(parameters.generate_figure, list)
                    and 'neighbors-wrt-cell-number' in parameters.generate_figure) \
                or generate_figure:
            monitoring.to_log_and_console("... generate neighbors histogram figure file", 1)
            _figures_neighbor_histogram(self, parameters)
            monitoring.to_log_and_console("... done", 1)

        if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'symmetry-axis') \
                or (isinstance(parameters.generate_figure, list) and 'symmetry-axis' in parameters.generate_figure) \
                or generate_figure:
            monitoring.to_log_and_console("... generate symmetry axis figure file", 1)
            _figure_symmetry_axis(self, parameters)
            monitoring.to_log_and_console("... done", 1)


################################################################################
#
# cell number wrt time without and with temporal registration
#
################################################################################

def _figures_temporal_registration(atlases, parameters):
    """
    Plot cell number curves without and with linear temporal alignment
    Parameters
    ----------
    atlases
    parameters

    Returns
    -------

    """
    proc = "_figures_temporal_registration"

    filename = 'figures_temporal_registration'
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

    ref_atlases = atlases.get_atlases()
    atlas_names = list(ref_atlases.keys())

    cells_per_time = {}
    temporal_coefficients = {}
    for n in atlas_names:
        lineage = ref_atlases[n].cell_lineage
        cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
        cells = sorted(cells)
        div = 10 ** ref_atlases[n].time_digits_for_cell_id
        cells_per_time[n] = {}
        for c in cells:
            t = int(c) // div
            cells_per_time[n][t] = cells_per_time[n].get(t, 0) + 1
        temporal_coefficients[n] = ref_atlases[n].temporal_alignment

    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    f.write("cells_per_time = " + str(cells_per_time) + "\n")
    f.write("temporal_coefficients = " + str(temporal_coefficients) + "\n")

    f.write("\n")
    f.write("fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True, figsize=(16, 8))\n")
    f.write("labels = []\n")
    f.write("for n in cells_per_time:\n")
    f.write("    labels += [n]\n")
    f.write("    x = list(cells_per_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    y = [cells_per_time[n][i] for i in x]\n")
    f.write("    ax1.plot(x, y)\n")
    f.write("ax1.set_title(\"cell number (without alignment)\", fontsize=15)\n")
    f.write("ax1.legend(labels, prop={'size': 10})\n")

    f.write("for n in cells_per_time:\n")
    f.write("    x = list(cells_per_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    t = [temporal_coefficients[n][0] * i + temporal_coefficients[n][1] for i in x]\n")
    f.write("    y = [cells_per_time[n][i] for i in x]\n")
    f.write("    ax2.plot(t, y)\n")
    f.write("ax2.set_title(\"cell number (with alignment)\", fontsize=15)\n")
    f.write("ax2.legend(labels, prop={'size': 10})\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('temporal_alignment_cell_number")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")
    f.close()


################################################################################
#
# cell-to-cell distance between successive cells in a branch with respect to distance from first cell
# (from first time point/division to last time point/division)
#
################################################################################

def _figures_embryo_volume(atlases, parameters):
    """

    Parameters
    ----------
    atlases
    parameters

    Returns
    -------

    """
    proc = "_figures_embryo_volume"

    filename = 'figures_embryo_volume'
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

    ref_atlases = atlases.get_atlases()
    atlas_names = list(ref_atlases.keys())

    volume_along_time = {}
    temporal_coefficients = {}
    volume_local_estimation = {}
    surface_along_time = {}
    voxelsize_correction = {}
    for n in atlas_names:
        volume_along_time[n] = {}
        surface_along_time[n] = {}
        voxelsize_correction[n] = {}
        volume = ref_atlases[n].cell_volume
        contact_surface = ref_atlases[n].cell_contact_surface
        div = 10 ** ref_atlases[n].time_digits_for_cell_id
        for c in volume:
            t = int(c) // div
            volume_along_time[n][t] = volume_along_time[n].get(t, 0) + volume[c]
        for c in contact_surface:
            t = int(c) // div
            for d in contact_surface[c]:
                if d % div == 1 or d % div == 0:
                    surface_along_time[n][t] = surface_along_time[n].get(t, 0) + contact_surface[c][d]
        for t in volume_along_time[n]:
            voxelsize_correction[n][t] = ref_atlases[n].get_voxelsize_correction(t)

        temporal_coefficients[n] = ref_atlases[n].temporal_alignment
        volume_local_estimation[n] = ref_atlases[n].volume_local_estimation
    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    f.write("volume_along_time = " + str(volume_along_time) + "\n")
    f.write("\n")
    f.write("surface_along_time = " + str(surface_along_time) + "\n")
    f.write("\n")
    f.write("voxelsize = " + str(voxelsize_correction) + "\n")
    f.write("\n")
    f.write("temporal_coefficients = " + str(temporal_coefficients) + "\n")
    f.write("volume_local_estimation = " + str(volume_local_estimation) + "\n")

    f.write("\n")
    f.write("fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(24, 8))\n")
    f.write("labels = []\n")
    f.write("for n in volume_along_time:\n")
    f.write("    labels += [n]\n")
    f.write("    x = list(volume_along_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    y = [volume_along_time[n][i] for i in x]\n")
    f.write("\n")
    f.write("    y2 = [volume_local_estimation[n][0] * i + volume_local_estimation[n][1] for i in x]\n")
    f.write("    p = ax1.plot(x, y, label=n)\n")
    f.write("    ax1.plot(x, y2, color=p[0].get_color())  \n")
    f.write("\n")
    f.write("    ploss = (y[-1] - y[0])/y[0]\n")
    f.write("    ploss2 = (y2[-1] - y2[0])/y2[0]\n")
    f.write("    print(str(n) + ' volume percentage loss from first time point')\n")
    f.write("    print('\\t to the last one, from measures = ' + str(100.0 * (y[-1] - y[0])/y[0]))\n")
    f.write("    print('\\t \\t from linear estimation = ' + str(100.0 * ((y2[-1] - y2[0])/y2[0])))\n")
    f.write("    print('\\t to the 50th one, from linear estimation = ' + str(100.0 * ((y2[50] - y2[0])/y2[0])))\n")
    f.write("ax1.set_title(\"embryo volume (raw)\", fontsize=15)\n")
    f.write("ax1.legend(prop={'size': 10})\n")

    f.write("\n")
    f.write("for n in volume_along_time:\n")
    f.write("    x = list(volume_along_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    t = [temporal_coefficients[n][0] * i + temporal_coefficients[n][1] for i in x]\n")
    f.write("    y = [volume_along_time[n][i] for i in x]\n")
    f.write("    ax2.plot(t, y)\n")
    f.write("ax2.set_title(\"embryo volume (+ temporal alignment)\", fontsize=15)\n")
    f.write("ax2.legend(labels, prop={'size': 10})\n")

    f.write("\n")
    f.write("for n in volume_along_time:\n")
    f.write("    x = list(volume_along_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    t = [temporal_coefficients[n][0] * i + temporal_coefficients[n][1] for i in x]\n")
    f.write("    y = [volume_along_time[n][i] * voxelsize[n][i] * voxelsize[n][i] * voxelsize[n][i] for i in x]\n")
    f.write("    ax3.plot(t, y)\n")
    f.write("ax3.set_title(\"embryo volume (+ constantness)\", fontsize=15)\n")
    f.write("ax3.legend(labels, prop={'size': 10})\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('temporal_alignment_embryo_volume")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    f.write("fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(24, 8))\n")
    f.write("labels = []\n")
    f.write("for n in surface_along_time:\n")
    f.write("    labels += [n]\n")
    f.write("    x = list(surface_along_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    y = [surface_along_time[n][i] for i in x]\n")
    f.write("    p = ax1.plot(x, y, label=n)\n")
    f.write("\n")
    f.write("ax1.set_title(\"embryo surface (raw)\", fontsize=15)\n")
    f.write("ax1.legend(prop={'size': 10})\n")

    f.write("\n")
    f.write("for n in surface_along_time:\n")
    f.write("    x = list(surface_along_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    t = [temporal_coefficients[n][0] * i + temporal_coefficients[n][1] for i in x]\n")
    f.write("    y = [surface_along_time[n][i] for i in x]\n")
    f.write("    ax2.plot(t, y)\n")
    f.write("ax2.set_title(\"embryo surface (+ temporal alignment)\", fontsize=15)\n")
    f.write("ax2.legend(labels, prop={'size': 10})\n")

    f.write("\n")
    f.write("for n in surface_along_time:\n")
    f.write("    x = list(surface_along_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    t = [temporal_coefficients[n][0] * i + temporal_coefficients[n][1] for i in x]\n")
    f.write("    y = [surface_along_time[n][i] * voxelsize[n][i] * voxelsize[n][i] for i in x]\n")
    f.write("    ax3.plot(t, y)\n")
    f.write("ax3.set_title(\"embryo surface (+ constantness)\", fontsize=15)\n")
    f.write("ax3.legend(labels, prop={'size': 10})\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('temporal_alignment_embryo_surface")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")
    f.close()


################################################################################
#
# cell neighbors number wrt total number of cells in the embryo
#
################################################################################

def _neighbor_histogram(neighbors, atlas, atlasname, threshold=0.05):
    contact = atlas.cell_contact_surface
    cellname = atlas.cell_name

    cell_number_threshold = 10
    #
    # nodespertime is a dictionary
    # nodespertime[t] = #nodes at time t
    #
    nodes = list(contact.keys())
    div = 10 ** atlas.time_digits_for_cell_id
    times = [n // div for n in nodes]
    nodespertime = Counter(times)

    for n in contact:
        surface = 0
        for k in contact[n]:
            surface += contact[n][k]
        frac = []
        for k in contact[n]:
            if k % div == 0 or k % div == 1:
                continue
            frac += [contact[n][k]/surface]
        #
        # frac: array of surface fraction (without the background)
        #
        ncells = len([x for x in frac if x >= threshold])
        #
        # print a message for large neighbor number
        #
        if ncells > cell_number_threshold:
            msg = "cell " + str(n)
            if n in cellname:
                msg += " (" + str(cellname[n]) + ")"
            msg += " of properties '" + str(atlasname) + "' has " + str(ncells) + " neighbors"
            msg += " (above " + str(threshold*100) + "%)"
            monitoring.to_log_and_console("\t " + msg)
        t = n // div
        neighbors[nodespertime[t]] = neighbors.get(nodespertime[t], []) + [ncells]
    return neighbors


def _figures_neighbor_histogram_period(f):

    f.write("\n")
    f.write("axtwin = ax.secondary_xaxis('top')\n")
    f.write("axtwin.tick_params('x', direction='inout', length=10, width=2)\n")
    f.write("axtwin.set_xticks([76, 110, 250, 332, 624])\n")

    f.write("\n")
    f.write("x = [2, 76, 76, 2]\n")
    f.write("y = [0, 0, 14, 14]\n")
    f.write("plt.fill(x, y, color='grey', alpha=0.4)\n")
    f.write("plt.text(15, 13.5, 'cleavage', ha='left', va='center', fontsize=12, wrap=True)\n")

    f.write("\n")
    f.write("x = [110, 250, 250, 110]\n")
    f.write("y = [0, 0, 14, 14]\n")
    f.write("plt.fill(x, y, color='grey', alpha=0.4)\n")
    f.write("plt.text(180, 13.5, 'gastrula', ha='center', va='center', fontsize=12, wrap=True)\n")

    f.write("\n")
    f.write("x = [332, 624, 624, 332]\n")
    f.write("y = [0, 0, 14, 14]\n")
    f.write("plt.fill(x, y, color='grey', alpha=0.4)\n")
    f.write("plt.text(480, 13.5, 'neurula', ha='center', va='center', fontsize=12, wrap=True)\n")


def _figures_neighbor_histogram(atlases, parameters):
    """
    Plot cell neighbor number with respect to total cell number
    Parameters
    ----------
    atlases
    parameters

    Returns
    -------

    """
    proc = "_figures_neighbor_histogram"

    ref_atlases = atlases.get_atlases()
    neighbors = {}
    for a in ref_atlases:
        neighbors = _neighbor_histogram(neighbors, ref_atlases[a], a)

    #
    # neighbors is a dictionary indexed by the number of cells of an embryo
    # neighbors[64] is an array containing the number of neighbors
    #
    filename = 'figures_neighbor_histogram'
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
    f.write("import scipy.stats as stats\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    f.write("neighbors = " + str(neighbors) + "\n")
    f.write("ncells = sorted(list(neighbors.keys()))\n")

    f.write("\n")
    f.write("neighbors_per_ncell = []\n")
    f.write("ticks = []\n")
    f.write("xp = []\n")
    f.write("mp = []\n")
    f.write("sp = []\n")
    f.write("for i in ncells:\n")
    f.write("    xp += [i]\n")
    f.write("    mp += [np.mean(neighbors[i])]\n")
    f.write("    sp += [np.std(neighbors[i])]\n")
    f.write("    neighbors_per_ncell.append(neighbors[i])\n")
    f.write("    if i % 20 == 0:\n")
    f.write("        ticks += [str(i)]\n")
    f.write("    else:\n")
    f.write("        ticks += ['']\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(16, 6.5))\n")

    f.write("\n")
    f.write("ax.set_xlabel('number of cells', fontsize=15)\n")
    f.write("ax.set_ylabel('number of neighbors', fontsize=15)\n")
    f.write("ax.boxplot(neighbors_per_ncell, positions=ncells)\n")
    f.write("ax.set_title(\"cell neighbors\", fontsize=15)\n")
    f.write("ax.set_xticklabels(ticks)\n")
    f.write("ax.set_xlim(left=0)\n")

    _figures_neighbor_histogram_period(f)

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('neighbor_histogram_boxplot")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(16, 6.5))\n")

    f.write("\n")
    f.write("ax.set_xlabel('number of cells', fontsize=15)\n")
    f.write("ax.set_ylabel('average number of neighbors (+/- std dev)', fontsize=15)\n")
    f.write("ax.errorbar(xp, mp, yerr=sp, color='red', fmt='-', ecolor='blue')\n")
    f.write("ax.set_title(\"cell neighbors\", fontsize=15)\n")
    f.write("ax.set_xlim(left=0)\n")

    _figures_neighbor_histogram_period(f)

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('neighbor_histogram_average")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")


########################################################################################
#
# direction distribution
#
########################################################################################

def _figure_distribution(d, colormap='plasma', figname=None):

    i_indices = [p['voxel'][0] for p in d]
    j_indices = [p['voxel'][1] for p in d]
    k_indices = [p['voxel'][2] for p in d]
    tab = np.zeros((max(i_indices)+1, max(j_indices)+1, max(k_indices)+1), dtype=bool)
    values = np.zeros(tab.shape)
    for p in d:
        tab[p['voxel'][0]][p['voxel'][1]][p['voxel'][2]] = True
        values[p['voxel'][0]][p['voxel'][1]][p['voxel'][2]] = p['value']
    # normalize values to 1 for visualization
    maxvalue = max([p['value'] for p in d])
    print("_distribution_plot maxvalue = " + str(maxvalue))
    cmap = plt.get_cmap(colormap)
    colors = cmap(values/maxvalue)

    ax = plt.figure().add_subplot(projection='3d')
    # ax.voxels(r, g, b, sphere,
    #          facecolors=colors,
    #          edgecolors=np.clip(2*colors - 0.5, 0, 1),  # brighter
    #          linewidth=0.5)
    ax.voxels(tab, facecolors=colors, edgecolors=np.clip(2*colors - 0.5, 0, 1), linewidth=0.5)
    ax.set(xlabel='x', ylabel='y', zlabel='z')
    if figname is not None:
        plt.savefig(figname)
    else:
        plt.show()

########################################################################################
#
# symmetry axis assessment
#
########################################################################################


def _print_candidates(candidates):
    for i, p in enumerate(candidates):
        msg = "#{:2d} [{:.2f} {:.2f} {:.2f}]".format(i, p['vector'][0], p['vector'][1], p['vector'][2])
        msg += " error: {:.2f}".format(p['error'])
        print("   " + msg)
        msg = "value: {:.3f} single_pairing: {:.3f}".format(p['value'], p['score_from_single_pairing'])
        msg += " multiple_pairing: {:.3f}".format(p['score_from_multiple_pairing'])
        msg += " symmetry_icp: {:.3f}".format(p['score_from_symmetry_icp'])
        print("      " + msg)


def _symmetry_axis_error_wrt_times(atlas, parameters):

    icp.monitoring.copy(monitoring)

    cells = list(atlas.cell_contact_surface.keys())
    cells = sorted(cells)
    cells_per_time = {}
    div = 10 ** atlas.time_digits_for_cell_id
    for c in cells:
        t = int(c) // div
        cells_per_time[t] = cells_per_time.get(t, 0) + 1

    unnamedcells = [c for c in atlas.cell_volume if c not in atlas.cell_name]
    unnamedtimes = set([int(c) // div for c in unnamedcells])
    if len(unnamedtimes) > 0:
        # print("   unnamed cells = " + str(unnamedcells))
        monitoring.to_log_and_console("    unnamed times = " + str(sorted(unnamedtimes)))
    alltimes = sorted(list(cells_per_time.keys()))

    times = []

    err_distribution = []
    err_single_pairing = []
    err_multiple_pairing = []
    err_symmetry_icp = []

    index_from_distribution_value = []
    index_from_single_pairing = []
    index_from_multiple_pairing = []
    index_from_symmetry_icp = []

    for t in alltimes:
        if t in unnamedtimes:
            continue
        times += [t]

        symdir = atlas.get_symmetry_axis_from_names(t)

        #
        # compute distribution, and evaluate with 'single_pairing'
        #
        parameters.maxima_evaluation = 'single_pairing'
        candidates = atlas.get_direction_distribution_candidates(t, parameters)

        monitoring.to_log_and_console("   ... processing time " + str(t) + "/" + str(alltimes[-1]) + " - #candidates = " + str(len(candidates)))
        # print("symdir = " + str(symdir))
        #
        # evaluate with 'multiple_pairing'
        #
        parameters.maxima_evaluation = 'multiple_pairing'
        candidates = atlas.get_direction_distribution_candidates(t, parameters)

        #
        # evaluate with 'symmetric_registration'
        #
        parameters.maxima_evaluation = 'symmetric_registration'
        candidates = atlas.get_direction_distribution_candidates(t, parameters)

        for p in candidates:
            err = math.acos(np.dot(p['vector'], symdir))
            if err > np.pi / 2:
                err = np.pi - err
            p['error'] = err * 180.0 / np.pi

            err = math.acos(np.dot(p['average_from_single_pairing'], symdir))
            if err > np.pi / 2:
                err = np.pi - err
            p['error_from_single_pairing'] = err * 180.0 / np.pi

            err = math.acos(np.dot(p['average_from_multiple_pairing'], symdir))
            if err > np.pi / 2:
                err = np.pi - err
            p['error_from_multiple_pairing'] = err * 180.0 / np.pi

            err = math.acos(np.dot(p['average_from_symmetry_icp'], symdir))
            if err > np.pi / 2:
                err = np.pi - err
            p['error_from_symmetry_icp'] = err * 180.0 / np.pi

        minimal_error = min([p['error'] for p in candidates])
        index = [p['error'] for p in candidates].index(minimal_error)
        # print("ERROR")
        # _print_candidates(candidates)
        # print("   minimal_error = " + str(minimal_error))
        # print("   index = " + str(index))

        err_distribution += [minimal_error]
        err_single_pairing += [candidates[index]['error_from_single_pairing']]
        err_multiple_pairing += [candidates[index]['error_from_multiple_pairing']]
        err_symmetry_icp += [candidates[index]['error_from_symmetry_icp']]

        candidates = sorted(candidates, reverse=True, key=lambda x: x['value'])
        index_from_distribution_value += [[p['error'] for p in candidates].index(minimal_error) + 1]
        # print("VALUE SORTED")
        # _print_candidates(candidates)
        # print("   index_from_distribution_value = " + str(index_from_distribution_value))

        candidates = sorted(candidates, reverse=True, key=lambda x: x['score_from_single_pairing'])
        index_from_single_pairing += [[p['error'] for p in candidates].index(minimal_error) + 1]
        # print("SINGLE PAIRING SORTED")
        # _print_candidates(candidates)
        # print("   index_from_single_pairing = " + str(index_from_single_pairing))

        candidates = sorted(candidates, reverse=True, key=lambda x: x['score_from_multiple_pairing'])
        index_from_multiple_pairing += [[p['error'] for p in candidates].index(minimal_error) + 1]
        # print("MULTIPLE PAIRING SORTED")
        # _print_candidates(candidates)
        # print("   index_from_multiple_pairing = " + str(index_from_multiple_pairing))

        candidates = sorted(candidates, reverse=True, key=lambda x: x['score_from_symmetry_icp'])
        index_from_symmetry_icp += [[p['error'] for p in candidates].index(minimal_error) + 1]
        # print("SYMMETRIC ICP SORTED")
        # _print_candidates(candidates)
        # print("   index_from_symmetry_icp = " + str(index_from_symmetry_icp))

        atlas.del_property('direction_distribution')

    ncells = [cells_per_time[t] for t in times]

    return times, ncells, err_distribution, err_single_pairing, err_multiple_pairing, err_symmetry_icp, \
           index_from_distribution_value, index_from_single_pairing, index_from_multiple_pairing, \
           index_from_symmetry_icp


def _figure_symmetry_axis(atlases, parameters):
    proc = "_figure_symmetry_axis"
    filename = 'figure_symmetry_axis'
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

    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    ref_atlases = atlases.get_atlases()
    for a in ref_atlases:
        monitoring.to_log_and_console("... processing " + str(a))
        na = a.replace('-', '_')
        output = _symmetry_axis_error_wrt_times(ref_atlases[a], parameters)
        times, ncells, edist, esing, emult, esymicp, idist, ising, imult, isymicp = output
        f.write("temporal_coefficients = " + str(ref_atlases[a].temporal_alignment) + "\n")
        f.write("times_" + na + " = " + str(times) + "\n")
        f.write("atimes_" + na + " = [temporal_coefficients[0] * i + temporal_coefficients[1]")
        f.write(" for i in times_" + na + "]\n")
        f.write("ncells_" + na + " = " + str(ncells) + "\n")
        f.write("edist_" + na + " = " + str(edist) + "\n")
        f.write("esing_" + na + " = " + str(esing) + "\n")
        f.write("emult_" + na + " = " + str(emult) + "\n")
        f.write("esymicp_" + na + " = " + str(esymicp) + "\n")
        f.write("idist_" + na + " = " + str(idist) + "\n")
        f.write("ising_" + na + " = " + str(ising) + "\n")
        f.write("imult_" + na + " = " + str(imult) + "\n")
        f.write("isymicp_" + na + " = " + str(isymicp) + "\n")
        f.write("\n")

    f.write("\n")
    f.write("fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(24, 6))\n")
    for a in ref_atlases:
        na = a.replace('-', '_')
        f.write("p = ax1.plot(ncells_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("p = ax2.plot(edist_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("p = ax3.plot(idist_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("\n")

    f.write("ax1.grid(True)\n")
    f.write("ax1.legend(prop={'size': 10})\n")
    f.write("ax1.set_xlabel('cell number')\n")
    f.write("ax1.set_ylabel('time')\n")
    f.write("ax1.set_title(\"cell number\", fontsize=15)\n")
    f.write("\n")
    f.write("ax2.grid(True)\n")
    f.write("ax2.set_xlabel('error (degrees)')\n")
    f.write("ax2.set_title(\"symmetry axis error of the closest vector\", fontsize=15)\n")
    f.write("\n")
    f.write("ax3.grid(True)\n")
    f.write("ax3.set_xlabel('rank')\n")
    f.write("ax3.set_title(\"rank of the closest vector\", fontsize=15)\n")
    f.write("\n")
    f.write("fig.suptitle('Symmetry axis candidate sorted by distribution value', fontsize=16)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('symmetry_axis_distribution")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    f.write("fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(24, 6))\n")
    for a in ref_atlases:
        na = a.replace('-', '_')
        f.write("p = ax1.plot(ncells_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("p = ax2.plot(esing_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("p = ax3.plot(ising_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("\n")

    f.write("ax1.grid(True)\n")
    f.write("ax1.legend(prop={'size': 10})\n")
    f.write("ax1.set_xlabel('cell number')\n")
    f.write("ax1.set_ylabel('time')\n")
    f.write("ax1.set_title(\"cell number\", fontsize=15)\n")
    f.write("\n")
    f.write("ax2.grid(True)\n")
    f.write("ax2.set_xlabel('error (degrees)')\n")
    f.write("ax2.set_title(\"symmetry axis error of estimated vector\", fontsize=15)\n")
    f.write("\n")
    f.write("ax3.grid(True)\n")
    f.write("ax3.set_xlabel('rank')\n")
    f.write("ax3.set_title(\"rank of the closest vector\", fontsize=15)\n")
    f.write("\n")
    f.write("fig.suptitle('Symmetry axis candidate sorted by single pairing score', fontsize=16)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('symmetry_axis_singlepairing")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    f.write("fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(24, 6))\n")
    for a in ref_atlases:
        na = a.replace('-', '_')
        f.write("p = ax1.plot(ncells_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("p = ax2.plot(emult_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("p = ax3.plot(imult_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("\n")

    f.write("ax1.grid(True)\n")
    f.write("ax1.legend(prop={'size': 10})\n")
    f.write("ax1.set_xlabel('cell number')\n")
    f.write("ax1.set_ylabel('time')\n")
    f.write("ax1.set_title(\"cell number\", fontsize=15)\n")
    f.write("\n")
    f.write("ax2.grid(True)\n")
    f.write("ax2.set_xlabel('error (degrees)')\n")
    f.write("ax2.set_title(\"symmetry axis error of estimated vector\", fontsize=15)\n")
    f.write("\n")
    f.write("ax3.grid(True)\n")
    f.write("ax3.set_xlabel('rank')\n")
    f.write("ax3.set_title(\"rank of the closest vector\", fontsize=15)\n")
    f.write("\n")
    f.write("fig.suptitle('Symmetry axis candidate sorted by multiple pairing score', fontsize=16)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('symmetry_axis_multiplepairing")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    f.write("fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(24, 6))\n")
    for a in ref_atlases:
        na = a.replace('-', '_')
        f.write("p = ax1.plot(ncells_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("p = ax2.plot(esymicp_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("p = ax3.plot(isymicp_" + na + ", atimes_" + na + ", label='" + na + "')\n")
        f.write("\n")

    f.write("ax1.grid(True)\n")
    f.write("ax1.legend(prop={'size': 10})\n")
    f.write("ax1.set_xlabel('cell number')\n")
    f.write("ax1.set_ylabel('time')\n")
    f.write("ax1.set_title(\"cell number\", fontsize=15)\n")
    f.write("\n")
    f.write("ax2.grid(True)\n")
    f.write("ax2.set_xlabel('error (degrees)')\n")
    f.write("ax2.set_title(\"symmetry axis error of estimated vector\", fontsize=15)\n")
    f.write("\n")
    f.write("ax3.grid(True)\n")
    f.write("ax3.set_xlabel('rank')\n")
    f.write("ax3.set_title(\"rank of the closest vector\", fontsize=15)\n")
    f.write("\n")
    f.write("fig.suptitle('Symmetry axis candidate sorted by symmetry registration', fontsize=16)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('symmetry_axis_symmetryregistration")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.close()
