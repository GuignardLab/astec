import copy
import math
import numpy as np
from sklearn.neighbors import NearestNeighbors

import astec.utils.common as common
monitoring = common.Monitoring()


########################################################################################
#
# transformations
#
########################################################################################

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

def icp(ref=None, flo=None, transformation_type="affine", estimation="lts", retained_fraction=0.75, verbose=False):
    """
    Iterative Closest Point procedure
    Parameters
    ----------
    ref : numpy.ndarray
        A numpy ndarray of 3D points, with ref.shape[0] = 3 and ref.shape[1] the number of
        reference points
    flo : numpy.ndarray
        A numpy ndarray of 3D points, with flo.shape[0] = 3 and flo.shape[1] the number of
        floating points
    transformation_type : str
        'rigid'
        'similitude'
        'affine'
    estimation : str
        'ls'
        'lts'
    retained_fraction : float
    verbose : bool
    Returns
    -------
    A 4x4 numpy.ndarray that allows to transform the floating points (in homogeneous
    coordinates) to superimpose them on reference points.

    """
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
            imat = _lts_transformation(sref, sflo, transformation_type=transformation_type,
                                       retained_fraction=retained_fraction, verbose=False)
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
