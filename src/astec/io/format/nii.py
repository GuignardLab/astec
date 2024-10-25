import numpy as np
import nibabel as nib

from astec.components.spatial_image import SpatialImage

__all__ = ["read_nii", "write_nii"]


def read_nii(filename):
    nii = nib.load(filename)

    #
    # dimensions
    #
    xdim = nii.header["dim"][1]
    ydim = nii.header["dim"][2]
    zdim = nii.header["dim"][3]
    if nii.header["dim"][4] > 1:
        vdim = nii.header["dim"][4]
    elif nii.header["dim"][5] > 1:
        vdim = nii.header["dim"][5]
    else:
        vdim = 1

    #
    # voxel size
    #
    res = [nii.header["pixdim"][1], nii.header["pixdim"][2], nii.header["pixdim"][3]]

    # read datas
    mat = np.array(nii.dataobj)
    if vdim != 1:
        mat = mat.reshape((vdim, xdim, ydim, zdim), order="F")
        mat = mat.transpose(1, 2, 3, 0)
    else:
        mat = mat.reshape((xdim, ydim, zdim), order="F")

    # create SpatialImage

    img = SpatialImage(mat, res, vdim)

    return img


def write_nii(filename, img):
    header = nib.Nifti1Header()
    # 2 millimeters
    # 3 micrometers
    # seconds
    header.set_xyzt_units(xyz=3, t=8)
    header.set_data_shape(img.shape)
    header.set_zooms(getattr(img, "voxelsize", (1, 1, 1)))
    header.set_data_dtype(img.dtype)
    img = nib.nifti1.Nifti1Image(img, None, header=header)
    nib.save(img, filename)
