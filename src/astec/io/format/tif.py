import tifffile
import numpy as np
import os, os.path, sys, time

from astec.components.spatial_image import SpatialImage

__all__ = []
__all__ += ["read_tif", "write_tif"]


def read_tif(filename):
    tif = tifffile.TiffFile(filename)
    _data = tif.asarray()
    _data = _data.T

    xtag = tif.pages[0].tags["XResolution"]
    ytag = tif.pages[0].tags["YResolution"]
    imagej_metadata = tif.imagej_metadata

    _vx = xtag.value[1] / xtag.value[0]
    _vy = ytag.value[1] / ytag.value[0]
    if imagej_metadata is not None and "spacing" in imagej_metadata:
        _vz = imagej_metadata["spacing"]
    else:
        # TODO make this more universally applicable
        _vx = 1
        _vy = 1
        _vz = 1

    tif.close()
    # -- dtypes are not really stored in a compatible way (">u2" instead of uint16)
    # but we can convert those --
    dt = np.dtype(_data.dtype.name)
    # -- Return a SpatialImage please! --
    im = SpatialImage(_data, dtype=dt)
    im.voxelsize = _vx, _vy, _vz

    return im


def write_tif(filename, obj):
    proc = "write_tif"
    if len(obj.shape) > 3:
        raise IOError(
            proc + ": vectorial images are currently unsupported by tif writer"
        )

    vsx, vsy, vsz = obj.voxelsize
    obj = obj.T

    tifffile.imwrite(
        filename,
        obj,
        imagej=True,
        resolution=(1.0 / vsx, 1.0 / vsy),
        metadata={"spacing": vsz, "unit": "um", "axes": "ZYX"},
    )
    return
