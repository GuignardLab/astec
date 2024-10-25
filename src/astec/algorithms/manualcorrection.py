import os
import sys
import time
import numpy as np
import shutil
from scipy import ndimage as nd
import copy

import astec.algorithms.mars as mars
import astec.algorithms.astec as astec
from astec.utils import common
from astec.components.spatial_image import SpatialImage
from astec.io.image import imread, imsave
import astec.utils.diagnosis as udiagnosis
import astec.utils.ioproperties as ioproperties
import astec.utils.reconstruction as reconstruction
import astec.wrapping.cpp_wrapping as cpp_wrapping

#
#
#
#
#

monitoring = common.Monitoring()


########################################################################################
#
# classes
# - computation parameters
#
########################################################################################


class ManualCorrectionParameters(udiagnosis.DiagnosisParameters, astec.AstecParameters):
    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix="mancor_"):
        if "doc" not in self.__dict__:
            self.doc = {}

        udiagnosis.DiagnosisParameters.__init__(self, prefix=prefix)
        astec.AstecParameters.__init__(self, prefix=prefix)

        doc = "\n"
        doc += "Manual correction overview:\n"
        doc += "===========================\n"
        doc += "1. Fuses cells/labels from the segmentation to correct over-segmentation errors\n"
        doc += "2. Splits cells/labels from the segmentation to correct under-segmentation errors.\n"
        doc += "   Division is done by exploring a rangle of values of 'hmin' until two seeds are\n"
        doc += "     found inside the cell to be divided.\n"
        doc += "   Division is propagated for segmented time series as in first step of 'astec'\n"
        doc += "   according a lineage is present in the input segmentation directory\n"
        doc += "     a. segmentation of previous time point is deformed onto the current time point.\n"
        doc += (
            "     b. cells are eroded to get two seeds inside the cell to be divided\n"
        )
        doc += "     c. segmentation (watershed-based) from the deformed seeds\n"
        doc += " Transform segmentation images from the SEG/SEG_'EXP_SEG_FROM' directory to "
        doc += " segmentation images in the SEG/SEG_'EXP_SEG_TO' directory"
        self.doc["manualcorrection_overview"] = doc

        doc = "\t Directory containing the manual correction files\n"
        self.doc["manualcorrection_dir"] = doc
        self.manualcorrection_dir = None

        doc = "\t File containing the labels to be fused.\n"
        doc += "\t Each line may either indicate a comment, a division or a fusion to be made\n"
        doc += "\t - lines beginning by '#' are ignored (comment)\n"
        doc += (
            "\t - lines with only numbers concern changes for the first time point:\n"
        )
        doc += "\t   - one single number: label of the cell to be splitted at the first time point\n"
        doc += "\t   - several numbers: labels of the cells to be fused\n"
        doc += "\t - lines beginning by 'timevalue:' concern changes for the given time point\n"
        doc += "\t   - 'timevalue:' + one single number: label of the cell to be splitted\n"
        doc += (
            "\t   - 'timevalue:' + several numbers: labels of the cells to be fused\n"
        )
        doc += "\t - lines beginning by 'timevalue-timevalue:' concern changes for the given time point range\n"
        doc += "\t   - 'timevalue-timevalue:' + several numbers: labels of the cells to be fused\n"
        doc += "\t      if a lineage is present, it is used to track the cell labels to be fused\n"
        doc += (
            "\t      else the same labels are fused at each time point of the range\n"
        )
        self.doc["manualcorrection_file"] = doc
        self.manualcorrection_file = None

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print("#")
        print("# ManualCorrectionParameters ")
        print("#")
        print("")

        common.PrefixedParameter.print_parameters(self)

        udiagnosis.DiagnosisParameters.print_parameters(self)
        astec.AstecParameters.print_parameters(self)

        for line in self.doc["manualcorrection_overview"].splitlines():
            print("# " + line)

        self.varprint(
            "manualcorrection_dir",
            self.manualcorrection_dir,
            self.doc["manualcorrection_dir"],
        )
        self.varprint(
            "manualcorrection_file",
            self.manualcorrection_file,
            self.doc["manualcorrection_file"],
        )

        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("#" + "\n")
        logfile.write("# ManualCorrectionParameters " + "\n")
        logfile.write("#" + "\n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        udiagnosis.DiagnosisParameters.write_parameters_in_file(self, logfile)
        astec.AstecParameters.write_parameters_in_file(self, logfile)

        for line in self.doc["manualcorrection_overview"].splitlines():
            logfile.write("# " + line + "\n")

        self.varwrite(
            logfile,
            "manualcorrection_dir",
            self.manualcorrection_dir,
            self.doc.get("manualcorrection_dir", None),
        )
        self.varwrite(
            logfile,
            "manualcorrection_file",
            self.manualcorrection_file,
            self.doc.get("manualcorrection_file", None),
        )

        logfile.write("\n")
        return

    def write_parameters(self, log_file_name):
        with open(log_file_name, "a") as logfile:
            self.write_parameters_in_file(logfile)
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameters):
        udiagnosis.DiagnosisParameters.update_from_parameters(self, parameters)
        astec.AstecParameters.update_from_parameters(self, parameters)

        self.manualcorrection_dir = self.read_parameter(
            parameters, "manualcorrection_dir", self.manualcorrection_dir
        )
        self.manualcorrection_file = self.read_parameter(
            parameters, "manualcorrection_file", self.manualcorrection_file
        )
        self.manualcorrection_file = self.read_parameter(
            parameters, "mapping_file", self.manualcorrection_file
        )

    def update_from_parameter_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            monitoring.to_log_and_console(
                "Error: '" + parameter_file + "' is not a valid file. Exiting."
            )
            sys.exit(1)

        parameters = common.load_source(parameter_file)
        self.update_from_parameters(parameters)


########################################################################################
#
# some internal procedures
#
########################################################################################


def _read_correction_file(
    filename, lineage=None, first_time_point=1, time_digits_for_cell_id=4
):
    """

    Parameters
    ----------
    filename: filename of the corrections.
    For fusion, the admitted syntax are
    - for fusion at one given time point (if 't:' is omitted, it is applied at the first time point)
        [t:] label label [label ...]
    - for fusion at an interval of time points:
        t-t: label label [label ...]
      If no lineage is given, the same labels are used for the whole interval. If there is a lineage,
      labels at next time points are issued from the lineage.
    For division, the admitted syntax is
        [t:] label
    lineage:
    time_digits_for_cell_id

    Returns
    -------
        2 dictionaries (fusion, division) indexed either by 'first_time_point' or by an integer designing a
        time point.
        - fusion[t] is also a dictionary, whose keys are cell labels. fusion[t][c] is an array of equivalent labels
          for label 'c'
        - division[t] is a list of cell labels

    """
    proc = "_read_correction_file"
    fusion = {}
    division = {}
    separator_time_interval = ["-", ":"]
    separator_time_label = [":", ","]

    if not os.path.isfile(filename):
        monitoring.to_log_and_console(
            "Error: '" + filename + "' is not a valid file. Exiting."
        )
        sys.exit(1)
    f = open(filename)
    i = 0
    for line in f:
        i += 1
        # remove leading and trailing whitespaces
        li = line.strip()
        # skip comment
        if li.startswith("#"):
            continue
        # empty line
        if len(li) == 0:
            continue
        li = li.split()

        # one component: single integer, division at first time point
        if len(li) == 1:
            if li[0].isdigit():
                division[first_time_point] = division.get(first_time_point, []) + [
                    int(li[0])
                ]
                continue
            msg = (
                "line #"
                + str(i)
                + ": '"
                + str(line)
                + "' should be a single integer?! Skip the line."
            )
            monitoring.to_log_and_console(proc + ": " + msg)
            continue

        # two or more components: only integers or "t:" + integers
        # last component(s) should be only integer
        error = False
        for i in range(1, len(li)):
            if not li[i].isdigit():
                msg = (
                    str(i + 1)
                    + "th term of line #"
                    + str(i)
                    + ": '"
                    + str(line)
                    + "' should be numbers."
                )
                msg += " Skip the line."
                monitoring.to_log_and_console(proc + ": " + msg)
                error = True
                break
        if error:
            continue

        # only cell labels: fusion
        if li[0].isdigit():
            labels = [int(lab) for lab in li]
            if min(labels) == max(labels):
                msg = (
                    "line #"
                    + str(i)
                    + ": '"
                    + str(line)
                    + "' has the same label '"
                    + str(li[0])
                )
                msg += "for fusion?! ... Skip the line."
                monitoring.to_log_and_console(proc + ": " + msg)
                continue
            if first_time_point not in fusion:
                fusion[first_time_point] = {}
            # use minimal label as destination label
            for lab in labels:
                if lab > min(labels):
                    fusion[first_time_point][lab] = fusion[first_time_point].get(
                        lab, []
                    ) + [min(labels)]
            continue

        # first term is not a number, try to recognize time (interval)
        # "t:" + integer -> division
        # "t:" + integers -> fusion
        # "t-t:" + integers -> fusion
        if li[0][-1] not in separator_time_label:
            msg = (
                "last character of first term of line #"
                + str(i)
                + ": '"
                + str(line)
                + "' should be in "
            )
            msg += str(separator_time_label) + "."
            msg += " Skip the line."
            monitoring.to_log_and_console(proc + ": " + msg)
            continue

        # "t:" + integer(s)
        if li[0][:-1].isdigit():
            ctime = int(li[0][:-1])
            # "t:" + integer -> division
            if len(li) == 2:
                division[ctime] = division.get(ctime, []) + [int(li[1])]
                continue
            # "t:" + integers -> fusion
            labels = [int(li[i]) for i in range(1, len(li))]
            if ctime not in fusion:
                fusion[ctime] = {}
            # use minimal label as destination label
            for lab in labels:
                if lab > min(labels):
                    fusion[ctime][lab] = fusion[ctime].get(lab, []) + [min(labels)]
            continue

        # "t-t:" + integers -> fusion
        # there is one time interval separator
        count_separator_time_interval = 0
        for c in li[0][:-1]:
            if c in separator_time_interval:
                count_separator_time_interval += 1
        if count_separator_time_interval != 1:
            msg = "weird first term of line #" + str(i) + ": '" + str(line) + "'."
            msg += "No or too many time interval separator. "
            msg += " Skip the line."
            monitoring.to_log_and_console(proc + ": " + msg)
            continue

        # the time separator cuts the line into two parts (should be true if there is only one separator)
        ti = []
        for s in separator_time_interval:
            if s in li[0][:-1]:
                ti = li[0][:-1].split(s)
        if len(ti) != 2:
            msg = "weird first term of line #" + str(i) + ": '" + str(line) + "'."
            msg += "Should be cut into two pieces."
            msg += " Skip the line."
            monitoring.to_log_and_console(proc + ": " + msg)
            continue
        for t in ti:
            if not t.isdigit():
                msg = "weird first term of line #" + str(i) + ": '" + str(line) + "'."
                msg += "Should be cut into two integers."
                msg += " Skip the line."
                monitoring.to_log_and_console(proc + ": " + msg)
                continue

        # get the time interval
        mintime = min([int(t) for t in ti])
        maxtime = max([int(t) for t in ti])

        # "t-t:" + integer -> division?
        if len(li) == 2:
            if mintime != maxtime:
                msg = "line #" + str(i) + ": '" + str(line) + "'."
                msg += "Can not do a division during a time interval."
                msg += " Skip the line."
                monitoring.to_log_and_console(proc + ": " + msg)
                continue
            division[mintime] = division.get(mintime, []) + [int(li[1])]
            continue

        # "t-t:" + integers -> fusion
        labels = [int(li[i]) for i in range(1, len(li))]
        # fuse labels for a time interval
        if lineage is None:
            for ctime in range(mintime, maxtime + 1):
                if ctime not in fusion:
                    fusion[ctime] = {}
                # use minimal label as destination label
                for lab in labels:
                    if lab > min(labels):
                        fusion[ctime][lab] = fusion[ctime].get(lab, []) + [min(labels)]
            continue

        # tracks labels in lineage
        for ctime in range(mintime, maxtime + 1):
            for lab in labels:
                if lab > min(labels):
                    fusion[ctime][lab] = fusion[ctime].get(lab, []) + [min(labels)]
            if ctime == maxtime:
                break
            prev_cell_ids = [10**time_digits_for_cell_id + lab for lab in labels]
            next_cell_ids = []
            for c in prev_cell_ids:
                if c not in lineage:
                    continue
                next_cell_ids += lineage[c]
            if len(next_cell_ids) == 0:
                msg = "line #" + str(i) + ": '" + str(line) + "'."
                msg += "There is no cells after time '" + str(ctime) + "'."
                monitoring.to_log_and_console(proc + ": " + msg)
                break
            labels = [c - (c // 10**time_digits_for_cell_id) for c in next_cell_ids]

    f.close()
    return fusion, division


########################################################################################
#
#
#
########################################################################################


def _diagnosis_volume_image(output_image, parameters):
    im = imread(output_image)
    voxelsize = im.voxelsize
    vol = voxelsize[0] * voxelsize[1] * voxelsize[2]

    cell_label = np.unique(im)
    cell_volume = nd.sum(np.ones_like(im), im, index=np.int16(cell_label))

    volumes = zip(cell_label, cell_volume)
    volumes = sorted(volumes, key=lambda x: x[1])

    monitoring.to_log_and_console("    Number of cells: " + str(len(cell_label)), 0)
    monitoring.to_log_and_console("    Maximal label: " + str(np.max(im)), 0)
    # monitoring.to_log_and_console('    Cell ids: ' + str(cell_label), 0)

    monitoring.to_log_and_console("    Sorted cell volumes: ", 0)
    monitoring.to_log_and_console("      Id :    voxels          (um^3)", 0)

    if parameters.items <= 0 or parameters.items >= len(volumes):
        for v in volumes:
            msg = "    {:>4d} : {:>9d} {:>15s}".format(
                v[0], int(v[1]), "({:.2f})".format(v[1] * vol)
            )
            monitoring.to_log_and_console(msg, 0)
    else:
        if int(parameters.items) > 0:
            for v in volumes[: parameters.items]:
                msg = "    {:>4d} : {:>9d} {:>15s}".format(
                    v[0], int(v[1]), "({:.2f})".format(v[1] * vol)
                )
                monitoring.to_log_and_console(msg, 0)
            monitoring.to_log_and_console("       ...", 0)
            for v in volumes[-parameters.items :]:
                msg = "    {:>4d} : {:>9d} {:>15s}".format(
                    v[0], int(v[1]), "({:.2f})".format(v[1] * vol)
                )
                monitoring.to_log_and_console(msg, 0)


########################################################################################
#
#
#
########################################################################################


def _image_cell_fusion(
    input_image,
    output_image,
    current_time,
    properties,
    fusion,
    time_digits_for_cell_id=4,
):
    """

    Parameters
    ----------
    input_image
    output_image
    current_time
    properties
    fusion
    time_digits_for_cell_id

    Returns
    -------

    """
    proc = "_image_cell_fusion"
    im = imread(input_image)
    voxelsize = im.get_voxelsize()
    datatype = im.dtype

    #
    #
    #
    immax = np.max(im)
    mapping = np.arange(immax + 1)
    for lab in fusion[current_time]:
        eqlabel = min(min(fusion[current_time][lab]), mapping[lab])
        labels = [lab] + fusion[current_time][lab] + [mapping[lab]]
        for i in range(immax + 1):
            if mapping[i] in labels:
                mapping[i] = eqlabel

    im = mapping[im]
    imsave(output_image, SpatialImage(im, voxelsize=voxelsize).astype(datatype))

    #
    #
    #
    if properties is None:
        return None
    property_keys = list(properties.keys())
    for k in property_keys:
        if k == "cell_volume":
            for i in range(immax + 1):
                if mapping[i] == i:
                    continue
                label = current_time * 10**time_digits_for_cell_id + i
                eqlabel = current_time * 10**time_digits_for_cell_id + i
                properties["cell_volume"][eqlabel] += properties["cell_volume"][label]
                del properties["cell_volume"][label]
        elif k == "cell_lineage":
            reverse_lineage = {
                v: k for k, values in properties["cell_lineage"].items() for v in values
            }
            for i in range(immax + 1):
                if mapping[i] == i:
                    continue
                label = current_time * 10**time_digits_for_cell_id + i
                eqlabel = current_time * 10**time_digits_for_cell_id + mapping[i]
                #
                # both labels are in reverse lineage, check whether they have the same parent cell
                #
                if label in reverse_lineage and eqlabel in reverse_lineage:
                    if reverse_lineage[label] != reverse_lineage[eqlabel]:
                        msg = (
                            "warning, cells "
                            + str(i)
                            + " and "
                            + str(mapping[i])
                            + " are fused at time "
                        )
                        msg += (
                            str(current_time)
                            + " but they have different parent cells.\n"
                        )
                        msg += "\t Do not change the lineage."
                        monitoring.to_log_and_console(proc + ": " + msg)
                    else:
                        properties["cell_lineage"][reverse_lineage[eqlabel]].remove(
                            label
                        )
                        if (
                            eqlabel in properties["cell_lineage"]
                            and label in properties["cell_lineage"]
                        ):
                            properties["cell_lineage"][eqlabel] += properties[
                                "cell_lineage"
                            ][label]
                            del properties["cell_lineage"][label]
                        elif (
                            eqlabel not in properties["cell_lineage"]
                            and label not in properties["cell_lineage"]
                        ):
                            pass
                        else:
                            msg = (
                                "weird, only one of "
                                + str(label)
                                + " and "
                                + str(eqlabel)
                                + " is in lineage."
                            )
                            msg += "\t Do not change it."
                            monitoring.to_log_and_console(proc + ": " + msg)
                #
                #
                #
                elif label not in reverse_lineage and eqlabel not in reverse_lineage:
                    if (
                        eqlabel in properties["cell_lineage"]
                        and label in properties["cell_lineage"]
                    ):
                        properties["cell_lineage"][eqlabel] += properties[
                            "cell_lineage"
                        ][label]
                        del properties["cell_lineage"][label]
                    elif (
                        eqlabel not in properties["cell_lineage"]
                        and label not in properties["cell_lineage"]
                    ):
                        pass
                    else:
                        msg = (
                            "weird, only one of "
                            + str(label)
                            + " and "
                            + str(eqlabel)
                            + " is in lineage."
                        )
                        msg += "\t Do not change it."
                        monitoring.to_log_and_console(proc + ": " + msg)
                    #
                    #
                    #
                    msg = (
                        "weird, only one of "
                        + str(label)
                        + " and "
                        + str(eqlabel)
                        + " is in reverse lineage."
                    )
                    msg += "\t Do not change it."
                    monitoring.to_log_and_console(proc + ": " + msg)
        else:
            del properties[k]

    return properties


########################################################################################
#
#
#
########################################################################################


def _get_reconstructed_image(previous_time, current_time, experiment, parameters):
    proc = "_get_reconstructed_image"

    input_dir = experiment.fusion_dir.get_directory(0)
    input_name = experiment.fusion_dir.get_image_name(current_time)
    input_image = common.find_file(
        input_dir,
        input_name,
        file_type="image",
        callfrom=proc,
        local_monitoring=monitoring,
    )

    #
    #
    #
    output_reconstructed_image = common.add_suffix(
        input_image,
        parameters.result_suffix,
        new_dirname=experiment.astec_dir.get_rec_directory(0),
        new_extension=experiment.default_image_suffix,
    )
    if (
        os.path.isfile(output_reconstructed_image)
        and monitoring.forceResultsToBeBuilt is False
    ):
        monitoring.to_log_and_console(
            "    .. reconstructed image is '"
            + str(output_reconstructed_image).split(os.path.sep)[-1]
            + "'",
            2,
        )
        return output_reconstructed_image

    #
    # try to recover the reconstructed image from mar/input dir
    #
    input_reconstructed_image = common.add_suffix(
        input_image,
        parameters.result_suffix,
        new_dirname=experiment.mars_dir.get_rec_directory(0),
        new_extension=experiment.default_image_suffix,
    )
    if (
        os.path.isfile(input_reconstructed_image)
        and monitoring.forceResultsToBeBuilt is False
    ):
        monitoring.to_log_and_console(
            "    .. reconstructed image is '"
            + str(input_reconstructed_image).split(os.path.sep)[-1]
            + "'",
            2,
        )
        shutil.copy2(input_reconstructed_image, output_reconstructed_image)
        return output_reconstructed_image

    #
    # build the reconstructed image
    #
    experiment.working_dir = experiment.astec_dir
    output_reconstructed_image = reconstruction.build_reconstructed_image(
        current_time, experiment, parameters=parameters, previous_time=previous_time
    )
    if output_reconstructed_image is None or not os.path.isfile(
        output_reconstructed_image
    ):
        monitoring.to_log_and_console(
            "    .. "
            + proc
            + ": no reconstructed image was found/built for time "
            + str(current_time)
        )
        return False
    return output_reconstructed_image


def _slices_dilation_iteration(slices, maximum):
    return tuple(
        [
            slice(max(0, s.start - 1), min(s.stop + 1, maximum[i]))
            for i, s in enumerate(slices)
        ]
    )


def _slices_dilation(slices, maximum, iterations=1):
    for i in range(iterations):
        slices = _slices_dilation_iteration(slices, maximum)
    return slices


def _get_bounding_box(
    input_image,
    experiment,
    division,
    tmp_prefix_name,
    dilation_iterations=20,
    label_width=5,
):
    """
    Compute bounding boxes (including a margin) of cell of interest (cells that will divide)
    and extract subimages around these cells.
    Parameters
    ----------
    input_image
    experiment
    division
    tmp_prefix_name
    dilation_iterations

    Returns
    -------
    A dictionary of bounding boxes indexed by cell labels.
    """
    segmentation = imread(input_image)
    voxelsize = segmentation.get_voxelsize()
    datatype = segmentation.dtype

    bounding_boxes = {}
    #
    # bounding box are of the form [(slice(191, 266, None), slice(413, 471, None), slice(626, 692, None))]
    # ie [(xmin, xmax), (ymin, ymax), (zmin, zmax)]
    # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.find_objects.html
    # segmentation[bounding_boxes[c][0]] yields a smaller array defining by the bounding box
    #
    for c in division:
        bb = nd.find_objects(segmentation == c)
        bounding_boxes[c] = _slices_dilation(
            bb[0], maximum=segmentation.shape, iterations=dilation_iterations
        )
        cellid = str("{:0{width}d}".format(c, width=label_width))
        cell_prefix_name = (
            tmp_prefix_name
            + "_cell"
            + cellid
            + "_seg."
            + experiment.default_image_suffix
        )
        if not os.path.isfile(cell_prefix_name):
            ext_segmentation = segmentation[bounding_boxes[c]]
            imsave(
                cell_prefix_name,
                SpatialImage(ext_segmentation, voxelsize=voxelsize).astype(datatype),
            )
            del ext_segmentation
    del segmentation
    return bounding_boxes


def _two_seeds_extraction(
    first_segmentation, cellid, image_for_seed, experiment, parameters
):
    #
    # inspired by _cell_based_h_minima() in algorithms/astec.py
    #
    """

    Parameters
    ----------
    first_segmentation
    cellid
    image_for_seed
    parameters

    Returns
    -------

    """
    proc = "_two_seeds_extraction"

    #
    # h-minima extraction with h = max value
    # the difference image is kept for further computation
    #
    h_max = parameters.watershed_seed_hmin_max_value
    wparam = mars.WatershedParameters(obj=parameters)
    wparam.seed_hmin = h_max
    h_min = h_max

    input_image = image_for_seed
    unmasked_seed_image = common.add_suffix(
        image_for_seed,
        "_unmasked_seed_h" + str("{:03d}".format(h_min)),
        new_dirname=experiment.astec_dir.get_tmp_directory(),
        new_extension=experiment.default_image_suffix,
    )
    seed_image = common.add_suffix(
        image_for_seed,
        "_seed_h" + str("{:03d}".format(h_min)),
        new_dirname=experiment.astec_dir.get_tmp_directory(),
        new_extension=experiment.default_image_suffix,
    )
    difference_image = common.add_suffix(
        image_for_seed,
        "_seed_diff_h" + str("{:03d}".format(h_min)),
        new_dirname=experiment.astec_dir.get_tmp_directory(),
        new_extension=experiment.default_image_suffix,
    )

    if (
        not os.path.isfile(seed_image)
        or not os.path.isfile(difference_image)
        or monitoring.forceResultsToBeBuilt is True
    ):
        #
        # computation of labeled regional minima
        # -> keeping the 'difference' image allows to speed up the further computation
        #    for smaller values of h
        #
        mars.build_seeds(
            input_image, difference_image, unmasked_seed_image, experiment, wparam
        )
        #
        # select only the 'seeds' that are totally included in cells
        #
        cpp_wrapping.mc_mask_seeds(unmasked_seed_image, first_segmentation, seed_image)

    #
    # make a binary image (cell_segmentation) with two labels 0 and cellid
    #
    im_segmentation = imread(first_segmentation)
    cell_segmentation = np.zeros_like(im_segmentation)
    cell_segmentation[im_segmentation == cellid] = cellid
    np_unique = np.unique(cell_segmentation)
    if len(np_unique) != 2:
        monitoring.to_log_and_console(
            "       .. weird, sub-image of cell "
            + str(cellid)
            + " contains "
            + str(len(np_unique))
            + " labels = "
            + str(np_unique),
            2,
        )
    del im_segmentation

    checking = True
    while checking:
        #
        # get the seeds inside the cell
        #
        im_seed = imread(seed_image)
        labels = list(np.unique(im_seed[cell_segmentation == cellid]))
        if 0 in labels:
            labels.remove(0)

        #
        # two seeds? we're done
        # create an image with 0, 2, and 3 labels
        #
        if len(labels) == 2:
            two_seeds = (
                2 * (im_seed == labels[0]) + 3 * (im_seed == labels[1])
            ).astype(im_seed.dtype)
            res_image = common.add_suffix(
                image_for_seed,
                "_two_seeds",
                new_dirname=experiment.astec_dir.get_tmp_directory(),
                new_extension=experiment.default_image_suffix,
            )
            imsave(res_image, two_seeds)
            del two_seeds
            del im_seed
            return 2
        #
        #
        #
        elif len(labels) > 2:
            monitoring.to_log_and_console(
                "       .. too many extrema/seeds ("
                + str(len(labels))
                + ") for cell "
                + str(cellid),
                2,
            )
            return len(labels)

        #
        # change h_min
        #
        h_min -= parameters.watershed_seed_hmin_delta_value
        if h_min < parameters.watershed_seed_hmin_min_value:
            monitoring.to_log_and_console(
                "       .. last extrema/seeds number was ("
                + str(len(labels))
                + ") for cell "
                + str(cellid),
                2,
            )
            return len(labels)

        wparam.seed_hmin = h_min

        input_image = difference_image
        unmasked_seed_image = common.add_suffix(
            image_for_seed,
            "_unmasked_seed_h" + str("{:03d}".format(h_min)),
            new_dirname=experiment.astec_dir.get_tmp_directory(),
            new_extension=experiment.default_image_suffix,
        )
        seed_image = common.add_suffix(
            image_for_seed,
            "_seed_h" + str("{:03d}".format(h_min)),
            new_dirname=experiment.astec_dir.get_tmp_directory(),
            new_extension=experiment.default_image_suffix,
        )
        difference_image = common.add_suffix(
            image_for_seed,
            "_seed_diff_h" + str("{:03d}".format(h_min)),
            new_dirname=experiment.astec_dir.get_tmp_directory(),
            new_extension=experiment.default_image_suffix,
        )

        if (
            not os.path.isfile(seed_image)
            or not os.path.isfile(difference_image)
            or monitoring.forceResultsToBeBuilt is True
        ):
            mars.build_seeds(
                input_image,
                difference_image,
                unmasked_seed_image,
                experiment,
                wparam,
                operation_type="max",
            )
            cpp_wrapping.mc_mask_seeds(
                unmasked_seed_image, first_segmentation, seed_image
            )


def _two_seeds_propagation(c, experiment, tmp_prefix_name, label_width=4):
    #
    # try to get two seeds
    #
    cellid = str("{:0{width}d}".format(c, width=label_width))
    cell_seg_name = (
        tmp_prefix_name + "_cell" + cellid + "_seg." + experiment.default_image_suffix
    )
    cell_def_seed_name = (
        tmp_prefix_name
        + "_cell"
        + cellid
        + "_deformed_seeds."
        + experiment.default_image_suffix
    )

    #
    #
    #
    imseg = imread(cell_seg_name)
    im_seed = imread(cell_def_seed_name)
    labels = list(np.unique(im_seed[imseg == c]))
    if 0 in labels:
        labels.remove(0)
    del imseg

    #
    # two seeds? we're done
    # create an image with 0, 2, and 3 labels
    #
    if len(labels) == 2:
        two_seeds = (2 * (im_seed == labels[0]) + 3 * (im_seed == labels[1])).astype(
            im_seed.dtype
        )
        res_image = (
            tmp_prefix_name
            + "_cell"
            + cellid
            + "_seeds_two_seeds."
            + experiment.default_image_suffix
        )
        imsave(res_image, two_seeds)
        del two_seeds
        del im_seed
        return labels

    return labels


def _two_seeds_watershed(c, experiment, parameters, tmp_prefix_name, label_width=4):
    #
    # try to get two seeds
    #
    cellid = str("{:0{width}d}".format(c, width=label_width))
    cell_seg_name = (
        tmp_prefix_name + "_cell" + cellid + "_seg." + experiment.default_image_suffix
    )
    #
    #
    #
    cell_seed_name = (
        tmp_prefix_name
        + "_cell"
        + cellid
        + "_seeds_two_seeds."
        + experiment.default_image_suffix
    )
    cell_memb_name = (
        tmp_prefix_name
        + "_cell"
        + cellid
        + "_membrane."
        + experiment.default_image_suffix
    )
    cell_new_memb_name = (
        tmp_prefix_name
        + "_cell"
        + cellid
        + "_new_membrane."
        + experiment.default_image_suffix
    )
    cell_new_seed_name = (
        tmp_prefix_name
        + "_cell"
        + cellid
        + "_new_seeds."
        + experiment.default_image_suffix
    )
    #
    # get cell border from cell segmentation
    # set membrane signal to maxmembrane: it ensures that the watershed will "stay" inside the cell
    # use also the eroded cell background as a seed
    #
    # seed labels are [2, 3]
    #
    imseg = imread(cell_seg_name)
    binseg = np.ones_like(imseg)
    binseg[imseg == c] = 0
    eroseg = nd.binary_erosion(binseg, iterations=2)
    cellborder = binseg - eroseg
    del imseg
    del binseg

    immemb = imread(cell_memb_name)
    maxmembrane = np.max(immemb) + 1
    immemb[cellborder == 1] = maxmembrane
    imsave(cell_new_memb_name, immemb)
    del immemb
    del cellborder

    imseed = imread(cell_seed_name)
    imseed[eroseg == 1] = 1
    imsave(cell_new_seed_name, imseed)
    del eroseg
    del imseed

    #
    # watershed from the two seeds
    #
    cell_wate_name = (
        tmp_prefix_name
        + "_cell"
        + cellid
        + "_watershed."
        + experiment.default_image_suffix
    )
    mars.watershed(cell_new_seed_name, cell_new_memb_name, cell_wate_name, parameters)

    return


def _update_image_properties(
    im_segmentation,
    tmp_prefix_name,
    current_time,
    bounding_boxes,
    c,
    newlabel,
    properties,
    experiment,
    new_propagated_division,
    previous_labels=None,
    label_width=4,
    time_digits_for_cell_id=4,
):
    #
    # watershed has 3 labels 1, 2 and 3
    #
    cellid = str("{:0{width}d}".format(c, width=label_width))
    cell_wate_name = (
        tmp_prefix_name
        + "_cell"
        + cellid
        + "_watershed."
        + experiment.default_image_suffix
    )
    #
    # update segmentation
    #
    newseg = imread(cell_wate_name)
    newseg[im_segmentation[bounding_boxes[c]] != c] = 0
    im_segmentation[bounding_boxes[c]][newseg == 2] = c
    im_segmentation[bounding_boxes[c]][newseg == 3] = newlabel

    #
    # update volumes
    #
    cell_label = np.unique(newseg)
    cell_volume = nd.sum(np.ones_like(newseg), newseg, index=np.int16(cell_label))
    newseg_volumes = dict(zip(cell_label, cell_volume))

    cellid1 = current_time * 10**time_digits_for_cell_id + c
    cellid2 = current_time * 10**time_digits_for_cell_id + newlabel
    if cellid2 in properties["cell_volume"]:
        monitoring.to_log_and_console(
            "    .. weird, " + str(cellid2) + " was already in volume dictionary"
        )
    properties["cell_volume"][cellid1] = int(newseg_volumes[2])
    properties["cell_volume"][cellid2] = int(newseg_volumes[3])

    #
    # update lineage
    #
    previous_time = current_time - experiment.delta_time_point

    if len(previous_labels) == 1:
        prevcell = previous_time * 10**time_digits_for_cell_id + previous_labels[0]
        properties["cell_lineage"][prevcell] = [cellid1, cellid2]
    elif len(previous_labels) == 2:
        prevcell = previous_time * 10**time_digits_for_cell_id + previous_labels[0]
        properties["cell_lineage"][prevcell] = [cellid1]
        prevcell = previous_time * 10**time_digits_for_cell_id + previous_labels[1]
        properties["cell_lineage"][prevcell] = [cellid2]
    else:
        msg = (
            "      .. weird, cell(s) "
            + str(previous_labels)
            + " divides into "
            + str([c, newlabel])
        )
        monitoring.to_log_and_console(msg)
        msg = "         lineage is not updated"
        monitoring.to_log_and_console(msg)
        return new_propagated_division

    if cellid1 in properties["cell_lineage"]:
        if len(properties["cell_lineage"][cellid1]) == 1:
            t = properties["cell_lineage"][cellid1][0] // 10**time_digits_for_cell_id
            label = (
                properties["cell_lineage"][cellid1][0]
                - t * 10**time_digits_for_cell_id
            )
            new_propagated_division[t] = new_propagated_division.get(t, []) + [label]
        else:
            msg = "      .. weird, cell " + str(cellid1) + " divides into "
            msg += str(properties["cell_lineage"][cellid1])
            monitoring.to_log_and_console(msg)
            msg = "         should only divides into one cell. Will not propagate division any further."
            monitoring.to_log_and_console(msg)

    return new_propagated_division


def _image_cell_division(
    input_image,
    output_image,
    current_time,
    properties,
    experiment,
    parameters,
    new_division,
    propagated_division={},
    time_digits_for_cell_id=4,
):
    proc = "_image_cell_division"

    new_propagated_division = {}
    #
    # set working dir to mars dir, try to recover reconstructed image, if any
    #
    previous_time = current_time - experiment.delta_time_point

    #
    # set and make temporary directory
    #
    experiment.astec_dir.set_tmp_directory(current_time)
    experiment.astec_dir.make_tmp_directory()

    if (
        parameters.seed_reconstruction.keep_reconstruction is False
        and parameters.membrane_reconstruction.keep_reconstruction is False
        and parameters.morphosnake_reconstruction.keep_reconstruction is False
    ):
        experiment.mars_dir.set_rec_directory_to_tmp()
        experiment.astec_dir.set_rec_directory_to_tmp()
    else:
        experiment.astec_dir.make_rec_directory()
    reconstruction.monitoring.copy(monitoring)

    #
    # get bounding boxes of cells of interest and extract subimages
    # => dictionary of bounding boxes indexed by cell labels.
    # tmp_prefix_name + "_cell" + cellid + "_seg." + experiment.default_image_suffix
    # save also segmentation subimages
    #
    e = common.get_image_extension(input_image)
    b = os.path.basename(input_image)
    prefix_name = "ext_" + b[0 : len(b) - len(e)]
    tmp_prefix_name = os.path.join(
        experiment.astec_dir.get_tmp_directory(), prefix_name
    )

    dividing_cells = []
    if current_time in new_division:
        dividing_cells += copy.deepcopy(new_division[current_time])
    if current_time in propagated_division:
        dividing_cells += copy.deepcopy(propagated_division[current_time])

    dilation_iterations = 20
    label_width = 5
    bounding_boxes = _get_bounding_box(
        input_image,
        experiment,
        dividing_cells,
        tmp_prefix_name,
        dilation_iterations=dilation_iterations,
        label_width=label_width,
    )
    #
    #
    #
    membrane_image = _get_reconstructed_image(
        previous_time, current_time, experiment, parameters.membrane_reconstruction
    )
    if membrane_image is None:
        monitoring.to_log_and_console(
            "    .. "
            + proc
            + ": no membrane image was found/built for time "
            + str(current_time),
            2,
        )
        return properties, new_propagated_division

    #
    # extract and save membranes sub-images
    #
    im_membrane = imread(membrane_image)
    voxelsize = im_membrane.get_voxelsize()
    datatype = im_membrane.dtype
    for c in dividing_cells:
        cellid = str("{:0{width}d}".format(c, width=label_width))
        cell_memb_name = (
            tmp_prefix_name
            + "_cell"
            + cellid
            + "_membrane."
            + experiment.default_image_suffix
        )
        if not os.path.isfile(cell_memb_name):
            ext_membrane = im_membrane[bounding_boxes[c]]
            imsave(
                cell_memb_name,
                SpatialImage(ext_membrane, voxelsize=voxelsize).astype(datatype),
            )
            del ext_membrane
    del im_membrane

    #
    # at this point, sub-images of
    # - segmentation
    # - membrane
    # have been saved
    #

    #
    # de novo division: compute seeds
    #
    if current_time in new_division:
        #
        # extract seeds intensity sub-images
        #
        if parameters.seed_reconstruction.is_equal(parameters.membrane_reconstruction):
            monitoring.to_log_and_console(
                "    .. seed image is identical to membrane image", 2
            )
            image_for_seed = membrane_image
        else:
            image_for_seed = _get_reconstructed_image(
                previous_time, current_time, experiment, parameters.seed_reconstruction
            )
        if image_for_seed is None:
            monitoring.to_log_and_console(
                "    .. "
                + proc
                + " no seed image was found/built for time "
                + str(current_time),
                2,
            )
            return properties, new_propagated_division

        if parameters.seed_reconstruction.is_equal(parameters.membrane_reconstruction):
            for c in new_division[current_time]:
                cellid = str("{:0{width}d}".format(c, width=label_width))
                cell_memb_name = (
                    tmp_prefix_name
                    + "_cell"
                    + cellid
                    + "_membrane."
                    + experiment.default_image_suffix
                )
                cell_seed_name = (
                    tmp_prefix_name
                    + "_cell"
                    + cellid
                    + "_seeds."
                    + experiment.default_image_suffix
                )
                if not os.path.isfile(cell_seed_name):
                    shutil.copy2(cell_memb_name, cell_seed_name)
        else:
            im_for_seeds = imread(image_for_seed)
            voxelsize = im_for_seeds.get_voxelsize()
            datatype = im_for_seeds.dtype
            for c in new_division[current_time]:
                cellid = str("{:0{width}d}".format(c, width=label_width))
                cell_seed_name = (
                    tmp_prefix_name
                    + "_cell"
                    + cellid
                    + "_seeds."
                    + experiment.default_image_suffix
                )
                if not os.path.isfile(cell_seed_name):
                    ext_for_seeds = im_for_seeds[bounding_boxes[c]]
                    imsave(
                        cell_seed_name,
                        SpatialImage(ext_for_seeds, voxelsize=voxelsize).astype(
                            datatype
                        ),
                    )
                    del ext_for_seeds
            del im_for_seeds

        #
        # transform input_image into division_image
        #
        if current_time in propagated_division:
            division_image = common.add_suffix(
                input_image,
                "_division",
                new_dirname=experiment.astec_dir.get_tmp_directory(0),
                new_extension=experiment.default_image_suffix,
            )
        else:
            division_image = output_image

        #
        # segmentation sub-image
        # tmp_prefix_name + "_cell" + cellid + "_seg." + experiment.default_image_suffix
        # intensity sub-image (for seed extraction)
        # tmp_prefix_name + "_cell" + cellid + "_seeds." + experiment.default_image_suffix
        # intensity sub-image (for watershed)
        # tmp_prefix_name + "_cell" + cellid + "_membrane." + experiment.default_image_suffix
        #
        #
        # two seeds sub-image (labels are 0, 1 and 2)
        # tmp_prefix_name + "_cell" + cellid + "_seeds_two_seeds." + experiment.default_image_suffix
        #

        im_segmentation = imread(input_image)
        newlabel = int(np.max(im_segmentation))

        for c in new_division[current_time]:
            newlabel += 1
            #
            # try to get two seeds + watershed
            #
            cellid = str("{:0{width}d}".format(c, width=label_width))
            cell_seg_name = (
                tmp_prefix_name
                + "_cell"
                + cellid
                + "_seg."
                + experiment.default_image_suffix
            )
            cell_seed_name = (
                tmp_prefix_name
                + "_cell"
                + cellid
                + "_seeds."
                + experiment.default_image_suffix
            )
            nseeds = _two_seeds_extraction(
                cell_seg_name, c, cell_seed_name, experiment, parameters
            )
            if nseeds != 2:
                msg = (
                    "    .. unable to find two seeds for cell "
                    + str(c)
                    + " at time "
                    + str(current_time)
                )
                monitoring.to_log_and_console(msg)
                continue

            _two_seeds_watershed(
                c, experiment, parameters, tmp_prefix_name, label_width=label_width
            )

            new_propagated_division = _update_image_properties(
                im_segmentation,
                tmp_prefix_name,
                current_time,
                bounding_boxes,
                c,
                newlabel,
                properties,
                experiment,
                new_propagated_division,
                previous_labels=[c],
                label_width=label_width,
                time_digits_for_cell_id=time_digits_for_cell_id,
            )

        #
        if current_time not in propagated_division:
            imsave(division_image, im_segmentation)
            return properties, new_propagated_division
        input_image = division_image

    #
    # division issued from a previously computed division
    # retrieve seeds from previous segmentation image
    #
    if current_time in propagated_division:
        #
        # build seeds by eroding previous segmentation and deforming it
        #
        # erosion iterations are set by default in voxel units
        # there is also a volume defined in voxel units
        #
        monitoring.to_log_and_console("    build seeds from previous segmentation", 2)

        previous_segmentation = experiment.get_segmentation_image(previous_time)
        if previous_segmentation is None:
            monitoring.to_log_and_console(
                "    .. "
                + proc
                + ": no segmentation image was found for time "
                + str(previous_time)
            )
            return False

        #
        # transform the previous segmentation image then erode the cells
        # adapted from astec.py
        #
        deformed_seeds = common.add_suffix(
            input_image,
            "_deformed_seeds_from_previous",
            new_dirname=experiment.astec_dir.get_tmp_directory(0),
            new_extension=experiment.default_image_suffix,
        )
        #
        # deformed_segmentation will be needed for morphosnake correction
        # may also be used in reconstruction.py
        #
        deformed_segmentation = common.add_suffix(
            input_image,
            "_deformed_segmentation_from_previous",
            new_dirname=experiment.astec_dir.get_tmp_directory(0),
            new_extension=experiment.default_image_suffix,
        )

        if (
            not os.path.isfile(deformed_segmentation)
            or monitoring.forceResultsToBeBuilt is True
        ):
            deformation = reconstruction.get_deformation_from_current_to_previous(
                current_time,
                experiment,
                parameters.membrane_reconstruction,
                previous_time,
            )
            if deformation is None:
                monitoring.to_log_and_console(
                    "    .. " + proc + ": error when getting deformation field"
                )
                return False

            cpp_wrapping.apply_transformation(
                previous_segmentation,
                deformed_segmentation,
                deformation,
                interpolation_mode="nearest",
                monitoring=monitoring,
            )

        astec.build_seeds_from_previous_segmentation(
            deformed_segmentation, deformed_seeds, parameters
        )

        #
        # save sub-images of the deformed seeds
        #
        im_for_seeds = imread(deformed_seeds)
        voxelsize = im_for_seeds.get_voxelsize()
        datatype = im_for_seeds.dtype
        for c in propagated_division[current_time]:
            cellid = str("{:0{width}d}".format(c, width=label_width))
            cell_seed_name = (
                tmp_prefix_name
                + "_cell"
                + cellid
                + "_deformed_seeds."
                + experiment.default_image_suffix
            )
            if not os.path.isfile(cell_seed_name):
                ext_for_seeds = im_for_seeds[bounding_boxes[c]]
                imsave(
                    cell_seed_name,
                    SpatialImage(ext_for_seeds, voxelsize=voxelsize).astype(datatype),
                )
                del ext_for_seeds
        del im_for_seeds

        im_segmentation = imread(input_image)
        newlabel = int(np.max(im_segmentation))

        for c in propagated_division[current_time]:
            newlabel += 1

            #
            # previous_labels are the labels in previous image that will correspond to [2, 3]
            #
            previous_labels = _two_seeds_propagation(
                c, experiment, tmp_prefix_name, label_width=label_width
            )
            if len(previous_labels) != 2:
                msg = (
                    "    .. unable to find two seeds for cell "
                    + str(c)
                    + " at time "
                    + str(current_time)
                )
                monitoring.to_log_and_console(msg)
                continue

            _two_seeds_watershed(
                c, experiment, parameters, tmp_prefix_name, label_width=label_width
            )

            new_propagated_division = _update_image_properties(
                im_segmentation,
                tmp_prefix_name,
                current_time,
                bounding_boxes,
                c,
                newlabel,
                properties,
                experiment,
                new_propagated_division,
                previous_labels=previous_labels,
                label_width=label_width,
                time_digits_for_cell_id=time_digits_for_cell_id,
            )

        imsave(output_image, im_segmentation)
        del im_segmentation

    return properties, new_propagated_division


########################################################################################
#
#
#
########################################################################################

#
#
#
#
#


def correction_process(
    current_time,
    properties,
    experiment,
    parameters,
    fusion,
    new_division,
    propagated_division={},
):
    proc = "correction_process"

    mars_dir = experiment.mars_dir.get_directory(0)
    astec_dir = experiment.astec_dir.get_directory(0)

    #
    # input image is in EXP_SEG_FROM ie MARS subdirectory
    # - it can be suffixed either by 'mars' or by 'seg'
    # output image will be in EXP_SEG_TO ie ASTEC subdirectory
    #
    mars_name = experiment.mars_dir.get_image_name(current_time)
    seg_name = experiment.astec_dir.get_image_name(current_time)
    input_image = None

    seg_image = common.find_file(
        mars_dir,
        mars_name,
        file_type="image",
        callfrom=proc,
        local_monitoring=None,
        verbose=False,
    )
    if seg_image is not None:
        input_image = os.path.join(mars_dir, seg_image)
    else:
        seg_image = common.find_file(
            mars_dir,
            seg_name,
            file_type="image",
            callfrom=proc,
            local_monitoring=None,
            verbose=False,
        )
        if seg_image is not None:
            input_image = os.path.join(mars_dir, seg_image)
    if input_image is None:
        monitoring.to_log_and_console(
            "    .. "
            + proc
            + ": no segmentation image was found for time "
            + str(current_time)
        )
        monitoring.to_log_and_console("    .. exiting.")
        sys.exit(1)

    #
    # output image
    #
    seg_image = common.find_file(
        astec_dir,
        seg_name,
        file_type="image",
        callfrom=proc,
        local_monitoring=None,
        verbose=False,
    )

    #
    # nothing to do
    #
    if seg_image is not None and monitoring.forceResultsToBeBuilt is False:
        monitoring.to_log_and_console(
            "    corrected image '" + str(seg_image) + "' exists", 2
        )
        return properties, {}

    #
    #
    #
    if seg_image is None:
        output_image = os.path.join(
            astec_dir, seg_name + "." + experiment.result_image_suffix
        )
    else:
        output_image = os.path.join(astec_dir, seg_image)

    #
    # nothing to do, copy the image
    #
    if (
        current_time not in fusion
        and current_time not in new_division
        and current_time not in propagated_division
    ):
        cpp_wrapping.copy(input_image, output_image)
        monitoring.to_log_and_console(
            "    no corrections to be done for time " + str(current_time), 2
        )
        return properties, {}

    #
    # something to do
    #
    new_propagated_division = {}
    monitoring.to_log_and_console(
        "... correction of '" + str(input_image).split(os.path.sep)[-1] + "'", 1
    )

    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()
    #
    # cell fusion
    #
    if current_time in fusion:
        if current_time in new_division or current_time in propagated_division:
            fusedcell_image = common.add_suffix(
                input_image,
                experiment.result_image_suffix + "_fusedcell",
                new_dirname=experiment.astec_dir.get_tmp_directory(0),
                new_extension=experiment.default_image_suffix,
            )
        else:
            fusedcell_image = output_image
        properties = _image_cell_fusion(
            input_image,
            fusedcell_image,
            current_time,
            properties,
            fusion,
            time_digits_for_cell_id=time_digits_for_cell_id,
        )
        input_image = fusedcell_image

    #
    # cell division
    #
    if current_time in new_division or current_time in propagated_division:
        properties, new_propagated_division = _image_cell_division(
            input_image,
            output_image,
            current_time,
            properties,
            experiment,
            parameters,
            new_division,
            propagated_division,
            time_digits_for_cell_id=time_digits_for_cell_id,
        )

    return properties, new_propagated_division


#
#
#
#
#


def correction_control(experiment, parameters):
    """

    :param experiment:
    :param parameters:
    :return:
    """

    proc = "correction_control"
    start_time = time.time()
    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(
            str(proc)
            + ": unexpected type for 'experiment' variable: "
            + str(type(experiment))
        )
        sys.exit(1)

    if not isinstance(parameters, ManualCorrectionParameters):
        monitoring.to_log_and_console(
            str(proc)
            + ": unexpected type for 'parameters' variable: "
            + str(type(parameters))
        )
        sys.exit(1)

    #
    # make sure that the result directory exists
    #

    experiment.astec_dir.make_directory()
    monitoring.to_log_and_console("", 1)

    mars_dir = experiment.mars_dir.get_directory(0)
    astec_dir = experiment.astec_dir.get_directory(0)

    #
    # set last time point if required
    #

    if (experiment.first_time_point is None or experiment.first_time_point < 0) and (
        experiment.last_time_point is None or experiment.last_time_point < 0
    ):
        monitoring.to_log_and_console(
            "... time interval does not seem to be defined in the parameter file"
        )
        monitoring.to_log_and_console("    set parameters 'begin' and 'end'")
        monitoring.to_log_and_console("\t Exiting")
        sys.exit(1)

    if experiment.first_time_point is not None and experiment.first_time_point >= 0:
        if experiment.last_time_point is None:
            experiment.last_time_point = experiment.first_time_point
        elif experiment.last_time_point < experiment.first_time_point:
            monitoring.to_log_and_console(
                "... weird time interval: 'begin' = "
                + str(experiment.first_time_point)
                + ", 'end' = "
                + str(experiment.last_time_point)
            )

    #
    # read lineage
    # use experiment.astec_dir.get_file_name("_lineage") to have a well-formed fille name
    #
    lineage = None
    properties = None
    lineage_tree_file = common.find_file(
        mars_dir,
        experiment.astec_dir.get_file_name("_lineage"),
        file_type="lineage",
        callfrom=proc,
        verbose=False,
    )
    if lineage_tree_file is not None and os.path.isfile(
        os.path.join(mars_dir, lineage_tree_file)
    ):
        lineage_tree_path = os.path.join(mars_dir, lineage_tree_file)
        properties = ioproperties.read_dictionary(lineage_tree_path)
        #
        # will only update volume and lineage
        #
        keylist = list(properties.keys())
        for k in keylist:
            if k != "cell_lineage" and k != "cell_volume":
                del properties[k]
        if "cell_lineage" in properties:
            lineage = properties["cell_lineage"]

    #
    #
    #
    fusion, division = _read_correction_file(
        parameters.manualcorrection_file,
        lineage=lineage,
        first_time_point=experiment.first_time_point,
        time_digits_for_cell_id=experiment.get_time_digits_for_cell_id(),
    )

    propagated_division = {}
    for current_time in range(
        experiment.first_time_point,
        experiment.last_time_point + 1,
        experiment.delta_time_point,
    ):
        properties, propagated_division = correction_process(
            current_time,
            properties,
            experiment,
            parameters,
            fusion,
            division,
            propagated_division,
        )
        #
        # save lineage here (if there is a lineage to be saved)
        #
        if properties is not None:
            lineage_tree_path = os.path.join(
                astec_dir,
                experiment.astec_dir.get_file_name("_lineage")
                + "."
                + experiment.result_lineage_suffix,
            )
            ioproperties.write_dictionary(lineage_tree_path, properties)
        #
        # cleaning
        #
        if monitoring.keepTemporaryFiles is False:
            experiment.astec_dir.rmtree_tmp_directory()

    end_time = time.time()
    monitoring.to_log_and_console(
        "    computation time = " + str(end_time - start_time) + " s", 1
    )
    monitoring.to_log_and_console("", 1)
