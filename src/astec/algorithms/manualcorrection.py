
import os
import sys
import time
import numpy as np
import shutil
from scipy import ndimage as nd

import astec.algorithms.astec as astec
from astec.utils import common
from astec.components.spatial_image import SpatialImage
from astec.io.image import imread, imsave
import astec.utils.diagnosis as udiagnosis
import astec.utils.ioproperties as ioproperties

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

        doc = "\n"
        doc += "Manual correction overview:\n"
        doc += "===========================\n"
        doc += "Fuses labels from the 'mars' segmentation to correct\n"
        doc += "over-segmentation errors.\n"
        doc += "\n"
        self.doc['manualcorrection_overview'] = doc

        doc = "\t Directory containing the manual correction files\n"
        self.doc['manualcorrection_dir'] = doc
        self.manualcorrection_dir = None

        doc = "\t File containing the labels to be fused.\n"
        doc += "\t The syntax of this file is:\n"
        doc += "\t - 1 line per label association, eg\n"
        doc += "\t   '8 7'\n"
        doc += "\t - background label has value 1\n"
        doc += "\t - the character '#' denotes commented lines \n"
        self.doc['manualcorrection_file'] = doc
        self.manualcorrection_file = None

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print('')
        print('#')
        print('# ManualCorrectionParameters ')
        print('#')
        print('')

        common.PrefixedParameter.print_parameters(self)

        udiagnosis.DiagnosisParameters.print_parameters(self)

        for line in self.doc['manualcorrection_overview'].splitlines():
            print('# ' + line)

        self.varprint('manualcorrection_dir', self.manualcorrection_dir, self.doc['manualcorrection_dir'])
        self.varprint('manualcorrection_file', self.manualcorrection_file, self.doc['manualcorrection_file'])

        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write('\n')
        logfile.write('#' + '\n')
        logfile.write('# ManualCorrectionParameters ' + '\n')
        logfile.write('#' + '\n')
        logfile.write('\n')

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        udiagnosis.DiagnosisParameters.write_parameters_in_file(self, logfile)

        for line in self.doc['manualcorrection_overview'].splitlines():
            logfile.write('# ' + line + '\n')

        self.varwrite(logfile, 'manualcorrection_dir', self.manualcorrection_dir,
                      self.doc.get('manualcorrection_dir', None))
        self.varwrite(logfile, 'manualcorrection_file', self.manualcorrection_file,
                      self.doc.get('manualcorrection_file', None))

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

        self.manualcorrection_dir = self.read_parameter(parameters, 'manualcorrection_dir', self.manualcorrection_dir)
        self.manualcorrection_file = self.read_parameter(parameters, 'manualcorrection_file',
                                                         self.manualcorrection_file)
        self.manualcorrection_file = self.read_parameter(parameters, 'mapping_file', self.manualcorrection_file)

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
# some internal procedures
#
########################################################################################

def _read_correction_file(filename, lineage=None, first_time_point=1, time_digits_for_cell_id=4):
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
    separator_time_interval = ['-', ':']
    separator_time_label = [':', ',']

    if not os.path.isfile(filename):
        monitoring.to_log_and_console("Error: '" + filename + "' is not a valid file. Exiting.")
        sys.exit(1)
    f = open(filename)
    i = 0
    for line in f:
        i += 1
        # remove leading and trailing whitespaces
        li = line.strip()
        # skip comment
        if li.startswith('#'):
            continue
        # empty line
        if len(li) == 0:
            continue
        li = li.split()

        # one component: single integer, division at first time point
        if len(li) == 1:
            if li[0].isdigit():
                division[first_time_point] = division.get(first_time_point, []) + [int(li[0])]
                continue
            msg = "line #" + str(i) + ": '" + str(line) + "' should be a single integer?! Skip the line."
            monitoring.to_log_and_console(proc + ": " + msg)
            continue

        # two or more components: only integers or "t:" + integers
        # last component(s) should be only integer
        error = False
        for i in range(1, len(li)):
            if not li[i].isdigit():
                msg = str(i+1) + "th term of line #" + str(i) + ": '" + str(line) + "' should be numbers."
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
                msg = "line #" + str(i) + ": '" + str(line) + "' has the same label '" + str(li[0])
                msg += "for fusion?! ... Skip the line."
                monitoring.to_log_and_console(proc + ": " + msg)
                continue
            if first_time_point not in fusion:
                fusion[first_time_point] = {}
            # use minimal label as destination label
            for lab in labels:
                if lab > min(labels):
                    fusion[first_time_point][lab] = fusion[first_time_point].get(lab, []) + [min(labels)]
            continue

        # first term is not a number, try to recognize time (interval)
        # "t:" + integer -> division
        # "t:" + integers -> fusion
        # "t-t:" + integers -> fusion
        if li[0][-1] not in separator_time_label:
            msg = "last character of first term of line #" + str(i) + ": '" + str(line) + "' should be in "
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
            for ctime in range(mintime, maxtime+1):
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
            prev_cell_ids = [10 ** time_digits_for_cell_id + lab for lab in labels]
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
            labels = [c - (c // 10 ** time_digits_for_cell_id) for c in next_cell_ids]

    f.close()
    print("division = " + str(division))
    print("fusion = " + str(fusion))
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

    monitoring.to_log_and_console('    Number of cells: ' + str(len(cell_label)), 0)
    monitoring.to_log_and_console('    Maximal label: ' + str(np.max(im)), 0)
    # monitoring.to_log_and_console('    Cell ids: ' + str(cell_label), 0)

    monitoring.to_log_and_console('    Sorted cell volumes: ', 0)
    monitoring.to_log_and_console('      Id :    voxels          (um^3)', 0)

    if parameters.items <= 0 or parameters.items >= len(volumes):
        for v in volumes:
            msg = '    {:>4d} : {:>9d} {:>15s}'.format(v[0], int(v[1]), '({:.2f})'.format(v[1]*vol))
            monitoring.to_log_and_console(msg, 0)
    else:
        if int(parameters.items) > 0:
            for v in volumes[:parameters.items]:
                msg = '    {:>4d} : {:>9d} {:>15s}'.format(v[0], int(v[1]), '({:.2f})'.format(v[1] * vol))
                monitoring.to_log_and_console(msg, 0)
            monitoring.to_log_and_console('       ...', 0)
            for v in volumes[-parameters.items:]:
                msg = '    {:>4d} : {:>9d} {:>15s}'.format(v[0], int(v[1]), '({:.2f})'.format(v[1] * vol))
                monitoring.to_log_and_console(msg, 0)

########################################################################################
#
#
#
########################################################################################


def _cell_fusion(input_image, output_image, fusion):

    im = imread(input_image)
    voxelsize = im.get_voxelsize()
    type = im.dtype

    #
    #
    #
    immax = np.max(im)
    mapping = np.arange(immax + 1)
    for lab in fusion:
        eqlabel = min(min(fusion[lab]), mapping[lab])
        labels = [lab] + fusion[lab] + [mapping[lab]]
        for i in range(immax + 1):
            if mapping[i] in labels:
                mapping[i] = eqlabel

    im = mapping[im]
    imsave(output_image, SpatialImage(im, voxelsize=voxelsize).astype(type))

#
#
#
#
#

def correction_process(time_value, fusion, division, lineage_tree_information, experiment, parameters):
    """

    Parameters
    ----------
    time_value
    fusion
    division
    lineage_tree_information
    experiment
    parameters

    Returns
    -------

    """

    proc = "correction_process"

    mars_dir = experiment.mars_dir.get_directory(0)
    astec_dir = experiment.astec_dir.get_directory(0)

    #
    # input image is in EXP_SEG_FROM ie MARS subdirectory
    # - it can be suffixed either by 'mars' or by 'seg'
    # output image will be in EXP_SEG_TO ie ASTEC subdirectory
    #
    mars_name = experiment.mars_dir.get_image_name(time_value)
    seg_name = experiment.astec_dir.get_image_name(time_value)
    input_image = None

    seg_image = common.find_file(mars_dir, mars_name, file_type='image', callfrom=proc, local_monitoring=None,
                                 verbose=False)
    if seg_image is not None:
        input_image = os.path.join(mars_dir, seg_image)
    else:
        seg_image = common.find_file(mars_dir, seg_name, file_type='image', callfrom=proc, local_monitoring=None,
                                     verbose=False)
        if seg_image is not None:
            input_image = os.path.join(mars_dir, seg_image)
    if input_image is None:
        monitoring.to_log_and_console("    .. " + proc + ": no segmentation image was found for time "
                                      + str(time_value))
        monitoring.to_log_and_console("    .. exiting.")
        sys.exit(1)

    #
    # output image
    #
    seg_image = common.find_file(astec_dir, seg_name, file_type='image', callfrom=proc, local_monitoring=None,
                                 verbose=False)

    #
    # nothing to do
    #
    if seg_image is not None and monitoring.forceResultsToBeBuilt is False:
        monitoring.to_log_and_console("    corrected image '" + str(seg_image) + "' exists", 2)
        return

    #
    #
    #
    if seg_image is None:
        output_image = os.path.join(astec_dir, seg_name + '.' + experiment.result_image_suffix)
    else:
        output_image = os.path.join(astec_dir, seg_image)

    #
    # nothing to do, copy the image
    #
    if time_value not in fusion and time_value not in division:
        shutil.copy2(input_image, output_image)
        monitoring.to_log_and_console("    no corrections to be done for time " + str(time_value), 2)
        return

    #
    # something to do
    #
    monitoring.to_log_and_console("... correction of '" + str(input_image).split(os.path.sep)[-1] + "'", 1)
    start_time = time.time()

    #
    # cell fusion
    #
    if time_value in fusion:
        if time_value in division:
            fusedcell_image = common.add_suffix(input_image, experiment.result_image_suffix + "_fusedcell",
                                                new_dirname=experiment.astec_dir.get_tmp_directory(0),
                                                new_extension=experiment.default_image_suffix)
        else:
            fusedcell_image = output_image

        _cell_fusion(input_image, fusedcell_image, fusion[time_value])



    #
    # get
    # - the list of cell label
    # - the list of cell volume
    #
    # build a dictionary and sort it (increasing order) wrt the volume
    _diagnosis_volume_image(output_image, parameters)

    end_time = time.time()
    monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
    monitoring.to_log_and_console('', 1)

    return

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

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, ManualCorrectionParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    #
    # make sure that the result directory exists
    #

    experiment.astec_dir.make_directory()
    monitoring.to_log_and_console('', 1)

    mars_dir = experiment.mars_dir.get_directory(0)
    astec_dir = experiment.astec_dir.get_directory(0)

    #
    # set last time point if required
    #

    if (experiment.first_time_point is None or experiment.first_time_point < 0) and \
            (experiment.last_time_point is None or experiment.last_time_point < 0):
        monitoring.to_log_and_console("... time interval does not seem to be defined in the parameter file")
        monitoring.to_log_and_console("    set parameters 'begin' and 'end'")
        monitoring.to_log_and_console("\t Exiting")
        sys.exit(1)

    if experiment.first_time_point is not None and experiment.first_time_point >= 0:
        if experiment.last_time_point is None:
            experiment.last_time_point = experiment.first_time_point
        elif experiment.last_time_point < experiment.first_time_point:
            monitoring.to_log_and_console("... weird time interval: 'begin' = " + str(experiment.first_time_point)
                                          + ", 'end' = " + str(experiment.last_time_point))

    #
    # read lineage
    #
    lineage = None
    lineage_tree_information = None
    lineage_tree_file = common.find_file(mars_dir, experiment.mars_dir.get_file_name("_lineage"),
                                         file_type='lineage', callfrom=proc, verbose=False)

    if lineage_tree_file is not None and os.path.isfile(os.path.join(mars_dir, lineage_tree_file)):
        lineage_tree_path = os.path.join(mars_dir, lineage_tree_file)
        lineage_tree_information = ioproperties.read_dictionary(lineage_tree_path)
        if 'cell_lineage' in lineage_tree_information:
            lineage = lineage_tree_information['cell_lineage']
    #
    #
    #
    fusion, division = _read_correction_file(parameters.manualcorrection_file, lineage=lineage,
                                             first_time_point=experiment.first_time_point,
                                             time_digits_for_cell_id=experiment.get_time_digits_for_cell_id())

    for time_value in range(experiment.first_time_point, experiment.last_time_point + 1, experiment.delta_time_point):

        correction_process(time_value, fusion, division, lineage_tree_information, experiment, parameters)

    #
    # save lineage here (if there is a lineage to be saved)
    #
    if lineage_tree_information is not None:
        keylist = lineage_tree_information.keys()
        for k in keylist:
            if k != 'cell_lineage' and k != 'cell_volume':
                del lineage_tree_information[k]
        lineage_tree_path = os.path.join(astec_dir, experiment.astec_dir.get_file_name("_lineage") + "."
                                         + experiment.result_lineage_suffix)
        ioproperties.write_dictionary(lineage_tree_path, lineage_tree_information)
