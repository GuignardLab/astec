
import os
import imp
import sys
import time
import shutil

import ACE
import commonTools
import nomenclature
import lineage
import reconstruction

from CommunFunctions.ImageHandling import imread

#
#
#
#
#

monitoring = commonTools.Monitoring()


########################################################################################
#
# classes
# - computation environment
# - computation parameters
#
########################################################################################


class AstecEnvironment(object):

    def __init__(self):

        #
        # fusion data paths
        #
        self.path_fuse_exp = None
        self.path_fuse_exp_files = None

        # segmentation data paths
        #
        self.path_seg_exp = None
        self.path_seg_exp_files = None
        self.path_seg_exp_lineage = None

        #
        #
        #
        self.path_reconstruction = None
        self.temporary_path = None

        #
        #
        #
        self.path_logdir = None
        self.path_history_file = None
        self.path_log_file = None

    def update_from_file(self, parameter_file, start_time):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print ("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        self.path_fuse_exp = nomenclature.replaceFlags(nomenclature.path_fuse_exp, parameters)
        self.path_fuse_exp_files = nomenclature.replaceFlags(nomenclature.path_fuse_exp_files, parameters)

        self.path_seg_exp = nomenclature.replaceFlags(nomenclature.path_seg_exp, parameters)
        self.path_seg_exp_files = nomenclature.replaceFlags(nomenclature.path_seg_exp_files, parameters)
        self.path_seg_exp_lineage = nomenclature.replaceFlags(nomenclature.path_seg_exp_lineage, parameters)

        self.path_logdir = nomenclature.replaceFlags(nomenclature.path_seg_logdir, parameters)
        self.path_history_file = nomenclature.replaceFlags(nomenclature.path_seg_historyfile, parameters)
        self.path_log_file = nomenclature.replaceFlags(nomenclature.path_seg_logfile, parameters, start_time)

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('AstecEnvironment\n')

            logfile.write('- path_fuse_exp = ' + str(self.path_fuse_exp) + '\n')
            logfile.write('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files) + '\n')

            logfile.write('- path_seg_exp = ' + str(self.path_seg_exp) + '\n')
            logfile.write('- path_seg_exp_files = ' + str(self.path_seg_exp_files) + '\n')
            logfile.write('- path_seg_exp_lineage = ' + str(self.path_seg_exp_lineage) + '\n')

            logfile.write('- path_reconstruction = ' + str(self.path_reconstruction) + '\n')
            logfile.write('- temporary_path = ' + str(self.temporary_path) + '\n')

            logfile.write('- path_logdir = ' + str(self.path_logdir) + '\n')
            logfile.write('- path_history_file = ' + str(self.path_history_file)+'\n')
            logfile.write('- path_log_file = ' + str(self.path_log_file)+'\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('AstecEnvironment')

        print('- path_fuse_exp = ' + str(self.path_fuse_exp))
        print('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files))

        print('- path_seg_exp = ' + str(self.path_seg_exp))
        print('- path_seg_exp_files = ' + str(self.path_seg_exp_files))
        print('- path_seg_exp_lineage = ' + str(self.path_seg_exp_lineage))

        print('- path_reconstruction = ' + str(self.path_reconstruction))
        print('- temporary_path = ' + str(self.temporary_path))

        print('- path_logdir = ' + str(self.path_logdir))
        print('- path_history_file = ' + str(self.path_history_file))
        print('- path_log_file = ' + str(self.path_log_file))
        print("")


#
#
#
#
#


class AstecParameters(object):

    def __init__(self):
        #
        #
        #
        self.intensity_transformation = 'Identity'
        self.intensity_enhancement = None

        #
        # membrane enhancement parameters
        #
        self.ace = ACE.AceParameters()

        #
        #
        #
        self.keep_reconstruction = True

        #
        # images suffixes/formats
        #
        self.result_image_suffix = 'inr'
        self.default_image_suffix = 'inr'

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('AstecParameters\n')

            logfile.write('- intensity_transformation = ' + str(self.intensity_transformation) + '\n')
            logfile.write('- intensity_enhancement = ' + str(self.intensity_enhancement) + '\n')

            self.ace.write_parameters(log_file_name)

            logfile.write('- keep_reconstruction = ' + str(self.keep_reconstruction) + '\n')

            logfile.write('- result_image_suffix = ' + str(self.result_image_suffix) + '\n')
            logfile.write('- default_image_suffix = '+str(self.default_image_suffix) + '\n')

            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('AstecParameters')

        print('- intensity_transformation = ' + str(self.intensity_transformation))
        print('- intensity_enhancement = ' + str(self.intensity_enhancement))

        self.ace.print_parameters()

        print('- keep_reconstruction = ' + str(self.keep_reconstruction))

        print('- result_image_suffix = ' + str(self.result_image_suffix))
        print('- default_image_suffix = ' + str(self.default_image_suffix))

        print("")

    def update_from_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        #
        # reconstruction method
        #

        if hasattr(parameters, 'intensity_transformation'):
            self.intensity_transformation = parameters.intensity_transformation
        if hasattr(parameters, 'astec_intensity_transformation'):
            self.intensity_transformation = parameters.astec_intensity_transformation

        if hasattr(parameters, 'intensity_enhancement'):
            self.intensity_enhancement = parameters.intensity_enhancement
        if hasattr(parameters, 'astec_intensity_enhancement'):
            self.intensity_enhancement = parameters.astec_intensity_enhancement

        #
        #
        #
        self.ace.update_from_file(parameter_file)

        #
        # images suffixes/formats
        #
        if hasattr(parameters, 'result_image_suffix'):
            if parameters.result_image_suffix is not None:
                self.result_image_suffix = parameters.result_image_suffix
        if hasattr(parameters, 'default_image_suffix'):
            if parameters.default_image_suffix is not None:
                self.default_image_suffix = parameters.default_image_suffix


########################################################################################
#
# some internal procedures
#
########################################################################################

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

def astec_process(previous_time, current_time, lineage_tree_information, experiment, environment, parameters):


    proc = "astec_process"
    default_width = 3

    acquisition_time = str('{:0{width}d}'.format(current_time, width=default_width))

    #
    # nothing to do if the segmentation image exists
    #

    astec_name = nomenclature.replaceTIME(environment.path_seg_exp_files, current_time)
    astec_image = commonTools.find_file(environment.path_seg_exp, astec_name, monitoring=None, verbose=False)

    if astec_image is not None:
        if monitoring.forceResultsToBeBuilt is False:
            monitoring.to_log_and_console('    astec image already existing', 2)
            return
        else:
            monitoring.to_log_and_console('    astec image already existing, but forced', 2)

    astec_image = os.path.join(environment.path_seg_exp, astec_name + '.' + parameters.result_image_suffix)

    #
    # build the membrane image to be segmented
    #

    reconstruction.monitoring.copy(monitoring)
    membrane_image = reconstruction.build_membrane_image(current_time, environment, parameters,
                                                         previous_time=previous_time)

    #
    #
    #

    return


#
# check whether a lineage file exists
# loops over the time points
#
#


def astec_control(experiment, environment, parameters):
    """

    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    proc = "astec_control"
    default_width = 3

    #
    # make sure that the result directory exists
    #

    if not os.path.isdir(environment.path_seg_exp):
        monitoring.to_log_and_console(proc + ": weird, '" + str(environment.path_seg_exp) + "' does not exists", 1)
        monitoring.to_log_and_console("\t Exiting")
        sys.exit(1)

    if not os.path.isdir(environment.path_logdir):
        os.makedirs(environment.path_logdir)

    monitoring.to_log_and_console('', 1)

    #
    # re-read the lineage file, if any
    # and check whether any time point should be re-computed
    #

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    lineage_tree_information = lineage.read_lineage_tree(environment.path_seg_exp_lineage)

    if len(lineage_tree_information) > 0 and  'lin_tree' in lineage_tree_information:
        monitoring.to_log_and_console("    .. test '" + str(environment.path_seg_exp_lineage) + "'", 1)
        cellat = {}
        for y in lineage_tree_information['lin_tree']:
            t = y/10**4
            if t not in cellat:
                cellat[t] = 1
            else:
                cellat[t] += 1

        restart = -1
        t = first_time_point
        while restart == -1 and t <= last_time_point:
            #
            # possible time point of segmentation, test if ok
            #
            time_value = t + experiment.delta_time_point
            acquisition_time = str('{:0{width}d}'.format(time_value, width=default_width))
            segmentation_file = environment.path_seg_exp_files.replace(nomenclature.FLAG_TIME, acquisition_time)
            if not os.path.isfile(os.path.join(environment.path_seg_exp, segmentation_file)):
                monitoring.to_log_and_console("       image '" + segmentation_file + "' not found", 1)
                restart=t
            else:
                if cellat[t] == 0:
                    monitoring.to_log_and_console("       lineage of image '" + segmentation_file + "' not found", 1)
                    restart=t
                else:
                    try :
                        segmentation_image = imread(segmentation_file)
                    except IOError:
                        monitoring.to_log_and_console("       error in image '" + segmentation_file + "'", 1)
                        restart=t
            #
            #
            #

            if restart == -1:
                monitoring.to_log_and_console("       time '" + str(t) + "' seems ok", 1)
            t+=1
        first_time_point = restart
        monitoring.to_log_and_console(proc + ": restart computation at time '" + str(first_time_point) + "'", 1)
    else:
        monitoring.to_log_and_console(proc + ": start computation at time '" + str(first_time_point) + "'", 1)

    #
    #
    #

    for current_time in range(first_time_point + experiment.delay_time_point + experiment.delta_time_point,
                            last_time_point + experiment.delay_time_point + 1, experiment.delta_time_point):

        acquisition_time = str('{:0{width}d}'.format(current_time, width=default_width))
        previous_time = current_time - experiment.delta_time_point

        #
        # start processing
        #

        monitoring.to_log_and_console('... astec processing of time ' + acquisition_time, 1)
        start_time = time.time()

        #
        # set and make temporary directory
        #

        environment.temporary_path = os.path.join(environment.path_seg_exp, "TEMP_$TIME")
        environment.temporary_path = environment.temporary_path.replace(nomenclature.FLAG_TIME, acquisition_time)

        if not os.path.isdir(environment.temporary_path):
            os.makedirs(environment.temporary_path)

        if parameters.keep_reconstruction is True:
            environment.path_reconstruction = os.path.join(environment.path_seg_exp, "RECONSTRUCTION")
            if not os.path.isdir(environment.path_reconstruction):
                os.makedirs(environment.path_reconstruction)
        else:
            environment.path_reconstruction = environment.temporary_path

        #
        # process
        #

        astec_process(previous_time, current_time, lineage_tree_information, experiment, environment, parameters)

        #
        # cleaning
        #

        if monitoring.keepTemporaryFiles is False:
            shutil.rmtree(environment.temporary_path)

        #
        # end processing for a time point
        #
        end_time = time.time()

        monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
        monitoring.to_log_and_console('', 1)

    return