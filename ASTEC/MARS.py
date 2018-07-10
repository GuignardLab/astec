
import os
import imp
import sys
import time
import shutil

import ACE
import commonTools
import nomenclature
import CommunFunctions.cpp_wrapping as cpp_wrapping

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


class MarsEnvironment(object):

    def __init__(self):

        #
        # fused data paths
        #
        self.path_fuse_exp = None
        self.path_fuse_exp_files = None

        #
        # mars data paths
        #
        self.path_seg_exp = None
        self.path_mars_exp_files = None

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
        self.path_mars_exp_files = nomenclature.replaceFlags(nomenclature.path_mars_exp_files, parameters)

        self.path_logdir = nomenclature.replaceFlags(nomenclature.path_seg_logdir, parameters)
        self.path_history_file = nomenclature.replaceFlags(nomenclature.path_seg_historyfile, parameters)
        self.path_log_file = nomenclature.replaceFlags(nomenclature.path_seg_logfile, parameters, start_time)

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('MarsEnvironment\n')

            logfile.write('- path_fuse_exp = ' + str(self.path_fuse_exp) + '\n')
            logfile.write('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files) + '\n')

            logfile.write('- path_seg_exp = ' + str(self.path_seg_exp) + '\n')
            logfile.write('- path_mars_exp_files = ' + str(self.path_mars_exp_files) + '\n')

            logfile.write('- path_reconstruction = ' + str(self.path_reconstruction) + '\n')
            logfile.write('- temporary_path = ' + str(self.temporary_path) + '\n')

            logfile.write('- path_logdir = ' + str(self.path_logdir) + '\n')
            logfile.write('- path_history_file = ' + str(self.path_history_file)+'\n')
            logfile.write('- path_log_file = ' + str(self.path_log_file)+'\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('MarsEnvironment')

        print('- path_fuse_exp = ' + str(self.path_fuse_exp))
        print('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files))

        print('- path_seg_exp = ' + str(self.path_seg_exp))
        print('- path_mars_exp_files = ' + str(self.path_mars_exp_files))

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


class MarsParameters(object):

    def __init__(self):
        #
        #
        #
        self.first_time_point = -1
        self.last_time_point = -1

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
        # watershed parameters
        #
        self.watershed_seed_sigma = 0.6
        self.watershed_membrane_sigma = 0.15
        self.watershed_seed_hmin = 4

        #
        # images suffixes/formats
        #
        self.result_image_suffix = 'inr'
        self.default_image_suffix = 'inr'

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('MarsParameters\n')

            logfile.write('- first_time_point = ' + str(self.first_time_point) + '\n')
            logfile.write('- last_time_point = ' + str(self.last_time_point) + '\n')

            logfile.write('- intensity_transformation = ' + str(self.intensity_transformation) + '\n')
            logfile.write('- intensity_enhancement = ' + str(self.intensity_enhancement) + '\n')

            self.ace.write_parameters(log_file_name)

            logfile.write('- keep_reconstruction = ' + str(self.keep_reconstruction) + '\n')

            logfile.write('- watershed_seed_sigma = ' + str(self.watershed_seed_sigma) + '\n')
            logfile.write('- watershed_membrane_sigma = ' + str(self.watershed_membrane_sigma) + '\n')
            logfile.write('- watershed_seed_hmin = ' + str(self.watershed_seed_hmin) + '\n')

            logfile.write('- result_image_suffix = ' + str(self.result_image_suffix) + '\n')
            logfile.write('- default_image_suffix = '+str(self.default_image_suffix) + '\n')

            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('MarsParameters')

        print('- first_time_point = ' + str(self.first_time_point))
        print('- last_time_point = ' + str(self.last_time_point))

        print('- intensity_transformation = ' + str(self.intensity_transformation))
        print('- intensity_enhancement = ' + str(self.intensity_enhancement))

        self.ace.print_parameters()

        print('- keep_reconstruction = ' + str(self.keep_reconstruction))

        print('- watershed_seed_sigma = ' + str(self.watershed_seed_sigma))
        print('- watershed_membrane_sigma = ' + str(self.watershed_membrane_sigma))
        print('- watershed_seed_hmin = ' + str(self.watershed_seed_hmin))

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
        #
        #
        if hasattr(parameters, 'mars_begin'):
            self.first_time_point = parameters.mars_begin
        if hasattr(parameters, 'mars_end'):
            self.last_time_point = parameters.mars_end

        #
        # reconstruction methods
        #

        if hasattr(parameters, 'mars_method'):
            if parameters.mars_method == 1:
                self.intensity_transformation = 'Identity'
                self.intensity_enhancement = None
            elif parameters.mars_method == 2:
                self.intensity_transformation = None
                self.intensity_enhancement = 'GACE'

        if hasattr(parameters, 'intensity_transformation'):
            self.intensity_transformation = parameters.intensity_transformation
        if hasattr(parameters, 'mars_intensity_transformation'):
            self.intensity_transformation = parameters.mars_intensity_transformation

        if hasattr(parameters, 'intensity_enhancement'):
            self.intensity_enhancement = parameters.intensity_enhancement
        if hasattr(parameters, 'mars_intensity_enhancement'):
            self.intensity_enhancement = parameters.mars_intensity_enhancement

        #
        #
        #
        self.ace.update_from_file(parameter_file)

        #
        #
        #
        if hasattr(parameters, 'mars_keep_reconstruction'):
            if parameters.mars_keep_reconstruction is not None:
                self.keep_reconstruction = parameters.mars_keep_reconstruction
        if hasattr(parameters, 'keep_reconstruction'):
            if parameters.keep_reconstruction is not None:
                self.keep_reconstruction = parameters.keep_reconstruction

        #
        # watershed parameters
        #
        if hasattr(parameters, 'mars_sigma1'):
            if parameters.mars_sigma1 is not None:
                self.watershed_seed_sigma = parameters.mars_sigma1
        if hasattr(parameters, 'watershed_seed_sigma'):
            if parameters.watershed_seed_sigma is not None:
                self.watershed_seed_sigma = parameters.watershed_seed_sigma

        if hasattr(parameters, 'mars_sigma2'):
            if parameters.mars_sigma2 is not None:
                self.watershed_membrane_sigma = parameters.mars_sigma2
        if hasattr(parameters, 'watershed_membrane_sigma'):
            if parameters.watershed_membrane_sigma is not None:
                self.watershed_membrane_sigma = parameters.watershed_membrane_sigma

        if hasattr(parameters, 'mars_h_min'):
            if parameters.mars_h_min is not None:
                self.watershed_seed_hmin = parameters.mars_h_min
        if hasattr(parameters, 'watershed_seed_hmin'):
            if parameters.watershed_seed_hmin is not None:
                self.watershed_seed_hmin = parameters.watershed_seed_hmin

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

def _build_membrane_image(input_image, environment, parameters):
    """

    :param input_image:
    :param environment:
    :param parameters:
    :return:
    """
    #
    # build input image(s) for segmentation
    # 0. build names
    # 1. do image enhancement if required
    # 2. transform input image if required
    # 3. mix the two results
    #
    # If any transformation is performed on the input image, the final result image (input of the watershed)
    # will be named (in the 'temporary_path' directory)
    # - 'input_image'_membrane
    # if fusion is required
    # - 'input_image'_enhanced
    #   is the membrance-enhanced image (by tensor voting method)
    # - 'input_image'_intensity
    #   is the intensity-transformed image (by normalization)
    #

    enhanced_image = None
    intensity_image = None
    membrane_image = None

    if parameters.intensity_enhancement is None or parameters.intensity_enhancement == 'None':
        #
        # no membrane enhancement image
        # only set intensity_image
        #
        if parameters.intensity_transformation is None or parameters.intensity_transformation == 'None' \
                or parameters.intensity_transformation == 'Identity':
            #
            # nothing to do
            #
            intensity_image = input_image
        elif parameters.intensity_transformation == 'Normalization_to_u8':
            intensity_image = commonTools.add_suffix(input_image, "_membrane",
                                                     new_dirname=environment.path_reconstruction,
                                                     new_extension=parameters.default_image_suffix)
        else:
            monitoring.to_log_and_console("    unknwon intensity transformation method: '"
                                          + str(parameters.intensity_transformation) + "'", 2)
    elif parameters.intensity_enhancement == 'GACE':
        #
        # only set enhanced_image
        # or set the 3 names: intensity_image, enhanced_image, and membrane_image
        #
        if parameters.intensity_transformation is None or parameters.intensity_transformation == 'None':
            enhanced_image = commonTools.add_suffix(input_image, "_membrane",
                                                    new_dirname=environment.path_reconstruction,
                                                    new_extension=parameters.default_image_suffix)
        elif parameters.intensity_transformation == 'Identity':
            intensity_image = input_image
            enhanced_image = commonTools.add_suffix(input_image, "_enhanced", new_dirname=environment.temporary_path,
                                                    new_extension=parameters.default_image_suffix)
            membrane_image = commonTools.add_suffix(input_image, "_membrane",
                                                    new_dirname=environment.path_reconstruction,
                                                     new_extension=parameters.default_image_suffix)
        elif parameters.intensity_transformation == 'Normalization_to_u8':
            intensity_image = commonTools.add_suffix(input_image, "_intensity", new_dirname=environment.temporary_path,
                                                     new_extension=parameters.default_image_suffix)
            enhanced_image = commonTools.add_suffix(input_image, "_enhanced", new_dirname=environment.temporary_path,
                                                    new_extension=parameters.default_image_suffix)
            membrane_image = commonTools.add_suffix(input_image, "_membrane",
                                                    new_dirname=environment.path_reconstruction,
                                                     new_extension=parameters.default_image_suffix)
        else:
            monitoring.to_log_and_console("    unknwon intensity transformation method: '"
                                          + str(parameters.intensity_transformation) + "'", 2)
    else:
        monitoring.to_log_and_console("    unknwon membrane enhancement method: '"
                                      + str(parameters.intensity_enhancement) + "'", 2)

    #
    #
    #

    if parameters.intensity_enhancement is None or parameters.intensity_enhancement == 'None':
        pass
    elif parameters.intensity_enhancement == 'GACE':
        monitoring.to_log_and_console("    .. enhance membranes of '"
                                      + str(input_image).split(os.path.sep)[-1] + "'", 2)
        if not os.path.isfile(enhanced_image):
            ACE.monitoring.copy(monitoring)
            ACE.global_membrane_enhancement(input_image, enhanced_image, temporary_path=environment.temporary_path,
                                            parameters=parameters.ace)
    else:
        monitoring.to_log_and_console("    unknwon membrane enhancement method: '"
                                  + str(parameters.intensity_enhancement) + "'", 2)

    arit_options = "-o 2"

    if parameters.intensity_transformation is None or parameters.intensity_transformation == 'None':
        pass
    elif parameters.intensity_transformation == 'Identity':
        pass
    elif parameters.intensity_transformation == 'Normalization_to_u8':
        monitoring.to_log_and_console("    .. intensity normalization '"
                                      + str(input_image).split(os.path.sep)[-1] + "'", 2)
        if not os.path.isfile(intensity_image):
            cpp_wrapping.inline_to_u8(input_image, intensity_image, min_percentile=0.01, max_percentile=0.99,
                                      monitoring=monitoring)
        arit_options = "-o 1"
    else:
        monitoring.to_log_and_console("    unknwon intensity transformation method: '"
                                      + str(parameters.intensity_transformation) + "'", 2)

    #
    # fuse the two images (if fusion is required)
    #

    if enhanced_image is None:
        return intensity_image
    else:
        if intensity_image is None:
            return enhanced_image
        else:
            #
            # mix images
            #
            monitoring.to_log_and_console("       fusion of intensity and enhancement", 2)
            arit_options += " -max"
            if not os.path.isfile(membrane_image):
                cpp_wrapping.arithmetic(intensity_image, enhanced_image, membrane_image, other_options=arit_options)
            return membrane_image

    #
    # should not reach this point
    #
    return None

#
#
#
#
#

def mars_process(input_image, mars_image, environment, parameters):
    """

    :param input_image:
    :param mars_image:
    :param environment:
    :param parameters:
    :return:
    """

    proc = "mars_process"

    if monitoring.debug > 2:
        print ""
        print proc + " was called with:"
        print "- input_image = " + str(input_image)
        print "- mars_image = " + str(mars_image)
        print ""

    #
    # nothing to do if the fused image exists
    #

    if os.path.isfile(os.path.join(environment.path_seg_exp, mars_image)):
        if monitoring.forceResultsToBeBuilt is False:
            monitoring.to_log_and_console('    mars image already existing', 2)
            return
        else:
            monitoring.to_log_and_console('    mars image already existing, but forced', 2)

    #
    # build the membrane image to be segmented
    #

    membrane_image = _build_membrane_image(input_image, environment, parameters)

    if membrane_image is None or not os.path.isfile(membrane_image):
        monitoring.to_log_and_console("       '" + str(membrane_image).split(os.path.sep)[-1]
                                      + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    #
    # computation of seed image
    # - [smoothing]
    # - regional minima
    # - hysteresis thresholding
    #

    monitoring.to_log_and_console("    .. seed extraction '" + str(membrane_image).split(os.path.sep)[-1] + "'", 2)

    if parameters.watershed_seed_sigma > 0.0:
        monitoring.to_log_and_console("       smoothing '" + str(membrane_image).split(os.path.sep)[-1] + "'", 2)
        seed_preimage = commonTools.add_suffix(input_image, "_seed_preimage", new_dirname=environment.temporary_path,
                                               new_extension=parameters.default_image_suffix)
        if not os.path.isfile(seed_preimage) or monitoring.forceResultsToBeBuilt is True:
            cpp_wrapping.linear_smoothing(membrane_image, seed_preimage, parameters.watershed_seed_sigma,
                                          real_scale=True, other_options="-o 2", monitoring=monitoring)
    else:
        seed_preimage = membrane_image

    if not os.path.isfile(seed_preimage):
        monitoring.to_log_and_console("       '" + str(seed_preimage).split(os.path.sep)[-1] + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    monitoring.to_log_and_console("       extract regional minima '"
                                  + str(seed_preimage).split(os.path.sep)[-1] + "'", 2)

    minima_image = commonTools.add_suffix(input_image, "_minima", new_dirname=environment.temporary_path,
                                          new_extension=parameters.default_image_suffix)
    if not os.path.isfile(minima_image) or monitoring.forceResultsToBeBuilt is True:
        cpp_wrapping.regional_minima(seed_preimage, minima_image, h_min=parameters.watershed_seed_hmin,
                                     monitoring=monitoring)

    if not os.path.isfile(minima_image):
        monitoring.to_log_and_console("       '" + str(minima_image).split(os.path.sep)[-1] + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    monitoring.to_log_and_console("       label regional minima '"
                                  + str(minima_image).split(os.path.sep)[-1] + "'", 2)

    seed_image = commonTools.add_suffix(input_image, "_seeds", new_dirname=environment.temporary_path,
                                        new_extension=parameters.default_image_suffix)
    if not os.path.isfile(seed_image) or monitoring.forceResultsToBeBuilt is True:
        cpp_wrapping.connected_components(minima_image, seed_image, high_threshold=parameters.watershed_seed_hmin,
                                          monitoring=monitoring)

    #
    # computation of height image for watershed
    #

    if parameters.watershed_membrane_sigma > 0.0:
        monitoring.to_log_and_console("    .. smoothing '" + str(membrane_image).split(os.path.sep)[-1] + "'", 2)
        height_image = commonTools.add_suffix(input_image, "_height", new_dirname=environment.temporary_path,
                                              new_extension=parameters.default_image_suffix)
        if not os.path.isfile(height_image) or monitoring.forceResultsToBeBuilt is True:
            cpp_wrapping.linear_smoothing(membrane_image, height_image, parameters.watershed_membrane_sigma,
                                          real_scale=True, other_options="-o 2", monitoring=monitoring)
    else:
        height_image = membrane_image

    #
    # watershed
    #
    if not os.path.isfile(seed_image):
        monitoring.to_log_and_console("       '" + str(seed_image).split(os.path.sep)[-1] + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)
    if not os.path.isfile(height_image):
        monitoring.to_log_and_console("       '" + str(height_image).split(os.path.sep)[-1] + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    monitoring.to_log_and_console("    .. watershed '" + str(height_image).split(os.path.sep)[-1] + "'", 2)

    if not os.path.isfile(mars_image) or monitoring.forceResultsToBeBuilt is True:
        cpp_wrapping.watershed(seed_image, height_image, mars_image, monitoring=monitoring)

    return


#
#
#
#
#


def mars_control(experiment, environment, parameters):
    """

    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    default_width = 3

    #
    # make sure that the result directory exists
    #

    if not os.path.isdir(environment.path_seg_exp):
        os.makedirs(environment.path_seg_exp)

    if not os.path.isdir(environment.path_logdir):
        os.makedirs(environment.path_logdir)

    monitoring.to_log_and_console('', 1)

    #
    #
    #

    if parameters.first_time_point < 0 or parameters.last_time_point < 0:
        monitoring.to_log_and_console("... time interval does not seem to be defined in the parameter file")
        monitoring.to_log_and_console("    set parameters 'begin' and 'end'")
        monitoring.to_log_and_console("\t Exiting")
        sys.exit(1)

    if parameters.first_time_point > parameters.last_time_point:
        monitoring.to_log_and_console("... weird time interval: 'begin' = " + str(parameters.first_time_point)
                                      + ", 'end' = " + str(parameters.last_time_point))

    for time_value in range(parameters.first_time_point, parameters.last_time_point + 1, experiment.deltaTimePoint):

        acquisition_time = str('{:0{width}d}'.format(time_value, width=default_width))

        #
        # start processing
        #
        monitoring.to_log_and_console('... mars processing of time ' + acquisition_time, 1)
        start_time = time.time()

        #
        #
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
        # mars image name
        #

        mars_image = nomenclature.replaceTIME(environment.path_mars_exp_files,
                                              time_value + experiment.delayTimePoint) \
                     + '.' + parameters.result_image_suffix
        #
        #
        #
        name = nomenclature.replaceTIME(environment.path_fuse_exp_files, time_value + experiment.delayTimePoint)
        image = commonTools.find_file(environment.path_fuse_exp, name, monitoring)
        if image is None:
            monitoring.to_log_and_console("    .. image '" + name + "' not found in '"
                                          + str(environment.path_fuse_exp) + "'", 2)
            monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
            continue

        #
        #
        #
        image = os.path.join(environment.path_fuse_exp, image)
        mars_image = os.path.join(environment.path_seg_exp, mars_image)
        mars_process(image, mars_image, environment, parameters)

        if monitoring.keepTemporaryFiles is False:
            shutil.rmtree(environment.temporary_path)

        #
        # end processing for a time point
        #
        end_time = time.time()
        monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
        monitoring.to_log_and_console('', 1)

    return
