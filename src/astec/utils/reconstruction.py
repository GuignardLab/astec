import os
import sys
import shutil

import astec.utils.common as common
import astec.utils.ace as ace
import astec.wrapping.cpp_wrapping as cpp_wrapping

monitoring = common.Monitoring()


########################################################################################
#
#
#
########################################################################################


class ReconstructionParameters(ace.AceParameters):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix=None, suffix=None):

        if "doc" not in self.__dict__:
            self.doc = {}

        doc = "\n"
        doc += "Reconstruction parameters overview:\n"
        doc += "===================================\n"
        doc += "The reconstruction step serves at building an image more\n"
        doc += "suitable for the considered processing.\n"
        doc += "It may be built from 3 images, each of them issued from\n"
        doc += "the input image, after a prenormalization step.\n"
        doc += "1. an intensity-transformed image\n"
        doc += "2. a membrane-enhanced image, through g(l)ace procedure\n"
        doc += "3. am outer-membrane image\n"
        doc += "\n"
        self.doc['reconstruction_overview'] = doc

        ace.AceParameters.__init__(self, prefix=prefix)

        #
        #
        #
        doc = "\t Possible values are 'identity', 'normalization_to_u8', \n"
        doc += "\t or 'normalization_to_u16'\n"
        doc += "\t Performs a global robust normalization of the input image, the\n"
        doc += "\t intensity value corresponding to the min percentile is set\n"
        doc += "\t to 0, while the intensity value corresponding to the max\n"
        doc += "\t percentile is set either to 255 (u8) or 2047 (u16). In-between\n"
        doc += "\t values are linearly interpolated.\n"
        doc += "\t Should be left to 'identity' for integer-encoded images.\n"
        doc += "\t It has been introduced for real-encoded images.\n"
        self.doc['intensity_prenormalization'] = doc
        self.intensity_prenormalization = 'Identity'

        doc = "\t Percentile of the image histogram used to determine the value\n"
        doc += "\t  to be set to 0 (prenormalization step)\n"
        self.doc['prenormalization_min_percentile'] = doc
        self.prenormalization_min_percentile = 0.00

        doc = "\t Percentile of the image histogram used to determine the value\n"
        doc += "\t  to be set to either 255 or 2047 (prenormalization step)\n"
        self.doc['prenormalization_max_percentile'] = doc
        self.prenormalization_max_percentile = 1.00

        doc = "\t Intensity-transformation modes are 'identity', 'normalization_to_u8',\n"
        doc += "\t 'normalization_to_u16', 'cell_normalization_to_u8' or None\n"
        doc += "\t - 'identity': intensities of the input image are left unchanged\n"
        doc += "\t - 'normalization_to_u8': intensities of the input image are\n"
        doc += "\t    globally normalized into [0,255] \n"
        doc += "\t - 'normalization_to_u16': intensities of the input image are\n"
        doc += "\t    globally normalized into [0,2047] \n"
        doc += "\t - 'cell_normalization_to_u8': intensities of the input image are\n"
        doc += "\t    normalized into [0,255] on a cell-based manner, those cells \n"
        doc += "\t    being the cells at t-1 deformed on t. Works only in the \n"
        doc += "\t    propagation segmentation stage. Fragile\n"
        self.doc['intensity_transformation'] = doc
        self.intensity_transformation = 'Identity'

        doc = "\t Possible values are 'gace', 'glace', or None\n"
        self.doc['intensity_enhancement'] = doc
        self.intensity_enhancement = None

        doc = "\t Possible value are True or False.\n"
        doc += "\t Fragile. Kept for test purposes.\n"
        self.doc['outer_contour_enhancement'] = doc
        self.outer_contour_enhancement = False

        doc = "\t Possible values are 'addition' or 'maximum'\n"
        doc += "\t Set the fusion mode when several images are to be fused\n"
        doc += "\t 'maximum' should be used when all images are encoded\n"
        doc += "\t on one byte.\n"
        self.doc['reconstruction_images_combination'] = doc
        self.reconstruction_images_combination = "maximum"

        doc = "\t Possible values are 'cell', 'cellborder', 'cellinterior'\n"
        doc += "\t Set where to compute the histogram to get the minimum \n"
        doc += "\t value for the 'cell_normalization_to_u8' procedure\n"
        self.doc['cell_normalization_min_method'] = doc
        self.cell_normalization_min_method = 'cellinterior'

        doc = "\t Possible values are 'cell', 'cellborder', 'cellinterior'\n"
        doc += "\t Set where to compute the histogram to get the maximum \n"
        doc += "\t value for the 'cell_normalization_to_u8' procedure\n"
        self.doc['cell_normalization_max_method'] = doc
        self.cell_normalization_max_method = 'cellborder'

        doc = "\t Percentile of the image histogram used to determine the value\n"
        doc += "\t  to be set to 0 (normalization step)\n"
        self.doc['normalization_min_percentile'] = doc
        self.normalization_min_percentile = 0.01

        doc = "\t Percentile of the image histogram used to determine the value\n"
        doc += "\t  to be set to 0 (normalization step)\n"
        self.doc['normalization_max_percentile'] = doc
        self.normalization_max_percentile = 0.99

        doc = "\t Sigma of the gaussian used to smooth the minimum and maximum\n"
        doc += "\t images obtained through the 'cell_normalization_to_u8'\n"
        doc += "\t procedure\n"
        doc += "\t \n"
        doc += "\t \n"
        self.doc['cell_normalization_sigma'] = doc
        self.cell_normalization_sigma = 5.0

        doc = "\t Sigma (in real units) of the smoothing gaussian applied\n"
        doc += "\t to the intensity-transformed image, prior its fusion\n"
        self.doc['intensity_sigma'] = doc
        self.intensity_sigma = 0.0

        #
        # registration parameters
        #
        self.registration = []

        self.registration.append(common.RegistrationParameters(prefix=[self._prefix, 'linear_registration_']))
        self.registration[0].pyramid_highest_level = 5
        self.registration[0].pyramid_lowest_level = 3
        self.registration[0].gaussian_pyramid = True
        self.registration[0].transformation_type = 'affine'
        self.registration[0].transformation_estimation_type = 'wlts'
        self.registration[0].lts_fraction = 0.55

        self.registration.append(common.RegistrationParameters(prefix=[self._prefix, 'nonlinear_registration_']))
        self.registration[1].pyramid_highest_level = 5
        self.registration[1].pyramid_lowest_level = 3
        self.registration[1].gaussian_pyramid = True
        self.registration[1].transformation_type = 'vectorfield'
        self.registration[1].elastic_sigma = 2.0
        self.registration[1].transformation_estimation_type = 'wlts'
        self.registration[1].lts_fraction = 1.0
        self.registration[1].fluid_sigma = 2.0
        #
        #
        #
        doc = "\t Possible value are True or False.\n"
        doc += "\t If True, reconstructed images are kept.\n"
        self.doc['keep_reconstruction'] = doc
        self.keep_reconstruction = True

        #
        # the processing flow is as follows
        #
        # input image -> pre-normalized image -> normalized image -> intensity image }
        #                                     -> enhanced image                      } -> reconstructed image
        #                                     -> outer contour enhancement image     }
        #
        self._default_suffix = suffix
        self._prenormalized_suffix = None
        self._normalized_suffix = None
        self._intensity_suffix = None
        self._enhanced_suffix = None
        self._outercontour_suffix = None
        self._result_suffix = None

    def _set_default_suffix(self, suffix):
        self._default_suffix = suffix

    def _get_default_suffix(self):
        return self._default_suffix

    def _set_prenormalized_suffix(self, suffix):
        self._prenormalized_suffix = suffix

    def _get_prenormalized_suffix(self):
        return self._prenormalized_suffix

    def _set_normalized_suffix(self, suffix):
        self._normalized_suffix = suffix

    def _get_normalized_suffix(self):
        return self._normalized_suffix

    def _set_intensity_suffix(self, suffix):
        self._intensity_suffix = suffix

    def _get_intensity_suffix(self):
        return self._intensity_suffix

    def _set_enhanced_suffix(self, suffix):
        self._enhanced_suffix = suffix

    def _get_enhanced_suffix(self):
        return self._enhanced_suffix

    def _set_outercontour_suffix(self, suffix):
        self._outercontour_suffix = suffix

    def _get_outercontour_suffix(self):
        return self._outercontour_suffix

    def _set_result_suffix(self, suffix):
        self._result_suffix = suffix

    def _get_result_suffix(self):
        return self._result_suffix

    default_suffix = property(_get_default_suffix, _set_default_suffix)
    prenormalized_suffix = property(_get_prenormalized_suffix, _set_prenormalized_suffix)
    normalized_suffix = property(_get_normalized_suffix, _set_normalized_suffix)
    intensity_suffix = property(_get_intensity_suffix, _set_intensity_suffix)
    enhanced_suffix = property(_get_enhanced_suffix, _set_enhanced_suffix)
    outercontour_suffix = property(_get_outercontour_suffix, _set_outercontour_suffix)
    result_suffix = property(_get_result_suffix, _set_result_suffix)

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('#')
        print('# ReconstructionParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        for line in self.doc['reconstruction_overview'].splitlines():
            print('# ' + line)

        ace.AceParameters.print_parameters(self)
        self.varprint('intensity_prenormalization', self.intensity_prenormalization,
                      self.doc['intensity_prenormalization'])
        self.varprint('prenormalization_min_percentile', self.prenormalization_min_percentile,
                      self.doc['prenormalization_min_percentile'])
        self.varprint('prenormalization_max_percentile', self.prenormalization_max_percentile,
                      self.doc['prenormalization_max_percentile'])

        self.varprint('intensity_transformation', self.intensity_transformation, self.doc['intensity_transformation'])
        self.varprint('intensity_enhancement', self.intensity_enhancement, self.doc['intensity_enhancement'])
        self.varprint('outer_contour_enhancement', self.outer_contour_enhancement,
                      self.doc['outer_contour_enhancement'])
        self.varprint('reconstruction_images_combination', self.reconstruction_images_combination,
                      self.doc['reconstruction_images_combination'])

        self.varprint('cell_normalization_min_method', self.cell_normalization_min_method,
                      self.doc['cell_normalization_min_method'])
        self.varprint('cell_normalization_max_method', self.cell_normalization_max_method,
                      self.doc['cell_normalization_max_method'])

        self.varprint('normalization_min_percentile', self.normalization_min_percentile,
                      self.doc['normalization_min_percentile'])
        self.varprint('normalization_max_percentile', self.normalization_max_percentile,
                      self.doc['normalization_max_percentile'])
        self.varprint('cell_normalization_sigma', self.cell_normalization_sigma,
                      self.doc['cell_normalization_sigma'])
        self.varprint('intensity_sigma', self.intensity_sigma, self.doc['intensity_sigma'])

        for p in self.registration:
            p.print_parameters()

        self.varprint('keep_reconstruction', self.keep_reconstruction, self.doc['keep_reconstruction'])

        self.varprint('_default_suffix', str(self._default_suffix))
        self.varprint('_prenormalized_suffix', str(self._prenormalized_suffix))
        self.varprint('_normalized_suffix', str(self._normalized_suffix))
        self.varprint('_intensity_suffix', str(self._intensity_suffix))
        self.varprint('_enhanced_suffix', str(self._enhanced_suffix))
        self.varprint('_outercontour_suffix', str(self._outercontour_suffix))
        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# ReconstructionParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        for line in self.doc['reconstruction_overview'].splitlines():
            logfile.write('# ' + line + '\n')

        ace.AceParameters.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'intensity_prenormalization', self.intensity_prenormalization,
                      self.doc['intensity_prenormalization'])
        self.varwrite(logfile, 'prenormalization_min_percentile', self.prenormalization_min_percentile,
                      self.doc['prenormalization_min_percentile'])
        self.varwrite(logfile, 'prenormalization_max_percentile', self.prenormalization_max_percentile,
                      self.doc['prenormalization_max_percentile'])

        self.varwrite(logfile, 'intensity_transformation', self.intensity_transformation,
                      self.doc['intensity_transformation'])
        self.varwrite(logfile, 'intensity_enhancement', self.intensity_enhancement,
                      self.doc['intensity_enhancement'])
        self.varwrite(logfile, 'outer_contour_enhancement', self.outer_contour_enhancement,
                      self.doc['outer_contour_enhancement'])
        self.varwrite(logfile, 'reconstruction_images_combination', self.reconstruction_images_combination,
                      self.doc['reconstruction_images_combination'])

        self.varwrite(logfile, 'cell_normalization_min_method', self.cell_normalization_min_method,
                      self.doc['cell_normalization_min_method'])
        self.varwrite(logfile, 'cell_normalization_max_method', self.cell_normalization_max_method,
                      self.doc['cell_normalization_max_method'])

        self.varwrite(logfile, 'normalization_min_percentile', self.normalization_min_percentile,
                      self.doc['normalization_min_percentile'])
        self.varwrite(logfile, 'normalization_max_percentile', self.normalization_max_percentile,
                      self.doc['normalization_max_percentile'])
        self.varwrite(logfile, 'cell_normalization_sigma', self.cell_normalization_sigma,
                      self.doc['cell_normalization_sigma'])
        self.varwrite(logfile, 'intensity_sigma', self.intensity_sigma, self.doc['intensity_sigma'])

        for p in self.registration:
            p.write_parameters_in_file(logfile)

        self.varwrite(logfile, 'keep_reconstruction', self.keep_reconstruction, self.doc['keep_reconstruction'])

        self.varwrite(logfile, '_default_suffix', str(self._default_suffix))
        self.varwrite(logfile, '_prenormalized_suffix', str(self._prenormalized_suffix))
        self.varwrite(logfile, '_normalized_suffix', str(self._normalized_suffix))
        self.varwrite(logfile, '_intensity_suffix', str(self._intensity_suffix))
        self.varwrite(logfile, '_enhanced_suffix', str(self._enhanced_suffix))
        self.varwrite(logfile, '_outercontour_suffix', str(self._outercontour_suffix))

        logfile.write("\n")
        return

    def write_parameters(self, log_filename=None):
        if log_filename is not None:
            local_log_filename = log_filename
        else:
            local_log_filename = monitoring.log_filename
        if local_log_filename is not None:
            with open(local_log_filename, 'a') as logfile:
                self.write_parameters_in_file(logfile)
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameters):
        ace.AceParameters.update_from_parameters(self, parameters)

        self.intensity_prenormalization = self.read_parameter(parameters, 'intensity_prenormalization',
                                                              self.intensity_prenormalization)
        self.prenormalization_min_percentile = self.read_parameter(parameters, 'prenormalization_min_percentile',
                                                                   self.prenormalization_min_percentile)
        self.prenormalization_max_percentile = self.read_parameter(parameters, 'prenormalization_max_percentile',
                                                                   self.prenormalization_max_percentile)
        #
        # reconstruction methods
        #

        self.intensity_transformation = self.read_parameter(parameters, 'intensity_transformation',
                                                            self.intensity_transformation)
        self.intensity_enhancement = self.read_parameter(parameters, 'intensity_enhancement',
                                                         self.intensity_enhancement)
        self.outer_contour_enhancement = self.read_parameter(parameters, 'outer_contour_enhancement',
                                                             self.outer_contour_enhancement)
        self.reconstruction_images_combination = self.read_parameter(parameters, 'reconstruction_images_combination',
                                                                     self.reconstruction_images_combination)

        #
        # cell-based parameters
        #
        self.cell_normalization_min_method = self.read_parameter(parameters, 'cell_normalization_min_method',
                                                                 self.cell_normalization_min_method)
        self.cell_normalization_max_method = self.read_parameter(parameters, 'cell_normalization_max_method',
                                                                 self.cell_normalization_max_method)
        self.normalization_min_percentile = self.read_parameter(parameters, 'normalization_min_percentile',
                                                                self.normalization_min_percentile)
        self.normalization_max_percentile = self.read_parameter(parameters, 'normalization_max_percentile',
                                                                self.normalization_max_percentile)
        self.cell_normalization_sigma = self.read_parameter(parameters, 'cell_normalization_sigma',
                                                            self.cell_normalization_sigma)
        self.intensity_sigma = self.read_parameter(parameters, 'intensity_sigma', self.intensity_sigma)

        #
        #
        #
        self.registration[0].update_from_parameters(parameters)
        self.registration[1].update_from_parameters(parameters)

        #
        #
        #
        self.keep_reconstruction = self.read_parameter(parameters, 'keep_reconstruction', self.keep_reconstruction)

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
    #
    #
    ############################################################

    def is_prenormalization_equal(self, p):
        if self.intensity_prenormalization != p.intensity_prenormalization:
            return False
        if self.intensity_prenormalization == p.intensity_prenormalization:
            if self.intensity_prenormalization.lower() == 'identity':
                return True
            if self.prenormalization_min_percentile == p.prenormalization_min_percentile \
                    and self.prenormalization_max_percentile == p.prenormalization_max_percentile:
                return True
            return False
        return False

    def is_normalization_equal(self, p):
        if self.intensity_transformation != p.intensity_transformation:
            return False
        if self.intensity_transformation == p.intensity_transformation:
            if self.intensity_transformation.lower() == 'identity':
                return True
            if self.intensity_transformation.lower() == 'normalization_to_u8' \
                    or self.intensity_transformation.lower() == 'normalization_to_u16':
                if self.prenormalization_min_percentile == p.prenormalization_min_percentile \
                        and self.prenormalization_max_percentile == p.prenormalization_max_percentile:
                    return True
                return False
            if self.intensity_transformation.lower() == 'cell_normalization_to_u8':
                if self.cell_normalization_min_method == p.cell_normalization_min_method \
                        and self.cell_normalization_max_method == p.cell_normalization_max_method \
                        and self.cell_normalization_sigma == p.cell_normalization_sigma \
                        and self.prenormalization_min_percentile == p.prenormalization_min_percentile \
                        and self.prenormalization_max_percentile == p.prenormalization_max_percentile:
                    return True
                return False
        return False

    def is_intensity_smoothing_equal(self, p):
        if self.intensity_sigma == p.intensity_sigma:
            return True
        return False

    def is_intensity_enhancement_equal(self, p):
        if self.intensity_enhancement == p.intensity_enhancement:
            if self.intensity_enhancement is None or self.intensity_enhancement.lower() == 'none':
                return True
            return ace.AceParameters.is_equal(self, p)
        return False

    def is_outer_contour_enhancement_equal(self, p):
        if self.outer_contour_enhancement == p.outer_contour_enhancement:
            return True
        return False

    def is_equal(self, p, debug=False):
        proc = "is_equal"
        if not self.is_prenormalization_equal(p):
            if debug:
                print(proc + ": not equal at 'prenormalization'")
            return False
        if not self.is_normalization_equal(p):
            if debug:
                print(proc + ": not equal at 'normalization'")
            return False
        if not self.is_intensity_smoothing_equal(p):
            if debug:
                print(proc + ": not equal at 'intensity smoothing'")
            return False
        if not self.is_intensity_enhancement_equal(p):
            if debug:
                print(proc + ": not equal at 'intensity enhancement'")
            return False
        if not self.is_outer_contour_enhancement_equal(p):
            if debug:
                print(proc + ": not equal at 'outer contour enhancement'")
            return False
        #
        # registration parameters are not tested
        #
        return True

    ############################################################
    #
    #
    #
    ############################################################
    def set_suffixes(self, reference=None):

        if reference is None:
            if self.prenormalized_suffix is None:
                self.prenormalized_suffix = self.default_suffix
            if self.normalized_suffix is None:
                self.normalized_suffix = self.default_suffix
            if self.intensity_suffix is None:
                self.intensity_suffix = self.default_suffix
            if self.enhanced_suffix is None:
                self.enhanced_suffix = self.default_suffix
            if self.outercontour_suffix is None:
                self.outercontour_suffix = self.default_suffix
            if self.result_suffix is None:
                self.result_suffix = self.default_suffix
            return

        if self.is_prenormalization_equal(reference):
            self.prenormalized_suffix = reference.prenormalized_suffix
        else:
            self.set_suffixes()
            return

        if self.is_normalization_equal(reference):
            self.normalized_suffix = reference.normalized_suffix
            if self.is_intensity_smoothing_equal(reference):
                self.intensity_suffix = reference.intensity_suffix

        if self.is_intensity_enhancement_equal(reference):
            self.enhanced_suffix = reference.enhanced_suffix

        if self.is_outer_contour_enhancement_equal(reference):
            self.outercontour_suffix = reference.outercontour_suffix

        if self.is_equal(reference):
            self.result_suffix = reference.result_suffix

        self.set_suffixes()
        return


#
#
#
#
#

def get_deformation_from_current_to_previous(current_time, experiment, parameters, previous_time):
    """

    :param current_time:
    :param experiment:
    :param parameters:
    :param previous_time:
    :return:
    """

    proc = 'get_deformation_from_current_to_previous'

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, ReconstructionParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    #
    # image to be registered
    # floating image = image at previous time
    # reference image = image at current time
    #

    current_name = experiment.fusion_dir.get_image_name(current_time)
    current_image = common.find_file(experiment.fusion_dir.get_directory(), current_name, file_type='image',
                                     callfrom=proc, local_monitoring=None, verbose=False)
    if current_image is None:
        monitoring.to_log_and_console("    .. " + proc + " no fused image was found for time " + str(current_time), 2)
        return None
    current_image = os.path.join(experiment.fusion_dir.get_directory(), current_image)

    previous_name = experiment.fusion_dir.get_image_name(previous_time)
    previous_image = common.find_file(experiment.fusion_dir.get_directory(), previous_name, file_type='image',
                                      callfrom=proc, local_monitoring=None, verbose=False)
    if previous_image is None:
        monitoring.to_log_and_console("    .. " + proc + " no fused image was found for time " + str(previous_time), 2)
        return None
    previous_image = os.path.join(experiment.fusion_dir.get_directory(), previous_image)

    #
    # auxiliary file names
    #
    affine_image = common.add_suffix(previous_image, "_affine", new_dirname=experiment.working_dir.get_tmp_directory(0),
                                     new_extension=experiment.default_image_suffix)
    affine_trsf = common.add_suffix(previous_image, "_affine", new_dirname=experiment.working_dir.get_tmp_directory(0),
                                    new_extension="trsf")
    vector_image = common.add_suffix(previous_image, "_vector", new_dirname=experiment.working_dir.get_tmp_directory(0),
                                     new_extension=experiment.default_image_suffix)
    vector_trsf = common.add_suffix(previous_image, "_vector", new_dirname=experiment.working_dir.get_tmp_directory(0),
                                    new_extension="trsf")

    if os.path.isfile(vector_trsf) and not monitoring.forceResultsToBeBuilt:
        return vector_trsf

    #
    # linear registration
    #
    common.blockmatching(current_image, previous_image, affine_image, affine_trsf, None, parameters.registration[0])

    if not os.path.isfile(affine_image) or not os.path.isfile(affine_trsf):
        monitoring.to_log_and_console("    .. " + proc + " linear registration failed for time " + str(current_time), 2)
        return None

    #
    # non-linear registration
    #
    common.blockmatching(current_image, previous_image, vector_image, vector_trsf, affine_trsf,
                         parameters.registration[1])
    if not os.path.isfile(vector_image) or not os.path.isfile(vector_trsf):
        monitoring.to_log_and_console("    .. " + proc + " non-linear registration failed for time " +
                                      str(current_time), 2)
        return None

    return vector_trsf


#
#
#
#
#


def get_previous_deformed_segmentation(current_time, experiment, parameters, previous_time=None):
    """

    :param current_time:
    :param experiment:
    :param parameters:
    :param previous_time:
    :return:
    """

    #
    #
    #

    proc = 'get_previous_deformed_segmentation'

    #
    # it will generate (in the temporary_path directory)
    #
    # files related to the deformation computation
    # - $EN_fuse_t$TIME_affine.inr
    # - $EN_fuse_t$TIME_affine.trsf
    # - $EN_fuse_t$TIME_vector.inr
    # - $EN_fuse_t$TIME_vector.trsf
    #
    # the deformed segmentation
    # - $EN_[mars,seg]_t$TIME_deformed.inr
    #

    prev_segimage = experiment.get_segmentation_image(previous_time)
    if prev_segimage is None:
        if previous_time is None:
            monitoring.to_log_and_console("    .. " + proc + ": was called with 'previous_time = None'", 2)
        else:
            monitoring.to_log_and_console("    .. " + proc + ": no segmentation image was found for time "
                                          + str(previous_time), 2)
        return None

    prev_def_segimage = common.add_suffix(prev_segimage, "_deformed",
                                          new_dirname=experiment.working_dir.get_tmp_directory(0),
                                          new_extension=experiment.default_image_suffix)

    if os.path.isfile(prev_def_segimage):
        return prev_def_segimage

    deformation = get_deformation_from_current_to_previous(current_time, experiment, parameters, previous_time)

    if deformation is None:
        monitoring.to_log_and_console("    .. " + proc + ": no deformation was found for time "
                                      + str(current_time) + " towards " + str(previous_time), 2)
        return None

    monitoring.to_log_and_console("    .. resampling of '" + str(prev_segimage).split(os.path.sep)[-1] + "'", 2)

    cpp_wrapping.apply_transformation(prev_segimage, prev_def_segimage, deformation,
                                      interpolation_mode='nearest', monitoring=monitoring)

    return prev_def_segimage


########################################################################################
#
#
#
########################################################################################


def build_reconstructed_image(current_time, experiment, parameters, previous_time=None):
    """

    :param current_time:
    :param experiment:
    :param parameters:
    :param previous_time:
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
    # return the name of the reconstructed image. It can be the input image.
    #

    proc = "build_reconstructed_image"
    _trace_ = True
    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, ReconstructionParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    #
    # get fused image
    #
    input_dir = experiment.fusion_dir.get_directory(0)
    input_name = experiment.fusion_dir.get_image_name(current_time)

    input_image = common.find_file(input_dir, input_name, file_type='image', callfrom=proc, local_monitoring=monitoring)

    if input_image is None:
        monitoring.to_log_and_console("    .. image '" + input_name + "' not found in '" + str(input_dir) + "'", 2)
        return None
    input_image = os.path.join(input_dir, input_image)

    reconstructed_image = common.add_suffix(input_image, parameters.result_suffix,
                                            new_dirname=experiment.working_dir.get_rec_directory(0),
                                            new_extension=experiment.default_image_suffix)
    if os.path.isfile(reconstructed_image) and monitoring.forceResultsToBeBuilt is False:
        monitoring.to_log_and_console("    .. reconstructed image is '" +
                                      str(reconstructed_image).split(os.path.sep)[-1] + "'", 2)
        return reconstructed_image

    #
    # the processing flow is as follows
    #
    # input image -> pre-normalized image -> normalized image -> intensity image }
    #                                     -> enhanced image                      } -> reconstructed image
    #                                     -> outer contour enhancement image     }
    #

    #
    # intensity pre-normalization
    # if not 'identity', result image is prefixed with "_prenormalized"
    # further processing will be done on this "pre-normalized" image
    #
    if parameters.intensity_prenormalization.lower() == 'identity':
        pass
    else:
        local_suffix = parameters.prenormalized_suffix + "_prenormalized"
        prenormalized_image = common.add_suffix(input_image, local_suffix,
                                                new_dirname=experiment.working_dir.get_tmp_directory(0),
                                                new_extension=experiment.default_image_suffix)
        if os.path.isfile(prenormalized_image) and monitoring.forceResultsToBeBuilt is False:
            monitoring.to_log_and_console("    .. prenormalization image is '" +
                                          str(prenormalized_image).split(os.path.sep)[-1] + "'", 2)
            input_image = prenormalized_image
        elif parameters.intensity_prenormalization.lower() == 'normalization_to_u8' \
                or parameters.intensity_prenormalization.lower() == 'global_normalization_to_u8':
            monitoring.to_log_and_console("    .. intensity global prenormalization of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            cpp_wrapping.global_intensity_normalization(input_image, prenormalized_image,
                                                        min_percentile=parameters.prenormalization_min_percentile,
                                                        max_percentile=parameters.prenormalization_max_percentile,
                                                        other_options=None, monitoring=monitoring)
            input_image = prenormalized_image
        elif parameters.intensity_prenormalization.lower() == 'normalization_to_u16' \
                or parameters.intensity_prenormalization.lower() == 'global_normalization_to_u16':
            monitoring.to_log_and_console("    .. intensity global prenormalization of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            cpp_wrapping.global_intensity_normalization(input_image, prenormalized_image,
                                                        min_percentile=parameters.prenormalization_min_percentile,
                                                        max_percentile=parameters.prenormalization_max_percentile,
                                                        other_options="-o 2", monitoring=monitoring)
            input_image = prenormalized_image
        else:
            monitoring.to_log_and_console("    unknown intensity prenormalization method: '"
                                          + str(parameters.intensity_prenormalization) + "'", 2)

    #
    # set image bytes
    #
    intensity_image_bytes = 0
    enhanced_image_bytes = 0
    contour_image_bytes = 0

    #
    # histogram-based intensity transformation
    #

    if parameters.intensity_transformation is None or parameters.intensity_transformation.lower() == 'none':
        normalized_image = None
    elif parameters.intensity_transformation.lower() == 'identity':
        normalized_image = input_image
    else:
        normalized_image = common.add_suffix(input_image, parameters.normalized_suffix + "_normalized",
                                             new_dirname=experiment.working_dir.get_tmp_directory(0),
                                             new_extension=experiment.default_image_suffix)
        if os.path.isfile(normalized_image) and monitoring.forceResultsToBeBuilt is False:
            monitoring.to_log_and_console("    .. normalization image is '" +
                                          str(normalized_image).split(os.path.sep)[-1] + "'", 2)
            #
            # set byte counts from parameters
            #
            if parameters.intensity_transformation.lower() == 'normalization_to_u8' \
                    or parameters.intensity_transformation.lower() == 'global_normalization_to_u8' \
                    or parameters.intensity_transformation.lower() == 'cell_normalization_to_u8':
                intensity_image_bytes = 1
            else:
                intensity_image_bytes = 2
        #
        # global normalization on 1 byte
        #
        elif parameters.intensity_transformation.lower() == 'normalization_to_u8' \
                or parameters.intensity_transformation.lower() == 'global_normalization_to_u8':
            intensity_image_bytes = 1
            monitoring.to_log_and_console("    .. intensity global normalization of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            cpp_wrapping.global_intensity_normalization(input_image, normalized_image,
                                                        min_percentile=parameters.normalization_min_percentile,
                                                        max_percentile=parameters.normalization_max_percentile,
                                                        other_options=None, monitoring=monitoring)
        #
        # cell normalization on 1 byte
        #
        elif parameters.intensity_transformation.lower() == 'cell_normalization_to_u8':
            intensity_image_bytes = 1
            monitoring.to_log_and_console("    .. intensity cell-based normalization of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            if previous_time is None:
                monitoring.to_log_and_console("       previous time point was not given", 2)
                monitoring.to_log_and_console("    .. " + proc + ": switch to 'global normalization' ", 1)
                cpp_wrapping.global_intensity_normalization(input_image, normalized_image,
                                                            min_percentile=parameters.normalization_min_percentile,
                                                            max_percentile=parameters.normalization_max_percentile,
                                                            other_options=None, monitoring=monitoring)
            else:
                previous_deformed_segmentation = get_previous_deformed_segmentation(current_time, experiment,
                                                                                    parameters, previous_time)
                cpp_wrapping.cell_intensity_normalization(input_image, previous_deformed_segmentation, normalized_image,
                                                          min_percentile=parameters.normalization_min_percentile,
                                                          max_percentile=parameters.normalization_max_percentile,
                                                          cell_normalization_min_method=parameters.cell_normalization_min_method,
                                                          cell_normalization_max_method=parameters.cell_normalization_max_method,
                                                          sigma=parameters.cell_normalization_sigma,
                                                          other_options=None, monitoring=monitoring)
        #
        # global normalization on 2 bytes
        #
        elif parameters.intensity_transformation.lower() == 'normalization_to_u16' \
                or parameters.intensity_transformation.lower() == 'global_normalization_to_u16':
            intensity_image_bytes = 2
            monitoring.to_log_and_console("    .. intensity global normalization of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            cpp_wrapping.global_intensity_normalization(input_image, normalized_image,
                                                        min_percentile=parameters.normalization_min_percentile,
                                                        max_percentile=parameters.normalization_max_percentile,
                                                        other_options="-o 2", monitoring=monitoring)
        else:
            monitoring.to_log_and_console("    unknown intensity transformation method: '"
                                          + str(parameters.intensity_transformation) + "'", 2)
            normalized_image = input_image

    #
    # Gaussian smoothing of intensity transformed image
    #
    if normalized_image is None:
        intensity_image = None
    elif parameters.intensity_sigma > 0.0:
        intensity_image = common.add_suffix(input_image, parameters.intensity_suffix + "_intensity",
                                            new_dirname=experiment.working_dir.get_tmp_directory(0),
                                            new_extension=experiment.default_image_suffix)
        if os.path.isfile(intensity_image) and monitoring.forceResultsToBeBuilt is False:
            monitoring.to_log_and_console("    .. smoothed intensity transformed image is '" +
                                          str(intensity_image).split(os.path.sep)[-1]
                                          + "'", 2)
            #
            # set byte counts from parameters
            #
            if parameters.intensity_transformation.lower() == 'normalization_to_u8' \
                    or parameters.intensity_transformation.lower() == 'global_normalization_to_u8' \
                    or parameters.intensity_transformation.lower() == 'cell_normalization_to_u8':
                intensity_image_bytes = 1
            else:
                intensity_image_bytes = 2
        else:
            monitoring.to_log_and_console("    .. smoothing " + str(normalized_image).split(os.path.sep)[-1]
                                          + "' with sigma = " + str(parameters.intensity_sigma), 2)
            other_options = "-o " + str(intensity_image_bytes)
            cpp_wrapping.linear_smoothing(normalized_image, intensity_image, parameters.intensity_sigma,
                                          real_scale=True, filter_type='deriche', border=10,
                                          other_options=other_options, monitoring=monitoring)
            if not os.path.isfile(intensity_image):
                monitoring.to_log_and_console("    .. " + proc + ": error when smoothing intensity transformed image")
                intensity_image = normalized_image
    else:
        intensity_image = normalized_image

    #
    # there is a membrane enhancement image
    # enhanced_image has to be set
    #
    if parameters.intensity_enhancement is None or parameters.intensity_enhancement.lower() == 'none':
        enhanced_image = None
    else:
        local_suffix = parameters.enhanced_suffix + "_enhanced"
        #
        # trick to deal with already computed images with different image extension
        #
        enhanced_name = common.add_suffix(str(input_image).split(os.path.sep)[-1], local_suffix)
        enhanced_image = common.find_file(experiment.working_dir.get_tmp_directory(0), enhanced_name, file_type='image',
                                          callfrom=proc, local_monitoring=None, verbose=False)
        enhanced_image_bytes = 1
        if enhanced_image is not None and monitoring.forceResultsToBeBuilt is False:
            enhanced_image = os.path.join(experiment.working_dir.get_tmp_directory(0), enhanced_image)
            monitoring.to_log_and_console("    .. intensity enhanced image is '" +
                                          str(enhanced_image).split(os.path.sep)[-1] + "'", 2)
        else:
            enhanced_image = common.add_suffix(input_image, local_suffix,
                                               new_dirname=experiment.working_dir.get_tmp_directory(0),
                                               new_extension=experiment.default_image_suffix)
            if parameters.intensity_enhancement.lower() == 'gace' or parameters.intensity_enhancement.lower() == 'ace':
                monitoring.to_log_and_console("    .. global enhancement of '"
                                              + str(input_image).split(os.path.sep)[-1] + "'", 2)
                ace.monitoring.copy(monitoring)
                ace.global_membrane_enhancement(input_image, enhanced_image, experiment,
                                                temporary_path=experiment.working_dir.get_tmp_directory(0),
                                                parameters=parameters)
            elif parameters.intensity_enhancement.lower() == 'glace':
                monitoring.to_log_and_console("    .. cell enhancement of '"
                                              + str(input_image).split(os.path.sep)[-1] + "'", 2)
                ace.monitoring.copy(monitoring)
                previous_deformed_segmentation = get_previous_deformed_segmentation(current_time, experiment,
                                                                                    parameters,
                                                                                    previous_time)
                if previous_deformed_segmentation is None:
                    monitoring.to_log_and_console("    .. " + proc + ": switch to 'gace' ", 1)
                    ace.global_membrane_enhancement(input_image, enhanced_image, experiment,
                                                    temporary_path=experiment.working_dir.get_tmp_directory(0),
                                                    parameters=parameters)
                else:
                    ace.cell_membrane_enhancement(input_image, previous_deformed_segmentation, enhanced_image,
                                                  experiment,
                                                  temporary_path=experiment.working_dir.get_tmp_directory(0),
                                                  parameters=parameters)
            else:
                enhanced_image = None
                monitoring.to_log_and_console("    unknown enhancement method: '"
                                              + str(parameters.intensity_enhancement) + "'", 2)

    #
    # outer contour image
    #
    if parameters.outer_contour_enhancement is False:
        contour_image = None
    else:
        local_suffix = parameters.outercontour_suffix + "_outercontour"
        contour_image = common.add_suffix(input_image, local_suffix,
                                          new_dirname=experiment.working_dir.get_tmp_directory(0),
                                          new_extension=experiment.default_image_suffix)
        contour_image_bytes = 1
        if os.path.isfile(contour_image) and monitoring.forceResultsToBeBuilt is False:
            monitoring.to_log_and_console("    .. outer contour image is '" + str(contour_image).split(os.path.sep)[-1]
                                          + "'", 2)
        else:
            astec_name = experiment.astec_dir.get_image_name(current_time)
            deformed_segmentation = common.add_suffix(astec_name, '_deformed_segmentation_from_previous',
                                                      new_dirname=experiment.astec_dir.get_tmp_directory(),
                                                      new_extension=experiment.default_image_suffix)
            previous_segmentation = experiment.get_segmentation_image(previous_time)
            if not os.path.isfile(deformed_segmentation) or monitoring.forceResultsToBeBuilt is True:
                deformation = get_deformation_from_current_to_previous(current_time, experiment, parameters,
                                                                       previous_time)
                cpp_wrapping.apply_transformation(previous_segmentation, deformed_segmentation, deformation,
                                                  interpolation_mode='nearest', monitoring=monitoring)
            cpp_wrapping.outer_contour(deformed_segmentation, contour_image, monitoring=monitoring)

    #
    # Fusion of intensity_image, enhanced_image, contour_image
    #
    if intensity_image is not None and enhanced_image is None and contour_image is None:
        if parameters.keep_reconstruction:
            experiment.working_dir.make_rec_directory()
            shutil.copy2(intensity_image, reconstructed_image)
            return reconstructed_image
        return intensity_image
    elif intensity_image is None and enhanced_image is not None and contour_image is None:
        if parameters.keep_reconstruction:
            experiment.working_dir.make_rec_directory()
            shutil.copy2(enhanced_image, reconstructed_image)
            return reconstructed_image
        return enhanced_image
    elif intensity_image is None and enhanced_image is None and contour_image is not None:
        monitoring.to_log_and_console("    weird, only the outer contour is used for reconstruction", 2)
        if parameters.keep_reconstruction:
            experiment.working_dir.make_rec_directory()
            shutil.copy2(contour_image, reconstructed_image)
            return reconstructed_image
        return contour_image

    #
    # here there are at least two images to be fused
    #

    experiment.working_dir.make_rec_directory()

    if parameters.reconstruction_images_combination.lower() == "maximum" \
            or parameters.reconstruction_images_combination.lower() == "max":
        arit_options = " -max"
    elif parameters.reconstruction_images_combination.lower() == "addition" \
            or parameters.reconstruction_images_combination.lower() == "add":
        arit_options = " -add -o 2"

    else:
        monitoring.to_log_and_console("    unknown image combination operation: '"
                                      + str(parameters.reconstruction_images_combination) + "'", 2)
        return None

    if intensity_image is not None:
        if enhanced_image is not None:
            if parameters.reconstruction_images_combination.lower() == "maximum" \
                    or parameters.reconstruction_images_combination.lower() == "max":
                if intensity_image_bytes == 1 and enhanced_image_bytes == 1:
                    arit_options += " -o 1"
                else:
                    arit_options += " -o 2"
            cpp_wrapping.arithmetic_operation(intensity_image, enhanced_image, reconstructed_image,
                                              other_options=arit_options, monitoring=monitoring)
            if contour_image is not None:
                if parameters.reconstruction_images_combination.lower() == "maximum" \
                        or parameters.reconstruction_images_combination.lower() == "max":
                    if intensity_image_bytes == 1 and enhanced_image_bytes == 1 and contour_image_bytes == 1:
                        arit_options += " -o 1"
                    else:
                        arit_options += " -o 2"
                cpp_wrapping.arithmetic_operation(contour_image, reconstructed_image, reconstructed_image,
                                                  other_options=arit_options, monitoring=monitoring)
            return reconstructed_image
        else:
            if contour_image is not None:
                if parameters.reconstruction_images_combination.lower() == "maximum" \
                        or parameters.reconstruction_images_combination.lower() == "max":
                    if intensity_image_bytes == 1 and contour_image_bytes == 1:
                        arit_options += " -o 1"
                    else:
                        arit_options += " -o 2"
                cpp_wrapping.arithmetic_operation(intensity_image, contour_image, reconstructed_image,
                                                  other_options=arit_options, monitoring=monitoring)
                return reconstructed_image
            else:
                monitoring.to_log_and_console("    weird. There should be a fusion to be done. Returns intensity", 2)
                return intensity_image
    elif enhanced_image is not None:
        if contour_image is not None:
            if parameters.reconstruction_images_combination.lower() == "maximum" \
                    or parameters.reconstruction_images_combination.lower() == "max":
                if enhanced_image_bytes == 1 and contour_image_bytes == 1:
                    arit_options += " -o 1"
                else:
                    arit_options += " -o 2"
            cpp_wrapping.arithmetic_operation(enhanced_image, contour_image, reconstructed_image,
                                              other_options=arit_options, monitoring=monitoring)
            return reconstructed_image
        else:
            monitoring.to_log_and_console("    weird. There should be a fusion to be done. Returns enhanced", 2)
            return enhanced_image
    elif contour_image is not None:
        monitoring.to_log_and_console("    weird. There should be a fusion to be done. Returns contours", 2)
        return contour_image

    #
    # should not reach this point
    #
    return None
