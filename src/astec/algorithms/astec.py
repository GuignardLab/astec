
import os
import sys
import shutil
import time
import multiprocessing
import numpy as np
from scipy import ndimage as nd
import copy
import pandas as pd
from itertools import combinations
from skimage.morphology import h_minima, ball
import pickle as pkl
import raster_geometry as rg

import astec.utils.common as common
import astec.utils.ace as ace
import astec.algorithms.mars as mars
import astec.utils.reconstruction as reconstruction
import astec.utils.ioproperties as ioproperties
import astec.utils.diagnosis as diagnosis
from astec.utils import morphsnakes
from astec.components.spatial_image import SpatialImage
from astec.io.image import imread, imsave, imcopy
from astec.wrapping import cpp_wrapping
from astec.utils import membranes

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

#
#
#
class MorphoSnakeParameters(common.PrefixedParameter):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix=None):
            
        common.PrefixedParameter.__init__(self, prefix=prefix)

        if "doc" not in self.__dict__:
            self.doc = {}

        #
        # number of dilation iterations to get the initialization from 'previous' cell
        #
        doc = "\t Defines the initial shape of the the morphosnake.\n"
        doc += "\t It is the previous (deformed) cell dilated by\n"
        doc += "\t  this amount of iterations\n"
        self.doc['dilation_iterations'] = doc
        self.dilation_iterations = 10
        #
        # maximal number of iterations of the morphosnake
        #
        doc = "\t Maximal number of iteration\n"
        self.doc['iterations'] = doc
        self.iterations = 200
        #
        # threshold on the voxel number to break
        #
        doc = "\t Threshold on the voxel count to stop the morphosnake\n"
        doc += "\t evolution (stability is reached)\n"
        self.doc['delta_voxel'] = doc
        self.delta_voxel = 10**3

        doc = "\t Possible values are 'gradient' or 'image'\n"
        doc += "\t 'image' should be more adapted for membrane images\n"
        self.doc['energy'] = doc
        self.energy = 'image'

        doc = "\t Weight for the morphosnake energy\n"
        self.doc['smoothing'] = doc
        self.smoothing = 3

        doc = "\t Weight for the morphosnake energy\n"
        self.doc['balloon'] = doc
        self.balloon = 1

        doc = "\t Number of processors for parallelism\n"
        self.doc['processors'] = doc
        self.processors = 10

        doc = "\t Possible values are True or False\n"
        doc += "\t If true, mimics historical behavior.\n"
        doc += "\t Kept for test purposes\n"
        self.doc['mimic_historical_astec'] = doc
        self.mimic_historical_astec = False

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('#')
        print('# MorphoSnakeParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        self.varprint('dilation_iterations', self.dilation_iterations, self.doc['dilation_iterations'])
        self.varprint('iterations', self.iterations, self.doc['iterations'])
        self.varprint('delta_voxel', self.delta_voxel, self.doc['delta_voxel'])
        self.varprint('energy', self.energy, self.doc['energy'])
        self.varprint('smoothing', self.smoothing, self.doc['smoothing'])
        self.varprint('balloon', self.balloon, self.doc['balloon'])
        self.varprint('processors', self.processors, self.doc['processors'])
        self.varprint('mimic_historical_astec', self.mimic_historical_astec, self.doc['mimic_historical_astec'])
        print("")
        return

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# MorphoSnakeParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'dilation_iterations', self.dilation_iterations, self.doc['dilation_iterations'])
        self.varwrite(logfile, 'iterations', self.iterations, self.doc['iterations'])
        self.varwrite(logfile, 'delta_voxel', self.delta_voxel, self.doc['delta_voxel'])
        self.varwrite(logfile, 'energy', self.energy, self.doc['energy'])
        self.varwrite(logfile, 'smoothing', self.smoothing, self.doc['smoothing'])
        self.varwrite(logfile, 'balloon', self.balloon, self.doc['balloon'])
        self.varwrite(logfile, 'processors', self.processors, self.doc['processors'])
        self.varwrite(logfile, 'mimic_historical_astec', self.mimic_historical_astec,
                      self.doc['mimic_historical_astec'])
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
        self.dilation_iterations = self.read_parameter(parameters, 'dilation_iterations', self.dilation_iterations)
        self.dilation_iterations = self.read_parameter(parameters, 'astec_MorphosnakeIterations', 
                                                       self.dilation_iterations)
        
        self.iterations = self.read_parameter(parameters, 'iterations', self.iterations)
        self.iterations = self.read_parameter(parameters, 'astec_NIterations', self.iterations)
        
        self.delta_voxel = self.read_parameter(parameters, 'delta_voxel', self.delta_voxel)
        self.delta_voxel = self.read_parameter(parameters, 'astec_DeltaVoxels', self.delta_voxel)
        
        self.energy = self.read_parameter(parameters, 'energy', self.energy)
        self.smoothing = self.read_parameter(parameters, 'smoothing', self.smoothing)
        self.balloon = self.read_parameter(parameters, 'balloon', self.balloon)
        
        self.processors = self.read_parameter(parameters, 'processors', self.processors)
        self.processors = self.read_parameter(parameters, 'astec_nb_proc', self.processors)
        
        self.mimic_historical_astec = self.read_parameter(parameters, 'mimic_historical_astec', 
                                                          self.mimic_historical_astec)
        
    def update_from_parameter_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = common.load_source(parameter_file)
        self.update_from_parameters(parameters)


#
#
#
#
#

class AstecParameters(mars.WatershedParameters, MorphoSnakeParameters):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix="astec_"):

        if "doc" not in self.__dict__:
            self.doc = {}

        doc = "\n"
        doc += "Astec parameters overview:\n"
        doc += "==========================\n"
        doc += "Astec is a segmentation propagation procedure, where\n"
        doc += "images at time t-1 are used to segment the image at time t.\n"
        doc += "1. A first segmentation is done by a seeded watershed, where\n"
        doc += "   the seeds are extracted from the deformed segmentation\n"
        doc += "   image at t-1\n"
        doc += "2. For each cell of this first segmentation, the number of\n"
        doc += "   h-minima, for a range of h values, is studied, to decide\n"
        doc += "   whether a division will occur or not\n"
        doc += "3. This allows to build a new set of seeds, and then, a second\n"
        doc += "   segmentation image.\n"
        doc += "4. A first round of corrections corrects cells that experience\n"
        doc += "   a decrease of volume\n"
        doc += "5. A second round of corrections corrects cells that lose\n"
        doc += "   with respect to the background (morphosnake corrections)\n"
        doc += "   The morphosnakes may over-correct the cells, so this\n"
        doc += "   correction is corrected afterwards\n"
        doc += "\n"
        self.doc['astec_overview'] = doc

        ############################################################
        #
        # initialisation
        #
        ############################################################
        mars.WatershedParameters.__init__(self, prefix=prefix)

        self.seed_reconstruction = reconstruction.ReconstructionParameters(prefix=[self._prefix, "seed_"],
                                                                           suffix="_seed")
        self.seed_reconstruction.intensity_sigma = 0.6
        self.membrane_reconstruction = reconstruction.ReconstructionParameters(prefix=[self._prefix, "membrane_"],
                                                                           suffix="_membrane")
        self.membrane_reconstruction.intensity_sigma = 0.15
        self.morphosnake_reconstruction = reconstruction.ReconstructionParameters(prefix=[self._prefix, "morphosnake_"],
                                                                           suffix="_morphosnake")
        self.morphosnake_reconstruction.intensity_sigma = 0.15

        MorphoSnakeParameters.__init__(self, prefix=prefix)

        # general parameters
        doc = "\t Voxel size of the in- and output images\n"
        self.doc['voxelsize'] = doc
        self.voxelsize = None

        #
        #
        doc = "\t Possible values are 'seeds_from_previous_segmentation', \n"
        doc += "\t 'seeds_selection_without_correction', or None\n"
        doc += "\t Allows to stop the propagation before all stages\n"
        doc += "\t are completed\n"
        self.doc['propagation_strategy'] = doc
        self.propagation_strategy = None

        #
        # erosion of cell from previous segmentation
        #
        # previous_seg_erosion_cell_iterations: maximum number of erosion iteration for cells
        #   if the cell disappears, less iterations are done
        # previous_seg_erosion_cell_min_size: minimal size of a cell to perform erosion
        #

        # self.previous_seg_method = "erode_then_deform"
        doc = "\t Possible values are 'deform_then_erode' or 'erode_then_deform'\n"
        doc += "\t \n"
        doc += "\t \n"
        doc += "\t \n"
        doc += "\t \n"
        self.doc['previous_seg_method'] = doc
        self.previous_seg_method = "deform_then_erode"

        doc = "\t Number of erosions to get any cell seed at t\n"
        doc += "\t (except for the background) from the deformed\n"
        doc += "\t segmentation at t\n"
        self.doc['previous_seg_erosion_cell_iterations'] = doc
        self.previous_seg_erosion_cell_iterations = 10

        doc = "\t Number of erosions to get the background seed at t\n"
        doc += "\t from the deformed segmentation at t\n"
        self.doc['previous_seg_erosion_background_iterations'] = doc
        self.previous_seg_erosion_background_iterations = 25

        # TODO make this a relative value rather than a hardcoded cut-off
        doc = "\t Minimal size of a cell from segmentation at t-1\n"
        doc += "\t to generate a cell.\n"
        doc += "\t Lineage of too small cell then ends at t-1\n"
        self.doc['previous_seg_erosion_cell_min_size'] = doc
        self.previous_seg_erosion_cell_min_size = 1000

        #
        doc = "\t Possible values are True or False\n"
        doc += "\t Default behavior.\n"
        self.doc['background_seed_from_hmin'] = doc
        self.background_seed_from_hmin = True

        doc = "\t Possible values are True or False\n"
        doc += "\t Kept for test purposes.\n"
        self.doc['background_seed_from_previous'] = doc
        self.background_seed_from_previous = False

        #
        # to decide whether there will be division
        #
        doc = "\t Proportion of dynamic range, will be compared to the length of the plateau\n"
        doc += "\t of h-values that result in two seeds in a given cell, to decide whether a cell division\n"
        doc += "\t will occur or not\n"
        self.doc['seed_selection_tau'] = doc
        self.seed_selection_tau = 0.1

        #
        # to decide whether there will be division
        #
        doc = "\t Proportion of dynamic range as the minimum step size to map out the length of the plateau\n"
        doc += "\t of h-values that result in two seeds in a given cell."
        self.doc['minimum_step'] = doc
        self.minimum_step = 0.01

        #
        # threshold
        # cells deformed from previous timepoint that does not have any seed
        # and whose volume (in voxels) is below this threshold are discarded
        # they will correspond to dead-end in the lineage
        #
        doc = "\t Volume threshold\n"
        doc += "\t Too small 'mother' cell without corresponding\n"
        doc += "\t 'daughter' cells will not have a lineage\n"
        self.doc['minimum_volume_unseeded_cell'] = doc
        self.minimum_volume_unseeded_cell = 100
        #
        #
        #determines whether we will try to find a cavity 
        #critical for early mouse embryos, strictly do not use if cell volumes naturally differ a lot
        doc = "\t Blastocyst cavity detection\n"
        doc += "\t - 'cell' with volume larger than cut_off is given the background label and consider the cavity\n"
        self.doc['blastocyst_cavity_detection'] = doc
        self.blastocyst_cavity_detection = False

        #the cut_off (x median volume) from which a cell will be detected at the blastocyst cavity
        doc = "\t Blastocyst cavity detection cut_off\n"
        doc += "\t - 'cell' with volume larger than cut_off is given the background label and consider the cavity\n"
        self.doc['cavity_vol_cut_off'] = doc
        self.cavity_vol_cut_off = 3
        #
        # magic values for the volume checking
        # - volume_minimal_value is in voxel units
        #
        doc = "\t Volume ratio tolerance\n"
        doc += "\t The volume ratio is computed by (Leo's)\n"
        doc += "\t 1 - volume(cell at t-1) / sum volume(cells at t)\n"
        doc += "\t Admissible cells verify\n"
        doc += "\t - volume_ratio_tolerance <= ratio <= volume_ratio_tolerance\n"
        self.doc['volume_ratio_tolerance'] = doc
        self.volume_ratio_tolerance = 0.1

        doc = "\t Volume ratio threshold\n"
        doc += "\t To determine whether the volume change is too small\n"
        self.doc['volume_ratio_threshold'] = doc
        self.volume_ratio_threshold = 0.5

        doc = "\t Volume threshold\n"
        doc += "\t To determine whether a 'daughter' is large enough\n"
        self.doc['volume_minimal_value'] = doc
        self.volume_minimal_value = 1000

        doc = "\t Volume decrease cut-off\n"
        doc += "\t Minimsl ratio to determine whether the volume change between mother and daughter cell is too large\n"
        self.doc['volume_decrease_cut_off'] = doc
        self.volume_decrease_cut_off = 0.7
        #
        # magic values for the membrane sanity checking
        # 
        #
        doc = "\t Membrane sanity check - stringency\n"
        doc += "\t Membranes with a volume ratio smaller than\n"
        doc += "\t a cut-off = median(true membranes) + std(true mebranes) * stringency factor\n"
        doc += "\t are scored as true membranes, other cells are merged into respective bigger cells\n"
        self.doc['membrane_sanity_check'] = doc
        self.membrane_sanity_check = False

        doc = "\t Membrane sanity check - stringency\n"
        doc += "\t Membranes with a volume ratio smaller than\n"
        doc += "\t a cut-off = median(true membranes) + std(true mebranes) * stringency factor\n"
        doc += "\t are scored as true membranes, other cells are merged into respective bigger cells\n"
        self.doc['stringency_factor'] = doc
        self.stringency_factor = 4

        doc = "\t Membrane sanity check - time frame\n"
        doc += "\t This parameter determines how many time points\n"
        doc += "\t are taken into account for calculating the median and std \n"
        doc += "\t for the cut-off value\n"
        self.doc['gt_time_frame'] = doc
        self.gt_time_frame = None

        #
        # outer morphosnake correction
        #
        doc = "\t Possible values are true or False\n"
        self.doc['morphosnake_correction'] = doc
        self.morphosnake_correction = True

        doc = "\t Opening size for the correction of the morphosnake\n"
        doc += "\t correction.\n"
        self.doc['outer_correction_radius_opening'] = doc
        self.outer_correction_radius_opening = 20

        # diagnosis
        doc = "\t Possible values are True or False\n"
        self.doc['lineage_diagnosis'] = doc
        self.lineage_diagnosis = False

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('#')
        print('# AstecParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        for line in self.doc['astec_overview'].splitlines():
            print('# ' + line)

        mars.WatershedParameters.print_parameters(self)
        self.seed_reconstruction.print_parameters()
        self.membrane_reconstruction.print_parameters()
        self.morphosnake_reconstruction.print_parameters()
        MorphoSnakeParameters.print_parameters(self)

        self.varprint('voxelsize', self.voxelsize, self.doc['voxelsize'])

        self.varprint('propagation_strategy', self.propagation_strategy, self.doc['propagation_strategy'])

        self.varprint('previous_seg_method', self.previous_seg_method, self.doc['previous_seg_method'])
        self.varprint('previous_seg_erosion_cell_iterations', self.previous_seg_erosion_cell_iterations,
                      self.doc['previous_seg_erosion_cell_iterations'])
        self.varprint('previous_seg_erosion_background_iterations', self.previous_seg_erosion_background_iterations,
                      self.doc['previous_seg_erosion_background_iterations'])
        self.varprint('previous_seg_erosion_cell_min_size', self.previous_seg_erosion_cell_min_size,
                      self.doc['previous_seg_erosion_cell_min_size'])

        self.varprint('background_seed_from_hmin', self.background_seed_from_hmin,
                      self.doc['background_seed_from_hmin'])
        self.varprint('background_seed_from_previous', self.background_seed_from_previous,
                      self.doc['background_seed_from_previous'])

        self.varprint('seed_selection_tau', self.seed_selection_tau, self.doc['seed_selection_tau'])
        self.varprint('minimum_step', self.minimum_step, self.doc['minimum_step'])

        self.varprint('minimum_volume_unseeded_cell', self.minimum_volume_unseeded_cell,
                      self.doc['minimum_volume_unseeded_cell'])
        
        self.varprint('blastocyst_cavity_detection', self.blastocyst_cavity_detection, self.doc['blastocyst_cavity_detection'])
        self.varprint('cavity_vol_cut_off', self.cavity_vol_cut_off, self.doc['cavity_vol_cut_off'])
        
        self.varprint('volume_ratio_tolerance', self.volume_ratio_tolerance, self.doc['volume_ratio_tolerance'])
        self.varprint('volume_ratio_threshold', self.volume_ratio_threshold, self.doc['volume_ratio_threshold'])
        self.varprint('volume_minimal_value', self.volume_minimal_value, self.doc['volume_minimal_value'])
        self.varprint('volume_decrease_cut_off', self.volume_decrease_cut_off, self.doc['volume_decrease_cut_off'])

        self.varprint('morphosnake_correction', self.morphosnake_correction, self.doc['morphosnake_correction'])
        self.varprint('outer_correction_radius_opening', self.outer_correction_radius_opening,
                      self.doc['outer_correction_radius_opening'])

        self.varprint('lineage_diagnosis', self.lineage_diagnosis, self.doc['lineage_diagnosis'])
        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# AstecParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        for line in self.doc['astec_overview'].splitlines():
            logfile.write('# ' + line + '\n')

        mars.WatershedParameters.write_parameters_in_file(self, logfile)
        self.seed_reconstruction.write_parameters_in_file(logfile)
        self.membrane_reconstruction.write_parameters_in_file(logfile)
        self.morphosnake_reconstruction.write_parameters_in_file(logfile)
        MorphoSnakeParameters.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'voxelsize', self.voxelsize, self.doc['voxelsize'])
        
        self.varwrite(logfile, 'propagation_strategy', self.propagation_strategy, self.doc['propagation_strategy'])

        self.varwrite(logfile, 'previous_seg_method', self.previous_seg_method, self.doc['previous_seg_method'])
        self.varwrite(logfile, 'previous_seg_erosion_cell_iterations', self.previous_seg_erosion_cell_iterations,
                      self.doc['previous_seg_erosion_cell_iterations'])
        self.varwrite(logfile, 'previous_seg_erosion_background_iterations',
                      self.previous_seg_erosion_background_iterations,
                      self.doc['previous_seg_erosion_background_iterations'])
        self.varwrite(logfile, 'previous_seg_erosion_cell_min_size', self.previous_seg_erosion_cell_min_size,
                      self.doc['previous_seg_erosion_cell_min_size'])

        self.varwrite(logfile, 'background_seed_from_hmin', self.background_seed_from_hmin,
                      self.doc['background_seed_from_hmin'])
        self.varwrite(logfile, 'background_seed_from_previous', self.background_seed_from_previous,
                      self.doc['background_seed_from_previous'])

        self.varwrite(logfile, 'seed_selection_tau', self.seed_selection_tau, self.doc['seed_selection_tau'])
        self.varwrite(logfile, 'minimum_step', self.minimum_step, self.doc['minimum_step'])

        self.varwrite(logfile, 'minimum_volume_unseeded_cell', self.minimum_volume_unseeded_cell,
                      self.doc['minimum_volume_unseeded_cell'])
        
        self.varwrite(logfile, 'blastocyst_cavity_detection', self.blastocyst_cavity_detection,
                      self.doc['blastocyst_cavity_detection'])
        self.varwrite(logfile, 'cavity_vol_cut_off', self.cavity_vol_cut_off,
                      self.doc['cavity_vol_cut_off'])
        
        self.varwrite(logfile, 'volume_ratio_tolerance', self.volume_ratio_tolerance,
                      self.doc['volume_ratio_tolerance'])
        self.varwrite(logfile, 'volume_ratio_threshold', self.volume_ratio_threshold,
                      self.doc['volume_ratio_threshold'])
        self.varwrite(logfile, 'volume_minimal_value', self.volume_minimal_value, self.doc['volume_minimal_value'])
        self.varwrite(logfile, 'volume_decrease_cut_off', self.volume_decrease_cut_off, self.doc['volume_decrease_cut_off'])

        self.varwrite(logfile, 'membrane_stringency_factor', self.stringency_factor, self.doc['stringency_factor'])

        self.varwrite(logfile, 'membrane_ground_truth_time_frame', self.gt_time_frame, self.doc['gt_time_frame'])

        self.varwrite(logfile, 'morphosnake_correction', self.morphosnake_correction,
                      self.doc['morphosnake_correction'])
        self.varwrite(logfile, 'outer_correction_radius_opening', self.outer_correction_radius_opening,
                      self.doc['outer_correction_radius_opening'])

        self.varwrite(logfile, 'lineage_diagnosis', self.lineage_diagnosis, self.doc['lineage_diagnosis'])

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
        mars.WatershedParameters.update_from_parameters(self, parameters)
        self.seed_reconstruction.update_from_parameters(parameters)
        self.membrane_reconstruction.update_from_parameters(parameters)
        self.morphosnake_reconstruction.update_from_parameters(parameters)
        MorphoSnakeParameters.update_from_parameters(self, parameters)

        self.voxelsize = self.read_parameter(parameters, 'voxelsize', self.voxelsize)
        
        self.propagation_strategy = self.read_parameter(parameters, 'propagation_strategy', self.propagation_strategy)

        #
        #
        #
        self.previous_seg_method = self.read_parameter(parameters, 'previous_seg_method', self.previous_seg_method)
        self.previous_seg_erosion_cell_iterations = self.read_parameter(parameters,
                                                                        'previous_seg_erosion_cell_iterations',
                                                                        self.previous_seg_erosion_cell_iterations)
        self.previous_seg_erosion_background_iterations = self.read_parameter(parameters,
                                                                              'previous_seg_erosion_background_iterations',
                                                                              self.previous_seg_erosion_background_iterations)
        self.previous_seg_erosion_cell_min_size = self.read_parameter(parameters, 'previous_seg_erosion_cell_min_size',
                                                                      self.previous_seg_erosion_cell_min_size)


        #
        self.background_seed_from_hmin = self.read_parameter(parameters, 'background_seed_from_hmin',
                                                             self.background_seed_from_hmin)
        self.background_seed_from_previous = self.read_parameter(parameters, 'background_seed_from_previous',
                                                                 self.background_seed_from_previous)

        #
        # seed selection
        # added seed_reconstruction.intensity_sigma, seed_intensity_enhancement and seed_outer_contour_enhancement

        self.seed_reconstruction.intensity_sigma = self.read_parameter(parameters, 'seed_intensity_sigma', self.seed_reconstruction.intensity_sigma)
        
        self.seed_reconstruction.intensity_enhancement = self.read_parameter(parameters, 'seed_intensity_enhancement', self.seed_reconstruction.intensity_enhancement)
        
        self.seed_reconstruction.outer_contour_enhancement = self.read_parameter(parameters, 'seed_outer_contour_enhancement', self.seed_reconstruction.outer_contour_enhancement)
        ###

        self.seed_selection_tau = self.read_parameter(parameters, 'seed_selection_tau', self.seed_selection_tau)
        self.seed_selection_tau = self.read_parameter(parameters, 'Thau', self.seed_selection_tau)
        #added minimum step
        self.minimum_step = self.read_parameter(parameters, 'minimum_step', self.minimum_step)
        self.minimum_step = self.read_parameter(parameters, 'Minimum_step', self.minimum_step)

        self.minimum_volume_unseeded_cell = self.read_parameter(parameters, 'minimum_volume_unseeded_cell',
                                                                self.minimum_volume_unseeded_cell)
        
        self.blastocyst_cavity_detection = self.read_parameter(parameters, 'blastocyst_cavity_detection',
                                                                self.blastocyst_cavity_detection)
        self.cavity_vol_cut_off = self.read_parameter(parameters, 'cavity_vol_cut_off',
                                                                self.cavity_vol_cut_off)

        self.volume_ratio_tolerance = self.read_parameter(parameters, 'volume_ratio_tolerance',
                                                          self.volume_ratio_tolerance)

        self.volume_ratio_tolerance = self.read_parameter(parameters, 'volume_ratio_tolerance',
                                                          self.volume_ratio_tolerance)
        self.volume_ratio_tolerance = self.read_parameter(parameters, 'VolumeRatioSmaller',
                                                          self.volume_ratio_tolerance)

        self.volume_ratio_threshold = self.read_parameter(parameters, 'VolumeRatioBigger', self.volume_ratio_threshold)
        self.volume_minimal_value = self.read_parameter(parameters, 'volume_minimal_value', self.volume_minimal_value)
        self.volume_minimal_value = self.read_parameter(parameters, 'MinVolume', self.volume_minimal_value)
        #added volume_decrease_cut_off
        self.volume_decrease_cut_off = self.read_parameter(parameters, 'volume_decrease_cut_off', self.volume_decrease_cut_off)
        

        
        #added parameters for membrane sanity check
        self.membrane_sanity_check = self.read_parameter(parameters, 'membrane_sanity_check', self.membrane_sanity_check)
        self.stringency_factor = self.read_parameter(parameters, 'membrane_stringency_factor', self.stringency_factor)
        self.gt_time_frame = self.read_parameter(parameters, 'membrane_ground_truth_time_frame', self.gt_time_frame)
        ###

        self.morphosnake_correction = self.read_parameter(parameters, 'morphosnake_correction',
                                                          self.morphosnake_correction)
        self.outer_correction_radius_opening = self.read_parameter(parameters, 'outer_correction_radius_opening',
                                                                   self.outer_correction_radius_opening)
        self.outer_correction_radius_opening = self.read_parameter(parameters, 'RadiusOpening',
                                                                   self.outer_correction_radius_opening)

        self.lineage_diagnosis = self.read_parameter(parameters, 'lineage_diagnosis', self.lineage_diagnosis)

    def update_from_parameter_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = common.load_source(parameter_file)
        self.update_from_parameters(parameters)


########################################################################################
#
# some internal procedures
#
########################################################################################

#
# create seeds from previous segmentation
# cells are eroded either with a maximum number of iterations (10 for 'true' cells,
# 25 for the background) or with less iterations if the object to be eroded
# disappears
# Note (GM 15/07/2018): it should be more efficient to use distance maps,
# and to threshold them
#

def _erode_cell(parameters):
    """

    :param parameters:
    :return:
    """
    #
    # Erodes the label i in the label image
    # tmp : binary SpatialImage(bounding box around cell, where True == the cells position)
    # max_size_cell : size max allow for a cell (here put at np.inf)
    # size_cell : size of the cell to erode
    # iterations : maximum number of iterations for normal cells
    # out_iterations : maximum number of iterations for exterior
    # bb : bounding box if tmp in the global image (necessary when computing in parallel)
    # i : label of the cell to erode
    #

    ###### small modification by Gesa: adapting erosion to anisotropic images by creating a scaled structure element
    tmp, iterations, bb, i, voxelsize = parameters
    nb_iter = iterations
    
    # here we make sure that seeds of which we skip the erosion still make it into the seed image
    if nb_iter == 0:
        eroded = tmp
    else:
        if len(set(voxelsize)) > 1:
            # create 3D kernel for erosion which takes the anisotropy of the image into account
            x_dim, y_dim, z_dim = voxelsize
            # use scaling factor to find the smallest possible ellipsoid
            semiaxes = np.array([1/z_dim, 1/y_dim, 1/x_dim])*np.max(voxelsize)
            shape = [int(c) for c in np.ceil(semiaxes*2+1)]
            structure = rg.ellipsoid(shape, semiaxes)
        else:
            structure = None
        
        eroded = nd.binary_erosion(tmp, iterations=nb_iter, structure = structure)
        #making sure that we never remove a seed completely
        while list(np.unique(eroded)) == [False] and nb_iter >= 2:
                nb_iter -= 1
                eroded = nd.binary_erosion(tmp, iterations=nb_iter, structure = structure)
        if list(np.unique(eroded)) == [False] and nb_iter == 1:
            #if one anisotropic iteration is too much, try at least one isotropic iteration
            if structure is not None:
                eroded = nd.binary_erosion(tmp, iterations=nb_iter, structure = None)
                if list(np.unique(eroded)) == [False]:
                    eroded = tmp
            else:
                eroded = tmp
                print("not eroding this seed due to loss of seed in a single iteration")

    return eroded, i, bb


def build_seeds_from_previous_segmentation(label_image, output_image, parameters, nprocessors=26):
    """
    Erodes all the labels in the segmented image seg
    :param label_image: image whose cells are to be eroded
    :param output_image:
    :param parameters:
    :param nprocessors: number maximum of processors allowed to be used
    :return:
    """

    proc = '_build_seeds_from_previous_segmentation'

    #
    # parameter type checking
    #

    if not isinstance(parameters, AstecParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    #
    #
    #

    seg = imread(label_image)

    bboxes = nd.find_objects(seg)
    labels = np.unique(seg)

    pool = multiprocessing.Pool(processes=nprocessors)
    mapping = []

    #
    # since the label_image is obtained through the deformation of the segmentation at previous
    # time point, it may contains 0 values
    ###### small modification by Gesa: adapting erosion to anisotropic images by creating a scaled structure element
    if parameters.voxelsize != None:
        voxelsize = parameters.voxelsize
    else:
        voxelsize = seg.voxelsize
    
    for i in labels[labels != 0]:
        tmp = seg[bboxes[i - 1]] == i
        size_cell = np.sum(tmp)
        if size_cell > parameters.previous_seg_erosion_cell_min_size:
            if i == 1:
                mapping.append((tmp, parameters.previous_seg_erosion_background_iterations, bboxes[i - 1], i, voxelsize))
            else:
                mapping.append((tmp, parameters.previous_seg_erosion_cell_iterations, bboxes[i - 1], i, voxelsize))
        else:
            #adaptation by Gesa: we need to also send seeds we skip through this function, otherwise we loose them in the following
            if i != 1:
                mapping.append((tmp, 0, bboxes[i - 1], i, voxelsize))
                monitoring.to_log_and_console('     .. skip cell ' + str(i) + ', volume (' + str(size_cell) + ') <= '
                                            + str(parameters.previous_seg_erosion_cell_min_size), 2)

    outputs = pool.map(_erode_cell, mapping)
    pool.close()
    pool.terminate()

    seeds = np.zeros_like(seg)
    for eroded, i, bb in outputs:
        seeds[bb][eroded] = i
    print(f"{np.unique(seeds)=}")
    if parameters.voxelsize != None:
        voxelsize = parameters.voxelsize
    else:    
        seeds.voxelsize = seg.voxelsize
    print(f"{voxelsize=}")
    imsave(output_image, seeds)

    return

#
#
#############
#############
#############

#       start of 3rd big nodification by Gesa (Feb 2024):
#       this modification aims at restructuring the way h-values are determined 

#############
#############
#############


#new function:
def calculate_h_by_measuring_h_range_for_two_seeds_in_all_cells(raw_intensity_image, 
                                                                intensity_seed_image, 
                                                                segmentation_image, 
                                                                path_to_previous_seeds, 
                                                                output_path, 
                                                                parameters):
    
    """
    Function to determine for each cell in the embryo whether it has divided between the current and previous time point.

    params:
    raw_intensity_image = current time raw intensity image
    int_image_for_h_min = current time intensity image after pre-processing
    segmentation_image = eroded and transformed segmentation from t-1
    output_path = path where the full seed image will be saved
    
    
    """
    #load images:
    intensity_image = imread(raw_intensity_image)
    intensity_seed_image = imread(intensity_seed_image)
    segmentation_image  = imread(segmentation_image)
    #if the segmentation_image has a background of 1, we need to set it to 0 for the following to work smoothly
    #find background label --> it should be either 0 or 1, whichever has the bigger volume
    vol_0, vol_1 = nd.sum(np.ones_like(segmentation_image), segmentation_image, index = [0, 1])
    if vol_1 > vol_0:
        segmentation_image[segmentation_image == 1] = 0
    #create dilated bounding boxes for all labels, increment is 1 voxel if not specified
    b_boxes = membranes.calculate_incremented_bounding_boxes(segmentation_image)
    #calculate the dynamic range (dr) as the range of signal intensity between the 1 and 99%iles of raw intesnity values
    dr = np.percentile(intensity_image, 99) - np.percentile(intensity_image, 1)
    abs_tau = int(np.round((parameters.seed_selection_tau/100)*dr)) #tau is a percentage >> divide by 100
    results_dict = {} #final result of selected value and wanted number of seeds - label: [h_value, n_seeds]
    correspondences = {}
    full_seed_image = np.zeros_like(segmentation_image)
    previous_seeds = None
    max_label = max(np.unique(segmentation_image))
    for label, box in b_boxes.items():
        if box is not None:
            # temporary result: number of seeds for each h value in a dictionary h_value: n_seeds
            h_n_seeds = {}
            # create a intensity and segmentation sub-images to cut computation time
            sub_im = intensity_seed_image[box]
            seg_sub_im = segmentation_image[box]
            # stepsize = half the length between 1 and dr, first h in the middle of 1 and dr
            stepsize = (dr - 1) * 0.5
            h = int(np.round((dr - stepsize)))
            h_min = h_minima(sub_im, h = h)
            labels, _ = nd.label(h_min) # underscore because value isnt needed
            num_seeds = count_seeds_in_cell_area(labels, seg_sub_im, label) 
            h_n_seeds[h] = num_seeds
            # this defines the accuracy with which we search, in this case 1/10th of tau, rounded to the next highest integer
            minimum_step = int(np.round(dr * parameters.minimum_step))
            if minimum_step < 1:
                minimum_step = 1
            while (num_seeds != 2):
                if num_seeds < 2: 
                    stepsize = stepsize / 2
                    if stepsize < minimum_step:
                        stepsize = minimum_step
                    h = int(np.round(h - stepsize))
                    if (h > 0) and (h not in h_n_seeds.keys()):
                        h_min = h_minima(sub_im, h = h)
                        labels, _ = nd.label(h_min)
                        num_seeds = count_seeds_in_cell_area(labels, seg_sub_im, label)
                        h_n_seeds[h] = num_seeds
                    else: 
                        break  
                elif num_seeds > 2:
                    stepsize = stepsize / 2
                    if stepsize < minimum_step:
                        stepsize = minimum_step
                    h = int(np.round(h + stepsize))
                    if (h < dr) and (h not in h_n_seeds.keys()):
                        h_min = h_minima(sub_im, h = h)
                        labels, _ = nd.label(h_min)
                        num_seeds = count_seeds_in_cell_area(labels, seg_sub_im, label)
                        h_n_seeds[h] = num_seeds
                    else:
                        break
                    
                    
            if num_seeds == 2:
                stepsize = minimum_step
                h_former = h
                while num_seeds == 2:
                # find the lower end of h values that yield 2 seeds
                    h_low = int(np.round(h_former - stepsize))
                    if h_low > 0:
                        h_min = h_minima(sub_im, h = h_low)
                        labels, _ = nd.label(h_min)
                        num_seeds = count_seeds_in_cell_area(labels, seg_sub_im, label)
                        h_former = h_low
                        h_n_seeds[h_low] = num_seeds
                    else:
                        break
                    
                h_former = h
                num_seeds = 2
                while num_seeds == 2:
                    # find the upper end of h values that yield 2 seeds
                    h_up = int(np.round(h_former + stepsize))
                    if h_up < dr:
                        h_min = h_minima(sub_im, h = h_up)
                        labels, _ = nd.label(h_min)
                        num_seeds = count_seeds_in_cell_area(labels, seg_sub_im, label)
                        h_n_seeds[h_up] = num_seeds
                        h_former = h_up
                    else:
                        break
                        
            # after all calculations find the appropriate h value that should be used (and determine the number of seeds)
            h_value, n_seeds = find_appropriate_number_of_seeds(h_n_seeds, abs_tau)
            results_dict[label] = [h_value, n_seeds]
            if n_seeds == 1:
                # at this point I start putting together a seed image for the watershed, so I dont have to loop through all cells again later
                #run h-min with chosen h in sub_im, pick seeds inside cell, add to output array
                final_h_min = h_minima(sub_im, h = h_value)
                labels, _ = nd.label(final_h_min)
                # remove all labels that are not exclusively inside the predicted cell area
                seeds_to_be_removed = set(np.unique(labels[seg_sub_im != label]))
                if len(seeds_to_be_removed) > 0:
                    for seed in seeds_to_be_removed:                    
                        labels[labels == seed] = 0
                # change all seeds in cell to label of mother cell so we have coherent labels
                full_seed_image[box][labels > 0] = label
                correspondences[label] = [label]
                    
            elif n_seeds == 2: # that means the cell has divided
                # run h-min with chosen h in sub_im, pick seeds inside cell, add to output array
                final_h_min = h_minima(sub_im, h = h_value)
                labels, _ = nd.label(final_h_min)
                # find new labels for the two daughter cells, to start a new lineage
                new_labels = [max_label + 1, max_label + 2]
                max_label = max_label + 2
                correspondences[label] = new_labels
                # remove all labels that are not exclusively inside the predicted cell area
                seeds_to_be_removed = set(np.unique(labels[seg_sub_im != label]))
                if len(seeds_to_be_removed) > 0:
                    for seed in seeds_to_be_removed:                    
                        labels[labels == seed] = 0
                # rename the seeds already in labels and then place them into the full seed image
                labels_in_cell = list(np.unique(labels[labels > 0]))
                if len(labels_in_cell) != 2:
                    monitoring.to_log_and_console(f"Something is wrong with cell {label}, it was called dividing but the current time point does not have 2 seeds with the predicted parameters.")
                # this works because we know that there will only be two labels (otherwise this cell wouldnt have been classified as diving)
                # TODO this does not seem to work - a dividing cell pair was asigned the same label after division
                print(f"{labels_in_cell=}")
                for i, old_label in enumerate(labels_in_cell):
                    print(f"{old_label=}", f"{new_labels[i]=}")
                    labels = np.where(labels == old_label, new_labels[i], labels)
                # then map the new labels onto the old for coherent labeling
                full_seed_image[box] = np.where(labels > 0, labels, full_seed_image[box])
    
            else: # in case there are no seeds, we retrieve the seed from previous and assume the cell has not divided
                monitoring.to_log_and_console("       Processing unseeded cell: retrieving seed from previous for cell " + str(label), 2)
                if previous_seeds is None:
                    previous_seeds = imread(path_to_previous_seeds)
                full_seed_image[box] = np.where(previous_seeds[box] == label, previous_seeds[box], full_seed_image[box])
                correspondences[label] = [label]

    # add a seed for the background --> use the background seed from previous (this seed has already been eroded)
    if previous_seeds is None:
        previous_seeds = imread(path_to_previous_seeds)   
    full_seed_image[previous_seeds == 1] = 1

    # save full_seed_image after all cells were processed in output_path
    imsave(output_path, SpatialImage(full_seed_image))
    return results_dict, correspondences           

def count_seeds_in_cell_area(labels, cell_area, label):
    # only count seeds fully enclosed in the cell area prediction from previous segmentation
    # find labels inside and outside of the cell, make sure to only count labels that are only insde and remove background
    seeds_inside = set(np.unique(labels[cell_area == label]))
    seeds_outside = set(np.unique(labels[cell_area != label]))
    fully_enclosed_seeds = seeds_inside.difference(seeds_outside)
    count = len(fully_enclosed_seeds)
    if 0 in fully_enclosed_seeds:
        count -= 1
    return count

def find_appropriate_number_of_seeds(h_n_seeds_dictionary, abs_tau):
    # for each label: calculate the range of h-values that result in two seeds and compare it to tau (absolute value)

    two_seeds_list = [h for h, n in h_n_seeds_dictionary.items() if n == 2]
    if len(two_seeds_list) > 0:
        h_bottom = min(two_seeds_list)
        h_top = max(two_seeds_list)
        h_diff = h_top - h_bottom
    else:
        h_diff = 0
    print(f"{h_diff=}")
    print(f"{abs_tau=}")
    if h_diff > abs_tau:
        n_seeds = 2
        h_value = h_top
    else:
        n_seeds = 1
        one_seed_list = [h for h, n in h_n_seeds_dictionary.items() if n == 1]
        if len(one_seed_list) > 0:
            h_value = max(one_seed_list)
        else:
            #check if any tested h results in more than 0 seeds, otherwise this cell gets a seed from previous
            h_n_seeds_dictionary_no_zeros = {key: value for key, value in h_n_seeds_dictionary.items() if value > 0}
            if len(h_n_seeds_dictionary_no_zeros) > 0:
                h_value = min(h_n_seeds_dictionary_no_zeros, key=h_n_seeds_dictionary_no_zeros.get)
            else:
                h_value = min(h_n_seeds_dictionary, key=h_n_seeds_dictionary.get)
                n_seeds = 0

    return h_value, n_seeds


def _compute_volumes(im):
    """

    :param im:
    :return:
    """
    proc = "_compute_volumes"
    if type(im) is str:
        readim = imread(im)
    elif type(im) is SpatialImage:
        readim = im
    else:
        monitoring.to_log_and_console(str(proc) + ": unhandled type for 'im': " + str(type(im)) + "'")
        return

    labels = np.unique(readim)
    volume = nd.sum(np.ones_like(readim), readim, index=np.int16(labels))
    if type(im) is str:
        del readim
    return dict(list(zip(labels, volume)))


def _update_volume_properties(lineage_tree_information, segmented_image, current_time, experiment):

    time_digits = experiment.get_time_digits_for_cell_id()
    volumes = _compute_volumes(segmented_image)
    volume_key = ioproperties.keydictionary['volume']['output_key']
    # uncomment next line to have a lineage file readable by old post-correction version
    # volume_key = 'volumes_information'
    if volume_key not in lineage_tree_information:
        lineage_tree_information[volume_key] = {}

    dtmp = {}
    for key, value in volumes.items():
        newkey = current_time * 10 ** time_digits + int(key)
        dtmp[newkey] = value
    lineage_tree_information[volume_key].update(dtmp)
    return lineage_tree_information


def _build_correspondences_from_segmentation(segmented_image):
    volumes = _compute_volumes(segmented_image)
    tmp = {}
    # if the volume exists, it means that this cell has a mother cell with the same label
    for key in volumes:
        tmp[key] = [int(key)]
    return tmp


def _update_lineage_properties(lineage_tree_information, correspondences, previous_time, current_time, experiment):

    time_digits = experiment.get_time_digits_for_cell_id()
    lineage_key = ioproperties.keydictionary['lineage']['output_key']
    # uncomment next line to have a lineage file readable by old post-correction version
    # lineage_key = 'lin_tree'
    if lineage_key not in lineage_tree_information:
        lineage_tree_information[lineage_key] = {}

    dtmp = {}
    for key, value in correspondences.items():
        newkey = previous_time * 10**time_digits + int(key)
        vtmp = []
        for i in value:
            vtmp.append(current_time * 10**time_digits + int(i))
        dtmp[newkey] = vtmp

    lineage_tree_information[lineage_key].update(dtmp)
    return lineage_tree_information


def cavity_correction(input_segmentation, seed_image_path, output_segmentation, correspondences, parameters):
    monitoring.to_log_and_console(f"    .. running cavity detection ..", 1)
    input_segmentation = imread(input_segmentation)
    volumes = _compute_volumes(input_segmentation)
    seed_image = imread(seed_image_path)
    median_volume = np.median([v for k, v in volumes.items() if k > 1])
    bg_label = min(volumes.keys())
    print(f"{volumes=}")
    print(f"{median_volume=}")
    print(f"{parameters.cavity_vol_cut_off * median_volume=}")
    cavity_labels = [k for k, v in volumes.items() if v >= (parameters.cavity_vol_cut_off * median_volume)]
    if cavity_labels:
        corr_rev = {value: key for key, value_list in correspondences.items() for value in value_list}
        for label in cavity_labels:
            if label > 1:
                input_segmentation[input_segmentation == label] = bg_label
                seed_image[seed_image == label] = bg_label
                monitoring.to_log_and_console(f"    .. detected cavity: changing cell {label} to background", 1)
                # update correspondences
                if label in set(correspondences.keys()):
                    del correspondences[label]
                else:
                    mother_key = corr_rev[label]
                    correspondences[mother_key].remove(label)
    else:
        monitoring.to_log_and_console(f"    .. did not detect cavity: no changes", 1)

    imsave(output_segmentation, SpatialImage(input_segmentation))
    imsave(seed_image_path, SpatialImage(seed_image))

    return output_segmentation


    ###############
    ###############
    ###############
    ###############
    # 2ND MODIFICATION FOR MEMBRANE SANITY CHECK STARTS HERE (2/2)

def new_membrane_sanity_check(segmentation_image, 
                              previous_segmentation, 
                              selected_seeds, 
                              dataframe_path, 
                              experiment, 
                              parameters, 
                              correspondences, 
                              current_time):

    """
    Function that looks at all new membranes in cells thate have divided and decides whether 
    they should be respected as true membranes, or whether neighbouring cells should be merged.
    This decision is based on the measure "volume_ratio" that describes the ratio of membrane volume before and
    after binary_closing which is supposed to assess their "noisyness". The cut-off for the volume ratio depends on the distribution of membranes in the
    ground truth image, which is used as entry point for the astec_astec pipeline.
    
    Args:
        segmentation_image (str): path to labels image of the current time point
        previous_segmentation (str): path to the image that the segmentation propagatio starts. All membranes in this images are
                                considered correct (ground truth). When loaded this is a 3D array.
        correspondences (dict): Dictionary that maps the labels of the segmentation of t-1 
                                to the labels of the current time point. If cells divided, two labels will be mapped to t-1.
        current_time (int): time point of the time-lapse that is being processed.
        stringency_factor (int): factor for how strict the cut-off for membrane sanity should be calculated. 
                                cut_off = median(volume_ratios of ground truth membranes) + stringency_factor * std(volume_ratios of ground truth membranes)
        gt_time_frame (int or None): number of time points that should be considered for 
                                    calculating the cut-off for the membrane sanity check
        
    returns:
        segmentation_image (ndarray): Merged or original segmentation image.  
        correspondences (dict): Dictionary that maps the labels of the segmentation of t-1 
                                to the labels of the current time point. If cells divided, two labels will be mapped to t-1.
                                This dictionary is updated if cells were fused.
                        
   
    """
    ################ TODO!
    # Make sure the dataframe stays correct even if I repeat a time frame in the middle and the df is already created
    # it seems to be working, but better double check again


    proc = "_membrane_sanity_check"

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, AstecParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)


    # load segmentation image and current seeds (they need to be updated for downstream volume decrease correction)
    curr_seg = imread(segmentation_image)
    updated_seeds = imread(selected_seeds)
    # find interfaces between all cells
    interfaces, mapper = membranes.extract_touching_surfaces(curr_seg)

    #reverse correspondences dictionary to make cell_id mapping easier
    corr_rev = {value: key for key, value_list in correspondences.items() for value in value_list}
    if parameters.voxelsize != None:
        voxelsize = parameters.voxelsize
    else:
        voxelsize = curr_seg.voxelsize
    # get labels of cells that have divided between the previous and current time point and make sure there are only combinations of two
    newly_div_cells = [value for value in correspondences.values() if len(value) > 1]
    newly_div_cell_pairs = []
    for pair in newly_div_cells:
        newly_div_cell_pairs = newly_div_cell_pairs + list(combinations(pair, 2))

    # derive newly created membrane ids from new cell pairs
    new_membrane_ids = [value for key, value in mapper.items() if
                            key in newly_div_cell_pairs]

    # calculate volume_ratios for all cells, subset for cells that have divided 
    volume_ratios_all_current, volumes_all_current = membranes.volume_ratio_after_closing(interfaces, mapper, voxelsize)
    volume_ratios_new = {key: value for key, value in volume_ratios_all_current.items() if key in new_membrane_ids}


    ################ adapt this to 0 or 1 as background or remove?  
    if 0 in new_membrane_ids:
        print("    !! weird case: background label in dividing cells")
        sys.exit(0)
    stable_membrane_ids = set(np.unique(interfaces)).difference(set(new_membrane_ids))
    if 0 in stable_membrane_ids:
        stable_membrane_ids.remove(0)
        
    # load or create dataframe that has metrics of all true membranes
    if not os.path.isfile(dataframe_path):
        
        # load previous image that was processed, since if there was no dataframe constructed yet, 
        # this will be the first image and therefore considered ground truth (gt)
        gt_image = imread(previous_segmentation)

        # calculate volumes for ground truth
        gt_interfaces, gt_mapper = membranes.extract_touching_surfaces(gt_image)
        gt_true_membranes_volumes, gt_volumes = membranes.volume_ratio_after_closing(gt_interfaces, gt_mapper, voxelsize)
        # create dataframe with ground truth metrics
        gt_volumes_df = pd.DataFrame(columns = ["mem_id", "cells"])
        gt_volumes_df["mem_id"] = gt_mapper.values()
        gt_volumes_df["cells"] = gt_mapper.keys()
        gt_volumes_df["time_point"] = current_time - 1
        gt_volumes_df["mem_lineage_id"] = gt_mapper.values()
        # create extra dataframes with volume ratios/ volumes and  membrane ids to make sure we don't mix up values
        add_vol_ratio_df = pd.DataFrame.from_dict(gt_true_membranes_volumes, orient = "index", columns = ["mem_volume_ratio"])
        gt_volumes_df = gt_volumes_df.join(add_vol_ratio_df, on = "mem_id")
        add_vol_df = pd.DataFrame.from_dict(gt_volumes, orient = "index", columns = ["volume"])
        gt_volumes_df = gt_volumes_df.join(add_vol_df, on = "mem_id")
        gt_volumes_df["membrane_classification"] = True
        gt_volumes_df["new_contact"] = False
        gt_volumes_df["cut_off"] = None
    # if the dataframe already exists, just load it        
    else:
        gt_volumes_df = pd.read_pickle(dataframe_path)

    # add stable membranes (stable_membrane_ids) from current time point to ground truth
    mem_id_add = list(stable_membrane_ids)
    former_uncertain_pairs = list(gt_volumes_df.loc[(gt_volumes_df["membrane_classification"] == "uncertain") & (gt_volumes_df["time_point"] == current_time-1)]["cells"])
    
    
    previous_cell_pairs = set(gt_volumes_df.loc[gt_volumes_df["time_point"] == current_time-1]["cells"])
    max_lineage_id = gt_volumes_df["mem_lineage_id"].astype(int).max()

    append_df = pd.DataFrame(index = list(range (0, len(mem_id_add))), 
                            columns = ["mem_id", "cells", "time_point", "mem_lineage_id", "mem_volume_ratio", 
                                        "volume", "membrane_classification", "new_contact", "cut_off"])
    index = 0
    for mem_id in mem_id_add:
        append_df["mem_id"][index] = mem_id
        cell_pair = [k for k, v in mapper.items() if v == mem_id][0]
        append_df["cells"][index] = cell_pair
        append_df["time_point"][index] = current_time
        append_df["mem_volume_ratio"][index] = [v for k, v in volume_ratios_all_current.items() if k == mem_id][0]
        append_df["volume"][index] = [v for k, v in volumes_all_current.items() if k == mem_id][0]

        # convert to id's from the previous time point
        pair_in_prev_ids = membranes.translate_cell_pair_to_previous(cell_pair, corr_rev)
        if pair_in_prev_ids in former_uncertain_pairs:
            append_df["membrane_classification"][index] = "uncertain"
        else:
            append_df["membrane_classification"][index] = True
        if pair_in_prev_ids in previous_cell_pairs:
            append_df["new_contact"][index] = False
            mem_lineage = gt_volumes_df.loc[(gt_volumes_df["time_point"] == current_time-1) & 
                                            (gt_volumes_df["cells"] == pair_in_prev_ids)]["mem_lineage_id"].values[0]
            append_df["mem_lineage_id"][index] = mem_lineage
        else:
            append_df["new_contact"][index] = True
            max_lineage_id += 1
            append_df["mem_lineage_id"][index] = max_lineage_id
        append_df["cut_off"][index] = None
        index += 1
    gt_volumes_df = pd.concat([gt_volumes_df, append_df])
    gt_volumes_df.reset_index(drop = True, inplace = True)

    # sliding window for n time points that are used for calculating the volume_ratio cut-off from ground truth
    
    if not parameters.gt_time_frame == None:    
        lower_time_limit = current_time - parameters.gt_time_frame
    
    if (not parameters.gt_time_frame == None) and (experiment.first_time_point < lower_time_limit):  
        considered_gt_membranes = gt_volumes_df.loc[gt_volumes_df["time_point"] > lower_time_limit]
    # if the above condition is not true, we want to consider all gt membranes in the current gt dataframe
    else:
        considered_gt_membranes = gt_volumes_df

    considered_gt_membranes_certain = considered_gt_membranes.loc[considered_gt_membranes["membrane_classification"] == True]

    # calculate interfaces and volumes of considered gt membranes and find cut-off value for volume ratio 
    median_gt = considered_gt_membranes_certain["mem_volume_ratio"].median()
    std_gt = considered_gt_membranes_certain["mem_volume_ratio"].std()
    cut_off = float(median_gt + parameters.stringency_factor * std_gt)
    cut_off_grey_zone = float(median_gt + ((parameters.stringency_factor/2) * std_gt))

    # start sanity check only if any cells were called dividing or uncertain
    uncertain_membrane_ids = list(gt_volumes_df.loc[(gt_volumes_df["membrane_classification"] == "uncertain") 
                                                    & (gt_volumes_df["time_point"] == current_time)]["mem_id"])

    merged_segmentation_name = common.add_suffix(segmentation_image, "_after_sanity_check_t" + str('{:03d}'.format(current_time)),
                                                new_dirname=experiment.astec_dir.get_tmp_directory(),
                                                new_extension=experiment.result_image_suffix)
    
    merged_seeds_name = common.add_suffix(selected_seeds, "_after_sanity_check_t" + str('{:03d}'.format(current_time)),
                                                new_dirname=experiment.astec_dir.get_tmp_directory(),
                                                new_extension=experiment.result_image_suffix)
    
    if (len(newly_div_cell_pairs) == 0) & (len(uncertain_membrane_ids) == 0):
        monitoring.to_log_and_console('      .. found no new or uncertain cell membranes: not running membrane sanity check', 2)
        # save dataframe in main directory
        gt_volumes_df.to_pickle(dataframe_path)
        # save image here with new name to keep track that it went through the membrane sanity check
        #imsave(merged_segmentation_name, SpatialImage(curr_seg, voxelsize=voxelsize).astype(np.uint16))
        #imsave(merged_seeds_name, SpatialImage(selected_seeds, voxelsize=voxelsize).astype(np.uint16))
        return segmentation_image, selected_seeds, correspondences
    
    else:
        monitoring.to_log_and_console('      .. found dividing cells or uncertain membranes: running membrane sanity check', 2)
        volume_ratios_uncertain = {key: value for key, value in volume_ratios_all_current.items() if key in uncertain_membrane_ids}
        
        # define which membranes are false and merge pairs of cells that are connected through that membrane
        new_for_fusion =  [mem_id for mem_id, volume in volume_ratios_new.items() if volume > cut_off]
        uncertain_for_fusion = [mem_id for mem_id, volume in volume_ratios_uncertain.items() if volume > cut_off]

        passed_new_membranes =  [mem_id for mem_id, volume in volume_ratios_new.items() if volume <= cut_off]
        passed_uncertain_membranes = [mem_id for mem_id, volume in volume_ratios_uncertain.items() if volume <= cut_off] 

        # find pairs of cells that should be fused
        new_false_pairs = [key for key, value in mapper.items() if value in new_for_fusion]
        uncertain_false_pairs =  [key for key, value in mapper.items() if value in uncertain_for_fusion]
        false_pairs_list = new_false_pairs + uncertain_false_pairs
        if len(false_pairs_list) > 0:
            monitoring.to_log_and_console('      .. fusing cell pairs:' + f"{false_pairs_list}", 2)
            # merge cells that are seperated by false membranes in the segmentation image
            merged_segmentation, updated_seeds, changed_cells = membranes.merge_labels_with_false_membranes(false_pairs_list, 
                                                                                                      curr_seg, 
                                                                                                      updated_seeds,
                                                                                                      correspondences)
            # output will first be saved in tmp and copied to main in astec_process as final result of membrane sanity check
            imsave(merged_segmentation_name, SpatialImage(merged_segmentation, voxelsize=voxelsize).astype(np.uint16))
            imsave(merged_seeds_name, SpatialImage(updated_seeds, voxelsize=voxelsize).astype(np.uint16))
            correspondences = membranes.update_correspondences_dictionary(correspondences,
                                                                            changed_cells
                                                                            )
        else:
            monitoring.to_log_and_console('      .. did not find false membranes, not fusing any cell pairs', 2)
            imsave(merged_segmentation_name, SpatialImage(curr_seg, voxelsize=voxelsize).astype(np.uint16))
            imsave(merged_seeds_name, SpatialImage(updated_seeds, voxelsize=voxelsize).astype(np.uint16))
        mem_id_add = passed_new_membranes
        append_df = pd.DataFrame(index = list(range (0, len(mem_id_add))), 
                                columns = ["mem_id", "cells", "time_point", "mem_lineage_id", 
                                            "mem_volume_ratio", "volume", "membrane_classification", "new_contact", "cut_off"])
        index = 0
        for mem_id in mem_id_add:
            append_df["mem_id"][index] = mem_id
            cell_pair = [k for k, v in mapper.items() if v == mem_id][0]
            append_df["cells"][index] = cell_pair
            append_df["time_point"][index] = current_time
            mem_volume_ratio =  [v for k, v in volume_ratios_all_current.items() if k == mem_id][0]
            append_df["mem_volume_ratio"][index] = mem_volume_ratio
            append_df["volume"][index] = [v for k, v in volumes_all_current.items() if k == mem_id][0]
            
            if mem_volume_ratio >= cut_off_grey_zone:
                append_df["membrane_classification"][index] = "uncertain"
            else:
                append_df["membrane_classification"][index] = True

            append_df["new_contact"][index] = True
            max_lineage_id += 1
            append_df["mem_lineage_id"][index] = max_lineage_id
                
            append_df["cut_off"][index] = (cut_off_grey_zone, cut_off)
            index += 1
        gt_volumes_df = pd.concat([gt_volumes_df, append_df])
        gt_volumes_df.reset_index(drop = True, inplace = True)
        
        # add new membranes that are classified as false --> dividing cells will be fused
        mem_id_add = new_for_fusion
        append_df = pd.DataFrame(index = list(range (0, len(mem_id_add))), 
                                columns = ["mem_id", "cells", "time_point", "mem_lineage_id", 
                                            "mem_volume_ratio", "volume", "membrane_classification", "new_contact", "cut_off"])
        index = 0
        for mem_id in mem_id_add:
            append_df["mem_id"][index] = mem_id
            cell_pair = [k for k, v in mapper.items() if v == mem_id][0]
            append_df["cells"][index] = cell_pair            
            append_df["time_point"][index] = current_time
            mem_volume_ratio =  [v for k, v in volume_ratios_all_current.items() if k == mem_id][0]
            append_df["mem_volume_ratio"][index] = mem_volume_ratio
            append_df["volume"][index] = [v for k, v in volumes_all_current.items() if k == mem_id][0]
            append_df["membrane_classification"][index] = False
            append_df["new_contact"][index] = True
            max_lineage_id += 1
            append_df["mem_lineage_id"][index] = max_lineage_id
            append_df["cut_off"][index] = (cut_off_grey_zone, cut_off)
            index += 1
        
        # correct entries for uncertain membranes (from before)
        for mem_id in uncertain_for_fusion:
            index = gt_volumes_df.loc[(gt_volumes_df["mem_id"] == mem_id) 
                                              & (gt_volumes_df["time_point"] == current_time)].index[0]
            gt_volumes_df["membrane_classification"][index] = False
            gt_volumes_df["cut_off"][index] = (cut_off_grey_zone, cut_off)
            pair_in_prev_ids = membranes.translate_cell_pair_to_previous(cell_pair, corr_rev)

        for mem_id in passed_uncertain_membranes:
            index = gt_volumes_df.loc[(gt_volumes_df["mem_id"] == mem_id) \
                                              & (gt_volumes_df["time_point"] == current_time)].index[0]
            if gt_volumes_df["mem_volume_ratio"][index] < cut_off_grey_zone:
                gt_volumes_df["membrane_classification"][index] = True
            gt_volumes_df["cut_off"][index] = (cut_off_grey_zone, cut_off)

        # concatenate and save dataframe in main directory
        gt_volumes_df = pd.concat([gt_volumes_df, append_df])
        gt_volumes_df.reset_index(drop = True, inplace = True)
        gt_volumes_df.to_pickle(dataframe_path)
        return merged_segmentation_name, merged_seeds_name, correspondences


########################################################################################
#
#
#
########################################################################################

#########
#########
#########
# Modification by Gesa Feb 24: new volume diagnosis

def _volume_diagnosis(prev_volumes, curr_volumes, correspondences, volume_decrease_cut_off):
    """
    Function that calculates the volume of daughter cells and compares them to the mother cell. 
    If cells lost more than 30% of volume they will be added to a list. In case of cell divisions, total daughter cell volume will be compared
    to the mother's volume.

    params:
    prev_volumes: dictionary mother_cell_id: volume
    prev_volumes: dictionary daughter_cell_id: volume
    correspondences: dictionary with information about lineage mother:[daughter(s)]
    
    """
    proc = "_volume_diagnosis"
    cells_with_volume_loss = {}    
    for mother_c, daughters_c in correspondences.items():
        # skip background
        if mother_c == 1:
            continue
        if mother_c in prev_volumes is False:
            monitoring.to_log_and_console('    ' + proc + ': no volume for cell ' + str(mother_c)
                                          + ' in previous segmentation', 2)
            continue
        
        else:
            if len(daughters_c) == 1:
                vol_ratio = 1/(prev_volumes[mother_c]/curr_volumes[daughters_c[0]])
                if vol_ratio < volume_decrease_cut_off:
                    cells_with_volume_loss[mother_c] = daughters_c
            # measure the volume of the mother and daughter
            # calc a ratio
            # compare ratio to cut-off
            elif len(daughters_c) == 2:
                vol_ratio = 1/(prev_volumes[mother_c]/(curr_volumes[daughters_c[0]]+curr_volumes[daughters_c[1]]))
                if vol_ratio < volume_decrease_cut_off:
                    cells_with_volume_loss[mother_c] = daughters_c
            else:
                monitoring.to_log_and_console('    ' + proc + ': cell ' + str(mother_c)
                                          + ' has an anomalie in the number of daughter cells (more than 2)', 2)
                
    return cells_with_volume_loss


########
########
########
######## modification 3.3 by Gesa (Feb 24): new volume decrease correction function

def _volume_decrease_correction(astec_name, 
                                intensity_seed_image, 
                                previous_segmentation, 
                                segmentation_from_previous, 
                                segmentation_from_selection, 
                                selected_seeds, 
                                membrane_image, 
                                correspondences, 
                                parameters, 
                                experiment):

    """
    Function that adds extra seeds to cell that have significantly decreased in volume compared to the previous time point
    In this case we recompute seeds with h=1 and add all seeds that are fully enclosed in the predicted cell area    
    
    params:
    astec_name: base name for the output images
    intensity_seed_image: input image for h_min search (finding seeds for watershed)
    previous_segmentation: final segmentation from t-1
    segmentation_from_previous: segmentation from watershed on eroded and deformed seeds from t-1
    segmentation_from_selection: segmentation from watershed on seeds determined by h-min in current time point
    selected_seeds: seeds used to generate segmentation from selection
    membrane_image: input image for watershed
    results_division: dictionary mother_cell_id: [h_value, number of daughters]
    correspondences: dictionary with information about lineage mother:[daughter(s)]
    parameters
    experiment
    """

    proc = "_volume_decrease_correction"

    # parameter type checking
    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, AstecParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)
    
    full_seed_image = imread(selected_seeds)
    intensity_seed_image = imread(intensity_seed_image)
    previous_segmentation  = imread(previous_segmentation)
    current_segmentation = imread(segmentation_from_selection)
    segmentation_from_previous = imread(segmentation_from_previous)

    # calculate volumes of cells from current and previous segmentation
    prev_volumes = _compute_volumes(previous_segmentation)
    curr_volumes = _compute_volumes(current_segmentation)

    # check whole embryo volume (although this has no consequences)
    if 0 in set(np.unique(previous_segmentation)):
        prev_embryo_volume = previous_segmentation.size - prev_volumes[0]
    else:
        prev_embryo_volume = previous_segmentation.size - prev_volumes[1]
    if 0 in set(np.unique(current_segmentation)):
        curr_embryo_volume = current_segmentation.size - curr_volumes[0]
    else:
        curr_embryo_volume = current_segmentation.size - curr_volumes[1]
    volume_ratio = 1.0 - prev_embryo_volume / curr_embryo_volume
    if -parameters.volume_ratio_tolerance <= volume_ratio <= parameters.volume_ratio_tolerance:
        pass
    else:
        if volume_ratio < 0:
            monitoring.to_log_and_console('      .. warning: embryo volume has strongly diminished', 2)
        else:
            monitoring.to_log_and_console('      .. warning: embryo volume has strongly increased', 2)
    
    cells_with_volume_loss = _volume_diagnosis(prev_volumes, curr_volumes, correspondences, parameters.volume_decrease_cut_off)
    
    if len(cells_with_volume_loss) > 0:
        monitoring.to_log_and_console('        process cell(s) with unusually large decrease in volume: ' + str(cells_with_volume_loss), 2)
    else:
        monitoring.to_log_and_console('        .. no correction to be done', 2)
        return segmentation_from_selection, selected_seeds, cells_with_volume_loss

    b_boxes = membranes.calculate_incremented_bounding_boxes(segmentation_from_previous, increment = 3)
    # find cells that have not divided and add all seeds inside the cell, keep original seed because it will be the biggest
    change_in_seeds = 0
    #go through cells by mother cell
    for cell in cells_with_volume_loss.keys():
        change_in_seeds += 1
        #in the incremented bounding boxes function, we already adpat the indexing of the boxes (dict cell:box)
        box = b_boxes[cell]
        cell_area = segmentation_from_previous[box]
        # find maximum number of seeds within the predicted area of the cell
        h_min = h_minima(intensity_seed_image[box], h = 1)
        labels, _ = nd.label(h_min)
        # remove seeds that are not fully enclosed in the cell and then add all other seeds to the image
        seeds_to_be_removed = set(np.unique(labels[cell_area != cell]))
        if len(seeds_to_be_removed) > 0:
            for seed in seeds_to_be_removed:                    
                labels[labels == seed] = 0
            # check whether the cell has divided
            siblings = correspondences[cell]
        if len(siblings) == 1:
            # we re-seed the cell with all possible seeds, 
            # if this does not result in sufficient volume gain, we use the seed from previous                
            full_seed_image[box][labels > 0] = cell
        elif len(siblings) == 2:
            # iteratively assign seeds to the two daughter cells, even if just one cell lost volume, to make sure the new watershed doesnt leak
            # run a watershed on the labels image just to see which daughter cell the seed will belong to, later we only copy the seed
            #mars.watershed(labels, intensity_seed_image[box], "/Users/gesaloof/Desktop/dev/woon_second_batch/stack8/test_tmp_watershed.tif", parameters)
            from skimage.segmentation import watershed as skwatershed
            tmp_watershed = skwatershed(intensity_seed_image[box], labels)
            from tifffile import imwrite
            #tmp_watershed = imread("/Users/gesaloof/Desktop/dev/woon_second_batch/stack8/test_tmp_watershed.tif")

            round = 0
            # finding the two offspring cells --> should I just append both cells in the cells_with_volume_loss dict?
            seeds_to_assign = set([x for x in np.unique(tmp_watershed) if x != 0])
            while seeds_to_assign:
                #print(f"{seeds_to_assign=}")
                # TODO compare volume of siblinbgs and bring them to a similar volume before aasigning seeds
                # erosion of the bigger one? Or skipping the bigger one in the first rounds of assignment
                input_seg_image = current_segmentation[box]
                # create masks of +1 voxel around the daughter cells to understand which original cells the labels in tmp_watershed are touching
                if round == 0:
                    mask_cell_a = input_seg_image == siblings[0]
                    mask_cell_b = input_seg_image == siblings[1]
                dil_mask_cell_a = nd.binary_dilation(mask_cell_a)
                dil_mask_cell_b = nd.binary_dilation(mask_cell_b)
                # creates sets of what touches daughter cell a and set of what touches b
                touching_cell_a = set(np.unique(tmp_watershed[dil_mask_cell_a]))
                touching_cell_b = set(np.unique(tmp_watershed[dil_mask_cell_b]))
                all_unique_to_a = touching_cell_a.difference(touching_cell_b)
                unique_to_a = all_unique_to_a.intersection(set(seeds_to_assign))
                all_unique_to_b = touching_cell_b.difference(touching_cell_a)
                unique_to_b = all_unique_to_b.intersection(set(seeds_to_assign))
                all_touching_both = touching_cell_a.intersection(touching_cell_b)
                touching_both = all_touching_both.intersection(set(seeds_to_assign))
                if (round == 0) and (len(unique_to_a) == 0) and (len(unique_to_b) == 0) and (len(touching_both) == 0):
                    monitoring.to_log_and_console(f"       ... found no extra seeds for divided cell {cell}, seeding remains unchanged despite volume loss", 2)
                    break
                # TODO make this cut-off and input parameter so the user can change it in case cells actually divide very asymmetrically
                # if the volumes of the two daughters are more than 10% different, we skip one cell until the volumes have become more similar, then we add seeds to both
                vol_a = np.sum(mask_cell_a)
                vol_b = np.sum(mask_cell_b)

                if vol_a > vol_b * 1.1 and len(unique_to_b) != 0:
                    skip_a = True
                else:
                    skip_a = False
                if vol_b > vol_a * 1.1 and len(unique_to_a) != 0:
                    skip_b = True
                else:
                    skip_b = False

                # assign clear cases, decide doubles based on surface volume
                if skip_a == False:
                    for unique_seed in unique_to_a:
                        full_seed_image[box][labels == unique_seed] = siblings[0]
                        #mask_cell_a[tmp_watershed == unique_seed] = True
                        mask_cell_a = np.where(tmp_watershed == unique_seed, True, mask_cell_a)
                        # remove from list of all labels
                        seeds_to_assign.remove(unique_seed)
                if skip_b == False:
                    for unique_seed in unique_to_b:
                        full_seed_image[box][labels == unique_seed] = siblings[1]
                        #mask_cell_b[tmp_watershed == unique_seed] = True
                        mask_cell_b = np.where(tmp_watershed == unique_seed, True, mask_cell_b)
                        # remove from list of all labels
                        seeds_to_assign.remove(unique_seed)

                # TODO re-think whether both skips need to be false?
                if len(touching_both) > 0 and skip_a == False and skip_b == False:
                    for double_seed in touching_both:
                        if double_seed in seeds_to_assign:
                            #measure touching surface volume, assign to daughter cell with bigger shared volume
                            vol_with_a = np.sum(np.logical_and(tmp_watershed == double_seed, dil_mask_cell_a))
                            vol_with_b = np.sum(np.logical_and(tmp_watershed == double_seed, dil_mask_cell_b))
                            if vol_with_a > vol_with_b:
                                full_seed_image[box][labels == double_seed] = siblings[0]
                                mask_cell_a[tmp_watershed == double_seed] = True
                            else:
                                full_seed_image[box][labels == double_seed] = siblings[1]
                                mask_cell_b[tmp_watershed == double_seed] = True 
                            seeds_to_assign.remove(double_seed)
                round += 1

    if change_in_seeds == 0:
        monitoring.to_log_and_console('    .. no changes in seeds, do not recompute segmentation', 2)
        return segmentation_from_selection, selected_seeds, cells_with_volume_loss

    else:
        # if seeds were added:
        # 1. save the image of corrected seeds
        # 2. redo a watershed
        
        monitoring.to_log_and_console('    .. ' + str(change_in_seeds) + ' changes in seeds, recompute segmentation', 2)

        corr_selected_seeds = common.add_suffix(astec_name, '_seeds_from_corrected_selection',
                                                new_dirname=experiment.astec_dir.get_tmp_directory(),
                                                new_extension=experiment.default_image_suffix)
        if parameters.voxelsize != None:
            voxelsize = parameters.voxelsize
        else:    
            voxelsize = full_seed_image.voxelsize
        imsave(corr_selected_seeds, SpatialImage(full_seed_image, voxelsize=voxelsize).astype(np.uint16))
        del full_seed_image

        segmentation_from_corr_selection = common.add_suffix(astec_name, '_watershed_from_corrected_selection',
                                                            new_dirname=experiment.astec_dir.get_tmp_directory(),
                                                            new_extension=experiment.default_image_suffix)
        mars.watershed(corr_selected_seeds, membrane_image, segmentation_from_corr_selection, parameters)
        return segmentation_from_corr_selection, corr_selected_seeds, cells_with_volume_loss

def _volume_decrease_correction_check(astec_name,
                                        membrane_image,
                                        previous_segmentation,
                                        cells_with_volume_loss, 
                                        correspondences,
                                        parameters,
                                        experiment, 
                                        segmentation_from_selection,
                                        path_to_previous_seeds, 
                                        selected_seeds):
    """
    Function to test whether the reseeding has resulted in a sufficiently increased cell volume. 
    If not, the seeds for this cell will be replaced by the eroded and deformed seed from the previous segmentation.
    This will be done only for cells that have not divided, because otherwise the division will be lost.

    """
    previous_segmentation  = imread(previous_segmentation)
    current_segmentation = imread(segmentation_from_selection)
    seeds_from_previous = imread(path_to_previous_seeds)
    current_seed_image = imread(selected_seeds)

    # calculate volumes and cell with volume loss
    prev_volumes = _compute_volumes(previous_segmentation)
    curr_volumes = _compute_volumes(current_segmentation)
    cells_with_volume_loss = _volume_diagnosis(prev_volumes, curr_volumes, correspondences, parameters.volume_decrease_cut_off)
    
    # if the cell is still too small and has not divided, insert the seed from previous
    for cell in cells_with_volume_loss:
        siblings = correspondences[cell]
        if len(siblings) == 1:
            current_seed_image[seeds_from_previous == cell] = cell
    
    if parameters.voxelsize != None:
        voxelsize = parameters.voxelsize
    else:    
        voxelsize = current_seed_image.voxelsize
    imsave(selected_seeds, SpatialImage(current_seed_image, voxelsize=voxelsize).astype(np.uint16))
    del current_seed_image

    segmentation_from_corr_selection = common.add_suffix(astec_name, '_watershed_from_corrected_selection',
                                                            new_dirname=experiment.astec_dir.get_tmp_directory(),
                                                            new_extension=experiment.default_image_suffix)
    mars.watershed(selected_seeds, membrane_image, segmentation_from_corr_selection, parameters)
    return segmentation_from_corr_selection, selected_seeds       



########################################################################################
#
# Morphosnake correction for cells which experience a volume decrease
# due to background infiltration
#
########################################################################################


#
# MorphoSnakes can be found through the github site
# https://github.com/pmneila/morphsnakes
#
# It seems that the syntax has changed since Leo's original work (branch python2)
# MorphoSnakes can also be found in the scikit-image distribution
#

def _morphosnakes(parameters_for_parallelism):

    mother_c, bb, subimage_name, subsegmentation_name, astec_parameters = parameters_for_parallelism
    proc = "_morphosnakes"
    write_images = False

    #
    # dilation of the mother cell from subsegmentation to get the initial curve
    #
    subsegmentation = imread(subsegmentation_name)
    initialization = nd.binary_dilation(subsegmentation == mother_c, iterations=astec_parameters.dilation_iterations)
    if write_images:
        initialization_name = common.add_suffix(subsegmentation_name, '_initialization')
        imsave(initialization_name, SpatialImage(initialization.astype(np.uint8)))

    energy = astec_parameters.energy
    if energy .lower() == 'gradient':
        pass
    elif energy.lower() == 'image':
        pass
    else:
        monitoring.to_log_and_console(str(proc) + ": unknown energy function, switch to image-based")
        energy = 'image'

    if energy.lower() == 'gradient':
        gradient_name = common.add_suffix(subimage_name, '_gradient')
        cpp_wrapping.gradient_norm(subimage_name, gradient_name)
        energy = imread(gradient_name)
        energy = 1. / np.sqrt(1 + 100 * energy)
    else:
        # energy.lower() == 'image':
        subimage = imread(subimage_name)
        subimage_min = float(subimage.min())
        subimage_max = float(subimage.max())
        energy = 1. / np.sqrt(1 + 100 * (subimage.astype(float) - subimage_min)/(subimage_max - subimage_min))
        del subimage

    if write_images:
        energy_name = common.add_suffix(subsegmentation_name, '_energy')
        imsave(energy_name, SpatialImage(energy.astype(np.float32)))

    macwe = morphsnakes.MorphGAC(energy, smoothing=astec_parameters.smoothing, threshold=1,
                                 balloon=astec_parameters.balloon)
    macwe.levelset = initialization
    before = np.ones_like(initialization)

    step = 1
    for i in range(0, astec_parameters.iterations, step):
        bbefore = copy.deepcopy(before)
        before = copy.deepcopy(macwe.levelset)
        macwe.run(step)
        # print(str(i) + " condition 1 " + str(np.sum(before != macwe.levelset)))
        # print(str(i) + " condition 2 " + str(np.sum(bbefore != macwe.levelset)))
        if write_images:
            if i > 0 and i % 10 == 0:
                result_name = common.add_suffix(subsegmentation_name, '_step' + str(i))
                imsave(result_name, SpatialImage(macwe.levelset.astype(np.uint8)))
        if np.sum(before != macwe.levelset) < astec_parameters.delta_voxel \
                or np.sum(bbefore != macwe.levelset) < astec_parameters.delta_voxel:
            break

    cell_out = macwe.levelset
    cell_out = cell_out.astype(bool)
    del energy

    if write_images:
        levelset_name = common.add_suffix(subsegmentation_name, '_levelset')
        imsave(levelset_name, SpatialImage(cell_out.astype(np.uint8)))

    return mother_c, bb, cell_out


def _historical_morphosnakes(parameters_for_parallelism):

    mother_c, bb, subimage_name, subsegmentation_name, astec_parameters = parameters_for_parallelism
    proc = "_historical_morphosnakes"
    write_images = False

    #
    # dilation of the mother cell from subsegmentation to get the initial curve
    #
    subsegmentation = imread(subsegmentation_name)
    initialization = nd.binary_erosion(subsegmentation != mother_c, iterations=astec_parameters.dilation_iterations,
                                       border_value=1)
    if write_images:
        initialization_name = common.add_suffix(subsegmentation_name, '_initialization')
        imsave(initialization_name, SpatialImage(initialization.astype(np.uint8)))

    gradient_name = common.add_suffix(subimage_name, '_gradient')
    # cpp_wrapping.gradient_norm(subimage_name, gradient_name,  other_options=" -cont 0 ")
    cpp_wrapping.obsolete_gradient_norm(subimage_name, gradient_name)
    energy = imread(gradient_name)
    energy = 1. / np.sqrt(1 + 100 * energy)

    if write_images:
        energy_name = common.add_suffix(subsegmentation_name, '_energy')
        imsave(energy_name, SpatialImage(energy.astype(np.float32)))

    macwe = morphsnakes.MorphGAC(energy, smoothing=astec_parameters.smoothing, threshold=1,
                                 balloon=astec_parameters.balloon)
    macwe.levelset = initialization
    before = np.ones_like(initialization)

    step = 1
    for i in range(0, astec_parameters.iterations, step):
        bbefore = copy.deepcopy(before)
        before = copy.deepcopy(macwe.levelset)
        macwe.step()
        # print(str(i) + " condition 1 " + str(np.sum(before != macwe.levelset)))
        # print(str(i) + " condition 2 " + str(np.sum(bbefore != macwe.levelset)))
        if write_images:
            if i > 0 and i % 10 == 0:
                result_name = common.add_suffix(subsegmentation_name, '_step' + str(i))
                imsave(result_name, SpatialImage(macwe.levelset.astype(np.uint8)))
        if np.sum(before != macwe.levelset) < astec_parameters.delta_voxel \
                or np.sum(bbefore != macwe.levelset) < astec_parameters.delta_voxel:
            break

    cell_out = macwe.levelset
    cell_tmp = nd.binary_fill_holes(cell_out)
    cell_out = (cell_out.astype(np.bool) ^ cell_tmp)
    del energy

    if write_images:
        levelset_name = common.add_suffix(subsegmentation_name, '_levelset')
        imsave(levelset_name, SpatialImage(cell_out.astype(np.uint8)))

    return mother_c, bb, cell_out


def _slices_dilation_iteration(slices, maximum):
    return tuple([slice(max(0, s.start-1), min(s.stop+1, maximum[i])) for i, s in enumerate(slices)])


def _slices_dilation(slices, maximum, iterations=1):
    for i in range(iterations):
        slices = _slices_dilation_iteration(slices, maximum)
    return slices



############
# modified this function to adapt it to the changed volume_diagnosis (Gesa, April '24)

def _outer_volume_decrease_correction(astec_name, previous_segmentation, deformed_segmentation,
                                      segmentation_from_selection, membrane_image, correspondences,
                                      experiment, parameters):
    """
    Correction of cells that experiment a large volume decrease due to the background.
    The correction is done with the morphosnake algorithm.
    :param astec_name: generic name for image file name construction
    :param previous_segmentation: watershed segmentation obtained with segmentation image at previous timepoint
    :param deformed_segmentation: segmentation image at previous timepoint deformed to superimpose current time point
    :param membrane_image:
    :param correspondences: is a dictionary that gives, for each 'parent' cell (in the segmentation built from previous
    time segmentation) (ie the key), the list of 'children' cells (in the segmentation built from selected seeds)
    :param experiment:
    :param parameters:
    :return:
    """

    proc = "_outer_volume_decrease_correction"

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, AstecParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    #
    # compute volumes
    #
    prev_seg = imread(previous_segmentation)
    curr_seg = imread(segmentation_from_selection)

    prev_volumes = _compute_volumes(prev_seg)
    curr_volumes = _compute_volumes(curr_seg)

    #
    # cell volume diagnosis
    # - large_volume_ratio: list of mother cells such that volume(mother) << SUM volume(childrens)
    # - small_volume_ratio: list of mother cells such that volume(mother) >> SUM volume(childrens)
    # - small_volume_daughter: list of [mother_c, daughter_c] such that one of the daughter has
    #                          a too small volume (with its mother not in one of the two above lists)
    # - all_daughter_label: all labels of daughter cells

    cells_with_volume_loss = _volume_diagnosis(prev_volumes, curr_volumes, correspondences, parameters)

    #
    # here we look at cells that experiment a large decrease of volume
    # ie vol(mother) > vol(daughter(s))
    # this is the step (1) of section 2.3.3.6 of L. Guignard thesis
    # [corresponds to the list to_look_at in historical astec code]
    #

    ############################################################
    #
    # BEGIN: cell with large decrease of volume
    # volume_ratio < 0  <=> volume(mother) > SUM volume(childrens)
    # try to add seeds
    #
    ############################################################
    if len(cells_with_volume_loss) > 0:
        monitoring.to_log_and_console('        process cell(s) with large decrease of volume (morphosnake)', 2)
    else:
        del prev_seg
        del curr_seg
        monitoring.to_log_and_console('        .. no correction to be done', 2)
        return segmentation_from_selection, correspondences, []

    #
    # find cells with a volume decrease due to the background
    # volume decrease is checked wrt the deformed segmentation from previous time point
    #
    del prev_seg
    prev_seg = imread(deformed_segmentation)
    bounding_boxes = nd.find_objects(prev_seg)

    exterior_correction = []
    for mother_c in cells_with_volume_loss:

        #
        # create a sub-image where cells 'daughter_c' has the 'True' value
        # and a background at 'False'
        # create a sub-image where the cell 'mother_c' has the 'True' value
        # and a background at 'False'
        #
        # get the labels of the missed part (mother cell - all daughter cells)
        # 'labels' is no more a nd-array, just an array
        #

        bb = bounding_boxes[mother_c-1]
        submask_daughter_c = np.zeros_like(curr_seg[bb])
        for daughter_c in correspondences[mother_c]:
            submask_daughter_c[curr_seg[bb] == daughter_c] = mother_c
        submask_daughter_c = (submask_daughter_c == mother_c)
        submask_mother_c = (prev_seg[bb] == mother_c)

        labels = curr_seg[bb][submask_mother_c & (submask_daughter_c == False)]
        # print("labels are " + str(labels))

        #
        # is the main label 1 (the background label)?
        #
        labels_size = {}
        max_size = 0
        label_max = 0
        for v in np.unique(labels):
            labels_size[v] = np.sum(labels == v)
            if max_size < labels_size[v]:
                max_size = labels_size[v]
                label_max = v
        #
        # the main label is 1, and was also in the bounding box of the 'mother' cell
        # in the deformed segmentation from previous time (historical behavior)
        # [we should have check that the background is adjacent to the mother cell]
        #
        if label_max == 1 and 1 in prev_seg[bb]:
            exterior_correction.append(mother_c)

    if len(exterior_correction) == 0:
        del prev_seg
        del curr_seg
        monitoring.to_log_and_console('        .. no cells with large background part', 2)
        return segmentation_from_selection, correspondences, []

    #
    #
    #
    monitoring.to_log_and_console('        .. (mother) cell(s) to be corrected (morphosnake step): '
                                  + str(exterior_correction), 2)

    #
    # parameters for morphosnake
    #
    greylevel_image = imread(membrane_image)

    mapping = []
    for mother_c in exterior_correction:
        #
        # cell will be dilated by parameters.dilation_iterations to get the initial curve
        # add a margin of 5 voxels
        #
        bb = _slices_dilation(bounding_boxes[mother_c-1], maximum=prev_seg.shape,
                              iterations=parameters.dilation_iterations+5)
        #
        # subimages
        #
        subsegmentation = os.path.join(experiment.astec_dir.get_tmp_directory(), "cellsegment_" + str(mother_c) + "."
                                       + experiment.default_image_suffix)
        subimage = os.path.join(experiment.astec_dir.get_tmp_directory(), "cellimage_" + str(mother_c) + "."
                                + experiment.default_image_suffix)
        imsave(subsegmentation, prev_seg[bb])
        imsave(subimage, greylevel_image[bb])
        parameters_for_parallelism = (mother_c, bb, subimage, subsegmentation, parameters)
        mapping.append(parameters_for_parallelism)

    del greylevel_image
    del prev_seg

    #
    #
    #

    pool = multiprocessing.Pool(processes=parameters.processors)
    if parameters.mimic_historical_astec:
        outputs = pool.map(_historical_morphosnakes, mapping)
    else:
        outputs = pool.map(_morphosnakes, mapping)
    pool.close()
    pool.terminate()

    #
    # do the corrections issued from the morphosnake
    # check whether they are effective
    # is there is a correction to be done (cell has gained over the background),
    # then fuse all "daughter" cell
    #
    # For reproductibility, sort the outputs wrt mother cell id
    #
    outputs.sort()
    effective_exterior_correction = []
    for mother_c, bb, cell_out in outputs:
        daughter_c = correspondences[mother_c][0]
        if np.sum(curr_seg[bb] == 1 & cell_out) > 0:
            effective_exterior_correction.append(daughter_c)
            curr_seg[bb][curr_seg[bb] == 1 & cell_out] = daughter_c
            if len(correspondences[mother_c]) > 1:
                monitoring.to_log_and_console('           cell ' + str(mother_c) + ' will no more divide -> '
                                              + str(daughter_c), 2)
                for d in correspondences[mother_c]:
                    curr_seg[bb][curr_seg[bb] == d] = daughter_c
            correspondences[mother_c] = [daughter_c]

    #
    #
    #
    if len(effective_exterior_correction) == 0:
        del curr_seg
        monitoring.to_log_and_console('        .. no effective correction', 2)
        return segmentation_from_selection, correspondences, []

    #
    # there has been some effective corrections
    # save the image
    #
    monitoring.to_log_and_console('        .. cell(s) that have been corrected (morphosnake step): '
                                  + str(effective_exterior_correction), 2)

    segmentation_after_morphosnakes = common.add_suffix(astec_name, '_morphosnakes',
                                                        new_dirname=experiment.astec_dir.get_tmp_directory(),
                                                        new_extension=experiment.default_image_suffix)
    imsave(segmentation_after_morphosnakes, curr_seg)
    del curr_seg
    #
    #
    #
    return segmentation_after_morphosnakes, correspondences, effective_exterior_correction


def _outer_morphosnake_correction(astec_name, input_segmentation, reference_segmentation, output_segmentation,
                                  effective_exterior_correction, experiment, parameters):
    if len(effective_exterior_correction) == 0:
        return

    #
    #
    #
    opening_name = common.add_suffix(astec_name, '_morphosnakes_opening',
                                     new_dirname=experiment.astec_dir.get_tmp_directory(),
                                     new_extension=experiment.default_image_suffix)
    options = "-opening -R " + str(parameters.outer_correction_radius_opening)
    cpp_wrapping.mathematical_morphology(input_segmentation, opening_name, other_options=options, monitoring=monitoring)

    #
    # remove part of cell that are outside the opening,
    # that was already ouside in the reference segmentation (to prevent over-correction)
    #
    segmentation = imread(input_segmentation)
    reference = imread(reference_segmentation)
    opening = imread(opening_name)
    for daughter_c in effective_exterior_correction:
        segmentation[((segmentation == daughter_c) & (reference == 1) & (opening == 1))] = 1

    del opening
    del reference
    imsave(output_segmentation, segmentation)
    del segmentation
    return


def astec_process(previous_time, current_time, lineage_tree_information, experiment, parameters):
    """

    :param previous_time:
    :param current_time:
    :param lineage_tree_information:
    :param experiment:
    :param parameters:
    :return:
    """

    proc = "astec_process"

    #
    # 1. retrieve the membrane image
    #    it can be the fused image or a calculated image
    # 2. compute the "deformed" segmentation from previous time
    #    a. erode the segmentation from previous time to get seeds
    #    b. deform the seeds
    #    c. segmentation (watershed-based) from the deformed seeds
    # 3. For each cell, compute the number of h-minima for a collection of h
    # 4. For each cell, select a number of h-minima
    #    typically, 1 if no division, or 2 if division
    # 5. Build a seed image from the selected (cell-based) h-minima
    # 6. segmentation (watershed-based) from the built seeds
    #

    #
    # parameter type checking
    #
    print('\n \n \n RUNNING MODIFIED VERSION (Gesa) \n \n \n')


    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, AstecParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    #
    # nothing to do if the segmentation image exists
    #
    astec_dir = experiment.astec_dir.get_directory()
    astec_name = experiment.astec_dir.get_image_name(current_time)
    astec_image = common.find_file(astec_dir, astec_name, file_type='image', callfrom=proc, local_monitoring=None,
                                   verbose=False)

    if astec_image is not None:
        if monitoring.forceResultsToBeBuilt is False:
            monitoring.to_log_and_console('    astec image already existing', 2)
            return
        else:
            monitoring.to_log_and_console('    astec image already existing, but forced', 2)

    astec_image = os.path.join(astec_dir, astec_name + "." + experiment.result_image_suffix)

    #
    # build or read the membrane image to be segmented
    # this is the intensity (elevation) image used for the watershed
    #

    reconstruction.monitoring.copy(monitoring)

    membrane_image = reconstruction.build_reconstructed_image(current_time, experiment,
                                                              parameters.membrane_reconstruction,
                                                              previous_time=previous_time)

    if membrane_image is None:
        monitoring.to_log_and_console("    .. " + proc + ": no membrane image was found/built for time "
                                      + str(current_time), 2)
        return False

    #
    # build seeds by eroding previous segmentation and deforming it
    #
    # erosion iterations are set by default in voxel units
    # there is also a volume defined in voxel units
    #
    monitoring.to_log_and_console('    build seeds from previous segmentation', 2)

    previous_segmentation = experiment.get_segmentation_image(previous_time)
    if previous_segmentation is None:
        monitoring.to_log_and_console("    .. " + proc + ": no segmentation image was found for time "
                                      + str(previous_time), 2)
        return False

    #
    # two ways of getting the seeds from the previous segmentation
    # 1. transform the previous segmentation image then erode the cells
    # 2. erode the cells of the previous segmentation image then transform it
    #    this is the original version
    #
    deformed_seeds = common.add_suffix(astec_name, '_deformed_seeds_from_previous',
                                       new_dirname=experiment.astec_dir.get_tmp_directory(),
                                       new_extension=experiment.default_image_suffix)
    #
    # deformed_segmentation will be needed for morphosnake correction
    # may also be used in reconstruction.py
    #
    deformed_segmentation = common.add_suffix(astec_name, '_deformed_segmentation_from_previous',
                                              new_dirname=experiment.astec_dir.get_tmp_directory(),
                                              new_extension=experiment.default_image_suffix)
    if parameters.voxelsize != None:
        voxelsize = (1, 1, 1)
    else:    
        voxelsize = None

    if not os.path.isfile(deformed_segmentation) or monitoring.forceResultsToBeBuilt is True:
        deformation = reconstruction.get_deformation_from_current_to_previous(current_time, experiment,
                                                                              parameters.membrane_reconstruction,
                                                                              previous_time)
        if deformation is None:
            monitoring.to_log_and_console("    .. " + proc + ": error when getting deformation field")
            return False
        cpp_wrapping.apply_transformation(previous_segmentation, deformed_segmentation, deformation, voxel_size=voxelsize,
                                          interpolation_mode='nearest', monitoring=monitoring)

    if parameters.previous_seg_method.lower() == "deform_then_erode":
        monitoring.to_log_and_console('      transform previous segmentation then erode it', 2)

        if not os.path.isfile(deformed_seeds) or monitoring.forceResultsToBeBuilt is True:
            build_seeds_from_previous_segmentation(deformed_segmentation, deformed_seeds, parameters)

    elif parameters.previous_seg_method.lower() == "erode_then_deform":
        monitoring.to_log_and_console('      erode previous segmentation then transform it', 2)
        eroded_seeds = common.add_suffix(astec_name, '_eroded_seeds_from_previous',
                                         new_dirname=experiment.astec_dir.get_tmp_directory(),
                                         new_extension=experiment.default_image_suffix)
        if not os.path.isfile(eroded_seeds) or monitoring.forceResultsToBeBuilt is True:
            build_seeds_from_previous_segmentation(previous_segmentation, eroded_seeds, parameters)

        if not os.path.isfile(deformed_seeds) or monitoring.forceResultsToBeBuilt is True:
            deformation = reconstruction.get_deformation_from_current_to_previous(current_time, experiment,
                                                                                  parameters.membrane_reconstruction,
                                                                                  previous_time)
            cpp_wrapping.apply_transformation(eroded_seeds, deformed_seeds, deformation, voxel_size=voxelsize,
                                              interpolation_mode='nearest', monitoring=monitoring)
    else:
        monitoring.to_log_and_console("    .. " + proc + ": unknown method '" + str(parameters.previous_seg_method)
                                      + "' to build seeds from previous segmentation")
        return False

    #
    # watershed segmentation with seeds extracted from previous segmentation
    # $\tilde{S}_{t+1}$ in Leo's PhD
    #
    monitoring.to_log_and_console('    watershed from previous segmentation', 2)

    if parameters.propagation_strategy == 'seeds_from_previous_segmentation':
        segmentation_from_previous = astec_image
    else:
        segmentation_from_previous = common.add_suffix(astec_name, '_watershed_from_previous',
                                                       new_dirname=experiment.astec_dir.get_tmp_directory(),
                                                       new_extension=experiment.default_image_suffix)

    #
    # original astec: there is no smoothing of the membrane image
    #
    if not os.path.isfile(segmentation_from_previous) or monitoring.forceResultsToBeBuilt is True:
        mars.watershed(deformed_seeds, membrane_image, segmentation_from_previous, parameters)

    #
    # if the propagation strategy is only to get seeds from the erosion of the previous cells
    # we're done. Update the properties
    #

    if parameters.propagation_strategy == 'seeds_from_previous_segmentation':
        #
        # update volumes and lineage
        #
        lineage_tree_information = _update_volume_properties(lineage_tree_information, astec_image, current_time,
                                                             experiment)

        #
        # update lineage. It is somehow just to check since cells are not supposed to disappeared
        # _erode_cell() performs erosions until the maximal number of iterations is reached
        # or juste before the cell disappears
        #
        correspondences = _build_correspondences_from_segmentation(astec_image)
        lineage_tree_information = _update_lineage_properties(lineage_tree_information, correspondences, previous_time,
                                                              current_time, experiment)
        return lineage_tree_information

    #
    # here parameters.propagation_strategy is not 'seeds_from_previous_segmentation'
    # we continue
    #

    # #
    # # bounding_boxes: bounding boxes for each cell from the watershed segmentation
    # # cells: list of cell labels
    # #
    # first_segmentation = imread(segmentation_from_previous)
    # cells = list(np.unique(first_segmentation))
    # cells.remove(1)
    #bounding_boxes = dict(list(zip(list(range(1, max(cells) + 1)), nd.find_objects(first_segmentation))))
    
    #del first_segmentation

    #
    # for each cell and a collection of h values,
    # - compute a seed image for each value of h
    #   seeds are masked by the cells of the 'previous' segmentation
    # - get a number of seeds per cell of the previous segmentation
    #   and the parameters [h, sigma] that gives the corresponding seed image
    #
    # n_seeds, parameter_seeds are dictionaries indexed by mother cell index
    # n_seeds[mother_c] is an array of the number of seeds
    # parameter_seeds[mother_c] is an array (of same length) that gives the parameters ([h, sigma]) used for the
    #   computation, h being decreasing
    #

    #
    # In the original astec version, seeds are computed from the original fusion image
    # thus a 2-bytes encoded image
    #
    if parameters.seed_reconstruction.is_equal(parameters.membrane_reconstruction):
        monitoring.to_log_and_console("    .. seed image is identical to membrane image", 2)
        image_for_seed = membrane_image
    else:
        image_for_seed = reconstruction.build_reconstructed_image(current_time, experiment,
                                                                  parameters.seed_reconstruction,
                                                                  previous_time=previous_time)
    # seed_image = experiment.fusion_dir.get_image_name(current_time)
    # seed_image = common.find_file(experiment.fusion_dir.get_directory(), seed_image, file_type='image',
    # callfrom=proc, local_monitoring=None, verbose=False)
    if image_for_seed is None:
        monitoring.to_log_and_console("    .. " + proc + " no fused image was found for time " + str(current_time), 2)
        return None

    #
    # seed images will be named after 'image_for_seed'
    # seeds are masked by 'segmentation_from_previous'
    # ie the segmentation obtained with the seeds from the previous time point
    #




    ################
    ################
    ################

    # modification 3.2 by Gesa (Feb 2024)

    input_dir = experiment.fusion_dir.get_directory(0)
    input_name = experiment.fusion_dir.get_image_name(current_time)

    input_image = common.find_file(input_dir, input_name, file_type='image', callfrom=proc, local_monitoring=monitoring)

    if input_image is None:
        monitoring.to_log_and_console("    .. image '" + input_name + "' not found in '" + str(input_dir) + "'", 2)
        return None
    raw_intensity_image = os.path.join(input_dir, input_image)

    path_to_previous_seeds = deformed_seeds

    selected_seeds = common.add_suffix(astec_name, '_seeds_from_selection',
                                        new_dirname=experiment.astec_dir.get_tmp_directory(),
                                        new_extension=experiment.default_image_suffix)
    
    # results_division = dictionary mother_cell_label: [h_value, n_seeds]
    if not os.path.isfile(selected_seeds) or monitoring.forceResultsToBeBuilt is True:
        results_division, correspondences = calculate_h_by_measuring_h_range_for_two_seeds_in_all_cells(raw_intensity_image, 
                                                                    image_for_seed, 
                                                                    segmentation_from_previous, 
                                                                    path_to_previous_seeds, 
                                                                    selected_seeds, 
                                                                    parameters)
        results_division_path = os.path.join(str(experiment.astec_dir.get_tmp_directory()), 'results_division.pkl')
        with open(results_division_path, 'bw') as f:
            pkl.dump(results_division, f)
            f.close() 

        correspondences_path = os.path.join(str(experiment.astec_dir.get_tmp_directory()), 'correspondences.pkl')
        with open(correspondences_path, 'bw') as f:
            pkl.dump(correspondences, f)
            f.close() 
        # printing out the list of cells without seeds and the list of cells that may divide
        dividing_cells = [mother for mother, result_list in results_division.items() if result_list[1] >= 2]
        if len(dividing_cells) > 0:
            monitoring.to_log_and_console('    .. cells at time ' + str(previous_time) + ' that will divide: ' + str(dividing_cells), 2)
        else:
            monitoring.to_log_and_console('    .. no dividing cell between time point ' + str(previous_time) + " and " + str(current_time), 2)
        
    
    else:
        monitoring.to_log_and_console('    .. selected seeds image already exists, proceed to watershed ', 2)
        results_division_path = str(experiment.astec_dir.get_tmp_directory())+'/results_division.pkl'
        with open(results_division_path, "rb") as f:
            results_division = pkl.load(f, encoding="latin1")
            f.close()
        correspondences_path = str(experiment.astec_dir.get_tmp_directory())+'/correspondences.pkl'
        with open(correspondences_path, "rb") as f:
            correspondences = pkl.load(f, encoding="latin1")
            f.close()

    monitoring.to_log_and_console('    watershed from selection of seeds', 2)

    if parameters.propagation_strategy == 'seeds_selection_without_correction':
        segmentation_from_selection = astec_image
    else:
        segmentation_from_selection = common.add_suffix(astec_name, '_watershed_from_selection',
                                                        new_dirname=experiment.astec_dir.get_tmp_directory(),
                                                        new_extension=experiment.default_image_suffix)

    if not os.path.isfile(segmentation_from_selection) or monitoring.forceResultsToBeBuilt is True:
        mars.watershed(selected_seeds, membrane_image, segmentation_from_selection, parameters)

    #
    # if the propagation strategy is to get this segmentation without corrections
    # we're done. Update the properties
    #

    if parameters.propagation_strategy == 'seeds_selection_without_correction':
        #
        # update volumes and lineage
        #
        lineage_tree_information = _update_volume_properties(lineage_tree_information, astec_image, current_time,
                                                             experiment)

        #
        # update lineage.
        #
        lineage_tree_information = _update_lineage_properties(lineage_tree_information, correspondences, previous_time,
                                                              current_time, experiment)
        return lineage_tree_information
    ###############
    ###############
    ###############
    ###############
    # Modification to detect cavity in mouse blastocysts (Gesa, 12.07.24)
    input_segmentation = segmentation_from_selection
        

    if parameters.blastocyst_cavity_detection == True:
        print(f"{parameters.blastocyst_cavity_detection=}")
        output_segmentation = common.add_suffix(astec_name, '_cavity_corrected_segmentation',
                                                        new_dirname=experiment.astec_dir.get_tmp_directory(),
                                                        new_extension=experiment.default_image_suffix)
        print("running cavity detection")
        cavity_corrected_segmentation = cavity_correction(input_segmentation, 
                                                          selected_seeds, 
                                                          output_segmentation,
                                                          correspondences, 
                                                          parameters)
        # copy results from temporary folder to main folder
        imcopy(cavity_corrected_segmentation, astec_image)

        lineage_tree_information = _update_volume_properties(lineage_tree_information, astec_image,
                                                          current_time, experiment)
        
        lineage_tree_information = _update_lineage_properties(lineage_tree_information, correspondences, previous_time,
                                                           current_time, experiment)
        input_segmentation = cavity_corrected_segmentation

    ###############
    ###############
    ###############
    ###############


    
    ###############
    ###############
    ###############
    ###############
    # 1ST MODIFICATION FOR MEMBRANE SANITY CHECK STARTS HERE (1/2)

    if parameters.membrane_sanity_check == True:

        #path where the dataframe containing membranes metrics is stored
        membranes_df_path = os.path.join(astec_dir, "volume_ratio_membranes.pkl")
        #call membrane sanity check and return path to merged image (in tmp folder) and updated correspondences dictionary
        merged_segmentation, selected_seeds, correspondences = new_membrane_sanity_check(input_segmentation, 
                                                                                         previous_segmentation, 
                                                                                         selected_seeds,
                                                                                         membranes_df_path, 
                                                                                         experiment, 
                                                                                         parameters, 
                                                                                         correspondences, 
                                                                                         current_time)
                        
        #set input for next step to be the merged output of this step
        input_segmentation = merged_segmentation

        # copy results from temporary folder to main folder
        imcopy(merged_segmentation, astec_image)
        
        lineage_tree_information = _update_volume_properties(lineage_tree_information, astec_image,
                                                          current_time, experiment)
        
        lineage_tree_information = _update_lineage_properties(lineage_tree_information, correspondences, previous_time,
                                                           current_time, experiment)


    
    # END OF FIRST MODIFICATION
    ###############
    ###############
    ###############
    ###############
    
    

    #
    # Here, we have a first segmentation
    # let it correct it if required


    monitoring.to_log_and_console('    volume decrease correction', 2)

    segmentation_from_selection, selected_seeds, cells_with_volume_loss = _volume_decrease_correction(astec_name, 
                                                                              image_for_seed, 
                                                                              previous_segmentation, 
                                                                              segmentation_from_previous,
                                                                              input_segmentation, 
                                                                              selected_seeds, 
                                                                              membrane_image, 
                                                                              correspondences, 
                                                                              parameters,
                                                                              experiment)

    # test cells that were corrected, if volume has increased sufficiently, else use seed from previous
    segmentation_from_selection, selected_seeds = _volume_decrease_correction_check(astec_name,
                                                                                    image_for_seed,
                                                                                    previous_segmentation,
                                                                                    cells_with_volume_loss, 
                                                                                    correspondences,
                                                                                    parameters,
                                                                                    experiment,
                                                                                    input_segmentation, 
                                                                                    path_to_previous_seeds, 
                                                                                    selected_seeds)    

    input_segmentation = segmentation_from_selection
    

    #
    # Morphosnakes
    #
    # keep the current segmentation as a reference
    # - morphosnakes are supposed to gain volume over the background
    # -> morphosnakes have to be corrected if there is a gain
    # -> the correction should not result in a loss
    #
    reference_segmentation = input_segmentation

    if parameters.morphosnake_correction is True:
        monitoring.to_log_and_console('    outer volume decrease correction (morphosnakes)', 2)
        #
        # reference segmentation for morphosnakes is
        # the segmentation computed from previous time point seeds
        # segmentation_from_previous
        #
        if parameters.morphosnake_reconstruction.is_equal(parameters.membrane_reconstruction):
            monitoring.to_log_and_console("    .. morphosnake image is identical to membrane image", 2)
            image_for_morphosnake = membrane_image
        elif parameters.morphosnake_reconstruction.is_equal(parameters.seed_reconstruction):
            monitoring.to_log_and_console("    .. morphosnake image is identical to seed image", 2)
            image_for_morphosnake = image_for_seed
        else:
            image_for_morphosnake = reconstruction.build_reconstructed_image(current_time, experiment,
                                                                             parameters.morphosnake_reconstruction,
                                                                             previous_time=previous_time)

        output = _outer_volume_decrease_correction(astec_name, previous_segmentation, deformed_segmentation,
                                                   input_segmentation, image_for_morphosnake, correspondences,
                                                   experiment, parameters)
        output_segmentation, correspondences, effective_exterior_correction = output
        input_segmentation = output_segmentation
    else:
        effective_exterior_correction = []

    #
    # Outer correction (to correct morphosnake)
    #
    if len(effective_exterior_correction) > 0:
        monitoring.to_log_and_console('    outer morphosnake correction', 2)
        output_segmentation = common.add_suffix(astec_name, '_morphosnakes_correction',
                                                new_dirname=experiment.astec_dir.get_tmp_directory(),
                                                new_extension=experiment.default_image_suffix)
        _outer_morphosnake_correction(astec_name, input_segmentation, reference_segmentation, output_segmentation,
                                      effective_exterior_correction, experiment, parameters)
        input_segmentation = output_segmentation

    #
    # copy the last segmentation image (in the auxiliary directory) as the result
    imcopy(input_segmentation, astec_image)

    # update volumes and lineage
    #
    lineage_tree_information = _update_volume_properties(lineage_tree_information, astec_image,
                                                         current_time, experiment) 
    lineage_tree_information = _update_lineage_properties(lineage_tree_information, correspondences, previous_time,
                                                          current_time, experiment)
    
    return lineage_tree_information


########################################################################################
#
#
#
########################################################################################

#
# check whether a lineage file exists
# loops over the time points
#

def _get_last_time_from_lineage(lineage_tree_information, first_time_point, delta_time_point=1,
                                time_digits_for_cell_id=4):
    if lineage_tree_information == {}:
        return first_time_point
    if len(lineage_tree_information) > 0 \
            and ioproperties.keydictionary['lineage']['output_key'] in lineage_tree_information:
        monitoring.to_log_and_console("    .. test existing lineage tree", 1)
        cellinlineage = {}
        div = 10**time_digits_for_cell_id
        #
        # time points for 'parent cell' are marked into cellinlineage
        #
        for c in lineage_tree_information[ioproperties.keydictionary['lineage']['output_key']]:
            t = int(c)/div
            if t not in cellinlineage:
                cellinlineage[t] = 1
            else:
                cellinlineage[t] += 1
        first_time = first_time_point
        #
        # if a time point is present in cellinlineage, it means that the 'daughter cell' image has been
        # segmented
        #
        while True:
            if first_time in cellinlineage:
                first_time += delta_time_point
            else:
                return first_time

    return first_time_point


def _get_last_time_from_images(experiment, first_time_point, delta_time_point=1):
    current_time = first_time_point + delta_time_point
    while True:
        segmentation_file = experiment.get_segmentation_image(current_time, verbose=False)
        if segmentation_file is None or not os.path.isfile(segmentation_file):
            return current_time - delta_time_point
        current_time += delta_time_point


def _clean_lineage(lineage_tree_information, first_time_point, time_digits_for_cell_id=4):
    #
    # remove information for time points after the first_time_point
    #
    proc = "_clean_lineage"
    if lineage_tree_information == {}:
        return
    mul = 10 ** time_digits_for_cell_id
    for key in lineage_tree_information:
        if key == ioproperties.keydictionary['lineage']['output_key']:
            tmp = copy.deepcopy(lineage_tree_information[ioproperties.keydictionary['lineage']['output_key']])
            for k in tmp:
                if int(k) > first_time_point * mul:
                    del lineage_tree_information[ioproperties.keydictionary['lineage']['output_key']][k]
        elif key == ioproperties.keydictionary['volume']['output_key']:
            tmp = copy.deepcopy(lineage_tree_information[ioproperties.keydictionary['volume']['output_key']])
            for k in tmp:
                if int(k) > (first_time_point + 1) * mul:
                    del lineage_tree_information[ioproperties.keydictionary['volume']['output_key']][k]
        else:
            monitoring.to_log_and_console(str(proc) + ": unhandled key '" + str(key) + "'")
    return


def _clean_images(experiment, first_time_point, last_time_point, delta_time_point=1):
    #
    # rename images after the first_time_point
    #
    for current_time in range(first_time_point, last_time_point + 1, experiment.delta_time_point):
        segmentation_file = experiment.get_segmentation_image(current_time, verbose=False)
        if segmentation_file is not None and os.path.isfile(segmentation_file):
            print("rename " + str(segmentation_file) + " into " + segmentation_file + ".bak")
            shutil.move(segmentation_file, segmentation_file + ".bak")
            current_time += delta_time_point
    return


def _fill_volumes(lineage_tree_information, first_time_point, experiment):
    proc = "_fill_volumes"
    mul = 10 ** experiment.get_time_digits_for_cell_id()
    if lineage_tree_information != {}:
        if ioproperties.keydictionary['volume']['output_key'] in lineage_tree_information:
            tmp = lineage_tree_information[ioproperties.keydictionary['volume']['output_key']]
            for k in tmp:
                if k / mul == first_time_point:
                    monitoring.to_log_and_console("    .. cell volumes found for time #" + str(first_time_point), 1)
                    return lineage_tree_information
    #
    # no key found for first time point
    #
    first_segmentation = experiment.get_segmentation_image(first_time_point)
    if first_segmentation is None:
        monitoring.to_log_and_console(".. " + proc + ": no segmentation image found for time '" + str(first_time_point)
                                      + "'", 1)
        monitoring.to_log_and_console("\t Exiting", 1)
        sys.exit(1)
    monitoring.to_log_and_console("    .. computes cell volumes for time #" + str(first_time_point), 1)
    return _update_volume_properties(lineage_tree_information, first_segmentation, first_time_point, experiment)


def astec_control(experiment, parameters):
    """

    :param experiment:
    :param parameters:
    :return:
    """

    proc = "astec_control"
    _trace_ = True

    #
    # parameter type checking
    #

    monitoring.to_log_and_console("", 1)

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, AstecParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    if experiment.first_time_point is None:
        monitoring.to_log_and_console(str(proc) + ": first time point was not set: ")
        sys.exit(1)

    if experiment.last_time_point is None:
        monitoring.to_log_and_console(str(proc) + ": last time point was not set: ")
        sys.exit(1)

    #
    # copy monitoring information
    #
    ace.monitoring.copy(monitoring)
    common.monitoring.copy(monitoring)
    mars.monitoring.copy(monitoring)
    ioproperties.monitoring.copy(monitoring)
    reconstruction.monitoring.copy(monitoring)

    #
    # make sure that the result directory exists
    #

    experiment.astec_dir.make_directory()
    segmentation_dir = experiment.astec_dir.get_directory()

    #
    # re-read the lineage file, if any
    # find the restart time point, it is the minimum among
    # - the one issued from the lineage
    # - the one issued from the segmentation images
    # - the one from the parameters
    #
    lineage_tree_file = common.find_file(segmentation_dir, experiment.astec_dir.get_file_name("_lineage"),
                                         file_type='lineage', callfrom=proc, verbose=False)

    if lineage_tree_file is not None and os.path.isfile(os.path.join(segmentation_dir, lineage_tree_file)):
        lineage_tree_path = os.path.join(segmentation_dir, lineage_tree_file)
        lineage_tree_information = ioproperties.read_dictionary(lineage_tree_path)
    else:
        lineage_tree_path = os.path.join(segmentation_dir, experiment.astec_dir.get_file_name("_lineage") + "."
                                         + experiment.result_lineage_suffix)
        lineage_tree_information = {}

    # print(str(lineage_tree_information))

    #
    # check what is the last time point in lineage and in images
    #
    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point
    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()

    last_time_lineage = _get_last_time_from_lineage(lineage_tree_information, first_time_point,
                                                    delta_time_point=experiment.delta_time_point,
                                                    time_digits_for_cell_id=time_digits_for_cell_id)
    last_time_images = _get_last_time_from_images(experiment, first_time_point,
                                                  delta_time_point=experiment.delta_time_point)
    
    #
    # restart is the next time point to be computed
    #
    restart = min(last_time_lineage, last_time_images) + experiment.delta_time_point
    if type(experiment.restart_time_point) == int:
        if experiment.restart_time_point > first_time_point:
            restart = min(restart, experiment.restart_time_point)
    if restart <= last_time_point:
        monitoring.to_log_and_console(".. " + proc + ": start computation at time #" + str(restart), 1)
    else:
        monitoring.to_log_and_console(".. " + proc + ": nothing to do", 1)
    #
    # do some cleaning
    # - the lineage tree: remove lineage information for time points from restart - delta_time_point
    # - the segmentation images
    #
    delta_time_point = experiment.delta_time_point
    if restart <= last_time_point:
        monitoring.to_log_and_console("    .. clean lineage and segmentation directory", 1)
        _clean_lineage(lineage_tree_information, restart - delta_time_point,
                       time_digits_for_cell_id=time_digits_for_cell_id)
        _clean_images(experiment, restart, last_time_point, delta_time_point=delta_time_point)

    #
    # compute volumes for the previous time if required
    # (well, it may exist, if we just restart the segmentation)
    #
    if restart <= last_time_point:
        lineage_tree_information = _fill_volumes(lineage_tree_information, restart - delta_time_point, experiment)

    #
    #
    #
    for current_time in range(restart, last_time_point + 1, experiment.delta_time_point):

        acquisition_time = experiment.get_time_index(current_time)
        previous_time = current_time - experiment.delta_time_point

        #
        # start processing
        #

        monitoring.to_log_and_console('... astec processing of time #' + acquisition_time, 1)
        start_time = time.time()

        #
        # set and make temporary directory
        #
        experiment.astec_dir.set_tmp_directory(current_time)
        experiment.astec_dir.make_tmp_directory()

        if parameters.seed_reconstruction.keep_reconstruction is False \
                and parameters.membrane_reconstruction.keep_reconstruction is False \
                and parameters.morphosnake_reconstruction.keep_reconstruction is False:
            experiment.astec_dir.set_rec_directory_to_tmp()

        #
        # process
        #

        ret = astec_process(previous_time, current_time, lineage_tree_information, experiment, parameters)
        if ret is False:
            monitoring.to_log_and_console('    an error occurs when processing time ' + acquisition_time, 1)
            return False
        else:
            lineage_tree_information = ret

        #
        # cleaning
        #

        if monitoring.keepTemporaryFiles is False:
            experiment.astec_dir.rmtree_tmp_directory()

        #
        # save lineage here
        # thus, we have intermediary pkl in case of future failure
        #

        ioproperties.write_dictionary(lineage_tree_path, lineage_tree_information)

        #
        # end processing for a time point
        #
        end_time = time.time()

        monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
        monitoring.to_log_and_console('', 1)

    if parameters.lineage_diagnosis:
        monitoring.to_log_and_console("    .. test lineage", 1)
        diagnosis_parameters = diagnosis.DiagnosisParameters()
        diagnosis.diagnosis(lineage_tree_information, ['lineage', 'volume'], diagnosis_parameters,
                            time_digits_for_cell_id=time_digits_for_cell_id)

    return

