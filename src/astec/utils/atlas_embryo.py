
import copy
import os
import sys
from collections import Counter

import numpy as np
import sklearn.linear_model as sklm

import astec.utils.common as common
import astec.utils.properties as properties
import astec.utils.diagnosis as udiagnosis
import astec.utils.ioproperties as ioproperties

monitoring = common.Monitoring()


###########################################################
#
#
#
############################################################

class AtlasParameters(udiagnosis.DiagnosisParameters):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix=None):

        if "doc" not in self.__dict__:
            self.doc = {}

        udiagnosis.DiagnosisParameters.__init__(self, prefix=[prefix, "diagnosis_"])

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

        doc = "\t Output directory where to write atlas-individualized output files,"
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
        self.cell_normalization = 'local'

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
    def __init__(self, atlas_properties=None, time_digits_for_cell_id=4):
        proc = "Atlas.init"

        self._properties = {'temporal_alignment': (1.0, 0.0), 'volume_local_estimation': (0.0, 1.0),
                            'target_volume': 6000000}

        if isinstance(atlas_properties, dict):
            if 'cell_lineage' in atlas_properties:
                self.lineage = atlas_properties['cell_lineage']
            else:
                monitoring.to_log_and_console(str(proc) + ": 'cell_lineage' was not in dictionary")
            if 'cell_name' in atlas_properties:
                self.cell_name = atlas_properties['cell_name']
            else:
                monitoring.to_log_and_console(str(proc) + ": 'cell_name' was not in dictionary")
            if 'cell_contact_surface' in atlas_properties:
                self.cell_contact_surface = atlas_properties['cell_contact_surface']
            else:
                monitoring.to_log_and_console(str(proc) + ": 'cell_contact_surface' was not in dictionary")
            if 'cell_volume' in atlas_properties:
                self.cell_volume = atlas_properties['cell_volume']
                self.volume_local_fitting(time_digits_for_cell_id=time_digits_for_cell_id)
            else:
                monitoring.to_log_and_console(str(proc) + ": 'cell_volume' was not in dictionary")

    ############################################################
    #
    # Properties issued from the property files
    #
    ############################################################

    @property
    def lineage(self):
        """
        The cell lineage, as in the property file
        Returns
        -------

        """
        if 'cell_lineage' in self._properties:
            return self._properties['cell_lineage']
        return None

    @lineage.setter
    def lineage(self, atlas_properties):
        self._properties['cell_lineage'] = copy.deepcopy(atlas_properties)
        return

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

    ############################################################
    #
    # Properties (end)
    #
    ############################################################

    #
    # temporal alignment
    #
    def temporally_align_with(self, reference, time_digits_for_cell_id=4):
        proc = "Atlas.temporally_align_with"
        if not isinstance(reference, Atlas):
            monitoring.to_log_and_console(str(proc) + ": 'reference' should be of 'Atlas' class")
            return
        a, b = properties.temporal_alignment(reference.lineage, self.lineage,
                                             time_digits_for_cell_id=time_digits_for_cell_id)
        self._properties['temporal_alignment'] = (a, b)
        return

    #
    #
    #

    def volume_local_fitting(self, time_digits_for_cell_id=4):
        div = 10 ** time_digits_for_cell_id
        volume = self.cell_volume
        volume_along_time = {}
        for c in volume:
            t = int(c) // div
            volume_along_time[t] = volume_along_time.get(t, 0) + volume[c]

        x = list(volume_along_time.keys())
        x = sorted(x)
        y = [volume_along_time[i] for i in x]

        xnp = np.array(x)[:, np.newaxis]
        ynp = np.array(y)[:, np.newaxis]
        ransac = sklm.RANSACRegressor()
        ransac.fit(xnp, ynp)
        self._properties['volume_local_estimation'] = (ransac.estimator_.coef_[0][0], ransac.estimator_.intercept_[0])

    #
    #
    #

    def get_voxelsize_correction(self, timepoint, target_volume=60000000):
        if target_volume != self.target_volume:
            self.target_volume = target_volume
        v_coefficients = self.volume_local_estimation
        t_volume = v_coefficients[0] * timepoint + v_coefficients[1]
        return np.cbrt(self.target_volume / t_volume)


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
            monitoring.to_log_and_console(proc + ": set reference to " + str(names[0]))
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

    def temporal_alignment(self, time_digits_for_cell_id=4):
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
            atlases[a].temporally_align_with(atlases[ref_atlas], time_digits_for_cell_id=time_digits_for_cell_id)
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
        a nested dictionary of neighborhoods, where the keys are ['cell name']['reference name']
        where 'cell name' is the cell name (Conklin), and 'reference name' is the file name,
        a neighborhood is a dictionary of contact surfaces indexed by cell names
        it only considers the first time point after the division

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
        self.temporal_alignment(time_digits_for_cell_id=time_digits_for_cell_id)
        atlases = self.get_atlases()
        for n in atlases:
            msg = "    - "
            msg += "linear time warping of '" + str(n) + "' wrt '" + str(self._ref_atlas) + "' is "
            msg += str(atlases[n].temporal_alignment)
            monitoring.to_log_and_console(msg, 1)
        monitoring.to_log_and_console("    done", 1)

    ############################################################
    #
    #
    #
    ############################################################

    def generate_figure(self, parameters, time_digits_for_cell_id=4):
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
            _figures_temporal_registration(self, parameters, time_digits_for_cell_id=time_digits_for_cell_id)
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
            monitoring.to_log_and_console("... generate neighbors histogram file", 1)
            _figures_neighbor_histogram(self, parameters, time_digits_for_cell_id=time_digits_for_cell_id)
            monitoring.to_log_and_console("... done", 1)


################################################################################
#
# cell number wrt time without and with temporal registration
#
################################################################################

def _figures_temporal_registration(atlases, parameters, time_digits_for_cell_id=4):
    """
    Plot cell number curves without and with linear temporal alignment
    Parameters
    ----------
    atlases
    parameters
    time_digits_for_cell_id

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
        lineage = ref_atlases[n].lineage
        cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
        cells = sorted(cells)
        div = 10 ** time_digits_for_cell_id
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

def _figures_embryo_volume(atlases, parameters, time_digits_for_cell_id=4):
    """

    Parameters
    ----------
    atlases
    parameters
    time_digits_for_cell_id

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
    div = 10 ** time_digits_for_cell_id

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

def _neighbor_histogram(neighbors, atlas, atlasname, threshold=0.05, time_digits_for_cell_id=4):
    contact = atlas.cell_contact_surface
    cellname = atlas.cell_name

    cell_number_threshold = 10
    #
    # nodespertime is a dictionary
    # nodespertime[t] = #nodes at time t
    #
    nodes = list(contact.keys())
    div = 10 ** time_digits_for_cell_id
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


def _figures_neighbor_histogram(atlases, parameters, time_digits_for_cell_id=4):
    """
    Plot cell neighbor number with respect to total cell number
    Parameters
    ----------
    atlases
    parameters
    time_digits_for_cell_id

    Returns
    -------

    """
    proc = "_figures_neighbor_histogram"

    ref_atlases = atlases.get_atlases()
    neighbors = {}
    for a in ref_atlases:
        neighbors = _neighbor_histogram(neighbors, ref_atlases[a], a, time_digits_for_cell_id=time_digits_for_cell_id)

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

    f.close()
