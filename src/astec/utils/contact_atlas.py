import os
import sys
import copy
import operator

import scipy as sp
import scipy.stats as stats
import scipy.cluster.hierarchy as sch
import numpy as np

import astec.utils.common as common
import astec.utils.ioproperties as ioproperties
import astec.utils.properties as properties
import astec.utils.ascidian_name as uname
import astec.utils.contact as ucontact
import astec.utils.diagnosis as udiagnosis

#
# the atlas here is defined as the collection of (daughter) cell neighborhoods
# right after the division
#

monitoring = common.Monitoring()


class AtlasParameters(udiagnosis.DiagnosisParameters):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix='atlas_'):

        if "doc" not in self.__dict__:
            self.doc = {}

        udiagnosis.DiagnosisParameters.__init__(self, prefix=[prefix, "diagnosis_"])

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

        doc = "\t Write out morphonet selection files."
        self.doc['write_selection'] = doc
        self.write_selection = False

        #
        # how to build an atlas
        #
        # add the symmetric neighborhood as neighborhood
        # consider the other half of the embryo as a single cell
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
        doc = "\t Minimum number of atlases required to assess naming confidence. If there is not enough atlases\n"
        doc += "\t in the database for the aimed division, naming is not assessed."
        self.doc['confidence_atlases_nmin'] = doc
        self.confidence_atlases_nmin = 2
        doc = "\t Percentage of atlases used to assessed naming confidence. If the percentage is less than\n"
        doc += "\t 'confidence_atlases_nmin', 'confidence_atlases_nmin' atlases are used."
        self.doc['confidence_atlases_percentage'] = doc
        self.confidence_atlases_percentage = 50

        #
        #
        #
        doc = "\t How to compare two division patterns (a division is considered here\n"
        doc += "\t as the concatenation of the contact surface vectors of the 2 daughter\n"
        doc += "\t cells). Choices are:\n"
        doc += "\t - 'distance': the distance type is given by 'cell_contact_distance'\n"
        doc += "\t   distances are normalized between 0 (perfect match) and 1 (complete mismatch)\n"
        self.doc['division_contact_similarity'] = doc
        self.division_contact_similarity = 'distance'

        #
        #
        #
        doc = "\t True or False. Performs some diagnosis when reading an additional property file \n"
        doc += "\t into the atlases. Incrementing the verboseness ('-v' in the command line) may give\n"
        doc += "\t more details."
        self.doc['diagnosis_properties'] = doc
        self.diagnosis_properties = False

        #
        #
        #
        doc = "\t If True, will propose some daughters switches in the atlases. For a given division,\n"
        doc += "\t a global score is computed as the sum of all pairwise division similarity.\n"
        doc += "\t A switch is proposed for an atlas if it allows to decrease this global score."
        self.doc['division_permutation_proposal'] = doc
        self.division_permutation_proposal = False

        #
        #
        #
        doc = "\t Cluster distance used to build dendrograms. choices are\n"
        doc += "\t - 'single'\n"
        doc += "\t - 'complete'\n"
        doc += "\t - 'average'\n"
        doc += "\t - 'weighted'\n"
        doc += "\t - 'centroid'\n"
        doc += "\t - 'median'\n"
        doc += "\t - 'ward'\n"
        doc += "\t see https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html.\n"
        self.doc['dendrogram_cluster_distance'] = doc
        self.dendrogram_cluster_distance = 'single'
        #
        #
        #
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
        #
        #
        self.cells_to_be_traced = None

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('#')
        print('# CellAtlasParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        udiagnosis.DiagnosisParameters.print_parameters(self)

        self.varprint('atlasFiles', self.atlasFiles)
        self.varprint('referenceAtlas', self.referenceAtlas)

        self.varprint('outputDir', self.outputDir)
        self.varprint('write_selection', self.write_selection)

        self.varprint('add_symmetric_neighborhood', self.add_symmetric_neighborhood)
        self.varprint('differentiate_other_half', self.differentiate_other_half)
        self.varprint('use_common_neighborhood', self.use_common_neighborhood)
        self.varprint('name_delay_from_division', self.name_delay_from_division)
        self.varprint('confidence_delay_from_division', self.confidence_delay_from_division)
        self.varprint('self.confidence_atlases_nmin', self.confidence_atlases_nmin)
        self.varprint('confidence_atlases_percentage', self.confidence_atlases_percentage)

        self.varprint('division_contact_similarity', self.division_contact_similarity)

        self.varprint('diagnosis_properties', self.diagnosis_properties)

        self.varprint('division_permutation_proposal', self.division_permutation_proposal)

        self.varprint('dendrogram_cluster_distance', self.dendrogram_cluster_distance)

        self.varprint('generate_figure', self.generate_figure)
        self.varprint('figurefile_suffix', self.figurefile_suffix)

        self.varprint('cells_to_be_traced', self.cells_to_be_traced)
        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# CellAtlasParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        udiagnosis.DiagnosisParameters.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'atlasFiles', self.atlasFiles, self.doc.get('atlasFiles', None))
        self.varwrite(logfile, 'referenceAtlas', self.referenceAtlas, self.doc.get('referenceAtlas', None))

        self.varwrite(logfile, 'outputDir', self.outputDir, self.doc.get('outputDir', None))
        self.varwrite(logfile, 'write_selection', self.write_selection, self.doc.get('write_selection', None))

        self.varwrite(logfile, 'add_symmetric_neighborhood', self.add_symmetric_neighborhood,
                      self.doc.get('add_symmetric_neighborhood', None))
        self.varwrite(logfile, 'differentiate_other_half', self.differentiate_other_half,
                      self.doc.get('differentiate_other_half', None))
        self.varwrite(logfile, 'use_common_neighborhood', self.use_common_neighborhood,
                      self.doc.get('use_common_neighborhood', None))
        self.varwrite(logfile, 'name_delay_from_division', self.name_delay_from_division,
                      self.doc.get('name_delay_from_division', None))
        self.varwrite(logfile, 'confidence_delay_from_division', self.confidence_delay_from_division,
                      self.doc.get('confidence_delay_from_division', None))
        self.varwrite(logfile, 'confidence_atlases_nmin', self.confidence_atlases_nmin,
                      self.doc.get('confidence_atlases_nmin', None))
        self.varwrite(logfile, 'confidence_atlases_percentage', self.confidence_atlases_percentage,
                      self.doc.get('confidence_atlases_percentage', None))

        self.varwrite(logfile, 'division_contact_similarity', self.division_contact_similarity,
                      self.doc.get('division_contact_similarity', None))

        self.varwrite(logfile, 'diagnosis_properties', self.diagnosis_properties,
                      self.doc.get('diagnosis_properties', None))

        self.varwrite(logfile, 'division_permutation_proposal', self.division_permutation_proposal,
                      self.doc.get('division_permutation_proposal', None))

        self.varwrite(logfile, 'dendrogram_cluster_distance', self.dendrogram_cluster_distance,
                      self.doc.get('dendrogram_cluster_distance', None))

        self.varwrite(logfile, 'generate_figure', self.generate_figure, self.doc.get('generate_figure', None))
        self.varwrite(logfile, 'figurefile_suffix', self.figurefile_suffix, self.doc.get('figurefile_suffix', None))

        self.varwrite(logfile, 'cells_to_be_traced', self.cells_to_be_traced, self.doc.get('cells_to_be_traced', None))

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
        self.write_selection = self.read_parameter(parameters, 'write_selection', self.write_selection)

        self.add_symmetric_neighborhood = self.read_parameter(parameters, 'add_symmetric_neighborhood',
                                                              self.add_symmetric_neighborhood)
        self.differentiate_other_half = self.read_parameter(parameters, 'differentiate_other_half',
                                                            self.differentiate_other_half)
        self.use_common_neighborhood = self.read_parameter(parameters, 'use_common_neighborhood',
                                                           self.use_common_neighborhood)
        self.name_delay_from_division = self.read_parameter(parameters, 'name_delay_from_division',
                                                            self.name_delay_from_division)
        self.name_delay_from_division = self.read_parameter(parameters, 'delay_from_division',
                                                            self.name_delay_from_division)
        self.confidence_delay_from_division = self.read_parameter(parameters, 'confidence_delay_from_division',
                                                                  self.confidence_delay_from_division)
        self.confidence_delay_from_division = self.read_parameter(parameters, 'delay_from_division',
                                                                  self.confidence_delay_from_division)
        self.confidence_atlases_nmin = self.read_parameter(parameters, 'confidence_atlases_nmin',
                                                           self.confidence_atlases_nmin)
        self.confidence_atlases_percentage = self.read_parameter(parameters, 'confidence_atlases_percentage',
                                                                 self.confidence_atlases_percentage)

        self.division_contact_similarity = self.read_parameter(parameters, 'cell_contact_distance',
                                                               self.division_contact_similarity)
        self.division_contact_similarity = self.read_parameter(parameters, 'division_contact_similarity',
                                                               self.division_contact_similarity)

        self.diagnosis_properties = self.read_parameter(parameters, 'diagnosis_properties', self.diagnosis_properties)
        self.diagnosis_properties = self.read_parameter(parameters, 'naming_diagnosis', self.diagnosis_properties)
        self.diagnosis_properties = self.read_parameter(parameters, 'diagnosis_naming', self.diagnosis_properties)

        self.division_permutation_proposal = self.read_parameter(parameters, 'division_permutation_proposal',
                                                                 self.division_permutation_proposal)
        self.division_permutation_proposal = self.read_parameter(parameters, 'daughter_switch_proposal',
                                                                 self.division_permutation_proposal)

        self.dendrogram_cluster_distance = self.read_parameter(parameters, 'dendrogram_cluster_distance',
                                                               self.dendrogram_cluster_distance)

        self.generate_figure = self.read_parameter(parameters, 'generate_figure', self.generate_figure)
        self.figurefile_suffix = self.read_parameter(parameters, 'figurefile_suffix', self.figurefile_suffix)

        self.cells_to_be_traced = self.read_parameter(parameters, 'cells_to_be_traced', self.cells_to_be_traced)

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
# utilities on cell names
#
########################################################################################

def get_symmetric_neighborhood(neighborhood):
    """
    Changes the names of the cells in the neighborhood to get the symmetric neighborhood
    :param neighborhood:
    :return:
    """
    symneighborhood = {}
    for n in neighborhood:
        if n == 'background':
            sn = 'background'
        elif n == 'other-half':
            sn = 'other-half'
        else:
            sn = uname.get_symmetric_name(n)
        symneighborhood[sn] = neighborhood[n]
    return symneighborhood


########################################################################################
#
#
#
########################################################################################

def _write_list(listtobeprinted, firstheader="", otherheader="", maxlength=112, verboseness=0):
    txt = ""
    n = 0
    for i, item in enumerate(listtobeprinted):
        if i == 0:
            txt = firstheader
            n = 0
        if len(txt) + len(str(item)) <= maxlength:
            if n >= 1:
                txt += ","
            txt += " " + str(item)
            n += 1
        else:
            monitoring.to_log_and_console(txt, verboseness=verboseness)
            txt = otherheader + " " + str(item)
            n = 1
        if i == len(listtobeprinted) - 1:
            monitoring.to_log_and_console(txt, verboseness=verboseness)


def _write_summary_pairwise_switches(atlases, summary):

    divisions = atlases.get_divisions()
    percents = []
    mother_names = list(summary.keys())

    for n in mother_names:
        percents.append(100.0 * float(len(summary[n]['disagreement'])) / float(summary[n]['tested_couples']))
    [sorted_percents, sorted_mothers] = list(zip(*sorted(zip(percents, mother_names), reverse=True)))

    majority = {}
    equality = {}
    for n in sorted_mothers:
        majority[n] = {}
        equality[n] = {}
        if len(summary[n]['disagreement']) == 0:
            continue
        msg = " - " + str(n) + " cell division into "
        msg += str(uname.get_daughter_names(n)) + " has " + str(len(summary[n]['disagreement']))
        if len(summary[n]['disagreement']) > 1:
            msg += " disagreements"
        else:
            msg += " disagreement"
        percent = 100.0 * float(len(summary[n]['disagreement'])) / float(summary[n]['tested_couples'])
        msg += " (" + "{:2.2f}%".format(percent) + ")"
        monitoring.to_log_and_console(msg)
        msg = "\t over " + str(summary[n]['tested_couples']) + " tested configurations "
        msg += "and over " + str(len(divisions[n]))
        msg += " references: "
        monitoring.to_log_and_console(msg)
        #
        # print references
        #
        _write_list(sorted(divisions[n]), firstheader="\t     ", otherheader="\t     ", maxlength=112, verboseness=0)
        #
        # count the cases where one atlas disagrees
        #
        for a in divisions[n]:
            s = 0
            for pair in summary[n]['disagreement']:
                if a in pair:
                    s += 1
            if 2*s > len(divisions[n]):
                majority[n][a] = s
            elif 2*s == len(divisions[n]):
                equality[n][a] = s
        #
        # print detailed disagreements
        #
        _write_list(sorted(summary[n]['disagreement']), firstheader="\t - disagreement list:", otherheader="\t     ",
                    maxlength=112, verboseness=3)

    nitems = 0
    for n in sorted_mothers:
        nitems += len(majority[n]) + len(equality[n])
    if nitems == 0:
        return
    monitoring.to_log_and_console("")
    monitoring.to_log_and_console(" --- atlases pairwise disagreements: summary ---")
    for n in sorted(mother_names):
        if len(majority[n]) > 0:
            msg = " - " + str(n) + " division, atlas that mostly disagrees:"
            akeys = sorted(list(majority[n].keys()))
            for i, a in enumerate(akeys):
                msg += " " + str(a) + " (" + str(majority[n][a]) + "/" + str(len(divisions[n])) + ")"
                if i < len(akeys) - 1:
                    msg += ","
            monitoring.to_log_and_console(msg)
    for n in sorted(mother_names):
        if len(equality[n]) > 0:
            msg = " - " + str(n) + " division, atlas that equally disagrees:"
            akeys = sorted(list(equality[n].keys()))
            for i, a in enumerate(akeys):
                msg += " " + str(a) + " (" + str(equality[n][a]) + "/" + str(len(divisions[n])) + ")"
                if i < len(akeys) - 1:
                    msg += ","
            monitoring.to_log_and_console(msg)


def _diagnosis_pairwise_switches(atlases, parameters):
    proc = "_diagnosis_pairwise_switches"

    divisions = atlases.get_divisions()
    ccs = not parameters.use_common_neighborhood
    neighborhoods = atlases.get_neighborhoods(delay_from_division=parameters.name_delay_from_division)
    similarity = atlases.get_division_contact_similarity()

    summary = {}

    for n in divisions:
        #
        # only one reference/atlas for mother cell 'n': nothing to do
        #
        if len(divisions[n]) <= 1:
            continue
        summary[n] = {}
        summary[n]['tested_couples'] = 0
        summary[n]['disagreement'] = []

        d = uname.get_daughter_names(n)
        for r1 in divisions[n]:
            for r2 in divisions[n]:
                if r2 <= r1:
                    continue
                summary[n]['tested_couples'] += 1
                #
                # test reference r1 versus r2
                # it is assumed that (r1, switched(r2)) is similar to (switched(r1), r2)
                # so only (r1, switched(r2)) is tested
                #
                switch_neighs = {d[0]: copy.deepcopy(neighborhoods[d[1]][r2]),
                                 d[1]: copy.deepcopy(neighborhoods[d[0]][r2])}
                if d[0] in switch_neighs[d[0]]:
                    switch_neighs[d[0]][d[1]] = switch_neighs[d[0]][d[0]]
                    del switch_neighs[d[0]][d[0]]
                if d[1] in switch_neighs[d[1]]:
                    switch_neighs[d[1]][d[0]] = switch_neighs[d[1]][d[1]]
                    del switch_neighs[d[1]][d[1]]
                #
                # same
                #
                same_dist = division_contact_generic_distance(atlases, neighborhoods[d[0]][r1], neighborhoods[d[1]][r1],
                                                              neighborhoods[d[0]][r2], neighborhoods[d[1]][r2],
                                                              similarity=similarity, change_contact_surfaces=ccs)
                swit_dist = division_contact_generic_distance(atlases, neighborhoods[d[0]][r1], neighborhoods[d[1]][r1],
                                                              switch_neighs[d[0]], switch_neighs[d[1]],
                                                              similarity=similarity, change_contact_surfaces=ccs)
                if same_dist < swit_dist:
                    continue
                summary[n]['disagreement'] += [(r1, r2)]

    divisions_with_disagreement = [n for n in summary if len(summary[n]['disagreement']) > 0]

    msg = "tested divisions = " + str(len(summary))
    monitoring.to_log_and_console(str(proc) + ": " + msg)
    msg = "divisions with pairwise disagreement =  " + str(len(divisions_with_disagreement))
    monitoring.to_log_and_console("\t " + msg)
    msg = "  A disagreement means that a division from a reference is closer to the\n"
    msg += "  switched division of an other reference than the division itself."
    monitoring.to_log_and_console(msg)
    monitoring.to_log_and_console("")
    if len(divisions_with_disagreement) > 0:
        _write_summary_pairwise_switches(atlases, summary)
    monitoring.to_log_and_console("")

    return summary


########################################################################################
#
#
#
########################################################################################

def _dpp_switch_contact_surfaces(neighbors, reference, daughters):
    """
    Switch contact surfaces for the two daughters and atlas 'reference'.
    Parameters
    ----------
    neighbors
    reference
    daughters

    Returns
    -------

    """
    #
    # contact surfaces of daughters[0] for atlas 'reference'
    # replace contact surface with daughters[1] with a contact surface with daughters[0]
    #
    neighs = {0: copy.deepcopy(neighbors[1][reference]), 1: copy.deepcopy(neighbors[0][reference])}
    if daughters[0] in neighs[0]:
        neighs[0][daughters[1]] = neighs[0][daughters[0]]
        del neighs[0][daughters[0]]
    if daughters[1] in neighs[1]:
        neighs[1][daughters[0]] = neighs[1][daughters[1]]
        del neighs[1][daughters[1]]

    neighbors[0][reference] = neighs[0]
    neighbors[1][reference] = neighs[1]
    return neighbors


def _dpp_global_generic_distance(neighbors, references, atlases, parameters, debug=False):
    """
    Compute a global score. The global score is the average of local similarities over all
    couples of references/atlases.

    Parameters
    ----------
    neighbors: neighborhoods for the two daughters (dictionary indexed by [0,1] then by the references
    references: set of references
    debug

    Returns
    -------

    """

    similarity = atlases.get_division_contact_similarity()
    ccs = not parameters.use_common_neighborhood
    score = 0
    n = 0
    distances = {}
    for r1 in references:
        if debug:
            distances[r1] = {}
        for r2 in references:
            if r2 <= r1:
                continue
            dist = division_contact_generic_distance(atlases, neighbors[0][r1], neighbors[1][r1], neighbors[0][r2],
                                                     neighbors[1][r2], similarity=similarity,
                                                     change_contact_surfaces=ccs)
            if debug:
                distances[r1][r2] = dist
            score += dist
            n += 1
    if debug:
        print("---- _dpp_global_generic_distance")
        refs1 = distances.keys()
        refs1 = sorted(refs1)
        for r1 in refs1:
            refs2 = distances[r1].keys()
            refs2 = sorted(refs2)
            for r2 in refs2:
                print("   - dist[" + str(r1) + ", " + str(r2) + "] = " + str(distances[r1][r2]))
    return score / n


def _dpp_test_one_division(atlases, mother, parameters):
    """
    Test whether any daughter switch (for a given reference) improve a global score
    Parameters
    ----------
    atlases
    mother
    parameters

    Returns
    -------

    """
    divisions = atlases.get_divisions()
    if len(divisions[mother]) <= 1:
        return {}, []

    daughters = uname.get_daughter_names(mother)
    neighborhoods = atlases.get_neighborhoods(delay_from_division=parameters.name_delay_from_division)
    neighbors = {0: copy.deepcopy(neighborhoods[daughters[0]]), 1: copy.deepcopy(neighborhoods[daughters[1]])}

    # score before any changes
    debug = False
    if debug:
        print("")
        print("===== test division " + str(mother) + " : " + str(divisions[mother]))
    score = _dpp_global_generic_distance(neighbors, divisions[mother], atlases, parameters, debug=debug)

    returned_scores = [(None, score)]
    corrections = []
    i = 1
    while True:
        newscore = {}
        for r in sorted(divisions[mother]):
            #
            # switch contact surfaces for the daughters in atlas 'r'
            #
            tmp = copy.deepcopy(neighbors)
            tmp = _dpp_switch_contact_surfaces(tmp, r, daughters)
            if debug:
                print("===== test switch " + str(mother) + " / " + str(r))
            # compute a new score, keep it if it better than the one before any changes
            newscore[r] = _dpp_global_generic_distance(tmp, divisions[mother], atlases, parameters, debug=debug)
            if debug:
                print("     new score = " + str(newscore[r]) + " - original score = " + str(score))
            if newscore[r] > score:
                del newscore[r]
        # no found correction at this iteration
        if len(newscore) == 0:
            return corrections, returned_scores
        # found several correction
        # 1. pick the only one (if only one is found)
        # 2. or pick the one with maximal score change
        elif len(newscore) == 1:
            ref = list(newscore.keys())[0]
        else:
            ref = min(newscore, key=lambda key: newscore[key])
        # first iteration, keep the value of the global score decrease
        if i == 1:
            for r in newscore:
                returned_scores += [(r, score - newscore[r])]
        corrections += [(ref, score - newscore[ref])]
        i += 1
        # if one correction has been found, apply it
        # and look for an other additional correction
        tmp = copy.deepcopy(neighbors)
        tmp = _dpp_switch_contact_surfaces(tmp, ref, daughters)
        neighbors[0] = copy.deepcopy(tmp[0])
        neighbors[1] = copy.deepcopy(tmp[1])
        score = newscore[ref]


def division_permutation_proposal(atlases, parameters):

    # neighborhoods is a dictionary of dictionaries
    # ['cell name']['reference name']
    # first key is a cell name (daughter cell)
    # second key is the reference from which the neighborhood has been extracted

    #
    # mother cell name dictionary indexed by stage
    # stage 6: 32 cells
    # stage 7: 64 cells
    #

    proc = "division_permutation_proposal"

    divisions = atlases.get_divisions()
    mothers = {}
    for n in divisions:
        stage = n.split('.')[0][1:]
        mothers[stage] = mothers.get(stage, []) + [n]
    for s in mothers:
        mothers[s] = sorted(mothers[s])

    stages = list(mothers.keys())
    stages.sort()
    corrections = {}
    selection = {}

    for s in stages:
        corrections[s] = {}
        # if int(s) != 7:
        #    continue
        for m in mothers[s]:
            # if m != 'a7.0002_':
            #     continue
            correction, returned_score = _dpp_test_one_division(atlases, m, parameters)
            #
            # correction is a dictionary indexed by the iteration index
            # each value is a tuple ('atlas name', score increment)
            #
            if len(correction) > 0:
                corrections[s][m] = correction
                selection[m] = returned_score

    #
    # build output selections
    #
    ref_atlases = atlases.get_atlases()
    output_selections = atlases.get_output_selections()

    for m in selection:
        if len(selection[m]) <= 1:
            continue
        (a, score) = selection[m][0]

        for i, (ref, ds) in enumerate(selection[m]):

            if i == 0:
                continue

            # check if the reference is registered
            # discard symmetrical neighborhood for warning
            if ref not in ref_atlases:
                if not (ref[:4] == 'sym-' and ref[4:] in ref_atlases):
                    monitoring.to_log_and_console(proc + ": weird, '" + str(ref) + "' is not in reference atlases.", 4)
                continue

            keyscore = "morphonet_float_" + str(ref) + "_distance_average_before_permutation_proposal"
            keydecre = "morphonet_float_" + str(ref) + "_distance_decrement_percentage_after_permutation_proposal"
            output_selections[keyscore] = output_selections.get(keyscore, {})
            output_selections[keydecre] = output_selections.get(keydecre, {})

            lineage = ref_atlases[ref]['cell_lineage']
            name = ref_atlases[ref]['cell_name']
            cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))

            for c in cells:
                if c not in name:
                    continue
                if name[c] == m:
                    output_selections[keyscore][c] = score
                    output_selections[keydecre][c] = ds / score

    #
    # reporting
    #
    monitoring.to_log_and_console("====== division permutation proposal =====")
    monitoring.to_log_and_console("------ cell-based view")
    corrections_by_atlas = {}
    for s in stages:
        if len(corrections[s]) == 0:
            continue
        msg = "  - generation " + str(s)
        monitoring.to_log_and_console(msg)
        mothers = list(corrections[s].keys())
        mothers.sort()
        for m in mothers:
            if len(corrections[s][m]) == 0:
                continue
            tmp = [c[0] for c in corrections[s][m]]
            tmp.sort()
            for a in tmp:
                corrections_by_atlas[a] = corrections_by_atlas.get(a, []) + [m]
            msg = "    - division of '" + str(m) + "': "
            for i, r in enumerate(tmp):
                msg += str(r)
                if i < len(tmp)-1:
                    msg += ", "
            monitoring.to_log_and_console(msg)
    if len(corrections_by_atlas) > 0:
        monitoring.to_log_and_console("------ atlas-based view")
        refs = list(corrections_by_atlas.keys())
        refs.sort()
        for r in refs:
            msg = "  - reference '" + str(r) + "': " + str(corrections_by_atlas[r])
            monitoring.to_log_and_console(msg)

    if len(selection) > 0:
        quadruplet = []
        for m in selection:
            (a, score) = selection[m][0]
            if len(selection[m]) <= 1:
                continue
            for i, (ref, ds) in enumerate(selection[m]):
                if i == 0:
                    continue
                quadruplet += [(m, ref, score, 100.0 * ds / score)]
        quadruplet = sorted(quadruplet, key=operator.itemgetter(1))
        quadruplet = sorted(quadruplet, key=operator.itemgetter(3), reverse=True)
        monitoring.to_log_and_console("------ average distance percentage decrease view")
        for q in quadruplet:
            msg = "  - division of '" + str(q[0]) + "' in '" + str(q[1]) + "': "
            msg += "{:2.2f}% decrease of {:1.2f} average distance".format(q[3], q[2])
            monitoring.to_log_and_console(msg)

    monitoring.to_log_and_console("==========================================")

########################################################################################
#
#
#
########################################################################################


def call_to_scipy_linkage(atlases, config, cluster_distance='single', change_contact_surfaces=True):
    """

    Parameters
    ----------
    atlases
    config: dictionary of dictionary of neighborhoods indexed by [reference] then by [0,1],
        where 0 stands for one daughter, and 1 for the other
    cluster_distance
    change_contact_surfaces

    Returns
    -------
    conddist: the squareform vector built from the distance matrice
       (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html)
    z: the hierarchical clustering encoded as a linkage matrix
       (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)
    labels: the list of atlas names

    """
    similarity = atlases.get_division_contact_similarity()

    labels = []
    #
    # build a square matrix of distances
    #
    dist = np.zeros((len(config), len(config)))
    for i, r in enumerate(config):
        labels += [r]
        for j, s in enumerate(config):
            if r == s:
                dist[i][i] = 0.0
                continue
            if r > s:
                continue
            # if r == 'switched-' + str(s) or s == 'switched-' + str(r):
            #    continue
            dist[i][j] = 100.0 * division_contact_generic_distance(atlases, config[r][0], config[r][1], config[s][0],
                                                                   config[s][1], similarity=similarity,
                                                                   change_contact_surfaces=change_contact_surfaces)
            dist[j][i] = dist[i][j]

    conddist = sp.spatial.distance.squareform(dist)
    z = sch.linkage(conddist, method=cluster_distance)

    return conddist, z, labels


def _diagnosis_linkage(atlases, parameters):
    proc = "_diagnosis_linkage"
    divisions = atlases.get_divisions()
    ccs = not parameters.use_common_neighborhood
    cluster_distance = parameters.dendrogram_cluster_distance

    ref_atlases = atlases.get_atlases()
    output_selections = atlases.get_output_selections()

    merge_values = {}
    lastmerge_values = {}

    swmerge_values = {}
    swlastmerge_values = {}

    division_lastmerge_values = {}

    for n in divisions:
        stage = n.split('.')[0][1:]
        if len(divisions[n]) <= 2:
            continue
        #
        #
        #
        config = atlases.extract_division_neighborhoods(n, delay_from_division=parameters.name_delay_from_division)
        swconfig = switched_division_neighborhoods(config, n)

        #
        # distance array for couples of atlases/references
        #
        conddist, z, labels = call_to_scipy_linkage(atlases, config, cluster_distance=cluster_distance,
                                                    change_contact_surfaces=ccs)

        merge_values[stage] = merge_values.get(stage, []) + list(z[:, 2])
        lastmerge_value = z[:, 2][-1]
        lastmerge_values[stage] = lastmerge_values.get(stage, []) + [lastmerge_value]

        #
        # set the lastmerge_value in morphonet selection
        #
        for r in divisions[n]:
            if r not in ref_atlases:
                if not (r[:4] == 'sym-' and r[4:] in ref_atlases):
                    monitoring.to_log_and_console(proc + ": weird, '" + str(r) + "' is not in reference atlases.", 4)
                continue
            keyselection = "morphonet_float_" + str(r) + "_last_dendrogram_value"
            output_selections[keyselection] = output_selections.get(keyselection, {})
            lineage = ref_atlases[r]['cell_lineage']
            name = ref_atlases[r]['cell_name']
            cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
            for c in cells:
                if c not in name:
                    continue
                if name[c] == n:
                    output_selections[keyselection][c] = round(lastmerge_value)

        #
        # distance array for couples of atlases/references plus the switched ones
        #
        swconddist, swz, swlabels = call_to_scipy_linkage(atlases, swconfig, cluster_distance=cluster_distance,
                                                          change_contact_surfaces=ccs)

        swmerge_values[stage] = swmerge_values.get(stage, []) + list(swz[:, 2])
        swlastmerge_value = swz[:, 2][-1]
        swlastmerge_values[stage] = swlastmerge_values.get(stage, []) + [swlastmerge_value]

        division_lastmerge_values[n] = [lastmerge_value, swlastmerge_value]

    nochanges = {n: v for n, v in division_lastmerge_values.items() if v[1] <= v[0]}
    mother_with_nochanges = [n for n, v in division_lastmerge_values.items() if v[1] <= v[0]]
    #
    # set the lastmerge_value in morphonet selection
    #
    for n in mother_with_nochanges:
        for r in divisions[n]:
            if r not in ref_atlases:
                if not (r[:4] == 'sym-' and r[4:] in ref_atlases):
                    monitoring.to_log_and_console(proc + ": weird, '" + str(r) + "' is not in reference atlases.", 4)
                continue
            keyselection = "morphonet_selection_" + str(r) + "_dendrogram_warning"
            output_selections[keyselection] = output_selections.get(keyselection, {})
            lineage = ref_atlases[r]['cell_lineage']
            name = ref_atlases[r]['cell_name']
            cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
            for c in cells:
                if c not in name:
                    continue
                if name[c] == n:
                    output_selections[keyselection][c] = 100

    monitoring.to_log_and_console("------ division with same dendrogram last values (without and with switch)")
    msg = str(len(nochanges)) + "/" + str(len(division_lastmerge_values)) + " divisions"
    monitoring.to_log_and_console("\t " + msg)

    monitoring.to_log_and_console("------ cell-based view")
    division_by_generation = {}
    for m in mother_with_nochanges:
        g = m.split('.')[0][1:]
        division_by_generation[g] = division_by_generation.get(g, []) + [m]
    for g in division_by_generation:
        monitoring.to_log_and_console("  - generation " + str(g))
        mothers = list(division_by_generation[g])
        mothers.sort()
        for m in mothers:
            msg = "    - division of '" + str(m) + "': " + str(division_lastmerge_values[m])
            monitoring.to_log_and_console(msg)

    monitoring.to_log_and_console("------ dendrogram last-value view")
    sorted_nochanges = sorted(nochanges.items(), key=lambda v: v[1][0], reverse=True)
    for s in sorted_nochanges:
        msg = "    - division of '" + str(s[0]) + "': " + str(s[1])
        monitoring.to_log_and_console(msg)
    monitoring.to_log_and_console("")


########################################################################################
#
#
#
########################################################################################

def _build_common_neighborhoods(neighborhoods):
    """
    Build a new neighborhood dictionary where both a cell and its sister have
    the same neighbors.
    Parameters
    ----------
    neighborhoods: dictionary of dictionaries
        ['cell name']['reference name']['neighboring cell']
        first key is a cell name (daughter cell)
        second key is the reference from which the neighborhood has been extracted
        a neighborhood itself is a dictionary indexed by the neighboring cell names

    Returns
    -------

    """
    proc = "_build_common_neighborhoods"

    common_neighborhoods = {}

    for cell in neighborhoods:
        if cell in common_neighborhoods:
            continue
        sister = uname.get_sister_name(cell)
        if sister in common_neighborhoods:
            msg = "weird, '" + str(sister) + "' is in neighborhoods while '" + str(cell) + "' is not"
            monitoring.to_log_and_console(proc + ": " + msg)
        new_neighborhoods = ucontact.build_same_contact_surfaces(neighborhoods, [cell, sister])
        for n in new_neighborhoods:
            common_neighborhoods[n] = copy.deepcopy(new_neighborhoods[n])
    return common_neighborhoods


def switched_division_neighborhoods(config, mother_name):
    daughters = uname.get_daughter_names(mother_name)
    swconfig = {}
    for r in config:
        swconfig[r] = {}
        swconfig[r][0] = copy.deepcopy(config[r][0])
        swconfig[r][1] = copy.deepcopy(config[r][1])
        sr = 'switched-' + str(r)
        swconfig[sr] = {}
        swconfig[sr][0] = copy.deepcopy(config[r][1])
        swconfig[sr][1] = copy.deepcopy(config[r][0])
        if daughters[1] in swconfig[sr][0] and swconfig[sr][0][daughters[1]] > 0:
            msg = "  weird, " + str(daughters[1]) + " was found in its neighborhood for reference " + str(r)
            monitoring.to_log_and_console("      " + msg)
        if daughters[0] in swconfig[sr][0]:
            swconfig[sr][0][daughters[1]] = swconfig[sr][0][daughters[0]]
            del swconfig[sr][0][daughters[0]]
        if daughters[0] in swconfig[sr][1] and swconfig[sr][1][daughters[0]] > 0:
            msg = "  weird, " + str(daughters[0]) + " was found in its neighborhood for reference " + str(r)
            monitoring.to_log_and_console("      " + msg)
        if daughters[1] in swconfig[sr][1]:
            swconfig[sr][1][daughters[0]] = swconfig[sr][1][daughters[1]]
            del swconfig[sr][1][daughters[1]]
    return swconfig


########################################################################################
#
#
#
########################################################################################

def _get_branch_length(cell, lineage):
    length = 0
    c = cell
    while c in lineage and len(lineage[c]) == 1:
        length += 1
        c = lineage[c][0]
    return length


class Atlases(object):
    def __init__(self, parameters=None):

        self.cell_contact_distance = 'l1_distance'
        self.division_contact_similarity = 'distance'

        # reference atlas for time alignment
        self._ref_atlas = None

        # partial copy of the read atlases, required when figures have to be generated
        self._atlases = {}

        self._default_delay = None

        # nested dictionary of neighborhoods, where the keys are [delay_from_division]['cell name']['reference name']
        # where 'cell name' is the cell name (Conklin), and 'reference name' is the file name,
        # a neighborhood is a dictionary of contact surfaces indexed by cell names
        # it only considers the first time point after the division
        self._neighborhoods = {}

        # nested dictionary of volumes, where the keys are [delay_from_division]['cell name']['reference name']
        self._volumes = {}

        self._use_common_neighborhood = False

        # dictionary indexed by 'cell name' where 'cell name' is a mother cell giving the list
        # of references/atlases available for the two daughters
        self._divisions = {}
        # dictionary indexed by [delay_from_division]['cell name']
        # values are the pair of daughter cells (the largest first)
        # this is the set of divisions where the same daughter is always larger than the other
        # (and there are at least 5 atlases)
        self._unequal_divisions = {}

        # dictionary index by atlas name
        # to keep trace of some output (kept as morphonet selection)
        self._output_selections = {}

        if parameters is not None:
            self.update_from_parameters(parameters)

    ############################################################
    #
    # getters
    #
    ############################################################

    def get_cell_contact_distance(self):
        return self.cell_contact_distance

    def get_division_contact_similarity(self):
        return self.division_contact_similarity

    def get_default_delay(self):
        return self._default_delay

    def set_default_delay(self, delay=0):
        self._default_delay = delay

    def get_neighborhoods(self, delay_from_division=None):
        if delay_from_division is None:
            delay = self.get_default_delay()
        else:
            delay = delay_from_division
        if delay not in self._neighborhoods:
            self._neighborhoods[delay] = {}
        return self._neighborhoods[delay]

    def set_neighborhoods(self, neighborhoods, delay_from_division=0):
        self._neighborhoods[delay_from_division] = neighborhoods

    def get_volumes(self, delay_from_division=None):
        if delay_from_division is None:
            delay = self.get_default_delay()
        else:
            delay = delay_from_division
        if delay not in self._volumes:
            self._volumes[delay] = {}
        return self._volumes[delay]

    def set_volumes(self, volumes, delay_from_division=0):
        self._volumes[delay_from_division] = volumes

    def get_atlases(self):
        return self._atlases

    def set_atlases(self, atlases):
        self._atlases = atlases

    def get_use_common_neighborhood(self):
        return self._use_common_neighborhood

    def get_divisions(self):
        return self._divisions

    def get_unequal_divisions(self, delay_from_division=None):
        if delay_from_division is None:
            delay = self.get_default_delay()
        else:
            delay = delay_from_division
        if delay not in self._unequal_divisions:
            self._unequal_divisions[delay] = {}
        return self._unequal_divisions[delay]

    def get_output_selections(self):
        return self._output_selections

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameters):
        proc = "update_from_parameters"

        if not isinstance(parameters, ucontact.ContactSurfaceParameters):
            monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                          + str(type(parameters)))

        if parameters.cell_contact_distance.lower() == 'l1_distance' \
                or parameters.cell_contact_distance.lower() == 'l1-distance' \
                or parameters.cell_contact_distance.lower() == 'l1_norm' \
                or parameters.cell_contact_distance.lower() == 'l1-norm':
            self.cell_contact_distance = 'l1_distance'
        elif parameters.cell_contact_distance.lower() == 'l2_distance' \
                or parameters.cell_contact_distance.lower() == 'l2-distance' \
                or parameters.cell_contact_distance.lower() == 'l2_norm' \
                or parameters.cell_contact_distance.lower() == 'l2-norm':
            self.cell_contact_distance = 'l2_distance'
        else:
            monitoring.to_log_and_console(str(proc) + ": unhandled cell contact distance: '" +
                                          str(parameters.cell_contact_distance) + "'")

        if not isinstance(parameters, AtlasParameters):
            monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                          + str(type(parameters)))

        if parameters.referenceAtlas is not None:
            name = parameters.referenceAtlas.split(os.path.sep)[-1]
            if name.endswith(".xml") or name.endswith(".pkl"):
                name = name[:-4]
            self._ref_atlas = name

        if parameters.division_contact_similarity.lower() == 'distance' \
                or parameters.division_contact_similarity.lower() == 'norm' \
                or parameters.division_contact_similarity.lower() == 'l1_distance' \
                or parameters.division_contact_similarity.lower() == 'l1-distance' \
                or parameters.division_contact_similarity.lower() == 'l1_norm' \
                or parameters.division_contact_similarity.lower() == 'l1-norm' \
                or parameters.division_contact_similarity.lower() == 'l2_distance' \
                or parameters.division_contact_similarity.lower() == 'l2-distance' \
                or parameters.division_contact_similarity.lower() == 'l2_norm' \
                or parameters.division_contact_similarity.lower() == 'l2-norm':
            self.division_contact_similarity = 'distance'
        else:
            monitoring.to_log_and_console(str(proc) + ": unhandled division contact similarity: '" +
                                          str(parameters.division_contact_similarity) + "'")

    ############################################################
    #
    #
    #
    ############################################################

    def build_divisions(self, delay_from_division=None):
        """
        Build a dictionary index by mother cell name. Each entry contains reference names
        for which the both daughters exist
        Returns
        -------

        """
        proc = "build_divisions"

        if self._divisions is not None:
            del self._divisions
            self._divisions = {}

        #
        # get all references per division/mother cells
        #
        neighborhoods = self.get_neighborhoods(delay_from_division=delay_from_division)
        cell_names = sorted(list(neighborhoods.keys()))
        references = {}
        for cell_name in cell_names:
            mother_name = uname.get_mother_name(cell_name)
            references[mother_name] = references.get(mother_name, set()).union(set(neighborhoods[cell_name].keys()))

        #
        # remove references that does not exist for one daughter
        #
        mother_names = sorted(references.keys())
        for n in mother_names:
            daughters = uname.get_daughter_names(n)
            #
            # check whether each reference has the two daughters
            #
            refs = list(references[n])
            for r in refs:
                if daughters[0] in neighborhoods and r in neighborhoods[daughters[0]] and \
                        daughters[1] in neighborhoods and r in neighborhoods[daughters[1]]:
                    self._divisions[n] = self._divisions.get(n, []) + [r]
                else:
                    msg = "    " + str(proc) + ": remove atlas '" + str(r) + "' for division '" + str(n) + "'"
                    monitoring.to_log_and_console(msg)

    def build_unequal_divisions(self, delay_from_division=None):
        """
        Build a dictionary index by mother cell name. Each entry contains reference names
        for which the both daughters exist
        Returns
        -------

        """
        proc = "build_unequal_divisions"
        minimal_references = 5

        delay = delay_from_division
        if delay is None:
            delay = self.get_default_delay()

        if self._unequal_divisions is None:
            self._unequal_divisions = {}
        if delay in self._unequal_divisions:
            del self._unequal_divisions[delay]
        self._unequal_divisions[delay] = {}

        divisions = self.get_divisions()

        volumes = self.get_volumes(delay_from_division=delay)

        for mother in divisions:
            d = uname.get_daughter_names(mother)
            vol01 = 0
            vol10 = 0
            for r in divisions[mother]:
                if volumes[d[0]][r] > volumes[d[1]][r]:
                    vol01 += 1
                elif volumes[d[0]][r] < volumes[d[1]][r]:
                    vol10 += 1
            if vol01 > minimal_references and vol10 == 0:
                self._unequal_divisions[delay][mother] = [d[0], d[1]]
            elif vol10 > minimal_references and vol01 == 0:
                self._unequal_divisions[delay][mother] = [d[1], d[0]]

        msg = "found " + str(len(self._unequal_divisions[delay])) + " unequal divisions "
        msg += "at delay = " + str(delay) + " "
        msg += "(more than " + str(minimal_references) + " atlases) in " + str(len(divisions)) + " divisions"
        monitoring.to_log_and_console("\t" + msg)
        return

    ############################################################
    #
    #
    #
    ############################################################

    def temporal_alignment(self, time_digits_for_cell_id=4):

        all_atlases = self.get_atlases()
        ref_lineage = all_atlases[self._ref_atlas]['cell_lineage']
        for n in all_atlases:
            if n == self._ref_atlas:
                all_atlases[n]['temporal_alignment'] = (1.0, 0.0)
                continue
            lineage = all_atlases[n]['cell_lineage']
            a, b = properties.temporal_alignment(ref_lineage, lineage, time_digits_for_cell_id=time_digits_for_cell_id)
            all_atlases[n]['temporal_alignment'] = (a, b)
        self.set_atlases(all_atlases)

    ############################################################
    #
    #
    #
    ############################################################

    def print_neighborhood(self, cellname):
        neighborhoods = self.get_neighborhoods()
        print("============================================================")
        print("Neighborhoods of " + str(cellname))
        if cellname not in neighborhoods:
            return
        refs = list(neighborhoods[cellname].keys())
        refs = sorted(refs)
        for r in refs:
            neighbors = list(neighborhoods[cellname][r].keys())
            neighbors = sorted(neighbors)
            print("    - " + str(r))
            for n in neighbors:
                print("        - " + str(n) + " : " + str(neighborhoods[cellname][r][n]))
        print("============================================================")

    def add_neighborhoods(self, prop, parameters, atlas_name, delay_from_division=0, time_digits_for_cell_id=4):
        """

        Parameters
        ----------
        prop: embryo properties (atlas to be added)
        parameters:
        atlas_name: file or atlas name (will indexed the added neighborhoods for each cell name)
        delay_from_division:
        time_digits_for_cell_id:

        Returns
        -------

        """
        proc = "add_neighborhoods"

        if not isinstance(parameters, AtlasParameters):
            monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                          + str(type(parameters)))
            sys.exit(1)

        #
        # build a nested dictionary of neighborhood, where the keys are
        # ['cell name']['reference name']
        # where 'reference name' is the name
        # of the reference lineage, and neighborhood a dictionary of contact surfaces indexed by cell names
        # only consider the first time point after the division
        #

        if 'cell_lineage' not in prop:
            monitoring.to_log_and_console(str(proc) + ": 'cell_lineage' was not in dictionary")
            return

        if 'cell_contact_surface' not in prop:
            monitoring.to_log_and_console(str(proc) + ": 'cell_contact_surface' was not in dictionary")
            return

        if 'cell_name' not in prop:
            monitoring.to_log_and_console(str(proc) + ": 'cell_name' was not in dictionary")
            return

        if 'cell_volume' not in prop:
            monitoring.to_log_and_console(str(proc) + ": 'cell_volume' was not in dictionary")
            return

        neighborhoods = self.get_neighborhoods(delay_from_division=delay_from_division)
        volumes = self.get_volumes(delay_from_division=delay_from_division)

        #
        # remove empty names
        # leading or trailing spaces
        #
        cells = list(prop['cell_name'].keys())
        for c in cells:
            if prop['cell_name'][c] == '':
                del prop['cell_name'][c]
                continue
            prop['cell_name'][c] = prop['cell_name'][c].strip()

        #
        lineage = prop['cell_lineage']
        name = prop['cell_name']
        contact = prop['cell_contact_surface']

        div = 10 ** time_digits_for_cell_id

        #
        # get the daughter cells just after division
        #
        reverse_lineage = {v: k for k, values in lineage.items() for v in values}
        daughters = [lineage[c][0] for c in lineage if len(lineage[c]) == 2]
        daughters += [lineage[c][1] for c in lineage if len(lineage[c]) == 2]

        ancestor_name = []
        missing_name = []
        missing_contact = []
        missing_neighbors = []

        for daugh in daughters:
            #
            # mother cell should be named
            #
            if reverse_lineage[daugh] not in name:
                continue
            #
            # check whether the cell is in dictionaries
            #
            if daugh not in name:
                if daugh not in missing_name:
                    missing_name.append(daugh)
                    monitoring.to_log_and_console(str(proc) + ": daughter cell #" + str(daugh)
                                                  + " was not found in 'cell_name' dictionary. Skip it", 6)
                continue

            if daugh not in contact:
                if daugh not in missing_contact:
                    missing_contact.append(daugh)
                    monitoring.to_log_and_console(str(proc) + ": daughter cell #" + str(daugh)
                                                  + " was not found in 'cell_contact_surface' dictionary. Skip it")
                continue

            #
            # check whether the mother name is the right one
            #
            if name[reverse_lineage[daugh]] != uname.get_mother_name(name[daugh]):
                msg = "weird, name of daughter cell #" + str(daugh) + " is " + str(name[daugh])
                msg += " while its mother #" + str(reverse_lineage[daugh]) + " is named "
                msg += str(name[reverse_lineage[daugh]]) + ". Skip it"
                monitoring.to_log_and_console(str(proc) + ": " + msg)
                continue

            #
            # get the daughter cell after some additional delay
            # positive delay: count from the division
            # negative delay: count from the end of the shortest branch for the two sisters
            #
            d = daugh
            local_delay_from_division = 0
            if delay_from_division >= 0:
                local_delay_from_division = delay_from_division
            elif delay_from_division < 0:
                length0 = _get_branch_length(d, lineage)
                sisters = copy.deepcopy(lineage[reverse_lineage[d]])
                sisters.remove(d)
                length1 = _get_branch_length(sisters[0], lineage)
                local_delay_from_division = min(length0, length1) + delay_from_division
                if local_delay_from_division < 0:
                    local_delay_from_division = 0

            for i in range(local_delay_from_division):
                if d not in lineage:
                    break
                if len(lineage[d]) > 1:
                    break
                nextd = lineage[d][0]
                if nextd not in name or nextd not in contact:
                    break
                if nextd not in name or nextd not in contact:
                    break
                d = nextd

            #
            # build the neighborhood
            #
            neighbor = {}
            neighbor_is_complete = True
            # half_id = '*' or '_'
            half_id = prop['cell_name'][d][-1]
            for c in contact[d]:
                n = int(c) % div
                if n == 1 or n == 0:
                    neighbor['background'] = neighbor.get('background', 0) + contact[d][c]
                else:
                    cname = c
                    if cname not in name:
                        while cname in reverse_lineage and cname not in name:
                            cname = reverse_lineage[cname]
                    if cname not in name:
                        neighbor_is_complete = False
                        if c not in missing_neighbors:
                            missing_neighbors.append(c)
                            msg = "cell #" + str(c) + " was not found in 'cell_name' dictionary."
                            monitoring.to_log_and_console(proc + ": " + msg)
                        continue
                    # c in name:
                    cell_name = prop['cell_name'][cname]
                    if cname != c:
                        if c not in ancestor_name:
                            ancestor_name += [c]
                            msg = "use name '" + str(cell_name) + "' of cell #" + str(cname) + " for cell #" + str(c)
                            monitoring.to_log_and_console(proc + ": " + msg)
                    if parameters.differentiate_other_half:
                        neighbor[cell_name] = neighbor.get(cell_name, 0) + contact[d][c]
                    else:
                        if cell_name[-1] == half_id:
                            neighbor[cell_name] = neighbor.get(cell_name, 0) + contact[d][c]
                        else:
                            neighbor['other-half'] = neighbor.get('other-half', 0) + contact[d][c]

            if not neighbor_is_complete:
                msg = ": neighborhood of " + str(prop['cell_name'][d]) + " is not complete. Skip it"
                monitoring.to_log_and_console(str(proc) + msg)
                continue

            #
            # add neighborhood and volume
            #
            # if neighbor_is_complete:
            if prop['cell_name'][d] not in neighborhoods:
                neighborhoods[prop['cell_name'][d]] = {}
            if atlas_name in neighborhoods[prop['cell_name'][d]]:
                msg = "weird, " + str(atlas_name) + " was already indexed for neighbors of cell " + \
                      str(prop['cell_name'][d])
                monitoring.to_log_and_console(str(proc) + ": " + msg)
            neighborhoods[prop['cell_name'][d]][atlas_name] = neighbor

            if prop['cell_name'][d] not in volumes:
                volumes[prop['cell_name'][d]] = {}
            if atlas_name in volumes[prop['cell_name'][d]]:
                msg = "weird, " + str(atlas_name) + " was already indexed for volume of cell " + \
                      str(prop['cell_name'][d])
                monitoring.to_log_and_console(str(proc) + ": " + msg)
            if d not in prop['cell_volume']:
                msg = "\t cell #" + str(d) + " was not found in 'cell_volume' dictionary."
                monitoring.to_log_and_console(str(proc) + ": " + msg)
            else:
                volumes[prop['cell_name'][d]][atlas_name] = prop['cell_volume'][d]
            #
            # add symmetric neighborhood if asked
            #
            if parameters.add_symmetric_neighborhood:
                sname = uname.get_symmetric_name(prop['cell_name'][d])
                sreference = 'sym-' + atlas_name
                sneighbor = get_symmetric_neighborhood(neighbor)
                if sname not in neighborhoods:
                    neighborhoods[sname] = {}
                if sreference in neighborhoods[sname]:
                    msg = "weird, " + str(sreference) + " was already indexed for cell " + str(sname)
                    monitoring.to_log_and_console(str(proc) + ": " + msg)
                neighborhoods[sname][sreference] = sneighbor

                if sname not in volumes:
                    volumes[sname] = {}
                if sreference in volumes[sname]:
                    msg = "weird, " + str(sreference) + " was already indexed for volume of cell " + str(sname)
                    monitoring.to_log_and_console(str(proc) + ": " + msg)
                if d not in prop['cell_volume']:
                    msg = "\t cell #" + str(d) + " was not found in 'cell_volume' dictionary."
                    monitoring.to_log_and_console(str(proc) + ": " + msg)
                else:
                    volumes[sname][sreference] = prop['cell_volume'][d]

        self.set_neighborhoods(neighborhoods, delay_from_division=delay_from_division)
        self.set_volumes(volumes, delay_from_division=delay_from_division)

        if len(missing_name) > 0:
            monitoring.to_log_and_console(
                str(proc) + ": daughter cells without names = " + str(len(missing_name)) + "/" +
                str(len(daughters)))

        if len(missing_contact) > 0:
            monitoring.to_log_and_console(str(proc) + ": daughter cells without contact surfaces  = " +
                                          str(len(missing_contact)) + "/" + str(len(daughters)))

        if len(missing_neighbors) > 0:
            monitoring.to_log_and_console(
                str(proc) + ": neighboring cells without names  = " + str(len(missing_neighbors)))

        monitoring.to_log_and_console("")

        return

    def add_atlas(self, name, prop, parameters):
        """

        Parameters
        ----------
        name
        prop: embryo properties (atlas to be added)
        parameters:

        Returns
        -------

        """
        proc = "add_atlas"

        if not isinstance(parameters, AtlasParameters):
            monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                          + str(type(parameters)))
            sys.exit(1)

        #
        # build a nested dictionary of neighborhood, where the keys are
        # ['cell name']['reference name']
        # where 'reference name' is the name
        # of the reference lineage, and neighborhood a dictionary of contact surfaces indexed by cell names
        # only consider the first time point after the division
        #

        if 'cell_lineage' not in prop:
            monitoring.to_log_and_console(str(proc) + ": 'cell_lineage' was not in dictionary")
            return

        if 'cell_name' not in prop:
            monitoring.to_log_and_console(str(proc) + ": 'cell_name' was not in dictionary")
            return

        if 'cell_contact_surface' not in prop:
            monitoring.to_log_and_console(str(proc) + ": 'cell_contact_surface' was not in dictionary")
            return

        if 'cell_volume' not in prop:
            monitoring.to_log_and_console(str(proc) + ": 'cell_volume' was not in dictionary")
            return

        if self._ref_atlas is None:
            self._ref_atlas = name

        atlases = self.get_atlases()
        atlases[name] = {}
        atlases[name]['cell_lineage'] = copy.deepcopy(prop['cell_lineage'])
        atlases[name]['cell_name'] = copy.deepcopy(prop['cell_name'])
        atlases[name]['cell_contact_surface'] = copy.deepcopy(prop['cell_contact_surface'])
        atlases[name]['cell_volume'] = copy.deepcopy(prop['cell_volume'])
        self.set_atlases(atlases)

        return

    def build_neighborhoods(self, atlasfiles, parameters, time_digits_for_cell_id=4):
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
        proc = "build_neighborhoods"

        if not isinstance(parameters, AtlasParameters):
            monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                          + str(type(parameters)))
            sys.exit(1)

        #
        # extract neighborhoods for required delays
        #
        delays = [parameters.name_delay_from_division]
        if parameters.confidence_delay_from_division is not None and \
                parameters.confidence_delay_from_division not in delays:
            delays += [parameters.confidence_delay_from_division]
        self.set_default_delay(parameters.name_delay_from_division)

        if isinstance(atlasfiles, str):
            prop = ioproperties.read_dictionary(atlasfiles, inputpropertiesdict={})
            name = atlasfiles.split(os.path.sep)[-1]
            if name.endswith(".xml") or name.endswith(".pkl"):
                name = name[:-4]
            if parameters.diagnosis_properties:
                udiagnosis.diagnosis(prop, ['name', 'contact'], parameters,
                                     time_digits_for_cell_id=time_digits_for_cell_id)
            for d in delays:
                self.add_neighborhoods(prop, parameters, atlas_name=name, delay_from_division=d,
                                       time_digits_for_cell_id=time_digits_for_cell_id)
            self.add_atlas(name, prop, parameters)
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
                if parameters.diagnosis_properties:
                    udiagnosis.diagnosis(prop, ['name', 'contact'], parameters,
                                         time_digits_for_cell_id=time_digits_for_cell_id)
                for d in delays:
                    self.add_neighborhoods(prop, parameters, atlas_name=name, delay_from_division=d,
                                           time_digits_for_cell_id=time_digits_for_cell_id)
                self.add_atlas(name, prop, parameters)
                del prop

        #
        # check neighborhood extraction
        #
        for d in delays:
            neighborhoods = self.get_neighborhoods(delay_from_division=d)
            if neighborhoods is None:
                monitoring.to_log_and_console(str(proc) + ": empty neighborhoods for delay_from_division = " + str(d))
                sys.exit(1)

        #
        # build common neighborhood reference if required
        #
        if parameters.use_common_neighborhood:
            monitoring.to_log_and_console("... build common neighborhoods", 1)
            for d in delays:
                neighborhoods = self.get_neighborhoods(delay_from_division=d)
                self.set_neighborhoods(_build_common_neighborhoods(neighborhoods), delay_from_division=d)
            self._use_common_neighborhood = True
            monitoring.to_log_and_console("    done", 1)

        #
        # temporal alignment (done from the cell number)
        #
        monitoring.to_log_and_console("... temporal alignment of lineages", 1)
        self.temporal_alignment(time_digits_for_cell_id=time_digits_for_cell_id)
        atlases = self.get_atlases()
        for n in atlases:
            msg = "    - "
            msg += "linear time warping of '" + str(n) + "' wrt '" + str(self._ref_atlas) + "' is "
            msg += str(atlases[n]['temporal_alignment'])
            monitoring.to_log_and_console(msg, 1)
        monitoring.to_log_and_console("    done", 1)

        #
        # dictionary indexed by mother cell names,
        # give list of references for which both daughter cells exist
        #
        monitoring.to_log_and_console("... build division list", 1)
        self.build_divisions(delay_from_division=parameters.name_delay_from_division)
        monitoring.to_log_and_console("    done", 1)

        #
        # dictionary indexed by mother cell names,
        # give list of daughter cells, the largest one being the first one
        #
        monitoring.to_log_and_console("... build unequal division list", 1)
        for d in delays:
            self.build_unequal_divisions(delay_from_division=d)
        monitoring.to_log_and_console("    done", 1)

        if parameters.diagnosis_properties:
            monitoring.to_log_and_console("")
            monitoring.to_log_and_console("============================================================")
            monitoring.to_log_and_console("===== diagnosis: atlases pairwise disagreements")
            _diagnosis_pairwise_switches(self, parameters)
            monitoring.to_log_and_console("===== diagnosis: dendrogram/linkage diagnosis")
            _diagnosis_linkage(self, parameters)
            monitoring.to_log_and_console("============================================================")
            monitoring.to_log_and_console("")

    ############################################################
    #
    #
    #
    ############################################################

    def extract_division_neighborhoods(self, mother_name, delay_from_division=None):

        delay = delay_from_division
        if delay is None:
            delay = self.get_default_delay()

        divisions = self.get_divisions()
        neighborhoods = self.get_neighborhoods(delay_from_division=delay)
        daughters = uname.get_daughter_names(mother_name)

        config = {}
        for r in divisions[mother_name]:
            config[r] = {}
            config[r][0] = copy.deepcopy(neighborhoods[daughters[0]][r])
            config[r][1] = copy.deepcopy(neighborhoods[daughters[1]][r])
        return config

    ############################################################
    #
    #
    #
    ############################################################

    def get_cell_distance(self, cell_name1, ref1, cell_name2, ref2, change_contact_surfaces=True):
        """
        Returns a value in [0, 1].
        Parameters
        ----------
        cell_name1
        ref1
        cell_name2
        ref2
        change_contact_surfaces

        Returns
        -------

        """
        proc = "get_cell_distance"
        neighborhoods = self.get_neighborhoods()

        if ref1 not in neighborhoods[cell_name1]:
            msg = proc + ": '" + str(ref1) + "' in not in neighborhoods of '" + str(cell_name1) + "'"
            monitoring.to_log_and_console(msg, 1)
            return -1
        if ref2 not in neighborhoods[cell_name2]:
            msg = proc + ": '" + str(ref2) + "' in not in neighborhoods of '" + str(cell_name2) + "'"
            monitoring.to_log_and_console(msg, 1)
            return -1
        return ucontact.cell_contact_distance(neighborhoods[cell_name1][ref1], neighborhoods[cell_name2][ref2],
                                              distance=self.cell_contact_distance,
                                              change_contact_surfaces=change_contact_surfaces)

    def get_division_similarity(self, cell_name, ref1, ref2, change_contact_surfaces=True):
        """
        Returns a value in [0, 1]. 0 means dissimilarity while 1 means perfect similarity.
        Parameters
        ----------
        cell_name: division/mother cell name
        ref1: first reference/atlas
        ref2: first reference/atlas
        change_contact_surfaces:

        Returns
        -------

        """
        proc = "get_division_similarity"
        neighborhoods = self.get_neighborhoods()
        d = uname.get_daughter_names(cell_name)

        if ref1 not in neighborhoods[d[0]]:
            msg = proc + ": '" + str(ref1) + "' in not in neighborhoods of '" + str(d[0]) + "'"
            monitoring.to_log_and_console(msg, 1)
            return -1
        if ref2 not in neighborhoods[d[0]]:
            msg = proc + ": '" + str(ref2) + "' in not in neighborhoods of '" + str(d[0]) + "'"
            monitoring.to_log_and_console(msg, 1)
            return -1
        if ref1 not in neighborhoods[d[1]]:
            msg = proc + ": '" + str(ref1) + "' in not in neighborhoods of '" + str(d[1]) + "'"
            monitoring.to_log_and_console(msg, 1)
            return -1
        if ref2 not in neighborhoods[d[1]]:
            msg = proc + ": '" + str(ref2) + "' in not in neighborhoods of '" + str(d[1]) + "'"
            monitoring.to_log_and_console(msg, 1)
            return -1

        similarity = self.get_division_contact_similarity()
        d = division_contact_generic_distance(self, neighborhoods[d[0]][ref1], neighborhoods[d[1]][ref1],
                                              neighborhoods[d[0]][ref2], neighborhoods[d[1]][ref2],
                                              similarity=similarity, change_contact_surfaces=change_contact_surfaces)
        if d >= 0.0:
            return 1.0 - d
        else:
            msg = proc + ": distance computation failed" + " " + str(d)
            monitoring.to_log_and_console(msg, 1)
            return -1


########################################################################################
#
#
#
########################################################################################

def division_contact_generic_distance(atlases, daughter00, daughter01, daughter10, daughter11, similarity='distance',
                                      change_contact_surfaces=True):
    """

    Parameters
    ----------
    atlases
    daughter00: daughter #0 of ref #0
    daughter01: daughter #1 of ref #0
    daughter10: daughter #0 of ref #1
    daughter11: daughter #1 of ref #1
    similarity: 'l1_distance', 'l2_distance'. The latter is kept for historical
        reasons, but only 'l1_distance' should be used.
    change_contact_surfaces: True or False

    Returns
    -------

    """
    proc = "division_contact_generic_distance"

    if similarity.lower() == 'distance' or similarity.lower() == 'norm' or \
            similarity.lower() == 'l1_distance' or similarity.lower() == 'l1-distance' or \
            similarity.lower() == 'l1_norm' or similarity.lower() == 'l1-norm' or \
            similarity.lower() == 'l2_distance' or similarity.lower() == 'l2-distance' or \
            similarity.lower() == 'l2_norm' or similarity.lower() == 'l2-norm':
        return ucontact.division_contact_distance(daughter00, daughter01, daughter10, daughter11,
                                                  distance=atlases.cell_contact_distance,
                                                  change_contact_surfaces=change_contact_surfaces)

    msg = proc + ": unhandled similarity '" + str(similarity) + "'"
    monitoring.to_log_and_console(msg, 1)
    return -1.0
