import os
import sys
import copy

import scipy.stats as stats
import numpy as np

import astec.utils.common as common
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

        self.outputAtlasFile = None
        self.atlasFiles = []

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
        doc += "\t If 'False', consider all the cells of the symmetric half-embryo.\n"
        doc += "\t as a single cell.\n"
        self.doc['differentiate_other_half'] = doc
        self.differentiate_other_half = True

        doc = "\t The same cell has different neighbors from an atlas to the other.\n"
        doc += "\t If 'True' build and keep an unique common neighborhood (set of\n"
        doc += "\t neighbors) for all atlases."
        self.doc['use_common_neighborhood'] = doc
        self.use_common_neighborhood = True

        doc = "\t Delay from the division to extract the neighborhooods.\n"
        doc += "\t 0 means right after the division.\n"
        self.doc['delay_from_division'] = doc
        self.delay_from_division = 0

        #
        #
        #
        doc = "\t How to compare two division patterns (a division is considered here\n"
        doc += "\t as the concatenation of the contact surface vectors of the 2 daughter\n"
        doc += "\t cells). Choices are:\n"
        doc += "\t - 'distance': the distance type is given by 'cell_contact_distance'\n"
        doc += "\t   distances are normalized between 0 (perfect match) and 1 (complete mismatch)\n"
        doc += "\t - 'probability': 1-(division probability) is used to keep the same meaning\n"
        doc += "\t   for the 0 and 1 extremal values. Probabilities are built with the distance\n"
        doc += "\t   'cell_contact_distance'\n"
        self.doc['division_contact_similarity'] = doc
        self.division_contact_similarity = 'distance'

        #
        #
        #
        doc = "\t Performs some diagnosis when reading an additional property file into the atlases\n"
        doc += "\t Incrementing the verboseness ('-v' in the command line) may give more details."
        self.doc['diagnosis_properties'] = doc
        self.diagnosis_properties = False

        #
        #
        #
        self.naming_improvement = False

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

        self.varprint('outputAtlasFile', self.outputAtlasFile)
        self.varprint('atlasFiles', self.atlasFiles)

        self.varprint('add_symmetric_neighborhood', self.add_symmetric_neighborhood)
        self.varprint('differentiate_other_half', self.differentiate_other_half)
        self.varprint('use_common_neighborhood', self.use_common_neighborhood)
        self.varprint('delay_from_division', self.delay_from_division)

        self.varprint('division_contact_similarity', self.division_contact_similarity)

        self.varprint('diagnosis_properties', self.diagnosis_properties)

        self.varprint('naming_improvement', self.naming_improvement)

        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# CellAtlasParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        udiagnosis.DiagnosisParameters.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'outputAtlasFile', self.outputAtlasFile, self.doc.get('outputAtlasFile', None))
        self.varwrite(logfile, 'atlasFiles', self.atlasFiles, self.doc.get('atlasFiles', None))

        self.varwrite(logfile, 'add_symmetric_neighborhood', self.add_symmetric_neighborhood,
                      self.doc.get('add_symmetric_neighborhood', None))
        self.varwrite(logfile, 'differentiate_other_half', self.differentiate_other_half,
                      self.doc.get('differentiate_other_half', None))
        self.varwrite(logfile, 'use_common_neighborhood', self.use_common_neighborhood,
                      self.doc.get('use_common_neighborhood', None))
        self.varwrite(logfile, 'delay_from_division', self.delay_from_division,
                      self.doc.get('delay_from_division', None))

        self.varwrite(logfile, 'division_contact_similarity', self.division_contact_similarity,
                      self.doc.get('division_contact_similarity', None))

        self.varwrite(logfile, 'diagnosis_properties', self.diagnosis_properties,
                      self.doc.get('diagnosis_properties', None))

        self.varwrite(logfile, 'naming_improvement', self.naming_improvement, self.doc.get('naming_improvement', None))

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

        self.outputAtlasFile = self.read_parameter(parameters, 'outputAtlasFile', self.outputAtlasFile)
        self.atlasFiles = self.read_parameter(parameters, 'atlasFiles', self.atlasFiles)
        self.atlasFiles = self.read_parameter(parameters, 'referenceFiles', self.atlasFiles)

        self.add_symmetric_neighborhood = self.read_parameter(parameters, 'add_symmetric_neighborhood',
                                                              self.add_symmetric_neighborhood)
        self.differentiate_other_half = self.read_parameter(parameters, 'differentiate_other_half',
                                                            self.differentiate_other_half)
        self.use_common_neighborhood = self.read_parameter(parameters, 'use_common_neighborhood',
                                                           self.use_common_neighborhood)
        self.delay_from_division = self.read_parameter(parameters, 'delay_from_division', self.delay_from_division)

        self.division_contact_similarity = self.read_parameter(parameters, 'cell_contact_distance',
                                                               self.division_contact_similarity)
        self.division_contact_similarity = self.read_parameter(parameters, 'division_contact_similarity',
                                                               self.division_contact_similarity)

        self.diagnosis_properties = self.read_parameter(parameters, 'diagnosis_properties', self.diagnosis_properties)
        self.diagnosis_properties = self.read_parameter(parameters, 'naming_diagnosis', self.diagnosis_properties)
        self.diagnosis_properties = self.read_parameter(parameters, 'diagnosis_naming', self.diagnosis_properties)

        self.naming_improvement = self.read_parameter(parameters, 'naming_improvement', self.naming_improvement)

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

def _write_summary_pairwise_switches(atlases, summary):

    divisions = atlases.get_divisions()
    percents = []
    mother_names = list(summary.keys())

    for n in mother_names:
        percents.append(100.0 * float(len(summary[n]['disagreement'])) / float(summary[n]['tested_couples']))
    [sorted_percents, sorted_mothers] = list(zip(*sorted(zip(percents, mother_names), reverse=True)))

    for n in sorted_mothers:
        if len(summary[n]['disagreement']) == 0:
            continue
        msg = " - " + str(n) + " cell division into "
        msg += str(uname.get_daughter_names(n)) + " has " + str(len(summary[n]['disagreement']))
        if len(summary[n]['disagreement']) > 1:
            msg += " disagreements"
        else:
            msg += " disagreement"
        percent = 100.0 * float(len(summary[n]['disagreement'])) / float(summary[n]['tested_couples'])
        msg += " (" + "{:2.2f}%".format(percent) + ") "
        msg += " over " + str(summary[n]['tested_couples']) + " tested configurations "
        monitoring.to_log_and_console(msg)
        msg = "\t over " + str(len(divisions[n]))
        msg += " references: " + str(divisions[n])
        monitoring.to_log_and_console(msg)
        msg = "\t " + str(summary[n]['disagreement'])
        monitoring.to_log_and_console(msg, 3)


def _diagnosis_pairwise_switches(atlases, parameters):
    proc = "_diagnosis_pairwise_switches"

    divisions = atlases.get_divisions()
    ccs = not parameters.use_common_neighborhood
    neighborhoods = atlases.get_neighborhoods()
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
    monitoring.to_log_and_console("")
    if len(divisions_with_disagreement) > 0:
        _write_summary_pairwise_switches(atlases, summary)

    return summary


########################################################################################
#
#
#
########################################################################################

def _di_switch_contact_surfaces(neighbors, reference, daughters):
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


def _di_global_generic_distance(neighbors, references, atlases, parameters, debug=False):
    """
    Compute a global score. The global score is the sum of local similarities over all couple of references/atlases.

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
    for r1 in references:
        for r2 in references:
            if r2 <= r1:
                continue
            dist = division_contact_generic_distance(atlases, neighbors[0][r1], neighbors[1][r1], neighbors[0][r2],
                                                     neighbors[1][r2], similarity=similarity,
                                                     change_contact_surfaces=ccs)
            if debug:
                print("   - dist[" + str(r1) + ", " + str(r2) + "] = " + str(dist))
            score += dist
    return score


def _di_test_one_division(atlases, mother, parameters):
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
        return {}

    daughters = uname.get_daughter_names(mother)
    neighborhoods = atlases.get_neighborhoods()
    neighbors = {0: copy.deepcopy(neighborhoods[daughters[0]]), 1: copy.deepcopy(neighborhoods[daughters[1]])}

    # score before any changes
    score = _di_global_generic_distance(neighbors, divisions[mother], atlases, parameters)

    corrections = []
    i = 1
    while True:
        newscore = {}
        for r in divisions[mother]:
            #
            # switch contact surfaces for the daughters in atlas 'r'
            #
            tmp = copy.deepcopy(neighbors)
            tmp = _di_switch_contact_surfaces(tmp, r, daughters)
            # compute a new score, keep it if it better than the one before any changes
            newscore[r] = _di_global_generic_distance(tmp, divisions[mother], atlases, parameters)
            if newscore[r] > score:
                del newscore[r]
        # no found correction at this iteration
        if len(newscore) == 0:
            return corrections
        # found several correction
        # 1. pick the only one (if only one is found)
        # 2. or pick the one with maximal score change
        elif len(newscore) == 1:
            ref = list(newscore.keys())[0]
        else:
            ref = min(newscore, key=lambda key: newscore[key])
        corrections += [(ref, score - newscore[ref])]
        i += 1
        # if one correction has been found, apply it
        # and look for an other additional correction
        tmp = copy.deepcopy(neighbors)
        tmp = _di_switch_contact_surfaces(tmp, ref, daughters)
        neighbors[0] = copy.deepcopy(tmp[0])
        neighbors[1] = copy.deepcopy(tmp[1])
        score = newscore[ref]


def division_improvement(atlases, parameters):

    # neighborhoods is a dictionary of dictionaries
    # ['cell name']['reference name']
    # first key is a cell name (daughter cell)
    # second key is the reference from which the neighborhood has been extracted

    #
    # mother cell name dictionary indexed by stage
    # stage 6: 32 cells
    # stage 7: 64 cells
    #
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

    for s in stages:
        corrections[s] = {}
        # if int(s) != 7:
        #    continue
        for m in mothers[s]:
            # if m != 'a7.0002_':
            #     continue
            correction = _di_test_one_division(atlases, m, parameters)
            #
            # correction is a dictionary indexed by the iteration index
            # each value is a tuple ('atlas name', score increment)
            #
            if len(correction) > 0:
                corrections[s][m] = correction

    monitoring.to_log_and_console("====== daughter switch proposal =====")
    monitoring.to_log_and_console("------ cell-based view")
    corrections_by_atlas = {}
    for s in stages:
        if len(corrections[s]) == 0:
            continue
        msg = "- generation " + str(s)
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
            msg = "  - division of '" + str(m) + "': "
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
            msg = "- reference '" + str(r) + "': " + str(corrections_by_atlas[r])
            monitoring.to_log_and_console(msg)
    monitoring.to_log_and_console("==========================================")


########################################################################################
#
#
#
########################################################################################

def _add_neighborhoods(previous_neighborhoods, prop, parameters, atlas_name, time_digits_for_cell_id=4):
    """

    Parameters
    ----------
    previous_neighborhoods: already built neighborhood atlas
    prop: embryo properties (atlas to be added)
    parameters:
    atlas_name: file or atlas name (will indexed the added neighborhoods for each cell name)
    time_digits_for_cell_id;

    Returns
    -------

    """
    proc = "_add_neighborhoods"

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

    missing_name = []
    missing_contact = []
    missing_neighbors = []

    for daugh in daughters:
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
                monitoring.to_log_and_console(str(proc) + ": cell #" + str(daugh)
                                              + " was not found in 'cell_contact_surface' dictionary. Skip it")
            continue

        #
        # get the daughter cell after some additional delay
        #
        d = daugh
        for i in range(parameters.delay_from_division):
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
        half_id = prop['cell_name'][d][-1]
        for c in contact[d]:
            n = int(c) % div
            if n == 1 or n == 0:
                neighbor['background'] = neighbor.get('background', 0) + contact[d][c]
            elif c in name:
                if parameters.differentiate_other_half:
                    neighbor[prop['cell_name'][c]] = contact[d][c]
                else:
                    if prop['cell_name'][c][-1] == half_id:
                        neighbor[prop['cell_name'][c]] = contact[d][c]
                    else:
                        neighbor['other-half'] = neighbor.get('other-half', 0) + contact[d][c]
            else:
                neighbor_is_complete = False
                if c not in missing_neighbors:
                    missing_neighbors.append(c)
                    msg = "\t cell #" + str(c) + " was not found in 'cell_name' dictionary."
                    monitoring.to_log_and_console(msg)
                continue

        if not neighbor_is_complete:
            msg = ": neighborhood of " + str(prop['cell_name'][d]) + " is not complete. Skip it"
            monitoring.to_log_and_console(str(proc) + msg)

        if neighbor_is_complete:
            if prop['cell_name'][d] not in previous_neighborhoods:
                previous_neighborhoods[prop['cell_name'][d]] = {}
            if atlas_name in previous_neighborhoods[prop['cell_name'][d]]:
                msg = "weird, " + str(atlas_name) + " was already indexed for cell " + str(prop['cell_name'][d])
                monitoring.to_log_and_console(str(proc) + ": " + msg)
            previous_neighborhoods[prop['cell_name'][d]][atlas_name] = neighbor
            #
            # add symmetric neighborhood if asked
            #
            if parameters.add_symmetric_neighborhood:
                sname = uname.get_symmetric_name(prop['cell_name'][d])
                sreference = 'sym-' + atlas_name
                sneighbor = get_symmetric_neighborhood(neighbor)
                if sname not in previous_neighborhoods:
                    previous_neighborhoods[sname] = {}
                if sreference in previous_neighborhoods[sname]:
                    msg = "weird, " + str(sreference) + " was already indexed for cell " + str(sname)
                    monitoring.to_log_and_console(str(proc) + ": " + msg)
                previous_neighborhoods[sname][sreference] = sneighbor

    if len(missing_name) > 0:
        monitoring.to_log_and_console(str(proc) + ": daughter cells without names = " + str(len(missing_name)) + "/" +
                                      str(len(daughters)))

    if len(missing_contact) > 0:
        monitoring.to_log_and_console(str(proc) + ": daughter cells without contact surfaces  = " +
                                      str(len(missing_contact)) + "/" + str(len(daughters)))

    if len(missing_neighbors) > 0:
        monitoring.to_log_and_console(str(proc) + ": neighboring cells without names  = " + str(len(missing_neighbors)))

    monitoring.to_log_and_console("")

    return previous_neighborhoods


########################################################################################
#
#
#
########################################################################################

def _build_common_neighborhoods(neighborhoods):

    # neighborhoods is a dictionary of dictionaries
    # ['cell name']['reference name']['neighboring cell']
    # first key is a cell name (daughter cell)
    # second key is the reference from which the neighborhood has been extracted
    # a neighborhood itself is a dictionary indexed by the neighboring cell names

    common_neighborhoods = {}

    for cell in neighborhoods:
        common_neighborhoods[cell] = ucontact.build_same_contact_surfaces(neighborhoods[cell])
    return common_neighborhoods


########################################################################################
#
#
#
########################################################################################

class Atlases(object):
    def __init__(self, parameters=None):

        self.cell_contact_distance = 'l1_distance'
        self.division_contact_similarity = 'l1_distance'

        # nested dictionary of neighborhoods, where the keys are ['cell name']['reference name']
        # where 'cell name' is the cell name (Conklin), and 'reference name' is the file name,
        # a neighborhood is a dictionary of contact surfaces indexed by cell names
        # it only considers the first time point after the division
        self._neighborhoods = {}

        # dictionary indexed by 'cell name' where 'cell name' is a mother cell giving the list
        # of references/atlases available for the two daughters
        self._divisions = {}

        self._probability_step = 0.01
        self._probability = None

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

    def get_neighborhoods(self):
        return self._neighborhoods

    def get_divisions(self):
        return self._divisions

    def get_probability_step(self):
        return self._probability_step

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
        elif parameters.division_contact_similarity.lower() == 'probability':
            self.division_contact_similarity = 'probability'
        else:
            monitoring.to_log_and_console(str(proc) + ": unhandled division contact similarity: '" +
                                          str(parameters.division_contact_similarity) + "'")

    ############################################################
    #
    #
    #
    ############################################################

    def build_divisions(self):
        proc = "build_divisions"

        if self._divisions is not None:
            del self._divisions
            self._divisions = {}

        #
        # get all references per division/mother cells
        #
        neighborhoods = self._neighborhoods
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
                if r in neighborhoods[daughters[0]] and r in neighborhoods[daughters[1]]:
                    self._divisions[n] = self._divisions.get(n, []) + [r]
                else:
                    monitoring.to_log_and_console(str(proc) + ": remove atlas '" + str(r) + "' for division '" + str(n)
                                                  + "'")

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

        if isinstance(parameters, ucontact.ContactSurfaceParameters):
            self.update_from_parameters(parameters)

        if isinstance(atlasfiles, str):
            prop = properties.read_dictionary(atlasfiles, inputpropertiesdict={})
            name = atlasfiles.split(os.path.sep)[-1]
            if name.endswith(".xml") or name.endswith(".pkl"):
                name = name[:-4]
            if parameters.diagnosis_properties:
                udiagnosis.diagnosis(prop, ['name', 'contact'], parameters,
                                     time_digits_for_cell_id=time_digits_for_cell_id)
            self._neighborhoods = _add_neighborhoods(self._neighborhoods, prop, parameters, atlas_name=name,
                                                     time_digits_for_cell_id=time_digits_for_cell_id)
            del prop
        elif isinstance(atlasfiles, list):
            for f in atlasfiles:
                prop = properties.read_dictionary(f, inputpropertiesdict={})
                name = f.split(os.path.sep)[-1]
                if name.endswith(".xml") or name.endswith(".pkl"):
                    name = name[:-4]
                if parameters.diagnosis_properties:
                    udiagnosis.diagnosis(prop, ['name', 'contact'], parameters,
                                         time_digits_for_cell_id=time_digits_for_cell_id)
                self._neighborhoods = _add_neighborhoods(self._neighborhoods, prop, parameters, atlas_name=name,
                                                         time_digits_for_cell_id=time_digits_for_cell_id)
                del prop

        if parameters.use_common_neighborhood:
            monitoring.to_log_and_console("... build common neighborhoods", 1)
            self._neighborhoods = _build_common_neighborhoods(self._neighborhoods)
            monitoring.to_log_and_console("    done", 1)

        monitoring.to_log_and_console("... build division list", 1)
        self.build_divisions()
        monitoring.to_log_and_console("    done", 1)

        if parameters.diagnosis_properties:
            monitoring.to_log_and_console("")
            monitoring.to_log_and_console("============================================================")
            monitoring.to_log_and_console("===== diagnosis: atlases pairwise disagreements")
            _diagnosis_pairwise_switches(self, parameters)
            monitoring.to_log_and_console("============================================================")
            monitoring.to_log_and_console("")

    def build_probabilities(self):

        monitoring.to_log_and_console("... build probabilities", 1)
        neighborhoods = self.get_neighborhoods()

        #
        # compute all score couples
        #
        cscores = []
        sscores = []

        for cell in neighborhoods:
            sister = uname.get_sister_name(cell)
            for ref in neighborhoods[cell]:
                if sister not in neighborhoods:
                    continue
                if ref not in neighborhoods[sister]:
                    continue
                #
                # cell score
                #
                for r in neighborhoods[cell]:
                    if r == ref:
                        continue
                    if r not in neighborhoods[sister]:
                        continue
                    cscores.append(ucontact.cell_contact_distance(neighborhoods[cell][ref], neighborhoods[cell][r],
                                                                  distance=self.cell_contact_distance))
                    sscores.append(ucontact.cell_contact_distance(neighborhoods[sister][ref], neighborhoods[sister][r],
                                                                  distance=self.cell_contact_distance))
        values = np.array([cscores, sscores])
        kernel = stats.gaussian_kde(values)

        step = self.get_probability_step()
        x, y = np.mgrid[0:1:step, 0:1:step]
        positions = np.vstack([x.ravel(), y.ravel()])
        z = np.reshape(kernel(positions).T, x.shape)
        scale = 100.0 / z.sum()
        self._probability = np.reshape([[(z[z <= z[j][i]].sum()) * scale for i in range(z.shape[1])]
                                        for j in range(z.shape[0])], z.shape)

        monitoring.to_log_and_console("    done", 1)

    def get_probability(self, a, b):
        proc = "get_probability"

        if self._probability is None:
            self.build_probabilities()

        i = int(b // self.get_probability_step())
        j = int(a // self.get_probability_step())
        if i >= self._probability.shape[0] or i >= self._probability.shape[1]:
            monitoring.to_log_and_console(str(proc) + ": too large second value in (" + str(a) + ", " + str(b) + "")
            return -1.0
        if j >= self._probability.shape[0] or j >= self._probability.shape[1]:
            monitoring.to_log_and_console(str(proc) + ": too large first value in (" + str(a) + ", " + str(b) + "")
            return -1.0

        return (self._probability[j][i] + self._probability[i][j]) / 2.0

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
            msg = proc + ": distance computation failed"
            monitoring.to_log_and_console(msg, 1)
            return -1


def division_contact_generic_distance(atlases, daughter00, daughter01, daughter10, daughter11, similarity='l1_distance',
                                      change_contact_surfaces=True):
    proc = "division_contact_generic_distance"

    if similarity.lower() == 'distance' or similarity.lower() == 'norm' or \
            similarity.lower() == 'l1_distance' or similarity.lower() == 'l1-distance' or \
            similarity.lower() == 'l1_norm' or similarity.lower() == 'l1-norm' or \
            similarity.lower() == 'l2_distance' or similarity.lower() == 'l2-distance' or \
            similarity.lower() == 'l2_norm' or similarity.lower() == 'l2-norm':
        return ucontact.division_contact_distance(daughter00, daughter01, daughter10, daughter11,
                                                  distance=atlases.cell_contact_distance,
                                                  change_contact_surfaces=change_contact_surfaces)
    elif similarity.lower() == 'probability':
        d0 = ucontact.cell_contact_distance(daughter00, daughter10, distance=atlases.cell_contact_distance,
                                            change_contact_surfaces=change_contact_surfaces)
        d1 = ucontact.cell_contact_distance(daughter01, daughter11, distance=atlases.cell_contact_distance,
                                            change_contact_surfaces=change_contact_surfaces)
        return 1.0 - atlases.get_probability(d0, d1)/100.0

    msg = proc + ": unhandled similarity '" + str(similarity) + "'"
    monitoring.to_log_and_console(msg, 1)
    return -1.0

########################################################################################
#
#
#
########################################################################################
