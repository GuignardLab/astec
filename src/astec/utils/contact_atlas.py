import os
import sys
import copy
import operator

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
        doc = "\t if 'True', add the symmetric neighborhood as additionnal exemplar.\n"
        doc += "\t Default is 'True'"
        self.doc['add_symmetric_neighborhood'] = doc
        self.add_symmetric_neighborhood = True

        doc = "\t if 'True', differentiate the cells of the symmetric half-embryo.\n"
        doc += "\t If 'False', consider all the cells of the symmetric half-embryo.\n"
        doc += "\t as a single cell.\n"
        doc += "\t Default is 'True'"
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
        doc = "\t Performs some diagnosis when reading an additional property file into the atlases"
        self.doc['diagnosis_properties'] = doc
        self.diagnosis_properties = False

        #
        # investigate whether daughters have to be switched for a division
        #
        self.naming_improvement = False
        self.cluster_fusion_probability = 75.0

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

        self.varprint('diagnosis_properties', self.diagnosis_properties)

        self.varprint('naming_improvement', self.naming_improvement)
        self.varprint('cluster_fusion_probability', self.cluster_fusion_probability)

        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# CellAtlasParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        udiagnosis.DiagnosisParameters.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'outputAtlasFile', self.outputAtlasFile)
        self.varwrite(logfile, 'atlasFiles', self.atlasFiles)

        self.varwrite(logfile, 'add_symmetric_neighborhood', self.add_symmetric_neighborhood)
        self.varwrite(logfile, 'differentiate_other_half', self.differentiate_other_half)
        self.varwrite(logfile, 'use_common_neighborhood', self.use_common_neighborhood)
        self.varwrite(logfile, 'delay_from_division', self.delay_from_division)

        self.varwrite(logfile, 'diagnosis_properties', self.diagnosis_properties)

        self.varwrite(logfile, 'naming_improvement', self.naming_improvement)
        self.varwrite(logfile, 'cluster_fusion_probability', self.cluster_fusion_probability)

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

        self.diagnosis_properties = self.read_parameter(parameters, 'diagnosis_properties', self.diagnosis_properties)
        self.diagnosis_properties = self.read_parameter(parameters, 'naming_diagnosis', self.diagnosis_properties)
        self.diagnosis_properties = self.read_parameter(parameters, 'diagnosis_naming', self.diagnosis_properties)

        self.naming_improvement = self.read_parameter(parameters, 'naming_improvement', self.naming_improvement)
        self.cluster_fusion_probability = self.read_parameter(parameters, 'cluster_fusion_probability',
                                                              self.cluster_fusion_probability)

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


#######################################################################################
#
#
#
########################################################################################

def _get_probability_consistency(atlases, parameters):
    proc = "_get_probability_consistency"
    has_written_something = False

    neighborhoods = atlases.get_neighborhoods()

    #
    # list of daughter cells
    #
    cell_names = sorted(list(neighborhoods.keys()))

    #
    # get the list of references per division
    #
    references = {}
    for cell_name in cell_names:
        mother_name = uname.get_mother_name(cell_name)
        references[mother_name] = references.get(mother_name, set()).union(set(neighborhoods[cell_name].keys()))

    #
    # get discrepancies
    #

    discrepancy = {}
    tested_couples = {}

    # cell_names is the list of daughter cell names of the neighborhood list
    for cell_name in cell_names:

        #
        # get cell name and sister name
        #
        sister_name = uname.get_sister_name(cell_name)
        if sister_name not in neighborhoods:
            msg = "weird, cell " + str(cell_name) + " is in the reference neighborhoods, while its sister "
            msg += str(sister_name) + " is not "
            monitoring.to_log_and_console(str(proc) + ": " + msg)
            has_written_something = True
            # cell_names.remove(cell_name)
            continue

        #
        # only one neighborhood, nothing to test
        #
        if len(neighborhoods[cell_name]) == 1:
            # cell_names.remove(cell_name)
            # cell_names.remove(sister_name)
            continue

        mother_name = uname.get_mother_name(cell_name)
        #
        # get two reference embryos
        #
        warnings = []
        for ref1 in neighborhoods[cell_name]:
            if ref1 not in neighborhoods[sister_name]:
                if ref1 in warnings:
                    continue
                msg = "weird, reference " + str(ref1) + " is in the neighborhoods of cell "
                msg += str(cell_name) + " but not of its sister " + str(sister_name)
                monitoring.to_log_and_console(str(proc) + ": " + msg)
                has_written_something = True
                warnings.append(ref1)
                continue
            for ref2 in neighborhoods[cell_name]:
                if ref2 == ref1:
                    continue
                if ref2 not in neighborhoods[sister_name]:
                    if ref2 in warnings:
                        continue
                    msg = "weird, reference " + str(ref2) + " is in the neighborhoods of cell "
                    msg += str(cell_name) + " but not of its sister " + str(sister_name)
                    monitoring.to_log_and_console(str(proc) + ": " + msg)
                    has_written_something = True
                    warnings.append(ref2)
                    continue

                same_choice_1 = ucontact.contact_distance(neighborhoods[cell_name][ref1],
                                                          neighborhoods[cell_name][ref2],
                                                          similarity=parameters.contact_similarity)
                same_choice_2 = ucontact.contact_distance(neighborhoods[sister_name][ref1],
                                                          neighborhoods[sister_name][ref2],
                                                          similarity=parameters.contact_similarity)
                sister_choice_1 = ucontact.contact_distance(neighborhoods[cell_name][ref1],
                                                            neighborhoods[sister_name][ref2],
                                                            similarity=parameters.contact_similarity)
                sister_choice_2 = ucontact.contact_distance(neighborhoods[sister_name][ref1],
                                                            neighborhoods[cell_name][ref2],
                                                            similarity=parameters.contact_similarity)

                prob_same_choice = atlases.get_probability(same_choice_1, same_choice_2)
                prob_sister_choice = atlases.get_probability(sister_choice_1, sister_choice_2)

                if mother_name not in tested_couples:
                    tested_couples[mother_name] = 1
                else:
                    tested_couples[mother_name] += 1
                if prob_sister_choice < prob_same_choice:
                    continue
                if mother_name not in discrepancy:
                    discrepancy[mother_name] = []
                discrepancy[mother_name].append((ref1, ref2))

    if has_written_something:
        monitoring.to_log_and_console("")

    #
    # for each division/mother cell
    # - reference is a dictionary that gives the set of references
    # - tested_couples is a dictionary that gives the number of tests done
    # - discrepancy is a dictionary (key = mother_name) that gives the list of couple (of references) of non-agreement
    #
    return references, tested_couples, discrepancy


def _check_probability_consistency(atlases, parameters):
    proc = "_check_probability_consistency"

    references, tested_couples, discrepancy = _get_probability_consistency(atlases, parameters)

    monitoring.to_log_and_console("")
    monitoring.to_log_and_console(str(proc))
    monitoring.to_log_and_console("-------------------------------------------")
    monitoring.to_log_and_console("--- reference probabilities consistency ---")

    #
    # if some divisions have some discrepancies, the following ones in the lineage
    # are likely to exhibit discrepancies also
    #
    mother_names = sorted(discrepancy.keys())
    if len(mother_names) > 0:
        _write_neighborhood_consistency("probability discrepancies", mother_names, references, tested_couples,
                                        discrepancy)

    msg = "tested divisions = " + str(len(tested_couples))
    monitoring.to_log_and_console(str(proc) + ": " + msg)
    msg = "divisions with probability discrepancies =  " + str(len(discrepancy))
    monitoring.to_log_and_console("\t " + msg)

    monitoring.to_log_and_console("-------------------------------------------")
    monitoring.to_log_and_console("")


def diagnosis_atlases_probability(atlases, parameters):

    proc = "diagnosis_atlases_probability"

    if not isinstance(atlases, Atlases):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'atlases' variable: "
                                      + str(type(atlases)))
        sys.exit(1)

    if not isinstance(parameters, AtlasParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    neighborhoods = atlases.get_neighborhoods()

    processed_mothers = []
    mother_refs = {}
    mother_maxproba = {}

    sorted_cells = sorted(neighborhoods.keys())
    for cell in sorted_cells:
        #
        # check whether 'mother' has both daughters in neighborhoods
        #
        mother = uname.get_mother_name(cell)
        if mother in processed_mothers:
            continue
        processed_mothers.append(mother)

        sister = uname.get_sister_name(cell)
        if sister not in neighborhoods:
            continue

        #
        # get couple of atlases (ref, r) having bother d and sister
        # and compute score. Keep the one yielding the best probability
        #

        for ref in neighborhoods[cell]:
            key = (mother, ref)
            if key in mother_refs:
                continue
            if ref not in neighborhoods[sister]:
                continue

            probas = []
            for r in neighborhoods[cell]:
                if r == ref:
                    continue
                if r not in neighborhoods[sister]:
                    continue
                cscore = ucontact.contact_distance(neighborhoods[cell][ref], neighborhoods[cell][r],
                                                   similarity=parameters.contact_similarity)
                sscore = ucontact.contact_distance(neighborhoods[sister][ref], neighborhoods[sister][r],
                                                   similarity=parameters.contact_similarity)
                prob = atlases.get_probability(cscore, sscore)
                probas += [(prob, r)]

            if len(probas) > 0:
                mother_refs[key] = sorted(probas, key=operator.itemgetter(0), reverse=True)
                mother_maxproba[key] = mother_refs[key][0]

    sorted_mother = sorted(mother_maxproba.items(), key=operator.itemgetter(1))

    print("\n")
    print("\n")
    print("\n")
    print("\n")
    print(str(sorted_mother))
    print("\n")
    print("\n")
    print("\n")
    print("\n")
    monitoring.to_log_and_console(str(proc) + ": details")
    for t in sorted_mother:
        msg = "  division of " + str(t[0][0]) + " in embryo " + str(t[0][1]) + " has maximal probability of "
        msg += "{:.2f} with embryo ".format(t[1][0]) + str(t[1][1])
        monitoring.to_log_and_console(msg)
    #
    #
    #


########################################################################################
#
#
#
########################################################################################

def get_neighborhood_consistency(neighborhoods, parameters):
    proc = "get_neighborhood_consistency"
    has_written_something = False

    #
    # list of daughter cells
    #
    cell_names = sorted(list(neighborhoods.keys()))

    #
    # get the list of references per division
    #
    references = {}
    for cell_name in cell_names:
        mother_name = uname.get_mother_name(cell_name)
        references[mother_name] = references.get(mother_name, set()).union(set(neighborhoods[cell_name].keys()))

    #
    # get discrepancies
    #

    discrepancy = {}
    tested_couples = {}

    # cell_names is the list of daughter cell names of the neighborhood list
    for cell_name in cell_names:

        #
        # get cell name and sister name
        #
        sister_name = uname.get_sister_name(cell_name)
        if sister_name not in neighborhoods:
            msg = "weird, cell " + str(cell_name) + " is in the reference neighborhoods, while its sister "
            msg += str(sister_name) + " is not "
            monitoring.to_log_and_console(str(proc) + ": " + msg)
            has_written_something = True
            # cell_names.remove(cell_name)
            continue

        #
        # only one neighborhood, nothing to test
        #
        if len(neighborhoods[cell_name]) == 1:
            # cell_names.remove(cell_name)
            # cell_names.remove(sister_name)
            continue

        mother_name = uname.get_mother_name(cell_name)
        #
        # get two reference embryos
        #
        warnings = []
        for ref1 in neighborhoods[cell_name]:
            # print("ref1 = " + str(ref1) + " - " + str(type(ref1)) + " - " + str(ref1 == ))
            for ref2 in neighborhoods[cell_name]:
                if ref2 == ref1:
                    continue
                if ref2 not in neighborhoods[sister_name]:
                    if ref2 in warnings:
                        continue
                    msg = "weird, reference " + str(ref2) + " is in the neighborhoods of cell "
                    msg += str(cell_name) + " but not of its sister " + str(sister_name)
                    monitoring.to_log_and_console(str(proc) + ": " + msg)
                    has_written_something = True
                    warnings.append(ref2)
                    continue

                same_choice = ucontact.contact_distance(neighborhoods[cell_name][ref1], neighborhoods[cell_name][ref2],
                                                        similarity=parameters.contact_similarity)
                sister_choice = ucontact.contact_distance(neighborhoods[cell_name][ref1],
                                                          neighborhoods[sister_name][ref2],
                                                          similarity=parameters.contact_similarity)

                if mother_name not in tested_couples:
                    tested_couples[mother_name] = 1
                else:
                    tested_couples[mother_name] += 1
                if same_choice < sister_choice:
                    continue
                if mother_name not in discrepancy:
                    discrepancy[mother_name] = []
                discrepancy[mother_name].append((ref1, ref2))

    if has_written_something:
        monitoring.to_log_and_console("")

    #
    # for each division/mother cell
    # - reference is a dictionary that gives the set of references
    # - tested_couples is a dictionary that gives the number of tests done
    # - discrepancy is a dictionary (key = mother_name) that gives the list of couple (of references) of non-agreement
    #
    return references, tested_couples, discrepancy


def _write_neighborhood_consistency(txt, mother_names, references, tested_couples, discrepancy):
    #
    # get discrepancy/inconsistency percentage per division
    #
    percents = []
    for mother_name in mother_names:
        percents.append(100.0 * float(len(discrepancy[mother_name])) / float(tested_couples[mother_name]))
    [sorted_percents, sorted_mothers] = list(zip(*sorted(zip(percents, mother_names), reverse=True)))

    msg = "\n*** " + str(txt) + " = " + str(len(mother_names))
    monitoring.to_log_and_console(msg)

    for mother_name in sorted_mothers:
        msg = " - " + str(mother_name) + " cell division into "
        msg += str(uname.get_daughter_names(mother_name)) + " has " + str(len(discrepancy[mother_name]))
        if len(discrepancy[mother_name]) > 1:
            msg += " discrepancies"
        else:
            msg += " discrepancy"
        percent = 100.0 * float(len(discrepancy[mother_name])) / float(tested_couples[mother_name])
        msg += " (" + "{:2.2f}%".format(percent) + ") "
        msg += " over " + str(tested_couples[mother_name]) + " tested configurations "
        monitoring.to_log_and_console(msg)
        msg = "\t over " + str(len(references[mother_name]))
        msg += " references: " + str(references[mother_name])
        monitoring.to_log_and_console(msg)
        msg = "\t " + str(discrepancy[mother_name])
        monitoring.to_log_and_console(msg, 3)


def _check_neighborhood_consistency(neighborhoods, parameters):
    proc = "_check_neighborhood_consistency"

    references, tested_couples, discrepancy = get_neighborhood_consistency(neighborhoods, parameters)

    monitoring.to_log_and_console("")
    monitoring.to_log_and_console(str(proc))
    monitoring.to_log_and_console("-------------------------------------------")
    monitoring.to_log_and_console("--- reference neighborhoods consistency ---")

    #
    # if some divisions have some discrepancies, the following ones in the lineage
    # are likely to exhibit discrepancies also
    #
    mother_names = sorted(discrepancy.keys())
    if len(mother_names) > 0:
        _write_neighborhood_consistency("neighborhood discrepancies", mother_names, references, tested_couples,
                                        discrepancy)

    msg = "tested divisions = " + str(len(tested_couples))
    monitoring.to_log_and_console(str(proc) + ": " + msg)
    msg = "divisions with discrepancies =  " + str(len(discrepancy))
    monitoring.to_log_and_console("\t " + msg)

    monitoring.to_log_and_console("-------------------------------------------")
    monitoring.to_log_and_console("")


def _check_leftright_consistency(neighborhoods, parameters):
    monitoring.to_log_and_console("")
    monitoring.to_log_and_console("-------------------------------------------")
    monitoring.to_log_and_console("--- left/right neighborhood consistency ---")

    #
    # list of daughter cells
    #
    cell_names = sorted(list(neighborhoods.keys()))

    #
    # get the list of references per division
    #
    references = {}
    for cell_name in cell_names:
        mother_name = uname.get_mother_name(cell_name)
        references[mother_name] = references.get(mother_name, set()).union(set(neighborhoods[cell_name].keys()))

    #
    #
    #
    mother_names = sorted(list(references.keys()))
    processed_mothers = []
    discrepancy = {}
    tested_cells = {}
    messages = {}
    for mother in mother_names:
        if mother in processed_mothers:
            continue
        symmother = uname.get_symmetric_name(mother)
        processed_mothers.append(mother)

        if symmother not in references:
            continue
        daughters = uname.get_daughter_names(mother)

        for reference in references[mother]:
            if reference not in references[symmother]:
                continue

            if reference not in tested_cells:
                tested_cells[reference] = 1
            else:
                tested_cells[reference] += 1

            for daughter in daughters:
                if daughter not in neighborhoods:
                    continue
                if reference not in neighborhoods[daughter]:
                    continue
                symdaughter = uname.get_symmetric_name(daughter)
                symsister = uname.get_sister_name(symdaughter)
                if symdaughter not in neighborhoods or symsister not in neighborhoods:
                    continue
                if reference not in neighborhoods[symdaughter] or reference not in neighborhoods[symsister]:
                    continue
                #
                #
                #
                symsameneigh = get_symmetric_neighborhood(neighborhoods[symdaughter][reference])
                symsisterneigh = get_symmetric_neighborhood(neighborhoods[symsister][reference])

                same_choice = ucontact.contact_distance(neighborhoods[daughter][reference], symsameneigh,
                                                        similarity=parameters.contact_similarity)
                sister_choice = ucontact.contact_distance(neighborhoods[daughter][reference], symsisterneigh,
                                                          similarity=parameters.contact_similarity)

                if same_choice < sister_choice:
                    continue

                if reference not in discrepancy:
                    discrepancy[reference] = {}
                discrepancy[reference][mother] = discrepancy[reference].get(mother, []) + [(daughter, symdaughter)]

                messages[mother] = messages.get(mother, []) + [(reference, daughter, symsister, symdaughter)]
                # msg = "   - '" + str(reference) + "': " + str(daughter) + " neighborhood is closest to "
                # msg += str(symsister) + " neighborhood than to " + str(symdaughter) + " one"
                # monitoring.to_log_and_console(msg, 3)

    mothers = sorted(messages.keys())
    for mother in mothers:
        msg = "   - '" + str(mother) + "' division"
        monitoring.to_log_and_console(msg, 3)
        for (reference, daughter, symsister, symdaughter) in messages[mother]:
            msg = "      - '" + str(reference) + "': " + str(daughter) + " neighborhood is closest to "
            msg += str(symsister) + " neighborhood than to " + str(symdaughter) + " one"
            monitoring.to_log_and_console(msg, 3)

    for reference in discrepancy:
        msg = "- '" + str(reference) + "' tested divisions = " + str(tested_cells[reference])
        monitoring.to_log_and_console(msg)

        #
        # get the mother cell for each discrepancy value
        #
        mother_by_discrepancy = {}
        processed_mothers = []
        for mother in discrepancy[reference]:
            if mother in processed_mothers:
                continue
            symmother = uname.get_symmetric_name(mother)
            processed_mothers += [mother, symmother]

            d = len(discrepancy[reference][mother])
            if symmother in discrepancy[reference]:
                d += len(discrepancy[reference][symmother])
            mother_by_discrepancy[d] = mother_by_discrepancy.get(d, []) + [mother]

        divisions_with_discrepancy = 0
        for d in mother_by_discrepancy:
            divisions_with_discrepancy += len(mother_by_discrepancy[d])
        msg = "\t divisions with left/right discrepancies =  " + str(divisions_with_discrepancy)
        monitoring.to_log_and_console(msg)

        for d in sorted(mother_by_discrepancy.keys()):
            if d == 1:
                msg = "\t - divisions with left/right " + str(d) + " discrepancy = "
            else:
                msg = "\t - divisions with left/right " + str(d) + " discrepancies = "
            msg += str(len(mother_by_discrepancy[d]))
            monitoring.to_log_and_console(msg)
            processed_mothers = []
            for mother in mother_by_discrepancy[d]:
                if mother in processed_mothers:
                    continue
                symmother = uname.get_symmetric_name(mother)
                processed_mothers += [mother, symmother]
                msg = "(" + str(mother) + ", " + str(symmother) + ") divisions: " + str(d)
                if d == 1:
                    msg += " discrepancy"
                else:
                    msg += " discrepancies"
                monitoring.to_log_and_console("\t   " + msg)

    monitoring.to_log_and_console("-------------------------------------------")
    monitoring.to_log_and_console("")
    return discrepancy


def diagnosis_neighborhood_consistency(neighborhoods, parameters):
    _check_leftright_consistency(neighborhoods, parameters)
    _check_neighborhood_consistency(neighborhoods, parameters)


#######################################################################################
#
# consistency is a measure between two reference embryos exhibiting the same division
# there is an inconsistency when the score computed between the same daughter cells
# of two reference embryos is lower than the score of this daughter cell (of one reference)
# and its sister cell (of the other reference). Thus, the consistency is a measure in
# [0, 4]: 0 means a total consistency while 4 designs a total inconsistency.
#
########################################################################################

def _si_global_score(neighbors, references, parameters, debug=False):
    """
    Compute a global score. The global score is the sume of local scores over all couple of references/atlases.
    A local score is the score (ie the scalar product) of association of each daughter (of a reference) with
    the same daughter (of the other reference).
    Parameters
    ----------
    neighbors: neighborhoods for the two daughters (dictionary indexed by [0,1] then by the references
    references: set of references
    debug

    Returns
    -------

    """
    score = 0
    for r1 in references:
        for r2 in references:
            if r2 <= r1:
                continue
            score00 = ucontact.contact_distance(neighbors[0][r1], neighbors[0][r2],
                                                similarity=parameters.contact_similarity)
            score11 = ucontact.contact_distance(neighbors[1][r1], neighbors[1][r2],
                                                similarity=parameters.contact_similarity)
            if debug:
                print("   - score[" + str(r1) + ", " + str(r2) + "] = " + str(score00))
                print("   - score[" + str(r1) + ", " + str(r2) + "] = " + str(score11))
            score += score00
            score += score11
    return score


def _si_switch_contact_surfaces(neighbors, reference, daughters):
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
    n = copy.deepcopy(neighbors[0][reference])
    for c in n:
        if c == daughters[1]:
            neighbors[0][reference][daughters[0]] = neighbors[0][reference][daughters[1]]
            del neighbors[0][reference][daughters[1]]
    #
    # contact surfaces of daughters[0] for atlas 'reference'
    # replace contact surface with daughters[1] with a contact surface with daughters[0]
    #
    n = copy.deepcopy(neighbors[1][reference])
    for c in n:
        if c == daughters[0]:
            neighbors[1][reference][daughters[1]] = neighbors[1][reference][daughters[0]]
            del neighbors[1][reference][daughters[0]]
    #
    # switch the contact surfaces of daughters[0] with the ones of daughters[1] for atlas 'reference'
    #
    tmp = copy.deepcopy(neighbors[0][reference])
    neighbors[0][reference] = neighbors[1][reference]
    neighbors[1][reference] = tmp
    return neighbors


def _si_test_one_division(neighborhoods, mother, parameters):
    """
    Test whether any daughter switch (for a given reference) improve the global score
    Parameters
    ----------
    neighborhoods
    mother

    Returns
    -------

    """
    daughters = uname.get_daughter_names(mother)
    neighbors = {0: copy.deepcopy(neighborhoods[daughters[0]]), 1: copy.deepcopy(neighborhoods[daughters[1]])}

    references = set.intersection(set(neighbors[0].keys()), set(neighbors[1].keys()))
    if len(references) <= 1:
        return {}

    # score before any changes
    score = _si_global_score(neighbors, references, parameters)

    corrections = {}
    i = 1
    while True:
        newscore = {}
        for r in references:
            #
            # switch contact surfaces for the daughters in atlas 'r'
            #
            tmp = copy.deepcopy(neighbors)
            tmp = _si_switch_contact_surfaces(tmp, r, daughters)
            # compute a new score, keep it if it better than the one before any changes
            newscore[r] = _si_global_score(tmp, references, parameters)
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
        corrections[i] = (ref, score - newscore[ref])
        i += 1
        # if one correction has been found, apply it
        # and look for an other additional correction
        tmp = copy.deepcopy(neighbors)
        tmp = _si_switch_contact_surfaces(tmp, ref, daughters)
        neighbors[0] = copy.deepcopy(tmp[0])
        neighbors[1] = copy.deepcopy(tmp[1])
        score = newscore[ref]


def global_score_improvement(neighborhoods, parameters):
    """
    Check whether the global score can be improved by switching the two daughter cells of a
    mother cell in a given reference. The global score score (for a mother cell) is the sum
    (over all possible reference embryo couple) of the local scores.
    Parameters
    ----------
    neighborhoods: dictionary of dictionaries. First level of keys is the cell name (name of a daughter), and the
                   second level of keys is the reference embryo name.
    parameters:
    Returns
    -------
    There is no returned values. It just prints out some change suggestions.
    """
    proc = "global_score_improvement"

    # neighborhoods is a dictionary of dictionaries
    # ['cell name']['reference name']
    # first key is a cell name (daughter cell)
    # second key is the reference from which the neighborhood has been extracted

    corrections = {}

    processed_mothers = {}
    sorted_cells = sorted(neighborhoods.keys())
    for d in sorted_cells:
        mother = uname.get_mother_name(d)
        if mother in processed_mothers:
            continue
        sister = uname.get_sister_name(d)
        if sister not in neighborhoods:
            msg = "weird, cell " + str(d) + " is in the reference neighborhoods, while its sister "
            msg += str(sister) + " is not "
            monitoring.to_log_and_console(str(proc) + ": " + msg)
            continue
        correction = _si_test_one_division(neighborhoods, mother, parameters)
        #
        # correction is a dictionary indexed by the iteration index
        # each value is a tuple ('atlas name', score increment)
        #
        if len(correction) > 0:
            corrections[mother] = correction

    tmp = {}
    for c in corrections:
        for i in corrections[c]:
            tmp[corrections[c][i][0]] = tmp.get(corrections[c][i][0], []) + [c]

    monitoring.to_log_and_console("====== global score improvement =====")
    # print(str(corrections))
    refs = sorted(tmp.keys())
    for r in refs:
        msg = str(r) + ": "
        cells = sorted(tmp[r])
        for i in range(len(cells)):
            msg += str(cells[i])
            if i < len(cells)-1:
                msg += ", "
        monitoring.to_log_and_console(msg)
    monitoring.to_log_and_console("====================================")


########################################################################################
#
#
#
########################################################################################

def _pi_ref_in_cluster(clusters, ref):
    proc = "_pi_ref_in_cluster"
    ret = None
    for c in clusters:
        if ref in clusters[c]:
            if ret is not None:
                msg = proc + ": weird, " + str(ref) + " is in cluster " + str(clusters[c])
                msg += "\t but was already in " + str(ret)
            ret = (c, "non-switched")
    if ref[0:9] == 'switched-':
        sref = ref[9:]
    else:
        sref = 'switched-' + ref
    for c in clusters:
        if sref in clusters[c]:
            if ret is not None:
                msg = proc + ": weird, switched" + str(sref) + " is in cluster " + str(clusters[c])
                msg += "\t but was already in " + str(ret)
            ret = (c, "switched")
    return ret


def _pi_delete_refs(tmp, ref0, ref1):

    if ref0[0:9] == 'switched-':
        r0 = ref0[9:]
    else:
        r0 = ref0
    if ref1[0:9] == 'switched-':
        r1 = ref1[9:]
    else:
        r1 = ref1

    if (r0, r1) in tmp:
        del tmp[(r0, r1)]
    if (r1, r0) in tmp:
        del tmp[(r1, r0)]

    if ('switched-' + r0, r1) in tmp:
        del tmp[('switched-' + r0, r1)]
    if (r1, 'switched-' + r0) in tmp:
        del tmp[(r1, 'switched-' + r0)]

    if (r0, 'switched-' + r1) in tmp:
        del tmp[(r0, 'switched-' + r1)]
    if ('switched-' + r1, r0) in tmp:
        del tmp[('switched-' + r1, r0)]

    if ('switched-' + r0, 'switched-' + r1) in tmp:
        del tmp[('switched-' + r0, 'switched-' + r1)]
    if ('switched-' + r1, 'switched-' + r0) in tmp:
        del tmp[('switched-' + r1, 'switched-' + r0)]

    return tmp


def _pi_test_one_division(atlases, mother, parameters, debug=False):
    proc = "_pi_test_one_division"

    if debug:
        monitoring.to_log_and_console("      " + proc + ": clustering of divisions of " + str(mother))

    daughters = uname.get_daughter_names(mother)
    neighborhoods = atlases.get_neighborhoods()

    #
    #
    #
    references = set.intersection(set(neighborhoods[daughters[0]].keys()), set(neighborhoods[daughters[1]].keys()))
    stages = []
    for r in references:
        for d in daughters:
            for n in neighborhoods[d][r]:
                if n == 'background' or n == 'other-half':
                    continue
                stages += [str(n).split('.')[0][1:]]
    if debug:
        msg = "      - maximal stage in neighbors is " + str(max(stages))
        monitoring.to_log_and_console("    " + msg)

    #
    # construction of the dictionary of neighborhood config[(switched-)ref][0,1]
    #
    config = {}
    for r in references:
        config[r] = {}
        config[r][0] = copy.deepcopy(neighborhoods[daughters[0]][r])
        config[r][1] = copy.deepcopy(neighborhoods[daughters[1]][r])
        sr = 'switched-' + str(r)
        config[sr] = {}
        config[sr][0] = copy.deepcopy(neighborhoods[daughters[1]][r])
        config[sr][1] = copy.deepcopy(neighborhoods[daughters[0]][r])
        if daughters[1] in config[sr][0]:
            if debug:
                msg = "  weird, " + str(daughters[1]) + " was found in its neighborhood for reference " + str(r)
                monitoring.to_log_and_console("      " + msg)
        if daughters[0] in config[sr][0]:
            config[sr][0][daughters[1]] = config[sr][0][daughters[0]]
            del config[sr][0][daughters[0]]
        if daughters[0] in config[sr][1]:
            if debug:
                msg = "  weird, " + str(daughters[0]) + " was found in its neighborhood for reference " + str(r)
                monitoring.to_log_and_console("      " + msg)
        if daughters[1] in config[sr][1]:
            config[sr][1][daughters[0]] = config[sr][1][daughters[1]]
            del config[sr][1][daughters[1]]

    #
    # computation of probabilities
    #
    prob = {}
    for r in config:
        for s in config:
            if r >= s:
                continue
            if r == 'switched-' + str(s) or s == 'switched-' + str(r):
                continue
            a = ucontact.contact_distance(config[r][0], config[s][0], similarity=parameters.contact_similarity)
            b = ucontact.contact_distance(config[r][1], config[s][1], similarity=parameters.contact_similarity)
            #
            # probabilities are supposed to be symmetrical but are not ...
            #
            prob[(r, s)] = (atlases.get_probability(a, b) + atlases.get_probability(b, a)) / 2.0
            prob[(s, r)] = prob[(r, s)]

    #
    #
    #
    tmp = copy.deepcopy(prob)
    clusters = {}
    ncluster = 0
    i = 0
    while len(tmp) > 0:
        i += 1
        #
        # sortedprob is an array of [tuple, value] where the 2-tuple is a couple of reference
        #
        sortedprob = [[c, tmp[c]] for c in tmp]
        sortedprob = sorted(sortedprob, key=operator.itemgetter(1), reverse=True)
        if debug:
            msg = "  #" + str(i) + " max prob = " + str(sortedprob[0][1]) + " for " + str(sortedprob[0][0])
            monitoring.to_log_and_console("      " + msg)

        #
        # sortedprob[0]       = first value after sort, ie the largest value
        # sortedprob[0][0]    = the corresponding tuple, ie the couple of references yielding the largest value
        # sortedprob[0][0][0] = the first reference
        # sortedprob[0][0][1] = the second reference
        #
        c0 = _pi_ref_in_cluster(clusters, sortedprob[0][0][0])
        c1 = _pi_ref_in_cluster(clusters, sortedprob[0][0][1])

        clusters_have_changed = False

        if c0 is None:
            # new cluster
            if c1 is None:
                clusters_have_changed = True
                if debug:
                    monitoring.to_log_and_console("      " + "  - new cluster " + str(ncluster))
                if sortedprob[0][0][0][0:9] == 'switched-' and sortedprob[0][0][1][0:9] == 'switched-':
                    clusters[str(ncluster)] = [sortedprob[0][0][0][9:], sortedprob[0][0][1][9:]]
                else:
                    clusters[str(ncluster)] = [sortedprob[0][0][0], sortedprob[0][0][1]]
                ncluster += 1
            # add sortedprob[0][0][0] to c1 cluster
            else:
                clusters_have_changed = True
                if debug:
                    monitoring.to_log_and_console("      " + "  - add to cluster " + str(c1[0]))
                if c1[1] == "non-switched":
                    clusters[c1[0]] += [sortedprob[0][0][0]]
                else:
                    if sortedprob[0][0][0][0:9] == 'switched-':
                        clusters[c1[0]] += [sortedprob[0][0][0][9:]]
                    else:
                        clusters[c1[0]] += ['switched-' + sortedprob[0][0][0]]
                for r in clusters[c1[0]]:
                    if r == sortedprob[0][0][0]:
                        continue
                    _pi_delete_refs(tmp, sortedprob[0][0][0], r)
        else:
            # add sortedprob[0][0][1] to c0 cluster
            if c1 is None:
                clusters_have_changed = True
                if debug:
                    monitoring.to_log_and_console("      " + "  - add to cluster " + str(c0[0]))
                if c0[1] == "non-switched":
                    clusters[c0[0]] += [sortedprob[0][0][1]]
                else:
                    if sortedprob[0][0][1][0:9] == 'switched-':
                        clusters[c0[0]] += [sortedprob[0][0][1][9:]]
                    else:
                        clusters[c0[0]] += ['switched-' + sortedprob[0][0][1]]
                for r in clusters[c0[0]]:
                    if r == sortedprob[0][0][1]:
                        continue
                    _pi_delete_refs(tmp, sortedprob[0][0][1], r)
            # bridge between 2 clusters
            else:
                if debug:
                    msg = "  - could fuse clusters " + str(c0[0]) + " and " + str(c1[0])
                    if c0[1] == c1[1]:
                        msg += " without switch"
                    else:
                        msg += ", one of them being switched"
                    monitoring.to_log_and_console("      " + msg)

        #
        # remove couples of probabilities
        #
        tmp = _pi_delete_refs(tmp, sortedprob[0][0][0], sortedprob[0][0][1])

        if debug and clusters_have_changed:
            monitoring.to_log_and_console("      " + "  - clusters = " + str(clusters))

    return clusters, prob, references


def _pi_intra_cluster_statistics(cluster, prob):
    p = []
    for r in cluster:
        for s in cluster:
            if r >= s:
                continue
            p += [prob[(r, s)]]
    return p


def _pi_inter_cluster_statistics(cluster, dluster, prob):
    proc = "_pi_inter_cluster_statistics"
    p = {"non-switched": [], "switched": []}
    for r in cluster:
        for s in dluster:
            p["non-switched"] += [prob[(r, s)]]
            if r[0:9] == 'switched-':
                sr = r[9:]
            else:
                sr = 'switched-' + r
            if s[0:9] == 'switched-':
                ss = s[9:]
            else:
                ss = 'switched-' + s
            p["switched"] += [prob[(sr, s)]]
            p["switched"] += [prob[(r, ss)]]
    return p


def _pi_cluster_fusion(clusters, stats, prob, parameters):
    if len(clusters) == 1:
        return clusters, stats

    tmpprob = {}
    representative = {}
    #
    # clusters keys are 0, 1, 2, ...
    #
    for c in clusters:
        representative[c] = str(c)
        for d in clusters:
            if c >= d:
                continue
            tmpprob[(c, d, "non-switched")] = max(stats[c][d]["non-switched"])
            tmpprob[(c, d, "switched")] = max(stats[c][d]["switched"])

    #
    #
    #
    tmp = copy.deepcopy(tmpprob)
    i = 0
    while len(tmp) > 0:

        sortedprob = [[c, tmp[c]] for c in tmp]
        sortedprob = sorted(sortedprob, key=operator.itemgetter(1), reverse=True)
        #
        # sortedprob[0][0] = (c, d, ["switched", "non-switched"] with d > c
        # c, d are 0, 1, 2, ...
        # representative keys are clusters keys, ie 0, 1, 2, ...
        # representative values are 0, 1, 2, ... potentially prefixed with 'switched-'
        #

        if sortedprob[0][1] < parameters.cluster_fusion_probability:
            break

        r = representative[sortedprob[0][0][0]]
        if r[0:9] == 'switched-':
            sr = r[9:]
        else:
            sr = 'switched-' + str(r)
        #
        # build equivalences
        #
        for c in clusters:
            #
            # look for representative = sortedprob[0][0][1] or 'switched-'sortedprob[0][0][1]
            #
            if representative[c] == str(sortedprob[0][0][1]):
                if sortedprob[0][0][2] == "switched":
                    representative[c] = sr
                else:
                    representative[c] = r
            elif representative[c] == 'switched-' + str(sortedprob[0][0][1]):
                if sortedprob[0][0][2] == "switched":
                    representative[c] = r
                else:
                    representative[c] = sr

        del tmp[sortedprob[0][0]]
        i += 1

    newclusters = {}
    ncluster = 0
    for c in representative:
        if representative[c] == c:
            newclusters[c] = clusters[c]
            ncluster += 1
        else:
            if representative[c][0:9] == 'switched-':
                r = representative[c][9:]
                for ref in clusters[c]:
                    if ref[0:9] == 'switched-':
                        newclusters[r] += [ref[9:]]
                    else:
                        newclusters[r] += ['switched-' + ref]
            else:
                newclusters[representative[c]] += clusters[c]

    newstats = {}
    for c in newclusters:
        newstats[c] = {}
        newstats[c][c] = _pi_intra_cluster_statistics(newclusters[c], prob)
        for d in newclusters:
            if c >= d:
                continue
            newstats[c][d] = _pi_inter_cluster_statistics(newclusters[c], newclusters[d], prob)

    return newclusters, newstats


def _pi_print_cluster_stats(clusters, stats, cluster_name="cluster", verbose=3):
    for c in clusters:
        msg = "    - " + str(cluster_name) + " #" + str(c) + " = " + str(clusters[c])
        monitoring.to_log_and_console(msg, verbose)
        if monitoring.verbose >= verbose:
            msg = "      intra-" + str(cluster_name) + " probabilities:"
            msg += " min = {:2.2f}".format(min(stats[c][c]))
            msg += " mean = {:2.2f}".format(np.mean(stats[c][c]))
            msg += " max = {:2.2f}".format(max(stats[c][c]))
            monitoring.to_log_and_console(msg, verbose)
    for c in clusters:
        for d in clusters:
            if c >= d:
                continue
            msg = "      - " + str(cluster_name) + " #" + str(c) + " <=> " + str(cluster_name) + " #" + str(d)
            monitoring.to_log_and_console(msg, verbose)
            msg = "        inter-" + str(cluster_name) + " probabilities (without switch):"
            msg += " min = {:2.2f}".format(min(stats[c][d]["non-switched"]))
            msg += " mean = {:2.2f}".format(np.mean(stats[c][d]["non-switched"]))
            msg += " max = {:2.2f}".format(max(stats[c][d]["non-switched"]))
            monitoring.to_log_and_console(msg, verbose)
            msg = "                                    (with switch):  "
            msg += " min = {:2.2f}".format(min(stats[c][d]["switched"]))
            msg += " mean = {:2.2f}".format(np.mean(stats[c][d]["switched"]))
            msg += " max = {:2.2f}".format(max(stats[c][d]["switched"]))
            monitoring.to_log_and_console(msg, verbose)


def global_probability_improvement(atlases, parameters):
    proc = "global_probability_improvement"

    # neighborhoods is a dictionary of dictionaries
    # ['cell name']['reference name']
    # first key is a cell name (daughter cell)
    # second key is the reference from which the neighborhood has been extracted

    #
    # mother cell name dictionary indexed by stage
    # stage 6: 32 cells
    # stage 7: 64 cells
    #
    mothers = {}
    neighborhoods = atlases.get_neighborhoods()
    for c in neighborhoods:
        mother = uname.get_mother_name(c)
        stage = mother.split('.')[0][1:]
        if stage not in mothers:
            mothers[stage] = []
        if mother not in mothers[stage]:
            mothers[stage] += [mother]

    for s in mothers:
        mothers[s] = sorted(mothers[s])

    stages = list(mothers.keys())
    stages.sort()

    monitoring.to_log_and_console("====== global probability improvement =====")
    for s in stages:
        # if int(s) != 7:
        #    continue
        monitoring.to_log_and_console("")
        msg = "- stage " + str(s)
        monitoring.to_log_and_console(msg)
        for m in mothers[s]:
            # if m != 'a7.0002_':
            #     continue

            clusters, prob, references = _pi_test_one_division(atlases, m, parameters, debug=False)
            stats = {}
            for c in clusters:
                stats[c] = {}
                stats[c][c] = _pi_intra_cluster_statistics(clusters[c], prob)
                for d in clusters:
                    if c >= d:
                        continue
                    stats[c][d] = _pi_inter_cluster_statistics(clusters[c], clusters[d], prob)

            msg = "  - division of " + str(m)
            monitoring.to_log_and_console(msg)

            _pi_print_cluster_stats(clusters, stats, cluster_name="cluster", verbose=4)

            newclusters, newstats = _pi_cluster_fusion(clusters, stats, prob, parameters)
            _pi_print_cluster_stats(newclusters, newstats, cluster_name="fused-cluster", verbose=3)

    monitoring.to_log_and_console("===========================================")


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
    def __init__(self):
        self._neighborhoods = {}
        self._probability_step = 0.01
        self._probability = None

    ############################################################
    #
    # getters
    #
    ############################################################
    def get_neighborhoods(self):
        return self._neighborhoods

    def get_probability_step(self):
        return self._probability_step

    ############################################################
    #
    #
    #
    ############################################################

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
        a neighborhood a dictionary of contact surfaces indexed by cell names
        it only considers the first time point after the division

        """
        proc = "build_neighborhoods"

        if not isinstance(parameters, AtlasParameters):
            monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                          + str(type(parameters)))
            sys.exit(1)

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

        # consistency is a measure between two reference embryos exhibiting the same division
        # there is an inconsistency when the score computed between the same daughter cells
        # of two reference embryos is lower than the score of this daughter cell (of one reference)
        # and its sister cell (of the other reference). Thus, the consistency is a measure in
        # [0, 4]: 0 means a total consistency while 4 designs a total inconsistency.

        if parameters.use_common_neighborhood:
            monitoring.to_log_and_console("... build common neighborhoods", 1)
            self._neighborhoods = _build_common_neighborhoods(self._neighborhoods)
            monitoring.to_log_and_console("    done", 1)

        monitoring.to_log_and_console("... build probabilities", 1)
        self.build_probabilities(parameters)
        monitoring.to_log_and_console("    done", 1)

        if parameters.diagnosis_properties:
            monitoring.to_log_and_console("============================================================")
            monitoring.to_log_and_console("===== diagnosis on consistency")
            # diagnosis_neighborhood_consistency(self._neighborhoods, parameters)
            monitoring.to_log_and_console("===== diagnosis on probabilities")
            _check_probability_consistency(self, parameters)
            # diagnosis_atlases_probability(self, parameters)
            monitoring.to_log_and_console("============================================================")

    def build_probabilities(self, parameters):
        proc = "build_probabilities"

        if not isinstance(parameters, AtlasParameters):
            monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                          + str(type(parameters)))
            sys.exit(1)

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
                    cscores.append(ucontact.contact_distance(neighborhoods[cell][ref], neighborhoods[cell][r],
                                                             similarity=parameters.contact_similarity))
                    sscores.append(ucontact.contact_distance(neighborhoods[sister][ref], neighborhoods[sister][r],
                                                             similarity=parameters.contact_similarity))
        values = np.array([cscores, sscores])
        kernel = stats.gaussian_kde(values)

        step = self.get_probability_step()
        x, y = np.mgrid[0:1:step, 0:1:step]
        positions = np.vstack([x.ravel(), y.ravel()])
        z = np.reshape(kernel(positions).T, x.shape)
        scale = 100.0 / z.sum()
        self._probability = np.reshape([[(z[z <= z[j][i]].sum()) * scale for i in range(z.shape[1])]
                                        for j in range(z.shape[0])], z.shape)

    def get_probability(self, a, b):
        proc = "get_probability"
        i = int(b // self.get_probability_step())
        j = int(a // self.get_probability_step())
        if i >= self._probability.shape[1]:
            monitoring.to_log_and_console(str(proc) + ": too large second value in (" + str(a) + ", " + str(b) + "")
            return -1.0
        if j >= self._probability.shape[0]:
            monitoring.to_log_and_console(str(proc) + ": too large first value in (" + str(a) + ", " + str(b) + "")
            return -1.0
        return self._probability[j][i]


########################################################################################
#
#
#
########################################################################################
