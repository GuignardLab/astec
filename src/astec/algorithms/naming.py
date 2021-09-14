import os
import sys
import copy

import numpy as np
import random

import astec.utils.common as common
import astec.utils.neighborhood as uneighborhood
import astec.algorithms.properties as properties


#
#
#
#
#

monitoring = common.Monitoring()

_instrumented_ = False

########################################################################################
#
# classes
# - computation parameters
#
########################################################################################


class NamingParameters(uneighborhood.NeighborhoodParameters):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self):

        uneighborhood.NeighborhoodParameters.__init__(self)

        if "doc" not in self.__dict__:
            self.doc = {}

        self.inputFile = []
        self.outputFile = None

        #
        # for test:
        # names will be deleted, and tried to be rebuilt
        self.testFile = None
        self.test_diagnosis = False

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('#')
        print('# NamingParameters')
        print('#')
        print("")

        uneighborhood.NeighborhoodParameters.print_parameters(self)

        self.varprint('inputFile', self.inputFile)
        self.varprint('outputFile', self.outputFile)
        self.varprint('testFile', self.testFile)
        self.varprint('test_diagnosis', self.test_diagnosis)

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# NamingParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        uneighborhood.NeighborhoodParameters.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'inputFile', self.inputFile)
        self.varwrite(logfile, 'outputFile', self.outputFile)
        self.varwrite(logfile, 'testFile', self.testFile)
        self.varwrite(logfile, 'test_diagnosis', self.test_diagnosis)

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
        uneighborhood.NeighborhoodParameters.update_from_parameters(self, parameters)
        self.inputFile = self.read_parameter(parameters, 'inputFile', self.inputFile)
        self.outputFile = self.read_parameter(parameters, 'outputFile', self.outputFile)
        self.testFile = self.read_parameter(parameters, 'testFile', self.testFile)
        self.test_diagnosis = self.read_parameter(parameters, 'test_diagnosis', self.test_diagnosis)

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
# procedures for test
#
########################################################################################

def _build_test_set(prop, time_digits_for_cell_id=4, ncells=64):
    #
    # from an already named embryo, delete names except at one time point
    #
    proc = "_build_test_set"

    #
    # copy dictionary
    #
    returned_prop = {'cell_lineage': copy.deepcopy(prop['cell_lineage']), 'cell_name': copy.deepcopy(prop['cell_name']),
                     'cell_contact_surface': copy.deepcopy(prop['cell_contact_surface'])}

    #
    # get cells to count cells per time point
    #
    lineage = returned_prop['cell_lineage']
    cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
    cells = sorted(cells)

    div = 10 ** time_digits_for_cell_id
    ncells_per_time = np.zeros((max(cells) // div) + 1)
    for c in cells:
        t = int(c) // div
        ncells_per_time[t] += 1

    #
    # draw one time point of the desired number of cells
    #
    indices = np.where(np.array(ncells_per_time == ncells))[0]
    if len(indices) == 0:
        monitoring.to_log_and_console(str(proc) + ": there is no time point with #cells=" + str(ncells))
        indices = np.where(np.array(ncells_per_time > ncells))[0]
        if len(indices) == 0:
            monitoring.to_log_and_console(str(proc) + ": there is no time point with #cells>" + str(ncells))
            for t in range(len(ncells_per_time)):
                if ncells_per_time[t] == 0:
                    continue
                msg = "\t time " + str(t) + ": " + str(ncells_per_time[t]) + " cells"
                monitoring.to_log_and_console(msg)
            return None
        else:
            draw = indices[0]
    else:
        draw = random.choice(indices)

    monitoring.to_log_and_console(str(proc) + ": rename from time point " + str(draw))

    cells = list(returned_prop['cell_name'].keys())
    for c in cells:
        if int(c) // div != draw:
            del returned_prop['cell_name'][c]

    return returned_prop


def _test_naming(prop, reference_prop, discrepancies):
    proc = "_test_naming"

    monitoring.to_log_and_console("")
    monitoring.to_log_and_console(str(proc))
    monitoring.to_log_and_console("------------")

    #
    # get the cell names with error
    #
    sisters_errors = {}
    other_errors = {}
    division_errors = {}
    key_list = sorted(prop['cell_name'].keys())
    for k in key_list:
        if k not in reference_prop['cell_name']:
            monitoring.to_log_and_console("\t weird, key " + str(k) + " is not in reference properties")
        elif prop['cell_name'][k] != reference_prop['cell_name'][k]:
            if reference_prop['cell_name'][k] == uneighborhood.get_sister_name(prop['cell_name'][k]):
                if prop['cell_name'][k] not in sisters_errors:
                    sisters_errors[prop['cell_name'][k]] = 1
                    msg = "cell_name[" + str(k) + "]=" + str(prop['cell_name'][k]) + " is named as its sister "
                    msg += str(reference_prop['cell_name'][k])
                    if uneighborhood.get_mother_name(prop['cell_name'][k]) in discrepancies:
                        msg += ", division is incoherent among references"
                    monitoring.to_log_and_console("\t " + msg)
                else:
                    sisters_errors[prop['cell_name'][k]] += 1
                indexed_name = uneighborhood.get_mother_name(prop['cell_name'][k])
                division_errors[indexed_name] = division_errors.get(indexed_name, 0) + 1
            else:
                if prop['cell_name'][k] not in other_errors:
                    other_errors[prop['cell_name'][k]] = 1
                    msg = "cell_name[" + str(k) + "]=" + str(prop['cell_name'][k]) + " is not equal to reference name "
                    msg += str(reference_prop['cell_name'][k])
                    msg += " for cell " + str(k)
                    if uneighborhood.get_mother_name(prop['cell_name'][k]) in discrepancies:
                        msg += ", division is incoherent among references"
                    monitoring.to_log_and_console("\t " + msg)
                else:
                    other_errors[prop['cell_name'][k]] += 1

    #
    # if an error occur, a bad choice has been made at the mother division (first_errors)
    # or at some ancestor (second_errors)
    #
    # first_errors = {}
    # second_errors = {}
    # first_error_mothers = []
    # names = division_errors.keys()
    # for n in names:
    #     if uneighborhood.get_mother_name(n) in names or uneighborhood.get_mother_name(n) in first_error_mothers:
    #         second_errors[n] = division_errors[n]
    #     else:
    #         first_errors[n] = division_errors[n]
    #         if uneighborhood.get_mother_name(n) not in first_error_mothers:
    #             first_error_mothers.append(uneighborhood.get_mother_name(n))

    name_missing = {}
    reference_name = {}
    key_list = sorted(reference_prop['cell_name'].keys())
    for k in key_list:
        indexed_name = reference_prop['cell_name'][k]
        reference_name[indexed_name] = reference_name.get(indexed_name, 0) + 1
        if k not in prop['cell_name']:
            if reference_prop['cell_name'][k] not in name_missing:
                name_missing[reference_prop['cell_name'][k]] = 1
                msg = "reference cell_name[" + str(k) + "]=" + str(reference_prop['cell_name'][k]) + " is not found"
                monitoring.to_log_and_console("\t " + msg)
            else:
                name_missing[reference_prop['cell_name'][k]] += 1

    division_count = 0
    for k in division_errors:
        division_count += division_errors[k]
    other_count = 0
    for k in other_errors:
        other_count += other_errors[k]
    missing_count = 0
    for k in name_missing:
        missing_count += name_missing[k]

    msg = "ground-truth cells = " + str(len(reference_prop['cell_name'])) + " --- "
    msg += "ground-truth names = " + str(len(reference_name)) + " --- "
    msg += "tested cells = " + str(len(prop['cell_name'])) + " --- "
    msg += "retrieved names = " + str(len(reference_name) - len(name_missing)) + "\n"
    msg += "\t missing cell names = " + str(missing_count) + " for " + str(len(name_missing)) + " names --- \n"
    msg += "\t division errors in lineage = " + str(division_count) + " for " + str(len(division_errors)) + " names \n"
    if len(division_errors) > 0:
        msg += "\t    division errors = " + str(sorted(division_errors.keys())) + "\n"
    msg += "\t other errors in lineage = " + str(other_count) + " for " + str(len(other_errors)) + " names --- \n"
    if len(other_errors) > 0:
        msg += "\t    other errors = " + str(sorted(other_errors.keys())) + "\n"
    monitoring.to_log_and_console("summary" + ": " + msg)

    return


########################################################################################
#
#
#
########################################################################################

def _build_scores(mother, daughters, ancestor_name, prop, neighborhoods, parameters, time_digits_for_cell_id=4):
    proc = "_build_scores"

    #
    # are daughter names indexed?
    #
    daughter_names = uneighborhood.get_daughter_names(prop['cell_name'][mother])
    for name in daughter_names:
        #
        # no reference for this name
        #
        if name not in neighborhoods:
            msg = ": no reference neighborhoods for name " + str(name)
            msg += ". Can not name cells " + str(daughters)
            msg += " from mother cell " + str(mother)
            msg += " named " + str(prop['cell_name'][mother])
            monitoring.to_log_and_console(str(proc) + msg, 4)
            return None

    div = 10 ** time_digits_for_cell_id

    score = {}

    #
    # daughters is an array of 2 cell ids
    #
    for d in daughters:
        score[d] = {}

        #
        # build contact surface as a dictionary of names
        # 1. background
        # 2. cell already named
        # 3. sister of d
        #    give prop['cell_contact_surface'][d][c] to the untested name
        # 4. daughter cell not named, then named after its mother
        #    there might be two cells with this name
        #
        contact = {}
        sister = None
        for c in prop['cell_contact_surface'][d]:
            if int(c) % div == 1 or int(c) % div == 0:
                contact['background'] = contact.get('background', 0) + prop['cell_contact_surface'][d][c]
            elif c in prop['cell_name']:
                contact[prop['cell_name'][c]] = prop['cell_contact_surface'][d][c]
            elif c in daughters:
                if c != d:
                    sister = c
            elif c in ancestor_name:
                contact[ancestor_name[c]] = contact.get(ancestor_name[c], 0) + prop['cell_contact_surface'][d][c]
            else:
                monitoring.to_log_and_console("\t cell  " + str(c) + " was not found in 'cell_name' dictionary")
                monitoring.to_log_and_console(str(proc) + ": neighborhood of cell " + str(d)
                                              + " is not complete. Skip it")
                return None

        #
        # compute score for each candidate name by comparison with reference neighborhood
        #
        for name in daughter_names:
            #
            # get sister name
            #
            if name == daughter_names[0]:
                sister_name = daughter_names[1]
            else:
                sister_name = daughter_names[0]
            score[d][name] = {}
            #
            # add contact for the sister
            #
            if sister is not None:
                contact[sister_name] = prop['cell_contact_surface'][d][sister]
            for reference_name in neighborhoods[name]:
                score[d][name][reference_name] = uneighborhood.get_score(contact, neighborhoods[name][reference_name],
                                                                         parameters.neighborhood_comparison)
            if sister is not None:
                del contact[sister_name]
    return score


def _old_analyse_scores(scores):
    proc = "_analyse_scores"
    #
    # scores is a dictionary of dictionary
    # scores[d in daughters][n in daughter_names(mother)] is an dictionary of scalar product
    # the length of the array is the occurrence of [n in daughter_names(mother)] in the
    # neighborhood dictionary
    #
    name = {}
    name_certainty = {}

    # cell ids
    # cell name candidates
    # reference names
    ids = list(scores.keys())
    candidates = list(scores[ids[0]].keys())
    references = set(scores[ids[0]][candidates[0]].keys()).intersection(set(scores[ids[0]][candidates[1]].keys()),
                                                                        set(scores[ids[1]][candidates[0]].keys()),
                                                                        set(scores[ids[1]][candidates[1]].keys()))
    if references != set(scores[ids[0]][candidates[0]].keys()) or \
            references != set(scores[ids[0]][candidates[1]].keys()) or \
            references != set(scores[ids[1]][candidates[0]].keys()) or \
            references != set(scores[ids[1]][candidates[1]].keys()):
        msg = "weird, the set of references is different for each score"
        monitoring.to_log_and_console(str(proc) + ": " + msg)

    new_scores = {}
    agreement00 = 0
    agreement01 = 0
    disagreement = 0
    for r in references:
        new_scores[r] = [scores[ids[0]][candidates[0]][r], scores[ids[0]][candidates[1]][r],
                         scores[ids[1]][candidates[0]][r], scores[ids[1]][candidates[1]][r]]
        # 00 > 01 and 11 > 10
        if new_scores[r][0] > new_scores[r][1] and new_scores[r][3] > new_scores[r][2]:
            agreement00 += 1
        # 01 > 00 and 10 > 11
        elif new_scores[r][1] > new_scores[r][0] and new_scores[r][2] > new_scores[r][3]:
            agreement01 += 1
        else:
            disagreement += 1

    if agreement00 > agreement01 and agreement00 > disagreement:
        name[ids[0]] = candidates[0]
        name[ids[1]] = candidates[1]
        name_certainty[ids[0]] = agreement00 / (agreement00 + agreement01 + disagreement)
        name_certainty[ids[1]] = name_certainty[ids[0]]
    elif agreement01 > agreement00 and agreement01 > disagreement:
        name[ids[0]] = candidates[1]
        name[ids[1]] = candidates[0]
        name_certainty[ids[0]] = agreement01 / (agreement00 + agreement01 + disagreement)
        name_certainty[ids[1]] = name_certainty[ids[0]]
    else:
        msg = "there is a disagreement for cells " + str(candidates)
        monitoring.to_log_and_console(str(proc) + ": " + msg)
        score00 = 0
        score01 = 0
        for r in references:
            score00 += (new_scores[r][0] + new_scores[r][3]) / 2.0
            score01 += (new_scores[r][1] + new_scores[r][2]) / 2.0
        if score00 > score01:
            name[ids[0]] = candidates[0]
            name[ids[1]] = candidates[1]
            name_certainty[ids[0]] = 0.0
            name_certainty[ids[1]] = 0.0
        else:
            name[ids[0]] = candidates[1]
            name[ids[1]] = candidates[0]
            name_certainty[ids[0]] = 0.0
            name_certainty[ids[1]] = 0.0

    return name, name_certainty


def _old_analyse_scores_2021_05_18(scores, debug=False):
    proc = "_analyse_scores"
    #
    # scores is a dictionary of dictionary of dictionary
    # scores[cell id][name][reference] is the scalar product obtained when
    # associating the cell 'cell id' with 'name' for 'reference' neighborhood
    #
    name = {}
    name_certainty = {}

    # cell ids
    # cell name candidates
    # reference names
    ids = list(scores.keys())
    candidates = list(scores[ids[0]].keys())
    #
    # selection des references qui ont les 2 voisinages
    #
    references = set(scores[ids[0]][candidates[0]].keys()).intersection(set(scores[ids[0]][candidates[1]].keys()),
                                                                        set(scores[ids[1]][candidates[0]].keys()),
                                                                        set(scores[ids[1]][candidates[1]].keys()))
    if references != set(scores[ids[0]][candidates[0]].keys()) or \
            references != set(scores[ids[0]][candidates[1]].keys()) or \
            references != set(scores[ids[1]][candidates[0]].keys()) or \
            references != set(scores[ids[1]][candidates[1]].keys()):
        msg = "weird, the set of references is different for each score for cells " + str(candidates)
        monitoring.to_log_and_console(str(proc) + ": " + msg)

    new_scores = {}

    agreement00 = 0
    agreement01 = 0
    disagreement00 = 0
    disagreement01 = 0
    disagreement = 0
    sum_agreement00 = 0.0
    sum_agreement01 = 0.0
    sum_disagreement00 = 0.0
    sum_disagreement01 = 0.0

    for r in references:
        new_scores[r] = [scores[ids[0]][candidates[0]][r], scores[ids[0]][candidates[1]][r],
                         scores[ids[1]][candidates[0]][r], scores[ids[1]][candidates[1]][r]]
        # 00 > 01 and 11 > 10
        if new_scores[r][0] > new_scores[r][1] and new_scores[r][3] > new_scores[r][2]:
            sum_agreement00 += (new_scores[r][0] + new_scores[r][3]) - (new_scores[r][1] + new_scores[r][2])
            agreement00 += 1
            if debug:
                print(str(r) + ": total agreement for 00-11")
        # 01 > 00 and 10 > 11
        elif new_scores[r][0] < new_scores[r][1] and new_scores[r][3] < new_scores[r][2]:
            sum_agreement01 += (new_scores[r][1] + new_scores[r][2]) - (new_scores[r][0] + new_scores[r][3])
            agreement01 += 1
            if debug:
                print(str(r) + ": total agreement for 01-10")
        elif new_scores[r][0] + new_scores[r][3] > new_scores[r][1] + new_scores[r][2]:
            sum_disagreement00 += (new_scores[r][0] + new_scores[r][3]) - (new_scores[r][1] + new_scores[r][2])
            disagreement00 += 1
            if debug:
                print(str(r) + ": mild agreement for 00-11")
        elif new_scores[r][0] + new_scores[r][3] < new_scores[r][1] + new_scores[r][2]:
            sum_disagreement01 += (new_scores[r][1] + new_scores[r][2]) - (new_scores[r][0] + new_scores[r][3])
            disagreement01 += 1
            if debug:
                print(str(r) + ": mild agreement for 01-10")
        else:
            disagreement += 1
            if debug:
                print(str(r) + ": no agreement")

    somme = agreement00 + agreement01 + disagreement00 + disagreement01 + disagreement
    score00 = int(100.0 * (agreement00 + 0.5 * disagreement00) / somme)
    score01 = int(100.0 * (agreement01 + 0.5 * disagreement01) / somme)

    if score00 > score01:
        name[ids[0]] = candidates[0]
        name[ids[1]] = candidates[1]
        name_certainty[ids[0]] = score00
        name_certainty[ids[1]] = score00
        if debug:
            print("\t concludes to 00-11")
    elif score00 < score01:
        name[ids[0]] = candidates[1]
        name[ids[1]] = candidates[0]
        name_certainty[ids[0]] = score01
        name_certainty[ids[1]] = score01
        if debug:
            print("\t concludes to 01-10")
    else:
        msg = "there is no clear agreement for cells " + str(candidates)
        monitoring.to_log_and_console(str(proc) + ": " + msg)
        if sum_agreement00 + 0.5 * sum_disagreement00 > sum_agreement01 + 0.5 * sum_disagreement01:
            name[ids[0]] = candidates[0]
            name[ids[1]] = candidates[1]
            name_certainty[ids[0]] = score00
            name_certainty[ids[1]] = score00
        elif sum_agreement00 + 0.5 * sum_disagreement00 < sum_agreement01 + 0.5 * sum_disagreement01:
            name[ids[0]] = candidates[1]
            name[ids[1]] = candidates[0]
            name_certainty[ids[0]] = score01
            name_certainty[ids[1]] = score01
        else:
            msg = "there is no agreement at all for cells " + str(candidates)
            monitoring.to_log_and_console(str(proc) + ": " + msg)
            name[ids[0]] = None
            name[ids[1]] = None
            name_certainty[ids[0]] = 0
            name_certainty[ids[1]] = 0

    return name, name_certainty


def _analyse_scores_2021_05_25(scores, debug=False):
    proc = "_analyse_scores"
    #
    # scores is a dictionary of dictionary of dictionary
    # scores[cell id][name][reference] is the scalar product obtained when
    # associating the cell 'cell id' with 'name' for 'reference' neighborhood
    #
    name = {}
    name_certainty = {}

    # cell ids
    # cell name candidates
    # reference names
    ids = list(scores.keys())
    candidates = list(scores[ids[0]].keys())
    #
    # selection des references qui ont les 2 voisinages
    #
    references = set(scores[ids[0]][candidates[0]].keys()).intersection(set(scores[ids[0]][candidates[1]].keys()),
                                                                        set(scores[ids[1]][candidates[0]].keys()),
                                                                        set(scores[ids[1]][candidates[1]].keys()))
    if references != set(scores[ids[0]][candidates[0]].keys()) or \
            references != set(scores[ids[0]][candidates[1]].keys()) or \
            references != set(scores[ids[1]][candidates[0]].keys()) or \
            references != set(scores[ids[1]][candidates[1]].keys()):
        msg = "weird, the set of references is different for each score for cells " + str(candidates)
        monitoring.to_log_and_console(str(proc) + ": " + msg)

    new_scores = {}

    agreement00 = 0
    agreement01 = 0
    sum_agreement00 = 0.0
    sum_agreement01 = 0.0

    if debug:
        print("scores = " + str(scores))

    for r in references:
        if debug:
            print("- reference " + str(r))
        new_scores[r] = [scores[ids[0]][candidates[0]][r], scores[ids[0]][candidates[1]][r],
                         scores[ids[1]][candidates[0]][r], scores[ids[1]][candidates[1]][r]]
        if debug:
            print("- reference " + str(r) + " --- scores = " + str(new_scores[r]))
        # 00 > 01
        if new_scores[r][0] > new_scores[r][1]:
            agreement00 += 1
            sum_agreement00 += new_scores[r][0] - new_scores[r][1]
            if debug:
                print("   00 > 01")
        elif new_scores[r][0] < new_scores[r][1]:
            agreement01 += 1
            sum_agreement01 += new_scores[r][1] - new_scores[r][0]
            if debug:
                print("   00 < 01")

        # 00 > 10
        if new_scores[r][0] > new_scores[r][2]:
            agreement00 += 1
            sum_agreement00 += new_scores[r][0] - new_scores[r][2]
            if debug:
                print("   00 > 10")
        elif new_scores[r][0] < new_scores[r][2]:
            agreement01 += 1
            sum_agreement01 += new_scores[r][2] - new_scores[r][0]
            if debug:
                print("   00 < 10")

        # 11 > 01
        if new_scores[r][3] > new_scores[r][1]:
            agreement00 += 1
            sum_agreement00 += new_scores[r][3] - new_scores[r][1]
            if debug:
                print("   11 > 01")
        elif new_scores[r][3] < new_scores[r][1]:
            agreement01 += 1
            sum_agreement01 += new_scores[r][1] - new_scores[r][3]
            if debug:
                print("   11 < 01")

        # 11 > 10
        if new_scores[r][3] > new_scores[r][2]:
            agreement00 += 1
            sum_agreement00 += new_scores[r][3] - new_scores[r][2]
            if debug:
                print("   11 > 10")
        elif new_scores[r][3] < new_scores[r][2]:
            agreement01 += 1
            sum_agreement01 += new_scores[r][2] - new_scores[r][3]
            if debug:
                print("   11 < 10")

    if debug:
        print("agreement00 = " + str(agreement00))
        print("agreement01 = " + str(agreement01))
        print("sum_agreement00 = " + str(sum_agreement00))
        print("sum_agreement01 = " + str(sum_agreement01))
    if agreement00 > agreement01:
        name[ids[0]] = candidates[0]
        name[ids[1]] = candidates[1]
        name_certainty[ids[0]] = int(100.0 * agreement00 / (agreement00 + agreement01))
        name_certainty[ids[1]] = int(100.0 * agreement00 / (agreement00 + agreement01))
    elif agreement01 > agreement00:
        name[ids[0]] = candidates[1]
        name[ids[1]] = candidates[0]
        name_certainty[ids[0]] = int(100.0 * agreement01 / (agreement00 + agreement01))
        name_certainty[ids[1]] = int(100.0 * agreement01 / (agreement00 + agreement01))
    else:
        msg = "there is no clear agreement for cells " + str(candidates)
        monitoring.to_log_and_console(str(proc) + ": " + msg)
        if sum_agreement00 > sum_agreement01:
            name[ids[0]] = candidates[0]
            name[ids[1]] = candidates[1]
            name_certainty[ids[0]] = int(100.0 * sum_agreement00 / (sum_agreement00 + sum_agreement01))
            name_certainty[ids[1]] = int(100.0 * sum_agreement00 / (sum_agreement00 + sum_agreement01))
        elif sum_agreement01 > sum_agreement00:
            name[ids[0]] = candidates[1]
            name[ids[1]] = candidates[0]
            name_certainty[ids[0]] = int(100.0 * sum_agreement01 / (sum_agreement00 + sum_agreement01))
            name_certainty[ids[1]] = int(100.0 * sum_agreement01 / (sum_agreement00 + sum_agreement01))
        else:
            msg = "there is no agreement at all for cells " + str(candidates)
            monitoring.to_log_and_console(str(proc) + ": " + msg)
            name[ids[0]] = None
            name[ids[1]] = None
            name_certainty[ids[0]] = 0
            name_certainty[ids[1]] = 0

    return name, name_certainty


def _analyse_scores(scores, debug=False):
    proc = "_analyse_scores"
    #
    # scores is a dictionary of dictionary of dictionary
    # scores[cell id][name][reference] is the scalar product obtained when
    # associating the cell 'cell id' with 'name' for 'reference' neighborhood
    #
    name = {}
    name_certainty = {}

    # cell ids
    # cell name candidates
    # reference names
    ids = list(scores.keys())
    candidates = list(scores[ids[0]].keys())
    #
    # selection des references qui ont les 2 voisinages
    #
    references = set(scores[ids[0]][candidates[0]].keys()).intersection(set(scores[ids[0]][candidates[1]].keys()),
                                                                        set(scores[ids[1]][candidates[0]].keys()),
                                                                        set(scores[ids[1]][candidates[1]].keys()))
    if references != set(scores[ids[0]][candidates[0]].keys()) or \
            references != set(scores[ids[0]][candidates[1]].keys()) or \
            references != set(scores[ids[1]][candidates[0]].keys()) or \
            references != set(scores[ids[1]][candidates[1]].keys()):
        msg = "weird, the set of references is different for each score for cells " + str(candidates)
        monitoring.to_log_and_console(str(proc) + ": " + msg)

    new_scores = {}

    sum_agreement00 = 0.0
    sum_agreement01 = 0.0

    if debug:
        print("scores = " + str(scores))

    for r in references:
        if debug:
            print("- reference " + str(r))
        new_scores[r] = [scores[ids[0]][candidates[0]][r], scores[ids[0]][candidates[1]][r],
                         scores[ids[1]][candidates[0]][r], scores[ids[1]][candidates[1]][r]]
        if debug:
            print("- reference " + str(r) + " --- scores = " + str(new_scores[r]))
        sum_agreement00 += new_scores[r][0] + new_scores[r][3]
        sum_agreement01 += new_scores[r][1] + new_scores[r][2]

    if debug:
        print("sum_agreement00 = " + str(sum_agreement00))
        print("sum_agreement01 = " + str(sum_agreement01))

    if sum_agreement00 < sum_agreement01:
        name[ids[0]] = candidates[0]
        name[ids[1]] = candidates[1]
        name_certainty[ids[0]] = int(100.0 * sum_agreement00 / (sum_agreement00 + sum_agreement01))
        name_certainty[ids[1]] = int(100.0 * sum_agreement00 / (sum_agreement00 + sum_agreement01))
    elif sum_agreement01 < sum_agreement00:
        name[ids[0]] = candidates[1]
        name[ids[1]] = candidates[0]
        name_certainty[ids[0]] = int(100.0 * sum_agreement01 / (sum_agreement00 + sum_agreement01))
        name_certainty[ids[1]] = int(100.0 * sum_agreement01 / (sum_agreement00 + sum_agreement01))
    else:
        msg = "there is no agreement at all for cells " + str(candidates)
        monitoring.to_log_and_console(str(proc) + ": " + msg)
        name[ids[0]] = None
        name[ids[1]] = None
        name_certainty[ids[0]] = 0
        name_certainty[ids[1]] = 0

    return name, name_certainty


########################################################################################
#
# naming procedure
#
########################################################################################

def _propagate_naming(prop, atlases, parameters, time_digits_for_cell_id=4):
    proc = "_propagate_naming"

    if not isinstance(atlases, uneighborhood.Atlases):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'atlases' variable: "
                                      + str(type(atlases)))
        sys.exit(1)

    if 'cell_lineage' not in prop:
        monitoring.to_log_and_console(str(proc) + ": 'cell_lineage' was not in dictionary")
        return None

    if 'cell_contact_surface' not in prop:
        monitoring.to_log_and_console(str(proc) + ": 'cell_contact_surface' was not in dictionary")
        return None

    if 'cell_name' not in prop:
        monitoring.to_log_and_console(str(proc) + ": 'cell_name' was not in dictionary")
        return None

    #
    #
    #
    lineage = prop['cell_lineage']

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
    # initialize 'selection_name_choice_certainty'
    #
    prop['selection_name_choice_certainty'] = {}
    for k in prop['cell_name']:
        prop['selection_name_choice_certainty'][k] = 100

    #
    #
    #
    reverse_lineage = {v: k for k, values in lineage.items() for v in values}
    cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))

    #
    # backward propagation
    #
    monitoring.to_log_and_console(str(proc) + ": backward propagation")
    cells = sorted(cells, reverse=True)
    for c in cells:
        if c not in prop['cell_name']:
            continue
        if c not in reverse_lineage:
            continue
        mother = reverse_lineage[c]
        if len(lineage[mother]) == 1:
            if mother in prop['cell_name']:
                if prop['cell_name'][mother] != prop['cell_name'][c]:
                    prop['selection_name_choice_certainty'][mother] = 0
                    msg = ": weird, cell " + str(mother) + " is named " + str(prop['cell_name'][mother])
                    msg += ", but should be named " + str(prop['cell_name'][c])
                    msg += " as its single daughter"
                    monitoring.to_log_and_console(str(proc) + msg)
            else:
                prop['cell_name'][mother] = prop['cell_name'][c]
                prop['selection_name_choice_certainty'][mother] = 100
        elif len(lineage[mother]) == 2:
            ancestor_name = uneighborhood.get_mother_name(prop['cell_name'][c])
            if mother in prop['cell_name']:
                if prop['cell_name'][mother] != ancestor_name:
                    prop['selection_name_choice_certainty'][mother] = 0
                    msg = ": weird, cell " + str(mother) + " is named " + str(prop['cell_name'][mother])
                    msg += ", but should be named " + str(ancestor_name)
                    msg += " since one of its daughter is named " + str(prop['cell_name'][c])
                    monitoring.to_log_and_console(str(proc) + msg)
            else:
                prop['cell_name'][mother] = ancestor_name
                prop['selection_name_choice_certainty'][mother] = 100
        else:
            msg = ": weird, cell " + str(mother) + " has " + str(len(lineage[mother])) + "daughter(s)"
            monitoring.to_log_and_console(str(proc) + msg)

    #
    # forward propagation
    #
    monitoring.to_log_and_console(str(proc) + ": forward propagation")

    cells = sorted(cells)
    div = 10 ** time_digits_for_cell_id
    cells_per_time = {}
    missing_name = {}
    for c in cells:
        t = int(c) // div
        #
        # get cells and cell names at each time point
        #
        if t not in cells_per_time:
            cells_per_time[t] = [c]
        else:
            cells_per_time[t].append(c)
        if c not in prop['cell_name']:
            if t not in missing_name:
                missing_name[t] = [c]
            else:
                missing_name[t].append(c)

    timepoints = sorted(missing_name.keys())
    ancestor_name = {}
    cell_not_named = []

    neighborhoods = atlases.get_neighborhoods()

    for t in timepoints:
        division_to_be_named = {}
        for c in missing_name[t]:

            # already named
            if c in prop['cell_name']:
                continue
            # root of a tree, can not be named
            if c not in reverse_lineage:
                if c not in cell_not_named:
                    cell_not_named.append(c)
                monitoring.to_log_and_console(str(proc) + ": weird, cell " + str(c) + " is root of a subtree")
                continue

            # get its mother
            mother = reverse_lineage[c]
            # mother not in lineage
            # to account for lineage errors
            if mother not in lineage:
                monitoring.to_log_and_console(str(proc) + ": weird, cell " + str(mother) + " is not in lineage")
                continue
            # mother not named, can name the cell either
            if mother in cell_not_named:
                if c not in cell_not_named:
                    cell_not_named.append(c)
                continue
            # mother not named, can name the cell either
            if mother not in prop['cell_name']:
                if mother not in cell_not_named:
                    cell_not_named.append(mother)
                if c not in cell_not_named:
                    cell_not_named.append(c)
                msg = "mother cell " + str(mother) + " is not named."
                msg += " Can not name cell " + str(c) + " either."
                monitoring.to_log_and_console(str(proc) + ": " + msg, 5)
                if mother in ancestor_name:
                    ancestor_name[c] = ancestor_name[mother]
                else:
                    msg = "weird, cell " + str(mother) + " is not named and have no ancestor"
                    monitoring.to_log_and_console(str(proc) + ": " + msg, 5)
                continue
            #
            # easy case
            # give name to cells that are only daughter
            #
            if len(lineage[mother]) == 1:
                prop['cell_name'][c] = prop['cell_name'][mother]
                prop['selection_name_choice_certainty'][c] = prop['selection_name_choice_certainty'][mother]
            #
            # in case of division:
            # 1. give name if the sister cell is named
            # 2. keep divisions to be solved
            #
            elif len(lineage[mother]) == 2:
                daughters = copy.deepcopy(lineage[mother])
                daughters.remove(c)
                daughter_names = uneighborhood.get_daughter_names(prop['cell_name'][mother])
                #
                # daughter cell is already named
                #
                if daughters[0] in prop['cell_name']:
                    if prop['cell_name'][daughters[0]] in daughter_names:
                        daughter_names.remove(prop['cell_name'][daughters[0]])
                        prop['cell_name'][c] = daughter_names[0]
                        prop['selection_name_choice_certainty'][c] = prop['selection_name_choice_certainty'][daughters[0]]
                    else:
                        msg = ": weird, cell " + str(daughters[0]) + " is named " + str(prop['cell_name'][daughters[0]])
                        msg += ", but should be named in " + str(daughter_names) + " since its mother cell "
                        msg += str(mother) + " is named " + str(prop['cell_name'][mother])
                        monitoring.to_log_and_console(str(proc) + msg)
                    continue
                #
                # both daughters are not named: ancestor_name keep trace of their mother name
                #
                division_to_be_named[mother] = [daughters[0], c]
                ancestor_name[daughters[0]] = prop['cell_name'][mother]
                ancestor_name[c] = prop['cell_name'][mother]
            else:
                if mother not in cell_not_named:
                    cell_not_named.append(mother)
                cell_not_named.append(c)
                msg = ": weird, cell " + str(mother) + " has " + str(len(lineage[mother])) + " daughter(s)."
                msg += " Its offspring will not be named."
                monitoring.to_log_and_console(str(proc) + msg)

        if division_to_be_named == {}:
            continue

        #
        # here we have a dictionary of divisions to be named (key = mother id, value = array of sister ids)
        #
        for mother, daughters in division_to_be_named.items():
            debug = False
            #
            # scores is a dictionary of dictionary
            # scores[cell id][name][ref name] is an array of scalar product
            # cell id = d in daughters
            # name = name in daughter_names(mother)
            # the length of the array is the occurrence of [n in daughter_names(mother)] in the
            # neighborhood dictionary
            #
            scores = _build_scores(mother, daughters, ancestor_name, prop, neighborhoods, parameters,
                                   time_digits_for_cell_id=time_digits_for_cell_id)
            if debug:
                print("scores = " + str(scores))
            if scores is None:
                for c in daughters:
                    if c not in cell_not_named:
                        cell_not_named.append(c)
                # msg = "\t error when building scores for daughters " + str(daughters) + " of cell " + str(mother)
                # monitoring.to_log_and_console(msg)
                msg = " Can not name cells " + str(daughters) + " and their offsprings."
                monitoring.to_log_and_console(str(proc) + ": " + msg)
                continue

            name, name_certainty = _analyse_scores(scores, debug=debug)

            for c in name:
                if name[c] is not None:
                    prop['cell_name'][c] = name[c]
                prop['selection_name_choice_certainty'][c] = name_certainty[c]
                if c in ancestor_name:
                    del ancestor_name[c]

    return prop


########################################################################################
#
#
#
########################################################################################

def naming_process(experiment, parameters):
    proc = "naming_process"

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, NamingParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()

    uneighborhood.monitoring.copy(monitoring)

    #
    # should we clean reference here?
    #
    atlases = uneighborhood.Atlases()
    atlases.build_neighborhoods(parameters.atlasFiles, parameters, time_digits_for_cell_id=time_digits_for_cell_id)

    if atlases.get_neighborhoods() is None or atlases.get_neighborhoods() == {}:
        monitoring.to_log_and_console(str(proc) + ": empty neighborhood atlas?! ... Exiting ")
        sys.exit(1)

    #
    # read input properties to be named
    #
    prop = {}
    reference_prop = {}
    discrepancies = {}

    if parameters.testFile is not None:
        reference_prop = properties.read_dictionary(parameters.testFile, inputpropertiesdict={})
        if parameters.test_diagnosis:
            monitoring.to_log_and_console("============================================================")
            monitoring.to_log_and_console("===== diagnosis on '" + str(parameters.testFile) + "'")
            uneighborhood.diagnosis_properties_naming(reference_prop)
            monitoring.to_log_and_console("============================================================")
        prop = _build_test_set(reference_prop, time_digits_for_cell_id=time_digits_for_cell_id, ncells=64)
        if prop is None:
            monitoring.to_log_and_console(str(proc) + ": error when building test set")
            sys.exit(1)
    elif parameters.inputFile is not None:
        prop = properties.read_dictionary(parameters.inputFile, inputpropertiesdict={})

    if prop == {}:
        monitoring.to_log_and_console(str(proc) + ": no properties?!")
        sys.exit(1)

    if 'cell_name' not in prop:
        monitoring.to_log_and_console(str(proc) + ": no 'cell_name' in input dictionary")
        sys.exit(1)

    # clean from empty names
    cells = list(prop['cell_name'].keys())
    for c in cells:
        if prop['cell_name'][c] == '':
            del prop['cell_name'][c]

    #
    # naming propagation
    #
    prop = _propagate_naming(prop, atlases, parameters, time_digits_for_cell_id=time_digits_for_cell_id)
    prop = properties.set_fate_from_names(prop, time_digits_for_cell_id=time_digits_for_cell_id)
    prop = properties.set_color_from_fate(prop)
    #
    #
    #
    if parameters.testFile is not None:
        _test_naming(prop, reference_prop, discrepancies)

    if isinstance(parameters.outputFile, str):
        properties.write_dictionary(parameters.outputFile, prop)
