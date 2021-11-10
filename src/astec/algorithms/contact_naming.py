import os
import sys
import copy

import numpy as np
import random

import astec.utils.common as common
import astec.utils.contact_atlas as ucontacta
import astec.utils.properties as properties
import astec.utils.ascidian_name as uname
import astec.utils.diagnosis as diagnosis

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


class NamingParameters(ucontacta.AtlasParameters):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix='naming_'):

        if "doc" not in self.__dict__:
            self.doc = {}

        ucontacta.AtlasParameters.__init__(self, prefix=[prefix, "atlas_"])

        self.inputFile = []
        self.outputFile = None

        #
        #
        #
        doc = "\t Method to name the daugthers after a division \n"
        doc += "\t - 'sum' \n"
        doc += "\t - 'min' \n"
        self.doc['selection_method'] = doc
        self.selection_method = 'sum'

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

        ucontacta.AtlasParameters.print_parameters(self)

        self.varprint('inputFile', self.inputFile, self.doc.get('inputFile', None))
        self.varprint('outputFile', self.outputFile, self.doc.get('outputFile', None))

        self.varprint('selection_method', self.selection_method, self.doc.get('selection_method', None))

        self.varprint('testFile', self.testFile, self.doc.get('testFile', None))
        self.varprint('test_diagnosis', self.test_diagnosis, self.doc.get('test_diagnosis', None))
        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# NamingParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        ucontacta.AtlasParameters.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'inputFile', self.inputFile, self.doc.get('inputFile', None))
        self.varwrite(logfile, 'outputFile', self.outputFile, self.doc.get('outputFile', None))

        self.varwrite(logfile, 'selection_method', self.selection_method, self.doc.get('selection_method', None))

        self.varwrite(logfile, 'testFile', self.testFile, self.doc.get('testFile', None))
        self.varwrite(logfile, 'test_diagnosis', self.test_diagnosis, self.doc.get('test_diagnosis', None))

        logfile.write("\n")

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
        ucontacta.AtlasParameters.update_from_parameters(self, parameters)
        self.inputFile = self.read_parameter(parameters, 'inputFile', self.inputFile)
        self.outputFile = self.read_parameter(parameters, 'outputFile', self.outputFile)

        self.selection_method = self.read_parameter(parameters, 'selection_method', self.selection_method)

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
            if reference_prop['cell_name'][k] == uname.get_sister_name(prop['cell_name'][k]):
                if prop['cell_name'][k] not in sisters_errors:
                    sisters_errors[prop['cell_name'][k]] = 1
                    msg = "cell_name[" + str(k) + "]=" + str(prop['cell_name'][k]) + " is named as its sister "
                    msg += str(reference_prop['cell_name'][k])
                    if uname.get_mother_name(prop['cell_name'][k]) in discrepancies:
                        msg += ", division is incoherent among references"
                    monitoring.to_log_and_console("\t " + msg)
                else:
                    sisters_errors[prop['cell_name'][k]] += 1
                indexed_name = uname.get_mother_name(prop['cell_name'][k])
                division_errors[indexed_name] = division_errors.get(indexed_name, 0) + 1
            else:
                if prop['cell_name'][k] not in other_errors:
                    other_errors[prop['cell_name'][k]] = 1
                    msg = "cell_name[" + str(k) + "]=" + str(prop['cell_name'][k]) + " is not equal to reference name "
                    msg += str(reference_prop['cell_name'][k])
                    msg += " for cell " + str(k)
                    if uname.get_mother_name(prop['cell_name'][k]) in discrepancies:
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
    #     if uname.get_mother_name(n) in names or uname.get_mother_name(n) in first_error_mothers:
    #         second_errors[n] = division_errors[n]
    #     else:
    #         first_errors[n] = division_errors[n]
    #         if uname.get_mother_name(n) not in first_error_mothers:
    #             first_error_mothers.append(uname.get_mother_name(n))

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

def _get_neighborhoods(mother, immediate_daughters, prop, parameters, time_digits_for_cell_id=4):
    """

    Parameters
    ----------
    mother: cell id of the mother cell
    immediate_daughters: cell ids of the daughter cells
    prop: property dictionary of the embryo to be named
    parameters

    Returns
    -------
    A dictionary 'contact' indexed by [0, 1]. The value 'contact[i]' is the neighborhood of 'daughter[i]',
    a neighborhood being a contact surface vector indexed by cell names.
    None if there exists a neighbors without any named ancestor (including itself)
    """
    proc = "_get_neighborhoods"

    lineage = prop['cell_lineage']
    reverse_lineage = {v: k for k, values in lineage.items() for v in values}
    daughters = immediate_daughters

    #
    # delay from division
    #
    for i in range(parameters.delay_from_division):
        if len(lineage(daughters[0])) == 1 and len(lineage(daughters[1])) == 1:
            daughters[0] = lineage(daughters[0])[0]
            daughters[1] = lineage(daughters[1])[0]
        else:
            msg = ": build neighborhoods with a delay of " + str(i)
            msg += " instead of " + str(parameters.delay_from_division)
            msg += " for cells " + str(immediate_daughters)
            monitoring.to_log_and_console(str(proc) + msg, 4)
            break

    contact = {}
    div = 10 ** time_digits_for_cell_id
    half_id = prop['cell_name'][mother][-1]

    for index, d in enumerate(daughters):
        #
        # build contact surface as a dictionary of cell names
        # 1. background
        # 2. cell already named
        #    it can be named 'other-half'
        # 3. sister of d
        #    the key 'sister' is introduced here
        # 4. daughter cell not named, then named after its closest named ancestor (if any)
        #
        contact[index] = {}

        for c in prop['cell_contact_surface'][d]:
            if int(c) % div == 1 or int(c) % div == 0:
                contact[index]['background'] = contact[index].get('background', 0) + prop['cell_contact_surface'][d][c]
            elif c in prop['cell_name']:
                if parameters.differentiate_other_half:
                    contact[index][prop['cell_name'][c]] = prop['cell_contact_surface'][d][c]
                else:
                    if prop['cell_name'][c][-1] == half_id:
                        contact[index][prop['cell_name'][c]] = prop['cell_contact_surface'][d][c]
                    else:
                        contact[index]['other-half'] = contact[index].get('other-half', 0) + \
                                                       prop['cell_contact_surface'][d][c]
            elif c in daughters:
                if c != d:
                    contact[index]['sister'] = prop['cell_contact_surface'][d][c]
                else:
                    msg = ": weird, cell " + str(c) + " is in its contact surfaces"
                    monitoring.to_log_and_console(str(proc) + msg)
            else:
                #
                # cell without name, find a named ancestor
                #
                cell = c
                while cell in reverse_lineage and cell not in prop['cell_name']:
                    cell = reverse_lineage[cell]
                if cell in prop['cell_name']:
                    contact[index][prop['cell_name'][cell]] = contact[index].get(prop['cell_name'][cell], 0) + \
                                                           prop['cell_contact_surface'][d][c]
                else:
                    msg = ": unable to find a named ancestor for cell " + str(c)
                    msg += ". Skip naming of cells " + str(immediate_daughters)
                    monitoring.to_log_and_console(str(proc) + msg)
                    return None
    return contact


def _compute_distances(mother, daughters, prop, atlases, parameters, time_digits_for_cell_id=4):
    """

    Parameters
    ----------
    mother: cell id of the mother cell
    daughters: cell ids of the daughter cells
    prop: property dictionary of the embryo to be named
    atlases:
    parameters:
    time_digits_for_cell_id

    Returns
    -------
    a dictionary of distances (in [0,1]) indexed by [d][reference_name] where
        - d:
          d = 0: test (contacts[0], contacts[1]) <-> (daughter_names[0], daughter_names[1])
          d = 1: test (contacts[1], contacts[0]) <-> (daughter_names[0], daughter_names[1])
          where contacts[0] = contact vector for daughters[0]
                contacts[1] = contact vector for daughters[1]
        - reference_name is the name of a reference atlas/embryo
    """
    proc = "_compute_distances"

    divisions = atlases.get_divisions()
    if prop['cell_name'][mother] not in divisions:
        msg = ": no reference neighborhoods for division of '" + str(prop['cell_name'][mother]) + "'"
        msg += ". Can not name cells " + str(daughters)
        msg += " from mother cell " + str(mother)
        monitoring.to_log_and_console(str(proc) + msg, 4)
        return None

    neighborhoods = atlases.get_neighborhoods()
    daughter_names = uname.get_daughter_names(prop['cell_name'][mother])

    contacts = _get_neighborhoods(mother, daughters, prop, parameters, time_digits_for_cell_id=time_digits_for_cell_id)
    #
    # contacts[0] = contact vector for daughters[0]
    # contacts[1] = contact vector for daughters[1]
    # contacts[i]['sister'] is the contact surface of the sister
    #
    if contacts is None:
        msg = ": can not extract contact vector for division of '" + str(prop['cell_name'][mother]) + "'"
        msg += ". Can not name cells " + str(daughters)
        msg += " from mother cell " + str(mother)
        monitoring.to_log_and_console(str(proc) + msg, 4)
        return None

    scores = {0: {}, 1: {}}
    for i in range(2):
        #
        # i = 0: test (contacts[0], contacts[1]) <-> (daughter_names[0], daughter_names[1])
        # i = 1: test (contacts[1], contacts[0]) <-> (daughter_names[0], daughter_names[1])
        #
        if i == 0:
            if 'sister' in contacts[0]:
                contacts[0][daughter_names[1]] = contacts[0]['sister']
                del contacts[0]['sister']
            if 'sister' in contacts[1]:
                contacts[1][daughter_names[0]] = contacts[1]['sister']
                del contacts[1]['sister']
        else:
            if daughter_names[1] in contacts[0]:
                contacts[0][daughter_names[0]] = contacts[0][daughter_names[1]]
                del contacts[0][daughter_names[1]]
            if daughter_names[0] in contacts[1]:
                contacts[1][daughter_names[1]] = contacts[1][daughter_names[0]]
                del contacts[1][daughter_names[0]]
        #
        #
        #
        for ref in divisions[prop['cell_name'][mother]]:
            if i == 0:
                ic0 = 0
                ic1 = 1
            else:
                ic0 = 1
                ic1 = 0
            scores[i][ref] = ucontacta.division_contact_generic_distance(atlases, neighborhoods[daughter_names[0]][ref],
                                                                         neighborhoods[daughter_names[1]][ref],
                                                                         contacts[ic0], contacts[ic1],
                                                                         similarity=atlases.division_contact_similarity,
                                                                         change_contact_surfaces=True)
    return scores


########################################################################################
#
#
#
########################################################################################

def _give_name(daughters, daughter_names, distance, parameters, debug=False):
    proc = "_give_name"

    if len(distance[0]) != len(distance[1]):
        monitoring.to_log_and_console(str(proc) + ": weird, atlases number are different")

    arr0 = [distance[0][a] for a in distance[0]]
    arr1 = [distance[1][a] for a in distance[1]]
    name = {}
    name_certainty = {}

    if parameters.selection_method.lower() == 'minimum' or parameters.selection_method.lower() == 'min':
        if min(arr0) < min(arr1):
            name[daughters[0]] = daughter_names[0]
            name[daughters[1]] = daughter_names[1]
            name_certainty[daughters[0]] = int(100.0 - 100.0 * min(arr0))
            name_certainty[daughters[1]] = name_certainty[daughters[0]]
        elif min(arr1) < min(arr0):
            name[daughters[1]] = daughter_names[0]
            name[daughters[0]] = daughter_names[1]
            name_certainty[daughters[0]] = int(100.0 - 100.0 * min(arr1))
            name_certainty[daughters[1]] = name_certainty[daughters[0]]
        return name, name_certainty

    if parameters.selection_method.lower() == 'sum':
        if sum(arr0) < sum(arr1):
            name[daughters[0]] = daughter_names[0]
            name[daughters[1]] = daughter_names[1]
            name_certainty[daughters[0]] = int(100.0 - 100.0 * sum(arr0) / len(arr0))
            name_certainty[daughters[1]] = name_certainty[daughters[0]]
        elif sum(arr1) < sum(arr0):
            name[daughters[1]] = daughter_names[0]
            name[daughters[0]] = daughter_names[1]
            name_certainty[daughters[0]] = int(100.0 - 100.0 * sum(arr1) / len(arr1))
            name_certainty[daughters[1]] = name_certainty[daughters[0]]
        return name, name_certainty

    monitoring.to_log_and_console(str(proc) + ": selection method '" + str(parameters.selection_method) +
                                  "' not handled yet")
    return name, name_certainty


########################################################################################
#
# naming procedure
#
########################################################################################

def _propagate_name_along_branch(prop, cell):
    proc = "_propagate_name_along_branch"
    lineage = prop['cell_lineage']
    if cell not in prop['cell_name']:
        return
    c = cell
    while c in lineage and len(lineage[c]) == 1:
        nc = lineage[c][0]
        if nc in prop['cell_name']:
            if prop['cell_name'][nc] != prop['cell_name'][c]:
                msg = ": weird, cell " + str(nc) + " is named " + str(prop['cell_name'][nc])
                msg += ", but should be named " + str(prop['cell_name'][c])
                msg += " as its mother"
                monitoring.to_log_and_console(str(proc) + msg)
        else:
            prop['cell_name'][nc] = prop['cell_name'][c]
        c = nc
    return


def _propagate_naming(prop, atlases, parameters, time_digits_for_cell_id=4):
    proc = "_propagate_naming"

    if not isinstance(atlases, ucontacta.Atlases):
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

    if parameters.delay_from_division > 0:
        monitoring.to_log_and_console(str(proc) + ": WARNING, delay_from_division > 0 is not handled yet")
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
            mother_name = uname.get_mother_name(prop['cell_name'][c])
            if mother in prop['cell_name']:
                if prop['cell_name'][mother] != mother_name:
                    prop['selection_name_choice_certainty'][mother] = 0
                    msg = ": weird, cell " + str(mother) + " is named " + str(prop['cell_name'][mother])
                    msg += ", but should be named " + str(mother_name)
                    msg += " since one of its daughter is named " + str(prop['cell_name'][c])
                    monitoring.to_log_and_console(str(proc) + msg)
            else:
                prop['cell_name'][mother] = mother_name
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
    cell_not_named = []

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
                continue
            #
            # easy case
            # give name to cells that are only daughter and propagate along branch
            #
            if len(lineage[mother]) == 1:
                prop['cell_name'][c] = prop['cell_name'][mother]
                prop['selection_name_choice_certainty'][c] = prop['selection_name_choice_certainty'][mother]
                _propagate_name_along_branch(prop, c)
            #
            # in case of division:
            # 1. give name if the sister cell is named
            # 2. keep divisions to be solved
            #
            elif len(lineage[mother]) == 2:
                daughters = copy.deepcopy(lineage[mother])
                daughters.remove(c)
                daughter_names = uname.get_daughter_names(prop['cell_name'][mother])
                #
                # daughter cell is already named
                #
                if daughters[0] in prop['cell_name']:
                    if prop['cell_name'][daughters[0]] in daughter_names:
                        daughter_names.remove(prop['cell_name'][daughters[0]])
                        prop['cell_name'][c] = daughter_names[0]
                        prop['selection_name_choice_certainty'][c] = \
                            prop['selection_name_choice_certainty'][daughters[0]]
                    else:
                        msg = ": weird, cell " + str(daughters[0]) + " is named " + str(prop['cell_name'][daughters[0]])
                        msg += ", but should be named in " + str(daughter_names) + " since its mother cell "
                        msg += str(mother) + " is named " + str(prop['cell_name'][mother])
                        monitoring.to_log_and_console(str(proc) + msg)
                    continue
                #
                # both daughters are not named
                #
                division_to_be_named[mother] = [daughters[0], c]
            else:
                if mother not in cell_not_named:
                    cell_not_named.append(mother)
                cell_not_named.append(c)
                msg = ": weird, cell " + str(mother) + " has " + str(len(lineage[mother])) + " daughter(s)."
                msg += " Its offspring '" + str(lineage[mother]) + "' will not be named."
                monitoring.to_log_and_console(str(proc) + msg)

        if division_to_be_named == {}:
            continue

        #
        # here we have a dictionary of divisions to be named (key = mother id, value = array of sister ids)
        #
        for mother, daughters in division_to_be_named.items():
            debug = False
            #
            # distance is a dictionary of dictionary
            # distance[d in 0,1][ref name]
            # d = 0: test (contacts[0], contacts[1]) <-> (daughter_names[0], daughter_names[1])
            # d = 1: test (contacts[1], contacts[0]) <-> (daughter_names[0], daughter_names[1])
            # ref name in embryos
            #
            distance = _compute_distances(mother, daughters, prop, atlases, parameters,
                                          time_digits_for_cell_id=time_digits_for_cell_id)
            if debug:
                print("distance = " + str(distance))
            if distance is None:
                for c in daughters:
                    if c not in cell_not_named:
                        cell_not_named.append(c)
                # msg = "\t error when building scores for daughters " + str(daughters) + " of cell " + str(mother)
                # monitoring.to_log_and_console(msg)
                msg = " Can not name cells " + str(daughters) + " and their offsprings."
                monitoring.to_log_and_console(str(proc) + ": " + msg)
                continue

            daughter_names = uname.get_daughter_names(prop['cell_name'][mother])
            name, name_certainty = _give_name(daughters, daughter_names, distance, parameters, debug=debug)

            for c in name:
                if name[c] is not None:
                    prop['cell_name'][c] = name[c]
                    _propagate_name_along_branch(prop, c)
                prop['selection_name_choice_certainty'][c] = name_certainty[c]

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

    ucontacta.monitoring.copy(monitoring)

    #
    # should we clean reference here?
    #
    atlases = ucontacta.Atlases(parameters=parameters)
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
            diagnosis.monitoring.copy(monitoring)
            monitoring.to_log_and_console("============================================================")
            monitoring.to_log_and_console("===== diagnosis on '" + str(parameters.testFile) + "'")
            diagnosis.diagnosis(reference_prop, features=['name'], parameters=parameters,
                                time_digits_for_cell_id=time_digits_for_cell_id)
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
