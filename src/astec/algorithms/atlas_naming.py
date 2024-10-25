import os
import sys
import copy

import numpy as np
import random
import statistics

import astec.utils.common as common
import astec.utils.atlas_division as uatlasd
import astec.utils.atlas_cell as uatlasc
import astec.utils.atlas_embryo as uatlase
import astec.utils.neighborhood_distance as uneighborhood
import astec.utils.properties as properties
import astec.utils.ioproperties as ioproperties
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


class NamingParameters(uatlasd.DivisionParameters):
    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix="naming_"):
        if "doc" not in self.__dict__:
            self.doc = {}

        uatlasd.DivisionParameters.__init__(self, prefix=[prefix, "atlas_"])

        doc = "\t Input property file to be named. Must contain lineage, volumes and contact surfaces\n"
        doc += (
            "\t as well as some input names (one time point should be entirely named)."
        )
        self.doc["inputFile"] = doc
        self.inputFile = []
        doc = "\t Output property file."
        self.doc["outputFile"] = doc
        self.outputFile = None

        #
        #
        #
        doc = (
            "\t Method to name the daughters after a division. Distances are computed\n"
        )
        doc += "\t between the couple of daughters as well as the couple of switched\n"
        doc += "\t daughters and all the atlas divisions.\n"
        doc += (
            "\t - 'mean': choose the couple of names that yield the minimal average \n"
        )
        doc += "\t   distance over all the atlases. \n"
        doc += "\t - 'minimum': choose the couple of names that yield a minimal distance.\n"
        doc += "\t   It comes to name after the closest atlas (for this division).\n"
        doc += "\t - 'sum': same as 'mean' \n"
        doc += "\t - 'majority': for each choice, order the distances in increasing order\n"
        doc += "\t   then compute the cumulated sum over the n first elements (so the first\n"
        doc += (
            "\t   is the minimum distance, and the last one is the average over all\n"
        )
        doc += "\t   atlases). Counts then the number of times when a choice is better than\n"
        doc += (
            "\t   the switched one, and choose the name couple with majority choices."
        )
        self.doc["selection_method"] = doc
        self.selection_method = "mean"

        doc = "\t Minimum number of atlases required to assess naming confidence. If there is not enough atlases\n"
        doc += "\t in the database for the aimed division, naming is not assessed."
        self.doc["confidence_atlases_nmin"] = doc
        self.confidence_atlases_nmin = 2
        doc = "\t Percentage of atlases used to assess naming confidence. If the percentage is less than\n"
        doc += (
            "\t 'confidence_atlases_nmin', 'confidence_atlases_nmin' atlases are used."
        )
        self.doc["confidence_atlases_percentage"] = doc
        self.confidence_atlases_percentage = 50

        #
        # for test:
        # names will be deleted, and tried to be rebuilt
        doc = "\t Input property file to be tested (must include cell names).\n"
        doc += "\t A 64-cells time point is searched and the embryo is renamed from\n"
        doc += (
            "\t this given time point. Comparison between new names and actual ones\n"
        )
        doc += "\t are reported.\n"
        doc += "\t If given, 'inputFile' is ignored.\n"
        self.doc["testFile"] = doc
        self.testFile = None
        doc = "\t If True, some diagnosis are conducted on the property file to\n"
        doc += "\t be named or to be tested.\n"
        self.doc["test_diagnosis"] = doc
        self.test_diagnosis = False

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print("#")
        print("# NamingParameters")
        print("#")
        print("")

        uatlasd.DivisionParameters.print_parameters(self)

        self.varprint("inputFile", self.inputFile, self.doc.get("inputFile", None))
        self.varprint("outputFile", self.outputFile, self.doc.get("outputFile", None))

        self.varprint(
            "selection_method",
            self.selection_method,
            self.doc.get("selection_method", None),
        )

        self.varprint("confidence_atlases_nmin", self.confidence_atlases_nmin)
        self.varprint(
            "confidence_atlases_percentage", self.confidence_atlases_percentage
        )

        self.varprint("testFile", self.testFile, self.doc.get("testFile", None))
        self.varprint(
            "test_diagnosis", self.test_diagnosis, self.doc.get("test_diagnosis", None)
        )
        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# NamingParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        uatlasd.DivisionParameters.write_parameters_in_file(self, logfile)

        self.varwrite(
            logfile, "inputFile", self.inputFile, self.doc.get("inputFile", None)
        )
        self.varwrite(
            logfile, "outputFile", self.outputFile, self.doc.get("outputFile", None)
        )

        self.varwrite(
            logfile,
            "selection_method",
            self.selection_method,
            self.doc.get("selection_method", None),
        )

        self.varwrite(
            logfile,
            "confidence_atlases_nmin",
            self.confidence_atlases_nmin,
            self.doc.get("confidence_atlases_nmin", None),
        )
        self.varwrite(
            logfile,
            "confidence_atlases_percentage",
            self.confidence_atlases_percentage,
            self.doc.get("confidence_atlases_percentage", None),
        )

        self.varwrite(
            logfile, "testFile", self.testFile, self.doc.get("testFile", None)
        )
        self.varwrite(
            logfile,
            "test_diagnosis",
            self.test_diagnosis,
            self.doc.get("test_diagnosis", None),
        )

        logfile.write("\n")

    def write_parameters(self, log_file_name):
        with open(log_file_name, "a") as logfile:
            self.write_parameters_in_file(logfile)
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameters):
        uatlasd.DivisionParameters.update_from_parameters(self, parameters)
        self.inputFile = self.read_parameter(parameters, "inputFile", self.inputFile)
        self.outputFile = self.read_parameter(parameters, "outputFile", self.outputFile)

        self.selection_method = self.read_parameter(
            parameters, "selection_method", self.selection_method
        )

        self.confidence_atlases_nmin = self.read_parameter(
            parameters, "confidence_atlases_nmin", self.confidence_atlases_nmin
        )
        self.confidence_atlases_percentage = self.read_parameter(
            parameters,
            "confidence_atlases_percentage",
            self.confidence_atlases_percentage,
        )

        self.testFile = self.read_parameter(parameters, "testFile", self.testFile)
        self.test_diagnosis = self.read_parameter(
            parameters, "test_diagnosis", self.test_diagnosis
        )

    def update_from_parameter_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            monitoring.to_log_and_console(
                "Error: '" + parameter_file + "' is not a valid file. Exiting."
            )
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
    returned_prop = {
        "cell_lineage": copy.deepcopy(prop["cell_lineage"]),
        "cell_name": copy.deepcopy(prop["cell_name"]),
        "cell_contact_surface": copy.deepcopy(prop["cell_contact_surface"]),
        "cell_volume": copy.deepcopy(prop["cell_volume"]),
    }

    #
    # get cells to count cells per time point
    #
    lineage = returned_prop["cell_lineage"]
    cells = list(
        set(lineage.keys()).union(
            set([v for values in list(lineage.values()) for v in values])
        )
    )
    cells = sorted(cells)

    div = 10**time_digits_for_cell_id
    ncells_per_time = np.zeros((max(cells) // div) + 1)
    for c in cells:
        t = int(c) // div
        ncells_per_time[t] += 1

    #
    # draw one time point of the desired number of cells
    #
    indices = np.where(np.array(ncells_per_time == ncells))[0]
    if len(indices) == 0:
        monitoring.to_log_and_console(
            str(proc) + ": there is no time point with #cells=" + str(ncells)
        )
        indices = np.where(np.array(ncells_per_time > ncells))[0]
        if len(indices) == 0:
            monitoring.to_log_and_console(
                str(proc) + ": there is no time point with #cells>" + str(ncells)
            )
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

    cells = list(returned_prop["cell_name"].keys())
    for c in cells:
        if int(c) // div != draw:
            del returned_prop["cell_name"][c]

    return returned_prop


def _test_naming(prop, reference_prop, discrepancies):
    proc = "_test_naming"

    monitoring.to_log_and_console("")
    monitoring.to_log_and_console("----- test naming vs reference -------")

    #
    # get the cell names with error
    #
    sisters_errors = {}
    other_errors = {}
    division_errors = {}
    key_list = sorted(prop["cell_name"].keys())
    for k in key_list:
        if k not in reference_prop["cell_name"]:
            monitoring.to_log_and_console(
                "\t weird, key " + str(k) + " is not in reference properties"
            )
        elif prop["cell_name"][k] != reference_prop["cell_name"][k]:
            #
            # named as its sister
            #
            if reference_prop["cell_name"][k] == uname.get_sister_name(
                prop["cell_name"][k]
            ):
                if prop["cell_name"][k] not in sisters_errors:
                    sisters_errors[prop["cell_name"][k]] = 1
                    msg = (
                        "cell_name["
                        + str(k)
                        + "]="
                        + str(prop["cell_name"][k])
                        + " is named as its sister "
                    )
                    msg += str(reference_prop["cell_name"][k])
                    if uname.get_mother_name(prop["cell_name"][k]) in discrepancies:
                        msg += ", division is incoherent among references"
                    monitoring.to_log_and_console("\t " + msg)
                else:
                    sisters_errors[prop["cell_name"][k]] += 1
                indexed_name = uname.get_mother_name(prop["cell_name"][k])
                division_errors[indexed_name] = division_errors.get(indexed_name, 0) + 1
            #
            # other errors
            #
            else:
                if prop["cell_name"][k] not in other_errors:
                    other_errors[prop["cell_name"][k]] = 1
                    msg = (
                        "cell_name["
                        + str(k)
                        + "]="
                        + str(prop["cell_name"][k])
                        + " is not equal to reference name "
                    )
                    msg += str(reference_prop["cell_name"][k])
                    msg += " for cell " + str(k)
                    if uname.get_mother_name(prop["cell_name"][k]) in discrepancies:
                        msg += ", division is incoherent among references"
                    monitoring.to_log_and_console("\t " + msg)
                else:
                    other_errors[prop["cell_name"][k]] += 1

    #
    # keep trace of the test result as a selection
    #
    keyleaveoneout = "morphonet_selection_leave_one_out_errors"
    lineage = prop["cell_lineage"]
    name = prop["cell_name"]
    prop[keyleaveoneout] = {}
    cells = list(
        set(lineage.keys()).union(
            set([v for values in list(lineage.values()) for v in values])
        )
    )

    for c in cells:
        if c not in name:
            continue
        if name[c] in other_errors:
            prop[keyleaveoneout][c] = 255
        elif name[c] in division_errors:
            prop[keyleaveoneout][c] = 100
        elif name[c] in sisters_errors:
            prop[keyleaveoneout][c] = 200

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
    key_list = sorted(reference_prop["cell_name"].keys())
    for k in key_list:
        indexed_name = reference_prop["cell_name"][k]
        reference_name[indexed_name] = reference_name.get(indexed_name, 0) + 1
        if k not in prop["cell_name"]:
            if reference_prop["cell_name"][k] not in name_missing:
                name_missing[reference_prop["cell_name"][k]] = 1
                msg = (
                    "reference cell_name["
                    + str(k)
                    + "]="
                    + str(reference_prop["cell_name"][k])
                    + " is not found"
                )
                monitoring.to_log_and_console("\t " + msg)
            else:
                name_missing[reference_prop["cell_name"][k]] += 1

    division_count = 0
    for k in division_errors:
        division_count += division_errors[k]
    other_count = 0
    for k in other_errors:
        other_count += other_errors[k]
    missing_count = 0
    for k in name_missing:
        missing_count += name_missing[k]

    msg = "ground-truth cells = " + str(len(reference_prop["cell_name"])) + " --- "
    msg += "ground-truth names = " + str(len(reference_name)) + " --- "
    msg += "tested cells = " + str(len(prop["cell_name"])) + " --- "
    msg += "retrieved names = " + str(len(reference_name) - len(name_missing)) + "\n"
    msg += (
        "\t missing cell names = "
        + str(missing_count)
        + " for "
        + str(len(name_missing))
        + " names --- \n"
    )
    msg += (
        "\t division errors in lineage = "
        + str(division_count)
        + " for "
        + str(len(division_errors))
        + " names \n"
    )
    if len(division_errors) > 0:
        msg += "\t    division errors = " + str(sorted(division_errors.keys())) + "\n"
    msg += (
        "\t other errors in lineage = "
        + str(other_count)
        + " for "
        + str(len(other_errors))
        + " names --- \n"
    )
    if len(other_errors) > 0:
        msg += "\t    other errors = " + str(sorted(other_errors.keys())) + "\n"
    monitoring.to_log_and_console("summary" + ": " + msg)

    monitoring.to_log_and_console("--------------------------------------")

    return prop


def _test_unequal_divisions(prop, atlases):
    proc = "_test_unequal_divisions"

    lineage = prop["cell_lineage"]
    name = prop["cell_name"]

    divisions = atlases.get_divisions()
    adivisions = atlases.get_unequal_divisions()
    volumes = atlases.get_volumes()

    mothercells = [c for c in lineage if len(lineage[c]) == 2]
    errors = []
    tests = []

    for m in mothercells:
        if m not in name:
            continue
        if name[m] not in adivisions:
            continue
        if lineage[m][0] not in name or lineage[m][1] not in name:
            continue
        #
        # here we've got an unequal named division
        #
        if (
            lineage[m][0] not in prop["cell_volume"]
            or lineage[m][1] not in prop["cell_volume"]
        ):
            msg = "weird, either " + str(lineage[m][0])
            msg += " or " + str(lineage[m][1]) + " is not in 'cell_volume' dictionary"
            monitoring.to_log_and_console(proc + ": " + msg)
            continue

        tests += [m]

        if prop["cell_volume"][lineage[m][0]] > prop["cell_volume"][lineage[m][1]]:
            if adivisions[name[m]][0] != name[lineage[m][0]]:
                errors += [m]
        elif prop["cell_volume"][lineage[m][1]] > prop["cell_volume"][lineage[m][0]]:
            if adivisions[name[m]][0] != name[lineage[m][1]]:
                errors += [m]

    if len(errors) > 0:
        error_names = [name[m] for m in errors]
        keyunequal = "morphonet_selection_unequal_division_errors"
        prop[keyunequal] = {}
        cells = list(
            set(lineage.keys()).union(
                set([v for values in list(lineage.values()) for v in values])
            )
        )
        for c in cells:
            if c not in name:
                continue
            if name[c] in error_names:
                prop[keyunequal][c] = 255

    monitoring.to_log_and_console("")
    monitoring.to_log_and_console("----- test unequal divisions -------")
    msg = "      test " + str(len(tests)) + " unequal divisions"
    msg += " over " + str(len(adivisions)) + " found in atlases "
    monitoring.to_log_and_console(msg)
    msg = "      found " + str(len(errors)) + " potential errors"
    monitoring.to_log_and_console(msg)
    if len(errors) > 0:
        for m in errors:
            msg = "    - division of cell " + str(m) + " (" + str(name[m]) + ") "
            monitoring.to_log_and_console(msg)
            if name[lineage[m][0]] == adivisions[name[m]][0]:
                msg = (
                    str(name[lineage[m][0]]) + " / " + str(name[lineage[m][1]]) + " = "
                )
                msg += (
                    str(prop["cell_volume"][lineage[m][0]])
                    + " / "
                    + str(prop["cell_volume"][lineage[m][1]])
                )
                msg += " = {:.4f}".format(
                    prop["cell_volume"][lineage[m][0]]
                    / prop["cell_volume"][lineage[m][1]]
                )
            else:
                msg = (
                    str(name[lineage[m][1]]) + " / " + str(name[lineage[m][0]]) + " = "
                )
                msg += (
                    str(prop["cell_volume"][lineage[m][1]])
                    + " / "
                    + str(prop["cell_volume"][lineage[m][0]])
                )
                msg += " = {:.4f}".format(
                    prop["cell_volume"][lineage[m][1]]
                    / prop["cell_volume"][lineage[m][0]]
                )
            monitoring.to_log_and_console("      " + msg)
            for r in divisions[name[m]]:
                msg = (
                    str(r)
                    + ": "
                    + str(adivisions[name[m]][0])
                    + " / "
                    + str(adivisions[name[m]][1])
                    + " = "
                )
                msg += (
                    str(volumes[adivisions[name[m]][0]][r])
                    + " / "
                    + str(volumes[adivisions[name[m]][1]][r])
                )
                msg += " = {:.4f}".format(
                    volumes[adivisions[name[m]][0]][r]
                    / volumes[adivisions[name[m]][1]][r]
                )
                monitoring.to_log_and_console("        - " + msg)
    monitoring.to_log_and_console("---------------------------------------")

    return prop


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


def _get_neighborhoods(
    mother,
    immediate_daughters,
    prop,
    embryo,
    parameters,
    delay_from_division=0,
    time_digits_for_cell_id=4,
):
    """

    Parameters
    ----------
    mother: cell id of the mother cell
    immediate_daughters: cell ids of the daughter cells
    prop: property dictionary of the embryo to be named
    embryo
    parameters:
    delay_from_division
    time_digits_for_cell_id:

    Returns
    -------
    A dictionary 'contact' indexed by [0, 1]. The value 'contact[i]' is the neighborhood of 'daughter[i]',
    a neighborhood being a contact surface vector indexed by cell names.
    None if there exists a neighbors without any named ancestor (including itself)
    """
    proc = "_get_neighborhoods"

    lineage = prop["cell_lineage"]
    reverse_lineage = {v: k for k, values in lineage.items() for v in values}
    #
    # since we change daughters, use a copy of immediate_daughters
    # else immediate_daughters will be changed
    #
    daughters = copy.deepcopy(immediate_daughters)

    #
    # delay from division
    #
    delay = 0
    if delay_from_division >= 0:
        delay = delay_from_division
    elif delay_from_division < 0:
        length0 = _get_branch_length(daughters[0], lineage)
        length1 = _get_branch_length(daughters[1], lineage)
        delay = min(length0, length1) + delay_from_division
        if delay < 0:
            delay = 0

    for i in range(delay):
        #
        # this test seems to be redundant with the next. Cells at the last time
        #
        if daughters[0] not in lineage:
            msg = ": build neighborhoods with a delay of " + str(i)
            msg += " instead of " + str(delay)
            msg += " for cells " + str(immediate_daughters)
            msg += " (" + str(daughters[0]) + " was not in lineage)"
            monitoring.to_log_and_console(str(proc) + msg, 3)
            break
        if daughters[1] not in lineage:
            msg = ": build neighborhoods with a delay of " + str(i)
            msg += " instead of " + str(delay)
            msg += " for cells " + str(immediate_daughters)
            msg += " (" + str(daughters[1]) + " was not in lineage)"
            monitoring.to_log_and_console(str(proc) + msg, 3)
            break
        if len(lineage[daughters[0]]) == 1 and len(lineage[daughters[1]]) == 1:
            daughters[0] = lineage[daughters[0]][0]
            daughters[1] = lineage[daughters[1]][0]
        else:
            msg = ": build neighborhoods with a delay of " + str(i)
            msg += " instead of " + str(delay_from_division)
            msg += " for cells " + str(immediate_daughters)
            monitoring.to_log_and_console(str(proc) + msg, 4)
            break

    contact = {}
    div = 10**time_digits_for_cell_id
    half_id = prop["cell_name"][mother][-1]

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

        for c in prop["cell_contact_surface"][d]:
            if int(c) % div == 1 or int(c) % div == 0:
                contact[index]["background"] = (
                    contact[index].get("background", 0)
                    + prop["cell_contact_surface"][d][c]
                )
            elif c in prop["cell_name"]:
                if parameters.differentiate_other_half:
                    contact[index][prop["cell_name"][c]] = prop["cell_contact_surface"][
                        d
                    ][c]
                else:
                    if prop["cell_name"][c][-1] == half_id:
                        contact[index][prop["cell_name"][c]] = prop[
                            "cell_contact_surface"
                        ][d][c]
                    else:
                        contact[index]["other-half"] = (
                            contact[index].get("other-half", 0)
                            + prop["cell_contact_surface"][d][c]
                        )
            elif c in daughters:
                if c != d:
                    contact[index]["sister"] = prop["cell_contact_surface"][d][c]
                else:
                    msg = ": weird, cell " + str(c) + " is in its contact surfaces"
                    monitoring.to_log_and_console(str(proc) + msg)
            else:
                #
                # cell without name, find a named ancestor
                #
                cell = c
                while cell in reverse_lineage and cell not in prop["cell_name"]:
                    cell = reverse_lineage[cell]
                if cell in prop["cell_name"]:
                    contact[index][prop["cell_name"][cell]] = (
                        contact[index].get(prop["cell_name"][cell], 0)
                        + prop["cell_contact_surface"][d][c]
                    )
                else:
                    msg = ": unable to find a named ancestor for cell " + str(c)
                    msg += ". Skip naming of cells " + str(immediate_daughters)
                    monitoring.to_log_and_console(str(proc) + msg)
                    return None

    #
    # normalize the neighborhoods
    #
    for i in contact:
        contact[i] = uatlasc.neighborhood_normalization(
            contact[i],
            embryo,
            int(daughters[0]) // div,
            cell_normalization=parameters.cell_normalization,
        )
    return contact


def _print_neighborhood(neighborhood, start=""):
    neighbors = sorted(list(neighborhood.keys()))
    n = 0
    txt = start + "{"
    if "background" in neighbors:
        txt += "'background': {:.3f}, ".format(neighborhood["background"])
        n += 1
        neighbors.remove("background")
    for i, neigh in enumerate(neighbors):
        txt += str(neigh) + ": {:.3f}".format(neighborhood[neigh])
        n += 1
        #  print(str(n) + " " + str(n % 4 == 0) + txt)
        if i < len(neighbors) - 1:
            txt += ", "
        elif i == len(neighbors) - 1:
            txt += "} "
        if n % 4 == 0 or i == len(neighbors) - 1:
            monitoring.to_log_and_console(txt)
            txt = "      "


def _compute_distances(
    mother,
    daughters,
    prop,
    embryo,
    atlases,
    parameters,
    delay_from_division=0,
    homogenize_neighborhood=False,
    time_digits_for_cell_id=4,
    debug=False,
):
    """

    Parameters
    ----------
    mother: cell id of the mother cell
    daughters: cell ids of the daughter cells
    prop: property dictionary of the embryo to be named
    embryo
    atlases:
    parameters:
    delay_from_division
    homogenize_neighborhood
    time_digits_for_cell_id
    debug

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
    if prop["cell_name"][mother] not in divisions:
        msg = (
            ": no reference neighborhoods for division of '"
            + str(prop["cell_name"][mother])
            + "'"
        )
        msg += ". Can not name cells " + str(daughters)
        msg += " from mother cell " + str(mother)
        monitoring.to_log_and_console(str(proc) + msg, 4)
        return None

    neighborhoods = atlases.get_cell_neighborhood(
        delay_from_division=delay_from_division
    )
    daughter_names = uname.get_daughter_names(prop["cell_name"][mother])

    #
    # returned the normalized neighborhood
    #
    contacts = _get_neighborhoods(
        mother,
        daughters,
        prop,
        embryo,
        parameters,
        delay_from_division=delay_from_division,
        time_digits_for_cell_id=time_digits_for_cell_id,
    )
    #
    # contacts[0] = contact vector for daughters[0]
    # contacts[1] = contact vector for daughters[1]
    # contacts[i]['sister'] is the contact surface of the sister
    #
    if contacts is None:
        msg = (
            ": can not extract contact vector for division of '"
            + str(prop["cell_name"][mother])
            + "'"
        )
        msg += ". Can not name cells " + str(daughters)
        msg += " from mother cell " + str(mother)
        monitoring.to_log_and_console(str(proc) + msg, 4)
        return None

    mother_generation = int(prop["cell_name"][mother].split(".")[0][1:])
    daughter_neighborhoods = {
        0: copy.deepcopy(neighborhoods[daughter_names[0]]),
        1: copy.deepcopy(neighborhoods[daughter_names[1]]),
    }
    if homogenize_neighborhood:
        new_neighborhoods = uneighborhood.build_same_contact_surfaces(
            daughter_neighborhoods,
            [0, 1],
            celltobeexcluded=daughter_names,
            maximal_generation=mother_generation,
        )
        for n in new_neighborhoods:
            daughter_neighborhoods[n] = copy.deepcopy(new_neighborhoods[n])

    if debug:
        monitoring.to_log_and_console("Distance computation:")
        msg = "- Neighborhood of cell '"
        if daughters[0] in prop["cell_name"]:
            msg += str(prop["cell_name"][daughters[0]])
        else:
            msg += str(daughters[0])
        monitoring.to_log_and_console(msg + "' in embryo to be tested is: ")
        _print_neighborhood(contacts[0], start="   ")
        msg = "- Neighborhood of cell '"
        if daughters[1] in prop["cell_name"]:
            msg += str(prop["cell_name"][daughters[1]])
        else:
            msg += str(daughters[1])
        monitoring.to_log_and_console(msg + "' in embryo to be tested is: ")
        _print_neighborhood(contacts[1], start="   ")

        monitoring.to_log_and_console(
            "- Neighborhoods of cell '" + str(daughter_names[0]) + "' in atlases are: "
        )
        for r in sorted(daughter_neighborhoods[0].keys()):
            _print_neighborhood(
                daughter_neighborhoods[0][r], start="   " + str(r) + ": "
            )
        monitoring.to_log_and_console(
            "- Neighborhoods of cell '" + str(daughter_names[1]) + "' in atlases are: "
        )
        for r in sorted(daughter_neighborhoods[1].keys()):
            _print_neighborhood(
                daughter_neighborhoods[1][r], start="   " + str(r) + ": "
            )

    innersurfaces = []
    if parameters.exclude_inner_surfaces:
        innersurfaces = daughter_names

    scores = {0: {}, 1: {}}
    for i in range(2):
        #
        # i = 0: test (contacts[0], contacts[1]) <-> (daughter_names[0], daughter_names[1])
        # i = 1: test (contacts[1], contacts[0]) <-> (daughter_names[0], daughter_names[1])
        #
        if i == 0:
            if "sister" in contacts[0]:
                contacts[0][daughter_names[1]] = contacts[0]["sister"]
                del contacts[0]["sister"]
            if "sister" in contacts[1]:
                contacts[1][daughter_names[0]] = contacts[1]["sister"]
                del contacts[1]["sister"]
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
        for ref in divisions[prop["cell_name"][mother]]:
            if i == 0:
                ic0 = 0
                ic1 = 1
            else:
                ic0 = 1
                ic1 = 0
            scores[i][ref] = uatlasd.division_distance(
                daughter_neighborhoods[0][ref],
                daughter_neighborhoods[1][ref],
                contacts[ic0],
                contacts[ic1],
                change_contact_surfaces=True,
                innersurfaces=innersurfaces,
            )

    if debug:
        for i in scores:
            msg = "- distance of ("
            if i == 0:
                if daughters[0] in prop["cell_name"]:
                    msg += str(prop["cell_name"][daughters[0]])
                else:
                    msg += str(daughters[0])
                msg += ", "
                if daughters[1] in prop["cell_name"]:
                    msg += str(prop["cell_name"][daughters[1]])
                else:
                    msg += str(daughters[1])
            elif i == 1:
                if daughters[1] in prop["cell_name"]:
                    msg += str(prop["cell_name"][daughters[1]])
                else:
                    msg += str(daughters[1])
                msg += ", "
                if daughters[0] in prop["cell_name"]:
                    msg += str(prop["cell_name"][daughters[0]])
                else:
                    msg += str(daughters[0])
            msg += (
                ") <-> (" + str(daughter_names[0]) + ", " + str(daughter_names[1]) + ")"
            )
            monitoring.to_log_and_console(msg)
            for r in sorted(scores[i].keys()):
                monitoring.to_log_and_console(
                    "   - atlas " + str(r) + ": " + str(scores[i][r])
                )
    return scores


########################################################################################
#
#
#
########################################################################################


def _get_ordered_means(distance):
    mean0 = []
    mean1 = []
    for i, v in enumerate(sorted(distance[0].values())):
        if i == 0:
            mean0 += [v]
        else:
            mean0 += [mean0[i - 1] + v]
    for i, v in enumerate(sorted(distance[1].values())):
        if i == 0:
            mean1 += [v]
        else:
            mean1 += [mean1[i - 1] + v]
    for i in range(len(mean1)):
        mean0[i] /= i + 1
        mean1[i] /= i + 1
    return mean0, mean1


def _get_name(daughters, daughter_names, distance, parameters, debug=False):
    """

    Parameters
    ----------
    daughters: cell ids of the daughter cells
    daughter_names: the 2 possible names
    distance: a dictionary of distances (in [0,1]) indexed by [d][reference_name] where
        - d:
          d = 0: test (contacts[0], contacts[1]) <-> (daughter_names[0], daughter_names[1])
          d = 1: test (contacts[1], contacts[0]) <-> (daughter_names[0], daughter_names[1])
          where contacts[0] = contact vector for daughters[0]
                contacts[1] = contact vector for daughters[1]
    parameters
    debug

    Returns
    -------
    a dictionary indexed by the cell ids of the daughter cells giving their names.

    """

    proc = "_get_name"

    if len(distance[0]) != len(distance[1]):
        monitoring.to_log_and_console(
            str(proc) + ": weird, atlases number are different"
        )

    # distance is a dictionary of dictionary
    # distance[d in 0,1][ref name]
    # d = 0: test (contacts[0], contacts[1]) <-> (daughter_names[0], daughter_names[1])
    # d = 1: test (contacts[1], contacts[0]) <-> (daughter_names[0], daughter_names[1])
    # ref name in embryos

    name = {}

    #
    # check whether the same decision will be made whatever the number of atlases
    #
    mean0, mean1 = _get_ordered_means(distance)
    count0less1 = sum([1 for i in zip(mean0, mean1) if i[0] <= i[1]])
    count1less0 = len(mean1) - count0less1

    if debug:
        msg = (
            "\t - "
            + str(count0less1)
            + " naming choices ("
            + str(daughters[0])
            + ", "
            + str(daughters[1])
            + ") = ("
        )
        msg += daughter_names[0] + ", " + daughter_names[1] + ") = ["
        for i, v in enumerate(mean0):
            msg += "{:1.4f}".format(v)
            if i < len(mean0) - 1:
                msg += ", "
        msg += "]"
        msg += " - sum = {:.4f}".format(np.sum(mean0))
        monitoring.to_log_and_console(msg)
        msg = (
            "\t - "
            + str(count1less0)
            + " naming choices ("
            + str(daughters[1])
            + ", "
            + str(daughters[0])
            + ") = ("
        )
        msg += daughter_names[0] + ", " + daughter_names[1] + ") = ["
        for i, v in enumerate(mean1):
            msg += "{:1.4f}".format(v)
            if i < len(mean1) - 1:
                msg += ", "
        msg += "]"
        msg += " - sum = {:.4f}".format(np.sum(mean1))
        monitoring.to_log_and_console(msg)

    if (
        parameters.selection_method.lower() == "minimum"
        or parameters.selection_method.lower() == "min"
    ):
        #
        # renvoie le tuple (ref, value) avec la valeur minimale
        #
        min0 = min(distance[0].items(), key=lambda x: x[1])
        min1 = min(distance[1].items(), key=lambda x: x[1])
        if min0[1] <= min1[1]:
            name[daughters[0]] = daughter_names[0]
            name[daughters[1]] = daughter_names[1]
            natlas = count0less1
        else:
            name[daughters[1]] = daughter_names[0]
            name[daughters[0]] = daughter_names[1]
            natlas = count1less0
        return name, natlas

    if (
        parameters.selection_method.lower() == "sum"
        or parameters.selection_method.lower() == "mean"
        or parameters.selection_method.lower() == "average"
    ):
        mean0 = statistics.mean([distance[0][a] for a in distance[0]])
        mean1 = statistics.mean([distance[1][a] for a in distance[1]])
        if mean0 <= mean1:
            name[daughters[0]] = daughter_names[0]
            name[daughters[1]] = daughter_names[1]
            natlas = count0less1
        else:
            name[daughters[1]] = daughter_names[0]
            name[daughters[0]] = daughter_names[1]
            natlas = count1less0
        return name, natlas

    if parameters.selection_method.lower() == "majority":
        if count0less1 > count1less0 or (
            count0less1 == count1less0 and np.sum(mean0) <= np.sum(mean1)
        ):
            name[daughters[0]] = daughter_names[0]
            name[daughters[1]] = daughter_names[1]
            natlas = count0less1
        else:
            name[daughters[1]] = daughter_names[0]
            name[daughters[0]] = daughter_names[1]
            natlas = count1less0
        return name, natlas

    monitoring.to_log_and_console(
        str(proc)
        + ": selection method '"
        + str(parameters.selection_method)
        + "' not handled yet"
    )
    return name, 0


########################################################################################
#
# naming procedure
#
########################################################################################


def _propagate_name_along_branch(prop, cell, keylist):
    proc = "_propagate_name_along_branch"
    lineage = prop["cell_lineage"]
    if cell not in prop["cell_name"]:
        return
    c = cell
    while c in lineage and len(lineage[c]) == 1:
        nc = lineage[c][0]
        if nc in prop["cell_name"]:
            if prop["cell_name"][nc] != prop["cell_name"][c]:
                msg = (
                    ": weird, cell "
                    + str(nc)
                    + " is named "
                    + str(prop["cell_name"][nc])
                )
                msg += ", but should be named " + str(prop["cell_name"][c])
                msg += " as its mother"
                monitoring.to_log_and_console(str(proc) + msg)
        else:
            prop["cell_name"][nc] = prop["cell_name"][c]
            for k in keylist:
                prop[k][nc] = prop[k][c]
        c = nc
    return


def _propagate_naming(prop, embryo, atlases, parameters, time_digits_for_cell_id=4):
    proc = "_propagate_naming"

    if "cell_lineage" not in prop:
        monitoring.to_log_and_console(
            str(proc) + ": 'cell_lineage' was not in dictionary"
        )
        return None

    if "cell_contact_surface" not in prop:
        monitoring.to_log_and_console(
            str(proc) + ": 'cell_contact_surface' was not in dictionary"
        )
        return None

    if "cell_name" not in prop:
        monitoring.to_log_and_console(str(proc) + ": 'cell_name' was not in dictionary")
        return None

    #
    #
    #
    lineage = prop["cell_lineage"]
    keyagreement = "morphonet_float_name_choice_agreement"
    prop[keyagreement] = {}

    #
    # remove empty names
    # leading or trailing spaces
    #
    cells = list(prop["cell_name"].keys())
    for c in cells:
        if prop["cell_name"][c] == "":
            del prop["cell_name"][c]
            continue
        prop["cell_name"][c] = prop["cell_name"][c].strip()

    #
    #
    #
    reverse_lineage = {v: k for k, values in lineage.items() for v in values}
    cells = list(
        set(lineage.keys()).union(
            set([v for values in list(lineage.values()) for v in values])
        )
    )

    #
    # backward propagation
    #
    monitoring.to_log_and_console(str(proc) + ": backward propagation")
    cells = sorted(cells, reverse=True)
    for c in cells:
        if c not in prop["cell_name"]:
            continue
        if c not in reverse_lineage:
            continue
        mother = reverse_lineage[c]
        if len(lineage[mother]) == 1:
            if mother in prop["cell_name"]:
                if prop["cell_name"][mother] != prop["cell_name"][c]:
                    msg = (
                        ": weird, cell "
                        + str(mother)
                        + " is named "
                        + str(prop["cell_name"][mother])
                    )
                    msg += ", but should be named " + str(prop["cell_name"][c])
                    msg += " as its single daughter"
                    monitoring.to_log_and_console(str(proc) + msg)
            else:
                prop["cell_name"][mother] = prop["cell_name"][c]
        elif len(lineage[mother]) == 2:
            mother_name = uname.get_mother_name(prop["cell_name"][c])
            if mother in prop["cell_name"]:
                if prop["cell_name"][mother] != mother_name:
                    msg = (
                        ": weird, cell "
                        + str(mother)
                        + " is named "
                        + str(prop["cell_name"][mother])
                    )
                    msg += ", but should be named " + str(mother_name)
                    msg += " since one of its daughter is named " + str(
                        prop["cell_name"][c]
                    )
                    monitoring.to_log_and_console(str(proc) + msg)
            else:
                prop["cell_name"][mother] = mother_name
        else:
            msg = (
                ": weird, cell "
                + str(mother)
                + " has "
                + str(len(lineage[mother]))
                + "daughter(s)"
            )
            monitoring.to_log_and_console(str(proc) + msg)

    #
    # forward propagation
    #
    monitoring.to_log_and_console(str(proc) + ": forward propagation")

    cells = sorted(cells)
    div = 10**time_digits_for_cell_id
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
        if c not in prop["cell_name"]:
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
            if c in prop["cell_name"]:
                continue
            # root of a tree, can not be named
            if c not in reverse_lineage:
                if c not in cell_not_named:
                    cell_not_named.append(c)
                monitoring.to_log_and_console(
                    str(proc) + ": weird, cell " + str(c) + " is root of a subtree"
                )
                continue

            # get its mother
            mother = reverse_lineage[c]
            # mother not in lineage
            # to account for lineage errors
            if mother not in lineage:
                monitoring.to_log_and_console(
                    str(proc) + ": weird, cell " + str(mother) + " is not in lineage"
                )
                continue
            # mother not named, can name the cell either
            if mother in cell_not_named:
                if c not in cell_not_named:
                    cell_not_named.append(c)
                continue
            # mother not named, can name the cell either
            if mother not in prop["cell_name"]:
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
                prop["cell_name"][c] = prop["cell_name"][mother]
                if mother in prop[keyagreement]:
                    prop[keyagreement][c] = prop[keyagreement][mother]
            #
            # in case of division:
            # 1. give name if the sister cell is named
            # 2. keep divisions to be solved
            #
            elif len(lineage[mother]) == 2:
                daughters = copy.deepcopy(lineage[mother])
                daughters.remove(c)
                daughter_names = uname.get_daughter_names(prop["cell_name"][mother])
                #
                # daughter cell is already named
                #
                if daughters[0] in prop["cell_name"]:
                    if prop["cell_name"][daughters[0]] in daughter_names:
                        daughter_names.remove(prop["cell_name"][daughters[0]])
                        prop["cell_name"][c] = daughter_names[0]
                    else:
                        msg = (
                            ": weird, cell "
                            + str(daughters[0])
                            + " is named "
                            + str(prop["cell_name"][daughters[0]])
                        )
                        msg += (
                            ", but should be named in "
                            + str(daughter_names)
                            + " since its mother cell "
                        )
                        msg += (
                            str(mother) + " is named " + str(prop["cell_name"][mother])
                        )
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
                msg = (
                    ": weird, cell "
                    + str(mother)
                    + " has "
                    + str(len(lineage[mother]))
                    + " daughter(s)."
                )
                msg += (
                    " Its offspring '" + str(lineage[mother]) + "' will not be named."
                )
                monitoring.to_log_and_console(str(proc) + msg)

        if division_to_be_named == {}:
            continue

        #
        # here we have a dictionary of divisions to be named (key = mother id, value = array of sister ids)
        #
        for mother, daughters in division_to_be_named.items():
            debug = False
            if (
                parameters.cells_to_be_traced is not None
                and prop["cell_name"][mother] in parameters.cells_to_be_traced
            ):
                debug = True
            #
            # distance is a dictionary of dictionary
            # distance[d in 0,1][ref name]
            # d = 0: test (contacts[0], contacts[1]) <-> (daughter_names[0], daughter_names[1])
            # d = 1: test (contacts[1], contacts[0]) <-> (daughter_names[0], daughter_names[1])
            # ref name in embryos
            #
            distance = _compute_distances(
                mother,
                daughters,
                prop,
                embryo,
                atlases,
                parameters,
                delay_from_division=parameters.name_delay_from_division,
                time_digits_for_cell_id=time_digits_for_cell_id,
                debug=debug,
            )
            if debug:
                monitoring.to_log_and_console("Naming distance = " + str(distance))
            if distance is None:
                for c in daughters:
                    if c not in cell_not_named:
                        cell_not_named.append(c)
                # msg = "\t error when building scores for daughters " + str(daughters) + " of cell " + str(mother)
                # monitoring.to_log_and_console(msg)
                msg = " Can not name cells " + str(daughters) + " and their offsprings."
                monitoring.to_log_and_console(str(proc) + ": " + msg)
                continue

            daughter_names = uname.get_daughter_names(prop["cell_name"][mother])
            name, natlas = _get_name(
                daughters, daughter_names, distance, parameters, debug=debug
            )

            for c in name:
                if name[c] is not None:
                    prop["cell_name"][c] = name[c]
            if 0 < natlas < len(distance[0]):
                for c in name:
                    prop[keyagreement][c] = natlas / len(distance[0])
                msg = "partial agreement to name division of cell " + str(
                    prop["cell_name"][mother]
                )
                msg += (
                    ": "
                    + str(natlas)
                    + " over "
                    + str(len(distance[0]))
                    + " means of ordered distances"
                )
                monitoring.to_log_and_console(str(proc) + ": " + msg)

    if len(prop[keyagreement]) == 0:
        del prop[keyagreement]

    return prop


########################################################################################
#
# evaluation procedure
#
########################################################################################


def _get_evaluation(given_names, daughter_names, distance, parameters, debug=False):
    proc = "_get_evaluation"

    # distance is a dictionary of dictionary
    # distance[d in 0,1][ref name]
    # d = 0: test (contacts[0], contacts[1]) <-> (daughter_names[0], daughter_names[1])
    # d = 1: test (contacts[1], contacts[0]) <-> (daughter_names[0], daughter_names[1])
    # ref name in embryos

    if given_names[0] == daughter_names[0] and given_names[1] == daughter_names[1]:
        dist = [(distance[0][a], distance[1][a]) for a in distance[0]]
    elif given_names[0] == daughter_names[1] and given_names[1] == daughter_names[0]:
        dist = [(distance[1][a], distance[0][a]) for a in distance[0]]
    else:
        msg = ": weird, names are not switched or non-switched ?!"
        monitoring.to_log_and_console(str(proc) + msg)
        return None, None

    sorted_dist = sorted(dist, key=lambda v: v[0])

    natlas = int(
        round(parameters.confidence_atlases_percentage * len(sorted_dist) / 100.0)
    )
    if parameters.confidence_atlases_nmin > 0:
        natlas = max(parameters.confidence_atlases_nmin, natlas)
    if natlas <= 0 or len(sorted_dist) < natlas:
        return None, None

    mean0 = 0.0
    mean1 = 0.0
    for i in range(natlas):
        mean0 += sorted_dist[i][0]
        mean1 += sorted_dist[i][1]
    mean0 /= natlas
    mean1 /= natlas
    return mean0, mean1 - mean0


def _evaluate_naming(prop, embryo, atlases, parameters, time_digits_for_cell_id=4):
    proc = "_evaluate_naming"

    if "cell_lineage" not in prop:
        monitoring.to_log_and_console(
            str(proc) + ": 'cell_lineage' was not in dictionary"
        )
        return None

    if "cell_contact_surface" not in prop:
        monitoring.to_log_and_console(
            str(proc) + ": 'cell_contact_surface' was not in dictionary"
        )
        return None

    if "cell_name" not in prop:
        monitoring.to_log_and_console(str(proc) + ": 'cell_name' was not in dictionary")
        return None

    #
    # initialize 'morphonet_float_name_choice_certainty'
    #
    keymindistance = "morphonet_float_name_choice_min_distance"
    hkeymindistance = "morphonet_float_name_choice_min_distance_homogeneize"
    prop[keymindistance] = {}
    prop[hkeymindistance] = {}

    #
    #
    #
    lineage = prop["cell_lineage"]
    reverse_lineage = {v: k for k, values in lineage.items() for v in values}
    cells = list(
        set(lineage.keys()).union(
            set([v for values in list(lineage.values()) for v in values])
        )
    )
    names = prop["cell_name"]

    #
    # forward propagation
    #
    monitoring.to_log_and_console(str(proc) + ": forward propagation")

    cells = sorted(cells)
    div = 10**time_digits_for_cell_id
    cells_per_time = {}
    for c in cells:
        t = int(c) // div
        #
        # get cells and cell names at each time point
        #
        if t not in cells_per_time:
            cells_per_time[t] = [c]
        else:
            cells_per_time[t].append(c)
    timepoints = sorted(cells_per_time.keys())

    processed_mothers = []

    for t in timepoints:
        for c in cells_per_time[t]:
            if c not in names:
                continue
            if c not in reverse_lineage:
                continue

            # get its mother
            mother = reverse_lineage[c]
            if mother in processed_mothers:
                continue
            processed_mothers += [mother]

            # mother not in lineage
            # to account for lineage errors
            if mother not in lineage:
                monitoring.to_log_and_console(
                    str(proc) + ": weird, cell " + str(mother) + " is not in lineage"
                )
                continue

            if len(lineage[mother]) == 1:
                if mother not in prop[keymindistance]:
                    continue
                prop[keymindistance][c] = prop[keymindistance][mother]
                if mother not in prop[hkeymindistance]:
                    continue
                prop[hkeymindistance][c] = prop[hkeymindistance][mother]
                continue

            if len(lineage[mother]) > 2:
                msg = (
                    ": weird, cell "
                    + str(mother)
                    + " divides in "
                    + str(len(lineage[mother]))
                    + " cells "
                )
                monitoring.to_log_and_console(str(proc) + msg)
                continue

            # here len(lineage[mother]) == 2, this is a division
            daughters = lineage[mother]
            if daughters[0] not in names or daughters[1] not in names:
                msg = ": weird, cells " + str(daughters) + " are not both named "
                monitoring.to_log_and_console(str(proc) + msg)
                continue

            #
            # distance is a dictionary of dictionary
            # distance[d in 0,1][ref name]
            # d = 0: test (contacts[0], contacts[1]) <-> (daughter_names[0], daughter_names[1])
            # d = 1: test (contacts[1], contacts[0]) <-> (daughter_names[0], daughter_names[1])
            # ref name in embryos
            #
            # certainty are probabilities (in [0, 1])
            #
            debug = (
                isinstance(parameters.cells_to_be_traced, list)
                and prop["cell_name"][mother] in parameters.cells_to_be_traced
            )

            daughter_names = uname.get_daughter_names(prop["cell_name"][mother])
            given_names = [names[daughters[0]], names[daughters[1]]]

            if debug:
                msg = "===== un-homogenized neighborhood for division of " + str(
                    prop["cell_name"][mother]
                )
                monitoring.to_log_and_console(msg)
            distance = _compute_distances(
                mother,
                daughters,
                prop,
                embryo,
                atlases,
                parameters,
                delay_from_division=parameters.confidence_delay_from_division,
                time_digits_for_cell_id=time_digits_for_cell_id,
                debug=debug,
            )
            if distance is None:
                continue
            mindist, difference = _get_evaluation(
                given_names, daughter_names, distance, parameters, debug=debug
            )

            mean0, mean1 = _get_ordered_means(distance)
            count0less1 = sum([1 for i in zip(mean0, mean1) if i[0] <= i[1]])
            if 0 < count0less1 < len(distance[0]):
                if given_names[0] == daughter_names[0]:
                    natlas = count0less1
                else:
                    natlas = len(distance[0]) - count0less1
                msg = "partial agreement to evaluate division of cell " + str(
                    prop["cell_name"][mother]
                )
                msg += (
                    ": "
                    + str(natlas)
                    + " over "
                    + str(len(distance[0]))
                    + " means of ordered distances"
                )
                msg += " (with un-homogenized neighborhood)"
                monitoring.to_log_and_console(str(proc) + ": " + msg)

            if difference is None:
                continue
            for d in daughters:
                prop[keymindistance][d] = mindist
            if debug:
                msg = "----- naming assessment with un-homogenized neighborhood for division of "
                msg += str(prop["cell_name"][mother])
                monitoring.to_log_and_console(msg)
                msg = "      difference = " + str(difference)
                monitoring.to_log_and_console(msg)
                msg = "      distance   = " + str(mindist)
                monitoring.to_log_and_console(msg)

            if debug:
                msg = "===== homogenized neighborhood for division of " + str(
                    prop["cell_name"][mother]
                )
                monitoring.to_log_and_console(msg)
            hdistance = _compute_distances(
                mother,
                daughters,
                prop,
                embryo,
                atlases,
                parameters,
                delay_from_division=parameters.confidence_delay_from_division,
                homogenize_neighborhood=True,
                time_digits_for_cell_id=time_digits_for_cell_id,
                debug=debug,
            )

            if hdistance is None:
                continue
            hmindist, hdifference = _get_evaluation(
                given_names, daughter_names, hdistance, parameters, debug=debug
            )

            mean0, mean1 = _get_ordered_means(hdistance)
            count0less1 = sum([1 for i in zip(mean0, mean1) if i[0] <= i[1]])
            if 0 < count0less1 < len(distance[0]):
                if given_names[0] == daughter_names[0]:
                    natlas = count0less1
                else:
                    natlas = len(distance[0]) - count0less1
                msg = "partial agreement to evaluate division of cell " + str(
                    prop["cell_name"][mother]
                )
                msg += (
                    ": "
                    + str(natlas)
                    + " over "
                    + str(len(distance[0]))
                    + " means of ordered distances"
                )
                msg += " (with homogenized neighborhood)"
                monitoring.to_log_and_console(str(proc) + ": " + msg)

            if hdifference is None:
                continue
            for d in daughters:
                prop[hkeymindistance][d] = hmindist
            if debug:
                msg = "----- naming assessment with homogenized neighborhood for division of "
                msg += str(prop["cell_name"][mother])
                monitoring.to_log_and_console(msg)
                msg = "      difference = " + str(hdifference)
                monitoring.to_log_and_console(msg)
                msg = "      distance   = " + str(hmindist)
                monitoring.to_log_and_console(msg)

    return prop


########################################################################################
#
#
#
########################################################################################


def _naming_diagnosis(prop, filename, parameters, time_digits_for_cell_id=4):
    diagnosis.monitoring.copy(monitoring)
    monitoring.to_log_and_console(
        "============================================================"
    )
    monitoring.to_log_and_console("===== diagnosis on '" + str(filename) + "'")
    diagnosis.diagnosis(
        prop,
        features=["name"],
        parameters=parameters,
        time_digits_for_cell_id=time_digits_for_cell_id,
    )
    monitoring.to_log_and_console(
        "============================================================"
    )


def naming_process(experiment, parameters):
    proc = "naming_process"

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(
            str(proc)
            + ": unexpected type for 'experiment' variable: "
            + str(type(experiment))
        )
        sys.exit(1)

    if not isinstance(parameters, NamingParameters):
        monitoring.to_log_and_console(
            str(proc)
            + ": unexpected type for 'parameters' variable: "
            + str(type(parameters))
        )
        sys.exit(1)

    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()

    uatlasd.monitoring.copy(monitoring)
    uatlasc.monitoring.copy(monitoring)

    #
    # should we clean reference here?
    #
    atlases = uatlasd.DivisionAtlases(parameters=parameters)
    atlases.add_atlases(
        parameters.atlasFiles,
        parameters,
        time_digits_for_cell_id=time_digits_for_cell_id,
    )
    atlases.build_division_atlases(parameters)

    #
    # read input properties to be named
    #
    prop = {}
    reference_prop = {}
    discrepancies = {}

    if parameters.testFile is not None:
        reference_prop = ioproperties.read_dictionary(
            parameters.testFile, inputpropertiesdict={}
        )
        if parameters.test_diagnosis:
            _naming_diagnosis(
                reference_prop,
                parameters.testFile,
                parameters,
                time_digits_for_cell_id=time_digits_for_cell_id,
            )
        prop = _build_test_set(
            reference_prop, time_digits_for_cell_id=time_digits_for_cell_id, ncells=64
        )
        if prop is None:
            monitoring.to_log_and_console(str(proc) + ": error when building test set")
            sys.exit(1)
    elif parameters.inputFile is not None:
        prop = ioproperties.read_dictionary(
            parameters.inputFile, inputpropertiesdict={}
        )
        if parameters.test_diagnosis:
            _naming_diagnosis(
                prop,
                parameters.inputFile,
                parameters,
                time_digits_for_cell_id=time_digits_for_cell_id,
            )

    if prop == {}:
        monitoring.to_log_and_console(str(proc) + ": no properties?!")
        sys.exit(1)

    if "cell_name" not in prop:
        monitoring.to_log_and_console(
            str(proc) + ": no 'cell_name' in input dictionary"
        )
        sys.exit(1)

    # clean from empty names
    cells = list(prop["cell_name"].keys())
    for c in cells:
        if prop["cell_name"][c] == "":
            del prop["cell_name"][c]

    #
    # build an atlas from the embryo to be named,
    # temporally align it with the reference embryo
    #
    embryo = uatlase.Atlas(prop)
    ref_atlases = atlases.get_atlases()
    ref_atlas = atlases.get_reference_atlas()
    embryo.temporally_align_with(
        ref_atlases[ref_atlas], time_digits_for_cell_id=time_digits_for_cell_id
    )

    #
    # naming propagation
    #
    prop = _propagate_naming(
        prop,
        embryo,
        atlases,
        parameters,
        time_digits_for_cell_id=time_digits_for_cell_id,
    )
    prop = _evaluate_naming(
        prop,
        embryo,
        atlases,
        parameters,
        time_digits_for_cell_id=time_digits_for_cell_id,
    )
    prop = properties.set_fate_from_names(prop)
    prop = properties.set_color_from_fate(prop)
    #
    #
    #
    prop = _test_unequal_divisions(prop, atlases)
    #
    #
    #
    if parameters.testFile is not None:
        prop = _test_naming(prop, reference_prop, discrepancies)

    if isinstance(parameters.outputFile, str):
        ioproperties.write_dictionary(parameters.outputFile, prop)

    return prop
