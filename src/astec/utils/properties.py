import os
import sys

import numpy as np
import math

from operator import itemgetter

from astec.utils import ioproperties
from astec.utils import common
import astec.utils.ascidian_name as uname
from astec.wrapping import cpp_wrapping

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


class CellPropertiesParameters(common.PrefixedParameter):
    def __init__(self, prefix=None):
        common.PrefixedParameter.__init__(self, prefix=prefix)

        if "doc" not in self.__dict__:
            self.doc = {}

        self.max_chunks_properties = None
        self.sigma_segmentation_smoothing = None

    def print_parameters(self):
        print("")
        print("#")
        print("# CellPropertiesParameters")
        print("#")
        print("")

        common.PrefixedParameter.print_parameters(self)

        self.varprint("max_chunks_properties", self.max_chunks_properties)
        self.varprint("sigma_segmentation_smoothing", self.sigma_segmentation_smoothing)

        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# CellPropertiesParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, "max_chunks_properties", self.max_chunks_properties)
        self.varwrite(
            logfile, "sigma_segmentation_smoothing", self.sigma_segmentation_smoothing
        )

        logfile.write("\n")
        return

    def write_parameters(self, log_file_name):
        with open(log_file_name, "a") as logfile:
            self.write_parameters_in_file(logfile)
        return

    def update_from_parameters(self, parameters):
        self.max_chunks_properties = self.read_parameter(
            parameters, "max_chunks_properties", self.max_chunks_properties
        )
        self.sigma_segmentation_smoothing = self.read_parameter(
            parameters,
            "sigma_segmentation_smoothing",
            self.sigma_segmentation_smoothing,
        )

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
# properties computation from a sequence
# calculation is done after sequence intra-registration
#
########################################################################################


def property_computation(experiment, parameters):
    """

    :param experiment:
    :param parameters:
    :return:
    """

    proc = "property_computation"

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(
            str(proc)
            + ": unexpected type for 'experiment' variable: "
            + str(type(experiment))
        )
        sys.exit(1)

    if not isinstance(parameters, CellPropertiesParameters):
        monitoring.to_log_and_console(
            str(proc)
            + ": unexpected type for 'parameters' variable: "
            + str(type(parameters))
        )
        sys.exit(1)
    #
    # get directory name where to find co-registered images of the sequence
    # as well as the common image suffix
    #

    intrareg_path = os.path.join(
        experiment.intrareg_dir.get_directory(), experiment.post_dir.get_sub_directory()
    )
    #
    # is there a post-segmentation directory in the intra-registration directory ?
    #
    if not os.path.isdir(intrareg_path):
        monitoring.to_log(proc + ": '" + str(intrareg_path) + "' does not exist")
        intrareg_path = os.path.join(
            experiment.intrareg_dir.get_directory(),
            experiment.astec_dir.get_sub_directory(),
        )
        #
        # if no, is there a segmentation directory in the intra-registration directory ?
        #
        if not os.path.isdir(intrareg_path):
            monitoring.to_log(proc + ": '" + str(intrareg_path) + "' does not exist")
            intrareg_path = experiment.intrareg_dir.get_directory()
            monitoring.to_log_and_console(
                proc
                + ": neither POST/ or SEG/ sub-directories in '"
                + str(intrareg_path)
                + "'",
                0,
            )
            monitoring.to_log_and_console("Exiting.", 0)
            return None
        else:
            working_dir = experiment.astec_dir
    else:
        working_dir = experiment.post_dir

    monitoring.to_log_and_console(
        "... will compute sequence properties from '" + str(intrareg_path) + "'", 0
    )

    #
    # build name format for (post-corrected) segmentation images
    #
    name_format = (
        experiment.intrareg_dir.get_file_prefix()
        + experiment.intrareg_dir.get_file_suffix()
        + working_dir.get_file_suffix()
        + experiment.intrareg_dir.get_time_prefix()
        + experiment.get_time_format()
    )

    suffix = common.get_file_suffix(
        experiment, intrareg_path, name_format, flag_time=experiment.get_time_format()
    )
    if suffix is None:
        monitoring.to_log_and_console(
            proc + ": no consistent naming was found in '" + str(intrareg_path) + "'", 1
        )
        monitoring.to_log_and_console("Exiting.", 0)
        sys.exit(1)
    name_format += "." + str(suffix)
    template_format = os.path.join(intrareg_path, name_format)

    #
    #
    #
    output_name = (
        experiment.intrareg_dir.get_file_prefix()
        + experiment.intrareg_dir.get_file_suffix()
        + working_dir.get_file_suffix()
        + "_lineage"
    )
    output_name = os.path.join(intrareg_path, output_name)

    if os.path.isfile(output_name + ".xml") and os.path.isfile(output_name + ".pkl"):
        if not monitoring.forceResultsToBeBuilt:
            monitoring.to_log_and_console("    xml file already existing", 2)
            return output_name + ".xml"
        else:
            monitoring.to_log_and_console(
                "    xml file already existing, but forced", 2
            )

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    cpp_wrapping.cell_properties(
        template_format,
        output_name + ".xml",
        first_time_point,
        last_time_point,
        diagnosis_file=output_name + ".txt",
        n_processors=parameters.max_chunks_properties,
        cell_based_sigma=parameters.sigma_segmentation_smoothing,
        monitoring=monitoring,
    )

    return output_name + ".xml"


########################################################################################
#
# comparison of two dictionaries
#
########################################################################################


def _intersection_lists(l1, l2, name1, name2):
    """

    :param e1:
    :param e2:
    :param name1:
    :param name2:
    :return:
    """
    intersection = sorted(list(set(l1).intersection(set(l2))))
    difference1 = sorted(list(set(l1).difference(set(l2))))
    difference2 = sorted(list(set(l2).difference(set(l1))))

    monitoring.to_log_and_console(
        "    ... " + str(len(list(l1))) + " key cells are in '" + str(name1) + "'"
    )
    monitoring.to_log_and_console(
        "    ... " + str(len(list(l2))) + " key cells are in '" + str(name2) + "'"
    )
    if len(difference1) > 0:
        monitoring.to_log_and_console(
            "    ... "
            + str(len(difference1))
            + " key cells are in '"
            + str(name1)
            + "' and are not in '"
            + str(name2)
            + "'",
            1,
        )
        s = repr(difference1)
        monitoring.to_log_and_console("        " + s, 1)

    if len(difference2) > 0:
        monitoring.to_log_and_console(
            "    ... "
            + str(len(difference2))
            + " key cells are not in '"
            + str(name1)
            + "' but are in '"
            + str(name2)
            + "'",
            1,
        )
        s = repr(difference2)
        monitoring.to_log_and_console("        " + s, 1)

    return intersection


def _intersection_cell_keys(e1, e2, name1, name2):
    return _intersection_lists(list(e1.keys()), list(e2.keys()), name1, name2)


def _compare_lineage(d1, d2, name1, name2, description):
    # e1 and e2 are respectively the lineages of name1 and name2
    e1 = d1[description]
    e2 = d2[description]
    msg = (
        "  === "
        + str(description)
        + " comparison between "
        + str(name1)
        + " and "
        + str(name2)
        + " === "
    )
    monitoring.to_log_and_console(msg, 1)

    intersection = _intersection_cell_keys(e1, e2, name1, name2)
    if len(intersection) > 0:
        n = 0
        for k in intersection:
            if sorted(e1[k]) != sorted(e2[k]):
                n += 1

        if n > 0:
            monitoring.to_log_and_console(
                "    ... " + str(n) + " cells have different lineages", 1
            )
            for k in intersection:
                if sorted(e1[k]) != sorted(e2[k]):
                    s = "cell " + str(k) + " has different lineages: "
                    monitoring.to_log_and_console("        " + s, 1)
                    s = str(e1[k]) + " in '" + str(name1) + "'"
                    monitoring.to_log_and_console("          " + s, 1)
                    s = str(e2[k]) + " in '" + str(name2) + "'"
                    monitoring.to_log_and_console("          " + s, 1)

    monitoring.to_log_and_console("", 1)
    msg = (
        "    = "
        + str(description)
        + " comparison between daughter cells in "
        + str(name1)
        + " and "
        + str(name2)
    )
    monitoring.to_log_and_console(msg, 1)
    #
    # reverse lineage
    #
    reverse_e1 = {}
    for k, values in e1.items():
        for v in values:
            reverse_e1[v] = reverse_e1.get(v, []) + [k]

    reverse_e2 = {}
    for k, values in e2.items():
        for v in values:
            reverse_e2[v] = reverse_e2.get(v, []) + [k]

    intersection = _intersection_cell_keys(
        reverse_e1, reverse_e2, "reversed " + name1, "reversed " + name2
    )
    if len(intersection) > 0:
        n = 0
        for k in intersection:
            if reverse_e1[k] != reverse_e2[k]:
                n += 1

        if n > 0:
            monitoring.to_log_and_console(
                "    ... " + str(n) + " cells have different mother cell", 1
            )
            for k in intersection:
                if reverse_e1[k] != reverse_e2[k]:
                    s = "cell " + str(k) + " has different mother cell: "
                    monitoring.to_log_and_console("        " + s, 1)
                    s = str(reverse_e1[k]) + " in '" + str(name1) + "'"
                    monitoring.to_log_and_console("          " + s, 1)
                    s = str(reverse_e2[k]) + " in '" + str(name2) + "'"
                    monitoring.to_log_and_console("          " + s, 1)

    monitoring.to_log_and_console("", 1)
    msg = (
        "    = "
        + str(description)
        + " comparison between all cells in "
        + str(name1)
        + " and "
        + str(name2)
    )
    monitoring.to_log_and_console(msg, 1)
    l1 = list(
        set(e1.keys()).union(set([v for values in list(e1.values()) for v in values]))
    )
    l2 = list(
        set(e2.keys()).union(set([v for values in list(e2.values()) for v in values]))
    )
    _intersection_lists(l1, l2, name1, name2)

    return


# def _compare_h_min(d1, d2, name1, name2, description):
#     return


def _compare_volume(d1, d2, name1, name2, description):
    """

    :param d1:
    :param d2:
    :param name1:
    :param name2:
    :param description:
    :return:
    """
    #     dictionary de int
    #     cell_volume.590002 = <type 'int'>
    #     590002: 236936

    e1 = d1[description]
    e2 = d2[description]

    msg = (
        "  === "
        + str(description)
        + " comparison between "
        + str(name1)
        + " and "
        + str(name2)
        + " === "
    )
    monitoring.to_log_and_console(msg, 1)

    intersection = _intersection_cell_keys(e1, e2, name1, name2)
    if len(intersection) > 0:
        n = 0
        for k in intersection:
            if e1[k] != e2[k]:
                n += 1
        if n > 0:
            monitoring.to_log_and_console(
                "    ... " + str(n) + " cells have different volumes", 1
            )
            for k in intersection:
                if e1[k] != e2[k]:
                    s = "cell " + str(k) + " has different volumes: "
                    monitoring.to_log_and_console("        " + s, 1)
                    s = str(e1[k]) + " in " + str(name1)
                    monitoring.to_log_and_console("          " + s, 1)
                    s = str(e2[k]) + " in " + str(name2)
                    monitoring.to_log_and_console("          " + s, 1)
    return


# def _compare_sigma(d1, d2, name1, name2, description):
#     return


# def _compare_label_in_time(d1, d2, name1, name2, description):
#     return


def _compare_barycenter(d1, d2, name1, name2, description):
    """

    :param d1:
    :param d2:
    :param name1:
    :param name2:
    :param description:
    :return:
    """

    # 'barycenter': 'cell_barycenter'
    #     dictionary de numpy.ndarray de numpy.float64
    #     cell_barycenter.590002 = <type 'numpy.ndarray'>
    #     590002: array([ 258.41037242,  226.74975943,  303.67167927])

    e1 = d1[description]
    e2 = d2[description]

    monitoring.to_log_and_console("  === " + str(description) + " comparison === ", 1)

    intersection = _intersection_cell_keys(e1, e2, name1, name2)

    residual = []

    if len(intersection) > 0:
        for k in intersection:
            residual.append([k, math.sqrt(sum((e1[k] - e2[k]) * (e1[k] - e2[k])))])

    residual = sorted(residual, key=itemgetter(1), reverse=True)
    monitoring.to_log_and_console(
        "    ... largest residual at cell #"
        + str(residual[0][0])
        + " is "
        + str(residual[0][1]),
        1,
    )
    s = str(e1[residual[0][0]]) + " <-> " + str(e2[residual[0][0]])
    monitoring.to_log_and_console("        " + s, 1)

    # for i in range(min(len(intersection),10)):
    #     print "#" + str(i) + ": " + str(residual[i])

    return


def _are_fates_equal(fate1, fate2):
    if isinstance(fate1, str):
        if isinstance(fate2, str):
            if fate1 != fate2:
                return False
        elif isinstance(fate2, list):
            if fate1 not in fate2 or len(set(fate2)) > 1:
                return False
        else:
            return False
    elif isinstance(fate1, list):
        if isinstance(fate2, str):
            if fate2 not in fate1 or len(set(fate1)) > 1:
                return False
        elif isinstance(fate2, list):
            if set(fate1) != set(fate2):
                return False
        else:
            return False
    else:
        return False
    return True


def _compare_fate(d1, d2, name1, name2, description):
    e1 = d1[description]
    e2 = d2[description]

    lineage1 = d1["cell_lineage"]
    reverse_lineage1 = {v: k for k, values in lineage1.items() for v in values}
    lineage2 = d2["cell_lineage"]
    reverse_lineage2 = {v: k for k, values in lineage2.items() for v in values}

    msg = (
        "  === "
        + str(description)
        + " comparison between "
        + str(name1)
        + " and "
        + str(name2)
        + " === "
    )
    monitoring.to_log_and_console(msg, 1)

    intersection = _intersection_cell_keys(e1, e2, name1, name2)
    if len(intersection) == 0:
        return

    n = 0
    for k in intersection:
        if not _are_fates_equal(e1[k], e2[k]):
            n += 1

    if n == 0:
        return
    monitoring.to_log_and_console(
        "    ... " + str(n) + " cells have different fates", 1
    )

    return


def _compare_all_cells(d1, d2, name1, name2, description):
    """

    :param d1:
    :param d2:
    :param name1:
    :param name2:
    :param description:
    :return:
    """

    # 'all-cells': 'all_cells'  # liste de toutes les cellules ?
    #     liste de numpy.int64
    #     all_cells = <type 'list'>

    e1 = d1[description]
    e2 = d2[description]

    monitoring.to_log_and_console("  === " + str(description) + " comparison === ", 1)

    difference1 = list(set(e1).difference(set(e2)))
    difference2 = list(set(e2).difference(set(e1)))

    if len(difference1) > 0:
        monitoring.to_log_and_console(
            "    ... cells that are in '"
            + str(name1)
            + "' and not in '"
            + str(name2)
            + "'",
            1,
        )
        s = repr(difference1)
        monitoring.to_log_and_console("        " + s, 1)

    if len(difference2) > 0:
        monitoring.to_log_and_console(
            "    ... cells that are not in '"
            + str(name1)
            + "' but in '"
            + str(name2)
            + "'",
            1,
        )
        s = repr(difference2)
        monitoring.to_log_and_console("        " + s, 1)

    return


def _compare_principal_value(d1, d2, name1, name2, description):
    """

    :param d1:
    :param d2:
    :param name1:
    :param name2:
    :param description:
    :return:
    """

    # 'principal-value': 'cell_principal_values'
    #     dictionary de liste de numpy.float64
    #     cell_principal_values.590002 = <type 'list'>
    #     590002: [1526.0489371146978, 230.60881177650205, 91.063513300019849]

    e1 = d1[description]
    e2 = d2[description]

    monitoring.to_log_and_console("  === " + str(description) + " comparison === ", 1)

    intersection = _intersection_cell_keys(e1, e2, name1, name2)

    residual = []

    if len(intersection) > 0:
        for k in intersection:
            residual.append([k, max(abs(np.array(e1[k]) - np.array(e2[k])))])

    residual = sorted(residual, key=itemgetter(1), reverse=True)
    monitoring.to_log_and_console(
        "    ... largest residual at cell #"
        + str(residual[0][0])
        + " is "
        + str(residual[0][1]),
        1,
    )
    s = str(e1[residual[0][0]]) + " <-> " + str(e2[residual[0][0]])
    monitoring.to_log_and_console("        " + s, 1)

    # for i in range(min(len(intersection),10)):
    #     print "#" + str(i) + ": " + str(residual[i])

    return


def _compare_name(d1, d2, name1, name2, description):
    e1 = d1[description]
    e2 = d2[description]

    lineage1 = d1["cell_lineage"]
    reverse_lineage1 = {v: k for k, values in lineage1.items() for v in values}
    lineage2 = d2["cell_lineage"]
    reverse_lineage2 = {v: k for k, values in lineage2.items() for v in values}

    msg = (
        "  === "
        + str(description)
        + " comparison between "
        + str(name1)
        + " and "
        + str(name2)
        + " === "
    )
    monitoring.to_log_and_console(msg, 1)

    intersection = _intersection_cell_keys(e1, e2, name1, name2)
    if len(intersection) == 0:
        return

    n = 0
    for k in intersection:
        if e1[k] != e2[k]:
            n += 1
    if n == 0:
        return
    monitoring.to_log_and_console(
        "    ... " + str(n) + " cells have different names", 1
    )

    cell_starting = []
    cell_with_diff_mothers = []
    cell_before_change = []

    for k in intersection:
        if e1[k] == e2[k]:
            continue
        if k not in reverse_lineage1 or k not in reverse_lineage2:
            cell_starting += [k]
            continue
        if reverse_lineage1[k] != reverse_lineage2[k]:
            cell_with_diff_mothers += [k]
            continue
        if e1[reverse_lineage1[k]] == e2[reverse_lineage2[k]]:
            cell_before_change += [reverse_lineage1[k]]
    cell_before_change = sorted(list(set(cell_before_change)))

    if len(cell_starting) > 0:
        s = str(len(cell_starting)) + " branches starting with different names"
        monitoring.to_log_and_console("    ... " + s, 1)
        for k in cell_starting:
            s = "cell " + str(k) + " (and cells of its branch) has different names: "
            monitoring.to_log_and_console("        " + s, 1)
            s = str(e1[k]) + " in '" + str(name1) + "'"
            monitoring.to_log_and_console("          " + s, 1)
            s = str(e2[k]) + " in '" + str(name2) + "'"
            monitoring.to_log_and_console("          " + s, 1)

    if len(cell_with_diff_mothers) > 0:
        s = (
            str(len(cell_with_diff_mothers))
            + " cells with different names and different mothers"
        )
        monitoring.to_log_and_console("    ... " + s, 1)
        for k in cell_with_diff_mothers:
            s = "cell " + str(k) + " has different names and mothers: "
            monitoring.to_log_and_console("        " + s, 1)
            s = (
                "named "
                + str(e1[k])
                + " and issued from "
                + str(reverse_lineage1[k])
                + " in '"
                + str(name1)
                + "'"
            )
            monitoring.to_log_and_console("          " + s, 1)
            s = (
                "named "
                + str(e2[k])
                + " and issued from "
                + str(reverse_lineage2[k])
                + " in '"
                + str(name2)
                + "'"
            )
            monitoring.to_log_and_console("          " + s, 1)

    if len(cell_before_change) > 0:
        s = str(len(cell_before_change)) + " cells with differently named daughters"
        monitoring.to_log_and_console("    ... " + s, 1)
        division_switches = []
        other_cases = []
        for c in cell_before_change:
            daughters = sorted(list(set(lineage1[c]).intersection(set(lineage2[c]))))
            if len(lineage1[c]) == 2 and len(lineage2[c]) == 2 and len(daughters) == 2:
                if (
                    e1[daughters[0]] == e2[daughters[1]]
                    and e1[daughters[1]] == e2[daughters[0]]
                ):
                    division_switches += [c]
                    continue
            other_cases += [c]
        if len(division_switches) > 0:
            s = str(len(division_switches)) + " divisions with switched names"
            monitoring.to_log_and_console("      - " + s, 1)
            for c in division_switches:
                s = (
                    "cell "
                    + str(c)
                    + " ("
                    + str(e1[c])
                    + ") has daughters with different names:"
                )
                monitoring.to_log_and_console("        " + s, 1)
                daughters = sorted(
                    list(set(lineage1[c]).intersection(set(lineage2[c])))
                )
                l1 = ""
                l2 = ""
                for i, d in enumerate(daughters):
                    l1 += str(e1[d])
                    l2 += str(e2[d])
                    if i < len(daughters) - 1:
                        l1 += ", "
                        l2 += ", "
                s = (
                    str(daughters)
                    + " are named ["
                    + str(l1)
                    + "] in '"
                    + str(name1)
                    + "'"
                )
                monitoring.to_log_and_console("          " + s, 1)
                s = (
                    " " * len(str(daughters))
                    + " and named ["
                    + str(l2)
                    + "] in '"
                    + str(name2)
                    + "'"
                )
                monitoring.to_log_and_console("          " + s, 1)

        if len(other_cases) > 0:
            s = str(len(other_cases)) + " other cases"
            monitoring.to_log_and_console("      - " + s, 1)
            for c in other_cases:
                s = (
                    "cell "
                    + str(c)
                    + " ("
                    + str(e1[c])
                    + ") has daughters with different names:"
                )
                monitoring.to_log_and_console("        " + s, 1)
                daughters = sorted(
                    list(set(lineage1[c]).intersection(set(lineage2[c])))
                )
                difference1 = sorted(
                    list(set(lineage1[c]).difference(set(lineage2[c])))
                )
                difference2 = sorted(
                    list(set(lineage2[c]).difference(set(lineage1[c])))
                )

                l1 = ""
                l2 = ""
                for i, d in enumerate(daughters):
                    l1 += str(e1[d])
                    l2 += str(e2[d])
                    if i < len(daughters) - 1:
                        l1 += ", "
                        l2 += ", "
                s = (
                    str(daughters)
                    + " are named ["
                    + str(l1)
                    + "] in '"
                    + str(name1)
                    + "'"
                )
                monitoring.to_log_and_console("          " + s, 1)
                s = (
                    " " * len(str(daughters))
                    + " and named ["
                    + str(l2)
                    + "] in '"
                    + str(name2)
                    + "'"
                )
                monitoring.to_log_and_console("          " + s, 1)

                for d in difference1:
                    spaces = ""
                    for i in range(len(str(d))):
                        spaces += " "
                    s = str(d) + " is named " + str(e1[d]) + " in '" + str(name1) + "'"
                    monitoring.to_log_and_console("          " + s, 1)
                    s = (
                        spaces
                        + "and is not a daughter of "
                        + str(c)
                        + " in '"
                        + str(name2)
                        + "'"
                    )
                    monitoring.to_log_and_console("          " + s, 1)
                for d in difference2:
                    spaces = ""
                    for i in range(len(str(d))):
                        spaces += " "
                    s = str(d) + " is named " + str(e2[d]) + " in '" + str(name2) + "'"
                    monitoring.to_log_and_console("          " + s, 1)
                    s = (
                        spaces
                        + "and is not a daughter of "
                        + str(c)
                        + " in '"
                        + str(name1)
                        + "'"
                    )
                    monitoring.to_log_and_console("          " + s, 1)
    return


def _compare_contact(d1, d2, name1, name2, description):
    """

    :param d1:
    :param d2:
    :param name1:
    :param name2:
    :param description:
    :return:
    """
    # 'contact': 'cell_contact_surface',
    #     dictionary de dictionary de int
    #     cell_contact_surface.590002.590019 = <type 'int'>

    e1 = d1[description]
    e2 = d2[description]

    msg = (
        "  === "
        + str(description)
        + " comparison between "
        + str(name1)
        + " and "
        + str(name2)
        + " === "
    )
    monitoring.to_log_and_console(msg, 1)

    intersection = _intersection_cell_keys(e1, e2, name1, name2)
    if len(intersection) > 0:
        n = 0
        for k in intersection:
            d = list(set(e1[k].keys()).symmetric_difference(set(e2[k].keys())))
            if len(d) > 0:
                n += 1
        if n > 0:
            monitoring.to_log_and_console(
                "    ... " + str(n) + " cells have different contact surfaces", 1
            )
            for k in intersection:
                d = list(set(e1[k].keys()).symmetric_difference(set(e2[k].keys())))
                if len(d) > 0:
                    difference1 = list(set(e1[k].keys()).difference(set(e2[k].keys())))
                    difference2 = list(set(e2[k].keys()).difference(set(e1[k].keys())))
                    s = "cell " + str(k) + " has different contact surfaces: "
                    monitoring.to_log_and_console("        " + s, 1)
                    if len(difference1) > 0:
                        s = (
                            str(difference1)
                            + " in '"
                            + str(name1)
                            + "' and not in '"
                            + str(name2)
                            + "'"
                        )
                        monitoring.to_log_and_console("          " + s, 1)
                    if len(difference2) > 0:
                        s = (
                            str(difference2)
                            + " in '"
                            + str(name2)
                            + "' and not in '"
                            + str(name1)
                            + "'"
                        )
                        monitoring.to_log_and_console("          " + s, 1)
    return


# def _compare_history(d1, d2, name1, name2, description):
#    return


def _compare_principal_vector(d1, d2, name1, name2, description):
    """

    :param d1:
    :param d2:
    :param name1:
    :param name2:
    :param description:
    :return:
    """
    # 'principal-vector': 'cell_principal_vectors'    # liste de numpy.ndarray
    #     dictionary de liste de numpy.ndarray de numpy.float64
    #     cell_principal_vectors.590002 = <type 'list'>
    #     590002: [array([ 0.17420991, -0.74923203,  0.63898534]),
    #         array([-0.24877611,  0.59437038,  0.7647446 ]),
    #         array([ 0.95276511,  0.29219037,  0.08284582])]

    e1 = d1[description]
    e2 = d2[description]

    monitoring.to_log_and_console("  === " + str(description) + " comparison === ", 1)

    intersection = _intersection_cell_keys(e1, e2, name1, name2)

    residual = []

    if len(intersection) > 0:
        for k in intersection:
            residual.append(
                [
                    k,
                    max(
                        math.sqrt(sum((e1[k][0] - e2[k][0]) * (e1[k][0] - e2[k][0]))),
                        math.sqrt(sum((e1[k][1] - e2[k][1]) * (e1[k][1] - e2[k][1]))),
                        math.sqrt(sum((e1[k][2] - e2[k][2]) * (e1[k][2] - e2[k][2]))),
                    ),
                ]
            )

    residual = sorted(residual, key=itemgetter(1), reverse=True)
    monitoring.to_log_and_console(
        "    ... largest residual at cell #"
        + str(residual[0][0])
        + " is "
        + str(residual[0][1]),
        1,
    )
    s = str(e1[residual[0][0]]) + "\n" + "        " + " <-> " + str(e2[residual[0][0]])
    monitoring.to_log_and_console("        " + s, 1)

    return


def comparison(d1, d2, features, name1, name2):
    """

    :param d1:
    :param d2:
    :param features:
    :param name1:
    :param name2:
    :return:
    """

    monitoring.to_log_and_console("\n", 1)
    monitoring.to_log_and_console(
        "... comparison between '" + str(name1) + "' and '" + str(name2) + "'", 1
    )

    #
    # 1. find common keys
    #

    unpairedkeys1 = []
    unpairedkeys2 = []
    pairedkeys = []
    unrecognizedkeys1 = []
    unrecognizedkeys2 = []

    for k1 in d1:
        #
        # loop on known dictionary
        #
        recognizedkey = False
        for k in ioproperties.keydictionary:
            if k1 in ioproperties.keydictionary[k]["input_keys"]:
                recognizedkey = True
                pairedkey = False
                #
                # got it, try to pair it
                #
                for k2 in d2:
                    if k2 in ioproperties.keydictionary[k]["input_keys"]:
                        pairedkey = True
                        pairedkeys.append([k1, k2, k])
                        break
                if pairedkey is False:
                    unpairedkeys1.append(k1)
                break

        if recognizedkey is False:
            unrecognizedkeys1.append(k1)

    #
    #
    #

    for k2 in d2:
        #
        # loop on known dictionary
        #
        recognizedkey = False
        for k in ioproperties.keydictionary:
            if k2 in ioproperties.keydictionary[k]["input_keys"]:
                recognizedkey = True
                pairedkey = False
                #
                # got it, try to pair it
                #
                for k1 in d1:
                    if k1 in ioproperties.keydictionary[k]["input_keys"]:
                        pairedkey = True
                        # pairedkeys.append([k1,k2])
                        break
                if pairedkey is False:
                    unpairedkeys2.append(k2)
                break

        if recognizedkey is False:
            unrecognizedkeys2.append(k2)

    #
    # first output, compare the dictionaries keys
    #

    print_summary = False

    if features is None or len(features) == 0:
        print_summary = True

    #
    #
    #

    if print_summary is True:
        monitoring.to_log_and_console("    found keys", 1)
        if len(pairedkeys) > 0:
            monitoring.to_log_and_console(
                "    ... common keys to '" + str(name1) + "' and '" + str(name2) + "'",
                1,
            )
            for p in pairedkeys:
                monitoring.to_log_and_console(
                    "        " + str(p[0]) + " <-> " + str(p[1]), 1
                )
        if len(unpairedkeys1) > 0:
            monitoring.to_log_and_console(
                "    ... keys in '" + str(name1) + "' and not in '" + str(name2) + "'",
                1,
            )
            for k in unpairedkeys1:
                monitoring.to_log_and_console("        " + str(k), 1)
        if len(unpairedkeys2) > 0:
            monitoring.to_log_and_console(
                "    ... keys not in '" + str(name1) + "' and in '" + str(name2) + "'",
                1,
            )
            for k in unpairedkeys2:
                monitoring.to_log_and_console("        " + str(k), 1)
        if len(unrecognizedkeys1) > 0:
            monitoring.to_log_and_console(
                "    ... keys in '" + str(name1) + "' not recognized", 1
            )
            for k in unrecognizedkeys1:
                monitoring.to_log_and_console("        " + str(k), 1)
        if len(unrecognizedkeys2) > 0:
            monitoring.to_log_and_console(
                "    ... keys in '" + str(name2) + "' not recognized", 1
            )
            for k in unrecognizedkeys2:
                monitoring.to_log_and_console("        " + str(k), 1)

    #
    # 2. perform a comparison key by key
    #

    if len(pairedkeys) == 0:
        monitoring.to_log_and_console(
            "... no common keys between '" + str(name1) + "' and '" + str(name2) + "'",
            1,
        )
        monitoring.to_log_and_console("    comparison is not possible", 1)
        return

    #
    # recall that the dictionary keys are the 'output_key' of the ioproperties.keydictionary
    #

    if features is None or len(features) == 0:
        comparison_keys = [k[2] for k in pairedkeys]
    else:
        comparison_keys = features

    monitoring.to_log_and_console("", 1)

    for f in comparison_keys:
        if f not in ioproperties.keydictionary:
            monitoring.to_log_and_console(
                "    unknown property '" + str(f) + "' for comparison", 1
            )
            continue

        outk = ioproperties.keydictionary[f]["output_key"]

        for i in range(len(pairedkeys)):
            if pairedkeys[i][0] == outk:
                if outk == ioproperties.keydictionary["lineage"]["output_key"]:
                    _compare_lineage(d1, d2, name1, name2, outk)
                elif outk == ioproperties.keydictionary["h_min"]["output_key"]:
                    pass
                    # monitoring.to_log_and_console("    comparison of '" + str(outk) + "' not implemented yet", 1)
                elif outk == ioproperties.keydictionary["volume"]["output_key"]:
                    _compare_volume(d1, d2, name1, name2, outk)
                elif outk == ioproperties.keydictionary["sigma"]["output_key"]:
                    pass
                    # monitoring.to_log_and_console("    comparison of '" + str(outk) + "' not implemented yet", 1)
                elif outk == ioproperties.keydictionary["label_in_time"]["output_key"]:
                    pass
                    # monitoring.to_log_and_console("    comparison of '" + str(outk) + "' not implemented yet", 1)
                elif outk == ioproperties.keydictionary["barycenter"]["output_key"]:
                    _compare_barycenter(d1, d2, name1, name2, outk)
                elif outk == ioproperties.keydictionary["fate"]["output_key"]:
                    _compare_fate(d1, d2, name1, name2, outk)
                elif outk == ioproperties.keydictionary["all-cells"]["output_key"]:
                    _compare_all_cells(d1, d2, name1, name2, outk)
                elif (
                    outk == ioproperties.keydictionary["principal-value"]["output_key"]
                ):
                    _compare_principal_value(d1, d2, name1, name2, outk)
                elif outk == ioproperties.keydictionary["name"]["output_key"]:
                    _compare_name(d1, d2, name1, name2, outk)
                elif outk == ioproperties.keydictionary["contact"]["output_key"]:
                    _compare_contact(d1, d2, name1, name2, outk)
                elif outk == ioproperties.keydictionary["history"]["output_key"]:
                    pass
                    # monitoring.to_log_and_console("    comparison of '" + str(outk) + "' not implemented yet", 1)
                elif (
                    outk == ioproperties.keydictionary["principal-vector"]["output_key"]
                ):
                    _compare_principal_vector(d1, d2, name1, name2, outk)
                else:
                    monitoring.to_log_and_console(
                        "    unknown key '" + str(outk) + "' for comparison", 1
                    )
                break

    return


########################################################################################
#
#
#
########################################################################################


def _find_fate(cell_fate, name):
    for n, v in cell_fate.items():
        if name in v[0]:
            return n
    return None


def _find_fate_with_fullname(cell_fate, name):
    for n, v in cell_fate.items():
        if name[:-1] in v[0]:
            return n
    return None


def _build_dict_indexed_by_name(lastgeneration=9):
    dict_name = {"a1.0001": None, "b1.0001": None}
    #
    # generate names for each generation
    #
    for g in range(1, lastgeneration + 1):
        names = list(dict_name.keys())
        for n in names:
            # pick cell from generation g
            if int(n.split(".")[0][1:]) != g:
                continue
            daughters = uname.get_daughter_names(n + "_")
            for d in daughters:
                dict_name[d[:-1]] = None
    return dict_name


def _build_dict_fate(cell_fate, lastgeneration=9):
    proc = "_build_dict_fate"
    dict_fate = _build_dict_indexed_by_name(lastgeneration=lastgeneration)

    name_by_generation = {}
    for name in dict_fate:
        g = int(name.split(".")[0][1:])
        name_by_generation[g] = name_by_generation.get(g, []) + [name]
    generations = name_by_generation.keys()
    generations = sorted(generations)

    #
    # give fate to names
    #
    for name in dict_fate:
        fate = _find_fate(cell_fate, name)
        if fate is not None:
            dict_fate[name] = fate

    #
    # forward propagation
    # fate[daughter] = fate[mother]
    #
    for g in generations:
        for name in name_by_generation[g]:
            if dict_fate[name] is None:
                continue
            daughters = uname.get_daughter_names(name + "_")
            for d in daughters:
                if d[:-1] not in dict_fate:
                    continue
                if dict_fate[d[:-1]] is not None:
                    continue
                dict_fate[d[:-1]] = dict_fate[name]

    #
    # backward propagation
    # fate[mother] = fate[daughter(s)]
    #
    generations = sorted(generations, reverse=True)
    for g in generations:
        for name in name_by_generation[g]:
            if dict_fate[name] is not None:
                continue
            fate = []
            daughters = uname.get_daughter_names(name + "_")
            for d in daughters:
                if d[:-1] not in dict_fate:
                    continue
                if dict_fate[d[:-1]] is None:
                    continue
                if isinstance(dict_fate[d[:-1]], str):
                    fate += [dict_fate[d[:-1]]]
                elif isinstance(dict_fate[d[:-1]], list):
                    fate += dict_fate[d[:-1]]
                else:
                    msg = (
                        ":type '"
                        + str(type(dict_fate[d[:-1]]))
                        + "' of '"
                        + str(d[:-1])
                        + "' not handled yet"
                    )
                    monitoring.to_log_and_console(str(proc) + msg)
            fate = list(set(fate))
            if len(fate) == 1:
                dict_fate[name] = fate[0]
            elif len(fate) > 1:
                dict_fate[name] = fate

    return dict_fate


def _set_fate_from_names_by_propagation(d, keyfate, cell_fate):
    """
    Obsolete.
    Parameters
    ----------
    d
    keyfate
    cell_fate

    Returns
    -------

    """
    proc = "_set_fate_from_names_by_propagation"
    #
    # give fate to cell that have a name
    #
    d[keyfate] = {}
    for c in d["cell_name"]:
        fate = _find_fate_with_fullname(cell_fate, d["cell_name"][c])
        if fate is not None:
            d[keyfate][c] = fate

    #
    # forward propagation
    # fate[daughter] = fate[mother]
    #
    lineage = d["cell_lineage"]
    cells = list(
        set(lineage.keys()).union(
            set([v for values in list(lineage.values()) for v in values])
        )
    )
    cells = sorted(cells)

    for mother in cells:
        if mother not in d[keyfate]:
            continue
        if mother not in lineage:
            continue
        for c in lineage[mother]:
            if c in d[keyfate]:
                continue
            d[keyfate][c] = d[keyfate][mother]

    #
    # backward propagation
    # fate[mother] = fate[daughter(s)]
    #
    cells = sorted(cells, reverse=True)
    for mother in cells:
        if mother in d[keyfate]:
            continue
        if mother not in lineage:
            continue
        fate = []
        for c in lineage[mother]:
            if c not in d[keyfate]:
                continue
            if isinstance(d[keyfate][c], str):
                fate += [d[keyfate][c]]
            elif isinstance(d[keyfate][c], list):
                fate += d[keyfate][c]
            else:
                msg = ":type '" + str(type(d[keyfate][c])) + "' of d['" + str(keyfate)
                msg += "'][" + str(c) + "]" + "not handled yet"
                monitoring.to_log_and_console(str(proc) + msg)
        fate = list(set(fate))
        if len(fate) == 1:
            d[keyfate][mother] = fate[0]
        elif len(fate) > 1:
            d[keyfate][mother] = fate

    return d


def _set_fate_from_names(d, keyfate, cell_fate):
    #
    # build a dictionary of fates indexed by name
    #
    dict_fate = _build_dict_fate(cell_fate, lastgeneration=9)
    d[keyfate] = {}
    for c in d["cell_name"]:
        if dict_fate[d["cell_name"][c][:-1]] is None:
            continue
        d[keyfate][c] = dict_fate[d["cell_name"][c][:-1]]

    #
    # forward propagation of fates (for cells with no name)
    #
    lineage = d["cell_lineage"]
    cells = list(
        set(lineage.keys()).union(
            set([v for values in list(lineage.values()) for v in values])
        )
    )
    cells = sorted(cells)

    for mother in cells:
        if mother not in d[keyfate]:
            continue
        if mother not in lineage:
            continue
        for c in lineage[mother]:
            if c in d[keyfate]:
                continue
            d[keyfate][c] = d[keyfate][mother]

    return d


def set_fate_from_names(d, fate=4):
    proc = "set_fate_from_names"

    cell_fate2 = {
        "Anterior Endoderm": (["a7.0001", "a7.0002", "a7.0005"], 1),
        "Posterior Endoderm": (["b7.0001", "b7.0002", "b9.0034"], 2),
        "germ line": (["b7.0006"], 3),
        "Mesoderm 1 Notochord": (["a7.0003", "a7.0007"], 4),
        "Mesoderm 2 Notochord": (["b8.0006"], 5),
        "Mesoderm Trunk Lateral Cell": (["a7.0006"], 6),
        "Mesoderm Trunk ventral Cell": (["b7.0005"], 7),
        "Mesoderm First Muscle": (["b7.0004", "b7.0008"], 8),
        "Mesoderm Second Muscle": (["a9.0031", "b9.0033"], 9),
        "Mesoderm Mesenchyme": (["b7.0007", "b8.0005"], 10),
        "Posterior ventral Neural Plate": (["a7.0004"], 11),
        "Anterior + Dorsal Neural Plate": (["a7.0009", "a7.0010", "b8.0019"], 12),
        "Lateral Neural Plate": (["a7.0013", "a8.0015", "a9.0032"], 13),
        "Trunk Epidermis": (
            ["a7.0011", "a7.0012", "a7.0014", "a7.0015", "a7.0016"],
            14,
        ),
        "Midline Tail Epidermis": (
            [
                "b8.0020",
                "b8.0018",
                "b9.0041",
                "b8.0027",
                "b9.0056",
                "b9.0062",
                "b9.0064",
            ],
            15,
        ),
        "Mediolateral Tail Epidermis": (
            [
                "b8.0024",
                "b9.0045",
                "b9.0042",
                "b9.0043",
                "b9.0049",
                "b9.0055",
                "b9.0061",
                "b9.0063",
            ],
            16,
        ),
        "Lateral Tail Epidermis": (
            ["b9.0044", "b8.0026", "b9.0050", "b9.0046", "b8.0029", "b8.0030"],
            17,
        ),
    }

    cell_fate3 = {
        "Head Endoderm": (
            ["a6.0001", "a7.0001", "a7.0002", "a7.0005", "b7.0001", "b8.0003"],
            1,
        ),
        "1st Endodermal Lineage": (["b8.0004"], 2),
        "2nd Endodermal Lineage": (["b9.0034"], 3),
        "1st Lineage, Notochord": (["a7.0003", "a7.0007"], 4),
        "2nd Lineage, Notochord": (["b8.0006"], 5),
        "Trunk Lateral Cell": (["a7.0006"], 6),
        "Trunk Ventral Cell": (["b7.0005"], 7),
        "Mesenchyme": (["b7.0007", "b8.0005"], 8),
        "1st Lineage, Tail Muscle": (["b7.0004", "b7.0008"], 9),
        "2nd Lineage, Tail Muscle": (["a9.0031", "b9.0033"], 10),
        "Anterior Dorsal Neural Plate": (["a7.0013"], 11),
        "Anterior Ventral Neural Plate": (["a6.0005", "a7.0009", "a7.0010"], 12),
        "Posterior Dorsal Neural Plate": (["b8.0019"], 13),
        "Posterior Lateral Neural Plate": (["a8.0015", "a9.0032"], 14),
        "Posterior Ventral Neural Plate": (["a7.0004"], 15),
        "Head Epidermis": (
            [
                "a6.0006",
                "a6.0008",
                "a7.0014",
                "a7.0015",
                "a7.0016",
                "a7.0011",
                "a7.0012",
            ],
            16,
        ),
        "Tail Epidermis": (
            [
                "b7.0011",
                "b7.0012",
                "b7.0013",
                "b7.0014",
                "b7.0015",
                "b7.0016" "b9.0044",
                "b8.0026",
                "b9.0050",
                "b9.0046",
                "b7.0015",
                "b8.0029",
                "b8.0030",
                "b8.0024",
                "b9.0045",
                "b9.0042",
                "b9.0043",
                "b9.0049",
                "b9.0055",
                "b9.0061",
                "b9.0063",
                "b8.0020",
                "b8.0018",
                "b9.0041",
                "b8.0027",
                "b9.0056",
                "b9.0062",
                "b9.0064",
            ],
            17,
        ),
        "Germ Line": (["b7.0006"], 18),
    }

    # new fate (fate4) for visualisation MorphoNet and Tulip tree

    cell_fate4 = {
        "Anterior Head Endoderm": (
            [
                "a6.0001",
                "a7.0001",
                "a7.0002",
                "a7.0005",
                "a8.0001",
                "a8.0002",
                "a8.0003",
                "a8.0004",
                "a8.0009",
                "a8.0010",
            ],
            1,
        ),
        "Posterior Head Endoderm": (["b7.0001", "b8.0003", "b8.0001", "b8.0002"], 2),
        "1st Endodermal Lineage": (["b8.0004"], 3),
        "2nd Endodermal Lineage": (["b9.0034"], 4),
        "1st Lineage, Notochord": (
            ["a7.0003", "a7.0007", "a8.0005", "a8.0006", "a8.0013", "a8.0014"],
            5,
        ),
        "2nd Lineage, Notochord": (["b8.0006"], 6),
        "Trunk Lateral Cell": (["a7.0006", "a8.0011", "a8.0012"], 7),
        "Trunk Ventral Cell": (["b7.0005", "b8.0009", "b8.0010"], 8),
        "Mesenchyme": (["b7.0007", "b8.0005", "b8.0013", "b8.0014"], 9),
        "1st Lineage, Tail Muscle": (
            ["b7.0004", "b7.0008", "b8.0007", "b8.0008", "b8.0015", "b8.0016"],
            10,
        ),
        "2nd Lineage, Tail Muscle": (["a9.0031", "b9.0033"], 11),
        "Anterior Dorsal Neural Plate": (["a7.0013", "a8.0025", "a8.0026"], 12),
        "Anterior Ventral Neural Plate": (["a6.0005", "a7.0009", "a7.0010"], 13),
        "Posterior Dorsal Neural Plate": (["b8.0019"], 14),
        "Posterior Lateral Neural Plate": (["a8.0015", "a9.0032"], 15),
        "Posterior Ventral Neural Plate": (
            [
                "a7.0004",
                "a8.0007",
                "a8.0008",
                "a9.0013",
                "a9.0014",
                "a9.0015",
                "a9.0016",
                "a10.0025",
                "a10.0026",
                "a10.0027",
                "a10.0028",
                "a10.0029",
                "a10.0030",
                "a10.0031",
                "a10.0032",
            ],
            16,
        ),
        "Head Epidermis": (
            [
                "a6.0006",
                "a6.0008",
                "a7.0014",
                "a7.0015",
                "a7.0016",
                "a7.0011",
                "a7.0012",
                "a8.0027",
                "a8.0028",
                "a8.0029",
                "a8.0030",
                "a8.0031",
                "a8.0032",
                "a8.0021",
                "a8.0022",
                "a8.0023",
                "a8.0024",
                "a10.0081",
                "a10.0082",
            ],
            17,
        ),
        "Lateral Tail Epidermis": (
            [
                "b9.0044",
                "b8.0026",
                "b9.0050",
                "b9.0047",
                "b7.0015",
                "b8.0029",
                "b8.0030",
            ],
            18,
        ),
        "Medio-Lateral Tail Epidermis": (
            [
                "b8.0023",
                "b9.0048",
                "b9.0042",
                "b9.0043",
                "b9.0049",
                "b9.0055",
                "b9.0061",
                "b9.0063",
                "b10.0097",
                "b10.0098",
            ],
            19,
        ),
        "Midline Tail Epidermis": (
            [
                "b8.0020",
                "b8.0018",
                "b9.0041",
                "b8.0027",
                "b9.0056",
                "b9.0062",
                "b9.0064",
                "b10.0081",
                "b10.0082",
            ],
            20,
        ),
        "Germ Line": (["b7.0006"], 21),
    }

    cell_fate = {}
    keyfate = None
    if fate == 2:
        cell_fate = cell_fate2
        keyfate = "cell_fate2"
    elif fate == 3:
        cell_fate = cell_fate3
        keyfate = "cell_fate3"
    elif fate == 4:
        cell_fate = cell_fate4
        keyfate = "cell_fate"
    else:
        monitoring.to_log_and_console(
            proc + ": fate index '" + str(fate) + "' not handled"
        )
        return

    # clean properties from previous fates
    for f in ["cell_fate", "cell_fate2", "cell_fat3"]:
        if f in d:
            del d[f]

    # d = _set_fate_from_names_by_propagation(d, keyfate, cell_fate)
    d = _set_fate_from_names(d, keyfate, cell_fate)

    return d


def _set_color_from_fate(d, colormap_version=2020):
    proc = "_set_color_from_fate"

    color_fate_2020 = {
        "1st Lineage, Notochord": 2,
        "Posterior Ventral Neural Plate": 19,
        "Anterior Ventral Neural Plate": 9,
        "Anterior Head Endoderm": 8,
        "Anterior Endoderm": 8,
        "Posterior Head Endoderm": 17,
        "Posterior Endoderm": 17,
        "Trunk Lateral Cell": 20,
        "Mesenchyme": 14,
        "1st Lineage, Tail Muscle": 3,
        "Trunk Ventral Cell": 21,
        "Germ Line": 10,
        "Lateral Tail Epidermis": 12,
        "Head Epidermis": 11,
        "Trunk Epidermis": 11,
        "Anterior Dorsal Neural Plate": 7,
        "Posterior Lateral Neural Plate": 18,
        "2nd Lineage, Notochord": 5,
        "Medio-Lateral Tail Epidermis": 13,
        "Midline Tail Epidermis": 15,
        "Posterior Dorsal Neural Plate": 16,
        "1st Endodermal Lineage": 1,
        "2nd Lineage, Tail Muscle": 6,
        "2nd Endodermal Lineage": 4,
    }

    color_fate_2009 = {
        "1st Lineage, Notochord": 78,
        "Posterior Ventral Neural Plate": 58,
        "Anterior Ventral Neural Plate": 123,
        "Anterior Head Endoderm": 1,
        "Anterior Endoderm": 1,
        "Posterior Head Endoderm": 27,
        "Posterior Endoderm": 27,
        "Trunk Lateral Cell": 62,
        "Mesenchyme": 63,
        "1st Lineage, Tail Muscle": 135,
        "Trunk Ventral Cell": 72,
        "Germ Line": 99,
        "Lateral Tail Epidermis": 61,
        "Head Epidermis": 76,
        "Trunk Epidermis": 76,
        "Anterior Dorsal Neural Plate": 81,
        "Posterior Lateral Neural Plate": 75,
        "2nd Lineage, Notochord": 199,
        "Medio-Lateral Tail Epidermis": 41,
        "Midline Tail Epidermis": 86,
        "Posterior Dorsal Neural Plate": 241,
        "1st Endodermal Lineage": 40,
        "2nd Lineage, Tail Muscle": 110,
        "2nd Endodermal Lineage": 44,
    }

    if colormap_version == 2020:
        colormap = color_fate_2020
        keycolormap = "morphonet_selection_tissuefate_guignard_2020"
    elif colormap_version == 2009:
        colormap = color_fate_2009
        keycolormap = "morphonet_selection_tissuefate_lemaire_2009"
    else:
        monitoring.to_log_and_console(
            proc + ": colormap version '" + str(colormap_version) + "' not handled"
        )
        return

    if keycolormap in d:
        del d[keycolormap]
    d[keycolormap] = {}

    for c in d["cell_fate"]:
        if isinstance(d["cell_fate"][c], str):
            if d["cell_fate"][c] not in colormap:
                monitoring.to_log_and_console(
                    proc
                    + ": fate '"
                    + str(d["cell_fate"][c])
                    + "' not handled in color map"
                )
                continue
            d[keycolormap][c] = [colormap[d["cell_fate"][c]]]
        elif isinstance(d["cell_fate"][c], list):
            for f in d["cell_fate"][c]:
                if not isinstance(f, str):
                    msg = (
                        ":type '" + str(f) + "' found in d['cell_fate'][" + str(c) + "]"
                    )
                    msg += "not handled yet"
                    monitoring.to_log_and_console(str(proc) + msg)
                    continue
                if f not in colormap:
                    monitoring.to_log_and_console(
                        proc + ": fate '" + str(f) + "' not handled in color map"
                    )
                    continue
                if c not in d[keycolormap]:
                    d[keycolormap][c] = [colormap[f]]
                else:
                    d[keycolormap][c].append(colormap[f])
        else:
            msg = (
                ":type '"
                + str(d["cell_fate"][c])
                + "' of d['cell_fate']["
                + str(c)
                + "]"
            )
            msg += "not handled yet"
            monitoring.to_log_and_console(str(proc) + msg)

    return d


def set_color_from_fate(d):
    d = _set_color_from_fate(d, colormap_version=2020)
    d = _set_color_from_fate(d, colormap_version=2009)
    return d


########################################################################################
#
#
#
########################################################################################


def _find_t(cells_per_time, n):
    if n in cells_per_time.values():
        times = [t for t in cells_per_time if cells_per_time[t] == n]
        return (min(times) + max(times)) / 2.0
    smaller_times = [t for t in cells_per_time if cells_per_time[t] < n]
    larger_times = [t for t in cells_per_time if cells_per_time[t] > n]
    ts = max(smaller_times)
    tl = min(larger_times)
    return ts + (tl - ts) * (n - cells_per_time[ts]) / (
        cells_per_time[tl] - cells_per_time[ts]
    )


def _temporal_alignement(ref_cells_per_time, cells_per_time):
    first = max(min(ref_cells_per_time.values()), min(cells_per_time.values()))
    last = min(max(ref_cells_per_time.values()), max(cells_per_time.values()))
    ref_times = []
    times = []
    for n in range(first, last + 1):
        if n not in ref_cells_per_time.values() and n not in cells_per_time.values():
            continue
        ref_times += [_find_t(ref_cells_per_time, n)]
        times += [_find_t(cells_per_time, n)]
    # Consider $\sum_i \left( a t_i +b - t'_i\right)^2$
    #
    # \begin{eqnarray*}
    # \frac{\partial}{\partial b} \sum_i \left( a t_i +b - t'_i\right)^2 = 0
    # & \Leftrightarrow &
    # \sum_i 2 \left( a t_i +b - t'_i\right) = 0 \\
    # & \Leftrightarrow & a \sum_i t_i + N b - \sum_i t'_i = 0 \\
    # & \Leftrightarrow & b = \frac{\sum_i t'_i - a \sum_i t_i}{N}
    # \end{eqnarray*}
    #
    # \begin{eqnarray*}
    # \frac{\partial}{\partial a} \sum_i \left( a t_i +b - t'_i\right)^2 = 0
    # & \Leftrightarrow &
    # \sum_i 2 \left( a t_i +b - t'_i\right) t_i = 0 \\
    # & \Leftrightarrow &
    # a \sum_i (t_i)^2 + b \sum_i t_i - \sum_i t_i t'_i = 0 \\
    # & \Leftrightarrow &
    # a \sum_i (t_i)^2 +
    # \frac{\left( \sum_i t_i \right) \left( \sum_i t'_i \right)}{N}
    # - a \frac{\left( \sum_i t_i \right)^2}{N} - \sum_i t_i t'_i = 0 \\
    # & \Leftrightarrow &
    # a \left( \sum_i (t_i)^2 - \frac{\left( \sum_i t_i \right)^2}{N} \right)
    # = \sum_i t_i t'_i - \frac{\left( \sum_i t_i \right) \left( \sum_i t'_i \right)}{N}
    # \end{eqnarray*}
    num = sum(np.multiply(times, ref_times)) - sum(times) * sum(ref_times) / len(times)
    den = sum(np.multiply(times, times)) - sum(times) * sum(times) / len(times)
    a = num / den
    b = (sum(ref_times) - a * sum(times)) / len(times)
    return a, b


def temporal_alignment(
    ref_lineage, lineage, ref_time_digits_for_cell_id=4, time_digits_for_cell_id=4
):
    """

    Parameters
    ----------
    ref_lineage: reference lineage
    lineage: lineage to be aligned with the reference lineage
    time_digits_for_cell_id: number of digits for the cell number in the unique cell identifier

    Returns
    -------
    a, b: coefficients of the linear warping. The number of cells at time 't' of the lineage is
        made comparable to the number of cells at time 'a * t + b' of the reference lineage.

    """
    ref_div = 10**ref_time_digits_for_cell_id
    div = 10**time_digits_for_cell_id

    cells = list(
        set(ref_lineage.keys()).union(
            set([v for values in list(ref_lineage.values()) for v in values])
        )
    )
    cells = sorted(cells)
    ref_cells_per_time = {}
    for c in cells:
        t = int(c) // ref_div
        ref_cells_per_time[t] = ref_cells_per_time.get(t, 0) + 1

    cells = list(
        set(lineage.keys()).union(
            set([v for values in list(lineage.values()) for v in values])
        )
    )
    cells = sorted(cells)
    cells_per_time = {}
    for c in cells:
        t = int(c) // div
        cells_per_time[t] = cells_per_time.get(t, 0) + 1

    return _temporal_alignement(ref_cells_per_time, cells_per_time)


########################################################################################
#
# utilities for debugging, etc.
#
########################################################################################


def print_keys(d, desc=None):
    monitoring.to_log_and_console("\n", 1)
    if desc is None:
        monitoring.to_log_and_console("... contents", 1)
    else:
        monitoring.to_log_and_console("... contents of '" + str(desc) + "'", 1)

    if type(d) is dict:
        if d == {}:
            monitoring.to_log_and_console("    " + "empty dictionary", 1)
        else:
            monitoring.to_log_and_console("    " + "keys are:", 1)
            for k in d:
                monitoring.to_log_and_console("    " + "- " + str(k), 1)
    else:
        monitoring.to_log_and_console("    " + "input is not a dictionary", 1)

    return


def print_type(d, t=None, desc=None):
    if desc is None:
        desc = ""
    if t is None:
        t = ""

    if type(d) is dict:
        print("type of " + desc + " is " + str(t) + str(type(d)))
        for k in d:
            print_type(d[k], t + str(type(d)) + ".", desc + "." + str(k))

    elif type(d) in (list, np.array, np.ndarray):
        print("type of " + desc + " is " + str(t) + str(type(d)))
        print_type(d[0], t + str(type(d)) + ".", desc + "[0]")

    else:
        print("type of " + desc + " is " + str(t) + str(type(d)))

    return
