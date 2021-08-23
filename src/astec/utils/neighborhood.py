import os
import sys
import copy
import collections
import math

import astec.utils.common as common
import astec.algorithms.properties as properties


#
# the atlas here is defined as the collection of (daughter) cell neighborhoods
# right after the division
#

monitoring = common.Monitoring()


class NeighborhoodParameters(common.PrefixedParameter):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self):

        common.PrefixedParameter.__init__(self)

        if "doc" not in self.__dict__:
            self.doc = {}

        self.outputAtlasFile = None
        self.atlasFiles = []

        #
        # how to build an atlas
        #
        # add the symmetric neighborhood as neighborhood
        # consider the other half of the embryo as a single cell
        #
        self.add_symmetric_neighborhood = True
        self.differentiate_other_half = True

        #
        #
        #
        self.neighborhood_diagnosis = False

        #
        #
        #
        self.neighborhood_improvement = False

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('#')
        print('# CellNeighborhoodParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        self.varprint('outputAtlasFile', self.outputAtlasFile)
        self.varprint('atlasFiles', self.atlasFiles)

        self.varprint('add_symmetric_neighborhood', self.add_symmetric_neighborhood)
        self.varprint('differentiate_other_half', self.differentiate_other_half)

        self.varprint('neighborhood_diagnosis', self.neighborhood_diagnosis)
        self.varprint('neighborhood_improvement', self.neighborhood_improvement)

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# CellNeighborhoodParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'outputAtlasFile', self.outputAtlasFile)
        self.varwrite(logfile, 'atlasFiles', self.atlasFiles)

        self.varwrite(logfile, 'add_symmetric_neighborhood', self.add_symmetric_neighborhood)
        self.varwrite(logfile, 'differentiate_other_half', self.differentiate_other_half)

        self.varwrite(logfile, 'neighborhood_diagnosis', self.neighborhood_diagnosis)
        self.varwrite(logfile, 'neighborhood_improvement', self.neighborhood_improvement)

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
        self.outputAtlasFile = self.read_parameter(parameters, 'outputAtlasFile', self.outputAtlasFile)
        self.atlasFiles = self.read_parameter(parameters, 'atlasFiles', self.atlasFiles)
        self.atlasFiles = self.read_parameter(parameters, 'referenceFiles', self.atlasFiles)

        self.add_symmetric_neighborhood = self.read_parameter(parameters, 'add_symmetric_neighborhood',
                                                              self.add_symmetric_neighborhood)
        self.differentiate_other_half = self.read_parameter(parameters, 'differentiate_other_half',
                                                            self.differentiate_other_half)

        self.neighborhood_diagnosis = self.read_parameter(parameters, 'neighborhood_diagnosis',
                                                          self.neighborhood_diagnosis)
        self.neighborhood_improvement = self.read_parameter(parameters, 'neighborhood_improvement',
                                                            self.neighborhood_improvement)

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

def get_daughter_names(name):
    #
    # build daughter names from parent name
    #
    # anterior or posterior character 'a' or 'b'
    # stage (round of division)
    # '.'
    # p value (cell index)
    # left or right character '_' or '*'
    #
    abvalue = name.split('.')[0][0]
    stage = name.split('.')[0][1:]
    p = name.split('.')[1][0:4]
    lrvalue = name.split('.')[1][4]
    #
    # build daughter names
    #
    daughters = [abvalue + str(int(stage) + 1) + "." + '{:0{width}d}'.format(2 * int(p) - 1, width=4) + lrvalue,
                 abvalue + str(int(stage) + 1) + "." + '{:0{width}d}'.format(2 * int(p), width=4) + lrvalue]
    # print("name = " + str(name) + " -> daughter names = " + str(daughters))
    return daughters


def get_mother_name(name):
    #
    # build daughter names from parent name
    #
    # anterior or posterior character 'a' or 'b'
    # stage (round of division)
    # '.'
    # p value (cell index)
    # left or right character '_' or '*'
    #
    abvalue = name.split('.')[0][0]
    stage = name.split('.')[0][1:]
    p = name.split('.')[1][0:4]
    lrvalue = name.split('.')[1][4]
    #
    # build parent names
    #
    parent = abvalue + str(int(stage)-1) + "."
    if int(p) % 2 == 1:
        parent += '{:0{width}d}'.format((int(p)+1) // 2, width=4)
    else:
        parent += '{:0{width}d}'.format(int(p) // 2, width=4)
    parent += lrvalue
    # print("name = " + str(name) + " -> parent name = " + str(parent))
    return parent


def get_sister_name(name):
    sister_names = get_daughter_names(get_mother_name(name))
    sister_names.remove(name)
    return sister_names[0]


def get_symmetric_name(name):
    symname = name[:-1]
    if name[-1] == '*':
        symname += '_'
    elif name[-1] == '_':
        symname += '*'
    else:
        return None
    return symname


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
            sn = get_symmetric_name(n)
        symneighborhood[sn] = neighborhood[n]
    return symneighborhood


#######################################################################################
#
# compute scores as a scalar product
#
########################################################################################

def _print_common_neighborhoods(neighborhood0, neighborhood1, title=None):
    #
    # used by get_score(), to display neighborhoods after being put in a common frame
    #
    msg = ""
    if title is not None and isinstance(title, str):
        msg += title + " = "
    msg += "{"
    key_list = sorted(list(set(neighborhood0.keys()).union(set(neighborhood1.keys()))))
    for k in key_list:
        msg += str(k) + ": "
        if k in neighborhood0.keys():
            msg += str(neighborhood0[k])
        else:
            msg += "NULL"
        msg += " <-> "
        if k in neighborhood1.keys():
            msg += str(neighborhood1[k])
        else:
            msg += "NULL"
        if k != key_list[-1]:
            msg += ",\n\t "
        else:
            msg += "}"
    monitoring.to_log_and_console(msg)


#
# build a common representation of two neighborhood
#
def _build_common_contact_surfaces(tmp0, tmp1, vect0, vect1):
    if tmp1 == {}:
        return tmp0, tmp1, vect0, vect1
    names1 = list(tmp1.keys())
    for n1 in names1:
        # the name could have been deleted at case 2
        if n1 not in tmp1:
            continue
        # case 1: the same name in both neighborhood
        if n1 in tmp0:
            vect0[n1] = tmp0[n1]
            vect1[n1] = tmp1[n1]
            del tmp0[n1]
            del tmp1[n1]
            continue
        # subcase 1: it's the background (or any generic name) and it is not in tmp0
        if n1 == 'background' or n1 == 'other-half':
            vect0[n1] = 0.0
            vect1[n1] = tmp1[n1]
            del tmp1[n1]
            continue
        # case 2: daughter name in tmp1, and mother name in tmp0
        # add daughter surfaces
        pname = get_mother_name(n1)
        if pname in tmp0:
            vect0[pname] = tmp0[pname]
            vect1[pname] = 0.0
            for n in get_daughter_names(pname):
                if n in tmp1:
                    vect1[pname] += tmp1[n]
                    del tmp1[n]
            del tmp0[pname]
            continue
        # case 3: mother name in tmp1, and daughter name(s) in tmp0
        dname = get_daughter_names(n1)
        if dname[0] in tmp0 or dname[1] in tmp0:
            vect0[n1] = 0.0
            vect1[n1] = tmp1[n1]
            for n in dname:
                if n in tmp0:
                    vect0[n1] += tmp0[n]
                    del tmp0[n]
            del tmp1[n1]
            continue
        # case 4: no association in tmp0
        vect0[n1] = 0.0
        vect1[n1] = tmp1[n1]
        del tmp1[n1]
        continue
    return tmp0, tmp1, vect0, vect1


#
#
#
def _scalar_product(vect0, vect1):
    n0 = 0.0
    n1 = 0.0
    ps = 0.0
    for k in vect0:
        ps += vect0[k] * vect1[k]
        n0 += vect0[k] * vect0[k]
        n1 += vect1[k] * vect1[k]
    return ps / (math.sqrt(n0 * n1))


def _l1_normalized_modulus(vect0, vect1):
    n0 = 0.0
    n1 = 0.0
    nm = 0.0
    for k in vect0:
        n0 += vect0[k]
    for k in vect1:
        n1 += vect1[k]
    for k in vect0:
        nm += abs(vect0[k]/n0 - vect1[k]/n1)
    return nm


#
#
#
def get_score(neigh0, neigh1, title=None):
    """
    Compute the similarity score of two neighborhoods, as the normalized scalar
    product of the vectors of contact surfaces.
    Parameters
    ----------
    neigh0: dictionary depicting a cell neighborhood. Each key is a neighbor, and the
            associated dictionary value give the contact surface
    neigh1: dictionary depicting a cell neighborhood
    title: if not None, print out the neighborhoods as well as the score. For debug purpose.

    Returns
    -------
    the similarity score.

    """
    tmp0 = copy.deepcopy(neigh0)
    tmp1 = copy.deepcopy(neigh1)
    # 1. parse names of tmp1
    tmp0, tmp1, vect0, vect1 = _build_common_contact_surfaces(tmp0, tmp1, {}, {})
    # 2. parse remaining names of tmp0
    tmp1, tmp0, vect1, vect0 = _build_common_contact_surfaces(tmp1, tmp0, vect1, vect0)
    #
    # scalar product in [0, 1] (0: full diagreement, 1: perfect agreement)
    # _l1_normalized_modulus in [0, 2] (0: perfect agreement, 2:full diagreement)
    #
    # score = _scalar_product(vect0, vect1)
    # score = 1.0 - _l1_normalized_modulus(vect0, vect1) / 2.0
    score = _scalar_product(vect0, vect1)
    if title is not None:
        _print_common_neighborhoods(vect0, vect1, title=title)
        monitoring.to_log_and_console("\t score = " + str(score) + "\n")
    # compute score as a scalar product
    return score


########################################################################################
#
# diagnosis on naming
# is redundant with the diagnosis made in properties.py
# (except that contact surfaces are assessed too)
#
########################################################################################

def _neighborhood_diagnosis(prop, time_digits_for_cell_id=4, verbose=1):
    """
    Diagnosis on names extracted from embryo properties
    Parameters
    ----------
    prop
    time_digits_for_cell_id
    verbose

    Returns
    -------

    """
    proc = "_neighborhood_diagnosis"
    if 'cell_lineage' not in prop:
        monitoring.to_log_and_console(str(proc) + ": 'cell_lineage' was not in dictionary")
        return None

    if 'cell_contact_surface' not in prop:
        monitoring.to_log_and_console(str(proc) + ": 'cell_contact_surface' was not in dictionary")
        return None

    if 'cell_name' not in prop:
        monitoring.to_log_and_console(str(proc) + ": 'cell_name' was not in dictionary")
        return None

    lineage = prop['cell_lineage']
    name = prop['cell_name']
    contact = prop['cell_contact_surface']

    reverse_lineage = {v: k for k, values in lineage.items() for v in values}

    div = 10 ** time_digits_for_cell_id

    cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
    cells = sorted(cells)

    cells_per_time = {}
    names_per_time = {}
    missing_name = {}
    missing_contact = {}
    error_name = {}
    for c in cells:
        t = int(c) // div
        n = int(c) % div
        #
        # get cells and cell names at each time point
        #
        cells_per_time[t] = cells_per_time.get(t, []) + [n]
        if c in name:
            names_per_time[t] = names_per_time.get(t, []) + [name[c]]

        #
        # check names
        #
        if c not in name:
            missing_name[t] = missing_name.get(t, []) + [n]
        elif c in reverse_lineage:
            mother = reverse_lineage[c]
            if mother not in name:
                if verbose >= 1:
                    msg = ": weird, cell " + str(c) + " has a name = " + str(name[c])
                    msg += ", but its mother cell " + str(mother) + " has no name"
                    monitoring.to_log_and_console(str(proc) + msg)
                error_name[t] = error_name.get(t, []) + [n]
            else:
                if len(lineage[mother]) == 1:
                    if name[mother] != name[c]:
                        if verbose >= 1:
                            msg = ": weird, cell " + str(c) + " has a name = " + str(name[c])
                            msg += " different than its mother cell " + str(mother) + " name = " + str(name[mother])
                            monitoring.to_log_and_console(str(proc) + msg)
                        error_name[t] = error_name.get(t, []) + [n]
                elif len(lineage[mother]) == 2:
                    daughter_names = get_daughter_names(name[mother])
                    if name[c] not in daughter_names:
                        if verbose >= 1:
                            msg = ": weird, name of cell " + str(c) + " is " + str(name[c])
                            msg += " but should be in " + str(daughter_names)
                            msg += " since its mother cell " + str(mother) + " is named " + str(name[mother])
                            monitoring.to_log_and_console(str(proc) + msg)
                        error_name[t] = error_name.get(t, []) + [n]
                    else:
                        siblings = copy.deepcopy(lineage[mother])
                        siblings.remove(c)
                        daughter_names.remove(name[c])
                        if siblings[0] not in name:
                            if verbose >= 1:
                                msg = ": weird, cell " + str(c) + " has no name "
                                msg += ", it should be " + str(daughter_names[0])
                                msg += " since its mother cell " + str(mother) + " is named " + str(name[mother])
                                msg += " and its sibling " + str(c) + " is named " + str(name[c])
                                monitoring.to_log_and_console(str(proc) + msg)
                            error_name[t] = error_name.get(t, []) + [n]
                        elif name[siblings[0]] == name[c]:
                            if verbose >= 1:
                                msg = ": weird, name of cell " + str(c) + ", " + str(name[c])
                                msg += ", is the same than its sibling " + str(siblings[0])
                                msg += ", their mother cell " + str(mother) + " is named " + str(name[mother])
                                monitoring.to_log_and_console(str(proc) + msg)
                            error_name[t] = error_name.get(t, []) + [n]
                        elif name[siblings[0]] != daughter_names[0]:
                            if verbose >= 1:
                                msg = ": weird, name of cell " + str(siblings[0]) + " is " + str(name[siblings[0]])
                                msg += " but should be " + str(daughter_names[0])
                                msg += " since its mother cell " + str(mother) + " is named " + str(name[mother])
                                msg += " and its sibling " + str(c) + " is named " + str(name[c])
                                monitoring.to_log_and_console(str(proc) + msg)
                            error_name[t] = error_name.get(t, []) + [n]
                else:
                    if verbose >= 1:
                        monitoring.to_log_and_console(str(proc) + ": weird, cell " + str(mother) + " has " +
                                                      str(len(lineage[mother])) + " daughter cells")

        #
        # check contact surfaces
        #
        if c not in contact:
            missing_contact[t] = missing_contact.get(t, []) + [n]

    #
    # interval without errors
    #
    first_time = min(cells_per_time.keys())
    last_time = max(cells_per_time.keys())
    if missing_name != {}:
        last_time = min(last_time, min(missing_name.keys())-1)
    if missing_contact != {}:
        last_time = min(last_time, min(missing_contact.keys())-1)
    if error_name != {}:
        last_time = min(last_time, min(error_name.keys())-1)

    #
    # report
    #
    if verbose >= 1:
        monitoring.to_log_and_console(str(proc) + ": details")
        monitoring.to_log_and_console("\t - first time in lineage = " + str(min(cells_per_time.keys())))
        monitoring.to_log_and_console("\t   last time in lineage = " + str(max(cells_per_time.keys())))
        monitoring.to_log_and_console("\t - interval without errors = [" + str(first_time) + ", " + str(last_time) +
                                      "]")
        monitoring.to_log_and_console("\t - cells in lineage = " + str(len(cells)))
        monitoring.to_log_and_console("\t   #cells at first time = " +
                                      str(len(cells_per_time[min(cells_per_time.keys())])))
        monitoring.to_log_and_console("\t   #cells at last time = " +
                                      str(len(cells_per_time[max(cells_per_time.keys())])))
        # compte le nombre de noms
        monitoring.to_log_and_console("\t - names in lineage = " +
                                      str(len(list(collections.Counter(list(name.keys())).keys()))))

        for t in names_per_time:
            repeats = {k: names_per_time[t].count(k) for k in set(names_per_time[t]) if names_per_time[t].count(k) > 1}
            if repeats != {}:
                monitoring.to_log_and_console("\t - there are " + str(len(repeats)) + " repeated names at time " +
                                              str(t))
                for n, p in repeats.items():
                    monitoring.to_log_and_console("\t   " + str(n) + " is repeated " + str(p) + " times ")

        if error_name != {}:
            monitoring.to_log_and_console("\t - first time with cells with name inconsistency = " +
                                          str(min(error_name.keys())))
            if verbose >= 2:
                for t in sorted(error_name.keys()):
                    monitoring.to_log_and_console(
                        "\t   cells with name inconsistency at time " + str(t) + " = " + str(error_name[t]))

        if missing_name != {}:
            monitoring.to_log_and_console("\t - first time with cells without name = " + str(min(missing_name.keys())))
            if verbose >= 2:
                for t in sorted(missing_name.keys()):
                    monitoring.to_log_and_console("\t   cells without name at time " + str(t) + " = " +
                                                  str(missing_name[t]))

        if missing_contact != {}:
            monitoring.to_log_and_console("\t - first time with cells without contact surfaces = "
                                          + str(min(missing_contact.keys())))
            if verbose >= 2:
                for t in sorted(missing_contact.keys()):
                    monitoring.to_log_and_console("\t   cells without contact surfaces at time " + str(t) + " = " +
                                                  str(missing_contact[t]))
        monitoring.to_log_and_console("")

    return [first_time, last_time]


def get_neighborhood_consistency(neighborhoods):
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
        mother_name = get_mother_name(cell_name)
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
        sister_name = get_sister_name(cell_name)
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
                title = None
                same_choice = get_score(neighborhoods[cell_name][ref1], neighborhoods[cell_name][ref2], title=title)
                sister_choice = get_score(neighborhoods[cell_name][ref1], neighborhoods[sister_name][ref2], title=title)

                mother_name = get_mother_name(cell_name)

                if mother_name not in tested_couples:
                    tested_couples[mother_name] = 1
                else:
                    tested_couples[mother_name] += 1
                if same_choice > sister_choice:
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
        msg += str(get_daughter_names(mother_name)) + " has " + str(len(discrepancy[mother_name]))
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


def _check_neighborhood_consistency(neighborhoods):
    proc = "_check_neighborhood_consistency"

    references, tested_couples, discrepancy = get_neighborhood_consistency(neighborhoods)

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


def _check_leftright_consistency(neighborhoods):
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
        mother_name = get_mother_name(cell_name)
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
        symmother = get_symmetric_name(mother)
        processed_mothers.append(mother)

        if symmother not in references:
            continue
        daughters = get_daughter_names(mother)

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
                symdaughter = get_symmetric_name(daughter)
                symsister = get_sister_name(symdaughter)
                if symdaughter not in neighborhoods or symsister not in neighborhoods:
                    continue
                if reference not in neighborhoods[symdaughter] or reference not in neighborhoods[symsister]:
                    continue
                #
                #
                #
                symsameneigh = get_symmetric_neighborhood(neighborhoods[symdaughter][reference])
                symsisterneigh = get_symmetric_neighborhood(neighborhoods[symsister][reference])

                same_choice = get_score(neighborhoods[daughter][reference], symsameneigh)
                sister_choice = get_score(neighborhoods[daughter][reference], symsisterneigh)

                if same_choice > sister_choice:
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
            symmother = get_symmetric_name(mother)
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
                symmother = get_symmetric_name(mother)
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


def neighborhood_consistency_diagnosis(neighborhoods):
    _check_leftright_consistency(neighborhoods)
    _check_neighborhood_consistency(neighborhoods)


#######################################################################################
#
# consistency is a measure between two reference embryos exhibiting the same division
# there is an inconsistency when the score computed between the same daughter cells
# of two reference embryos is lower than the score of this daughter cell (of one reference)
# and its sister cell (of the other reference). Thus, the consistency is a measure in
# [0, 4]: 0 means a total consistency while 4 designs a total inconsistency.
#
########################################################################################

def _si_global_score(neighbors, references, debug=False):
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
            if debug:
                print("   - score[" + str(r1) + ", " + str(r2) + "] = " + str(get_score(neighbors[0][r1],
                                                                                        neighbors[0][r2])))
                print("   - score[" + str(r1) + ", " + str(r2) + "] = " + str(get_score(neighbors[1][r1],
                                                                                        neighbors[1][r2])))
            score += get_score(neighbors[0][r1], neighbors[0][r2])
            score += get_score(neighbors[1][r1], neighbors[1][r2])
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


def _si_test_one_division(neighborhoods, mother):
    """
    Test whether any daughter switch (for a given reference) improve the global score
    Parameters
    ----------
    neighborhoods
    mother

    Returns
    -------

    """
    daughters = get_daughter_names(mother)
    neighbors = {0: copy.deepcopy(neighborhoods[daughters[0]]), 1: copy.deepcopy(neighborhoods[daughters[1]])}

    references = set.intersection(set(neighbors[0].keys()), set(neighbors[1].keys()))
    if len(references) <= 1:
        return {}

    # score before any changes
    score = _si_global_score(neighbors, references)

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
            newscore[r] = _si_global_score(tmp, references)
            if newscore[r] < score:
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
            ref = max(newscore, key=lambda key: newscore[key])
        corrections[i] = (ref, newscore[ref] - score)
        i += 1
        # if one correction has been found, apply it
        # and look for an other additional correction
        tmp = copy.deepcopy(neighbors)
        tmp = _si_switch_contact_surfaces(tmp, ref, daughters)
        neighbors[0] = copy.deepcopy(tmp[0])
        neighbors[1] = copy.deepcopy(tmp[1])
        score = newscore[ref]


def global_score_improvement(neighborhoods):
    """
    Check whether the global score can be improved by switching the two daughter cells of a
    mother cell in a given reference. The global score score (for a mother cell) is the sum
    (over all possible reference embryo couple) of the local scores.
    Parameters
    ----------
    neighborhoods: dictionary of dictionaries. First level of keys is the cell name (name of a daughter), and the
                   second level of keys is the reference embryo name.
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
        mother = get_mother_name(d)
        if mother in processed_mothers:
            continue
        sister = get_sister_name(d)
        if sister not in neighborhoods:
            msg = "weird, cell " + str(d) + " is in the reference neighborhoods, while its sister "
            msg += str(sister) + " is not "
            monitoring.to_log_and_console(str(proc) + ": " + msg)
            continue
        correction = _si_test_one_division(neighborhoods, mother)
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

    if not isinstance(parameters, NeighborhoodParameters):
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

    for d in daughters:
        if reverse_lineage[d] not in name:
            continue
        #
        # check whether the cell is in dictionaries
        #
        if d not in name:
            if d not in missing_name:
                missing_name.append(d)
                monitoring.to_log_and_console(str(proc) + ": daughter cell #" + str(d)
                                              + " was not found in 'cell_name' dictionary. Skip it", 6)
            continue

        if d not in contact:
            if d not in missing_contact:
                missing_contact.append(d)
                monitoring.to_log_and_console(str(proc) + ": cell #" + str(d)
                                              + " was not found in 'cell_contact_surface' dictionary. Skip it")
            continue

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
                sname = get_symmetric_name(prop['cell_name'][d])
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

def build_neighborhoods(atlasfiles, parameters, time_digits_for_cell_id=4):
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

    if not isinstance(parameters, NeighborhoodParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    neighborhoods = {}
    if isinstance(atlasfiles, str):
        prop = properties.read_dictionary(atlasfiles, inputpropertiesdict={})
        if parameters.neighborhood_diagnosis:
            _neighborhood_diagnosis(prop, time_digits_for_cell_id=time_digits_for_cell_id)
        name = atlasfiles.split(os.path.sep)[-1]
        if name.endswith(".xml") or name.endswith(".pkl"):
            name = name[:-4]
        neighborhoods = _add_neighborhoods(neighborhoods, prop, parameters, atlas_name=name,
                                           time_digits_for_cell_id=time_digits_for_cell_id)
        del prop
    elif isinstance(atlasfiles, list):
        for f in atlasfiles:
            prop = properties.read_dictionary(f, inputpropertiesdict={})
            if parameters.neighborhood_diagnosis:
                _neighborhood_diagnosis(prop, time_digits_for_cell_id=time_digits_for_cell_id)
            name = f.split(os.path.sep)[-1]
            if name.endswith(".xml") or name.endswith(".pkl"):
                name = name[:-4]
            neighborhoods = _add_neighborhoods(neighborhoods, prop, parameters, atlas_name=name,
                                               time_digits_for_cell_id=time_digits_for_cell_id)
            del prop

    # consistency is a measure between two reference embryos exhibiting the same division
    # there is an inconsistency when the score computed between the same daughter cells
    # of two reference embryos is lower than the score of this daughter cell (of one reference)
    # and its sister cell (of the other reference). Thus, the consistency is a measure in
    # [0, 4]: 0 means a total consistency while 4 designs a total inconsistency.

    if parameters.neighborhood_diagnosis:
        neighborhood_consistency_diagnosis(neighborhoods)

    return neighborhoods


########################################################################################
#
#
#
########################################################################################



