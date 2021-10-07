

from operator import itemgetter
import collections
import copy
import statistics
import os
import sys

from astec.utils import common
import astec.utils.ascidian_name as uname
import astec.utils.contact as ucontact
import astec.utils.properties as properties


#
#
#
#
#

monitoring = common.Monitoring()

########################################################################################
#
#
#
########################################################################################


class DiagnosisParameters(ucontact.ContactSurfaceParameters):

    def __init__(self, prefix='diagnosis_'):

        if "doc" not in self.__dict__:
            self.doc = {}

        ucontact.ContactSurfaceParameters.__init__(self, prefix=prefix)

        #
        # parametres pour le diagnostic sur le volume
        #
        self.minimal_volume = 0
        self.maximal_volume_variation = 40.0
        self.maximal_volume_derivative = 30.0

        #
        # nombre d'items a imprimer
        #
        self.items = 10
        
        #
        # diagnostic pour le lignage
        #
        self.minimal_length = 10
        
        #
        # parametres pour le diagnostic sur les surfaces de contact
        #
        self.maximal_contact_similarity = 0.1
        self.generate_figure = False
        self.figurefile_suffix = ""

    def print_parameters(self):
        print("")
        print('#')
        print('# DiagnosisParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        ucontact.ContactSurfaceParameters.print_parameters(self)

        self.varprint('minimal_volume', self.minimal_volume, self.get('minimal_volume', None))
        self.varprint('maximal_volume_variation', self.maximal_volume_variation,
                      self.get('maximal_volume_variation', None))
        self.varprint('maximal_volume_derivative', self.maximal_volume_derivative,
                      self.get('maximal_volume_derivative', None))

        self.varprint('items', self.items, self.get('items', None))

        self.varprint('minimal_length', self.minimal_length, self.get('minimal_length', None))

        self.varprint('maximal_contact_similarity', self.maximal_contact_similarity,
                      self.get('maximal_contact_similarity', None))
        self.varprint('generate_figure', self.generate_figure, self.get('generate_figure', None))
        self.varprint('figurefile_suffix', self.figurefile_suffix, self.get('figurefile_suffix', None))

        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# DiagnosisParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        ucontact.ContactSurfaceParameters.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'minimal_volume', self.minimal_volume, self.doc.get('minimal_volume', None))
        self.varwrite(logfile, 'maximal_volume_variation', self.maximal_volume_variation,
                      self.doc.get('maximal_volume_variation', None))
        self.varwrite(logfile, 'maximal_volume_derivative', self.maximal_volume_derivative,
                      self.doc.get('maximal_volume_derivative', None))

        self.varwrite(logfile, 'items', self.items, self.doc.get('items', None))

        self.varwrite(logfile, 'minimal_length', self.minimal_length, self.doc.get('minimal_length', None))

        self.varwrite(logfile, 'maximal_contact_similarity', self.maximal_contact_similarity,
                      self.doc.get('maximal_contact_similarity', None))
        self.varwrite(logfile, 'generate_figure', self.generate_figure, self.doc.get('generate_figure', None))
        self.varwrite(logfile, 'figurefile_suffix', self.figurefile_suffix, self.doc.get('figurefile_suffix', None))

        logfile.write("\n")
        return

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            self.write_parameters_in_file(logfile)
        return

    def update_from_parameters(self, parameters):

        ucontact.ContactSurfaceParameters.update_from_parameters(self, parameters)

        self.minimal_volume = self.read_parameter(parameters, 'minimal_volume', self.minimal_volume)
        self.maximal_volume_variation = self.read_parameter(parameters, 'maximal_volume_variation',
                                                            self.maximal_volume_variation)
        self.maximal_volume_derivative = self.read_parameter(parameters, 'maximal_volume_derivative',
                                                             self.maximal_volume_derivative)

        self.items = self.read_parameter(parameters, 'items', self.items)

        self.minimal_length = self.read_parameter(parameters, 'minimal_length', self.minimal_length)

        self.maximal_contact_similarity = self.read_parameter(parameters, 'maximal_contact_similarity',
                                                              self.maximal_contact_similarity)
        self.generate_figure = self.read_parameter(parameters, 'generate_figure', self.generate_figure)
        self.generate_figure = self.read_parameter(parameters, 'generate_figures', self.generate_figure)
        self.figurefile_suffix = self.read_parameter(parameters, 'figurefile_suffix', self.figurefile_suffix)

    def update_from_parameter_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = common.load_source(parameter_file)
        self.update_from_parameters(parameters)

    def update_from_args(self, args):
        if args.diagnosis_minimal_volume is not None:
            self.minimal_volume = args.diagnosis_minimal_volume
        if args.diagnosis_items is not None:
            self.items = args.diagnosis_items


########################################################################################
#
#
#
########################################################################################

def _get_time_interval_from_lineage(direct_lineage, time_digits_for_cell_id=4):

    nodes = list(set(direct_lineage.keys()).union(set([v for values in list(direct_lineage.values()) for v in values])))
    first_time = min(nodes) // 10 ** time_digits_for_cell_id
    last_time = max(nodes) // 10 ** time_digits_for_cell_id

    return first_time, last_time


def _get_time_interval_from_properties(d, time_digits_for_cell_id=4):

    proc = "_get_time_interval_from_properties"
    keylineage = None

    for k in d:
        if k in properties.keydictionary['lineage']['input_keys']:
            keylineage = k

    if keylineage is None:
        monitoring.to_log_and_console(str(proc) + ": no lineage in input dictionary", 1)
        return None, None

    return _get_time_interval_from_lineage(d[keylineage], time_digits_for_cell_id=time_digits_for_cell_id)


########################################################################################
#
#
#
########################################################################################

def _decode_cell_id(s, time_digits_for_cell_id=4):
    div = 10**time_digits_for_cell_id
    return "cell #{:4d}".format(int(s) % div) + " of image #{:4d}".format(int(s) // div)


def _cell_id(c, time_digits_for_cell_id=4):
    t = c // 10 ** time_digits_for_cell_id
    c -= t * 10 ** time_digits_for_cell_id
    return c


def _print_list(prop, tab, time_digits_for_cell_id=4, verbose=2):

    keyname = None
    keynameset = set(prop.keys()).intersection(properties.keydictionary['name']['input_keys'])
    if len(keynameset) > 0:
        keyname = list(keynameset)[0]

    for c in tab:
        t = c // 10**time_digits_for_cell_id
        cell_id = c % 10**time_digits_for_cell_id
        msg = "    - cell #" + str(cell_id) + " of time " + str(t)
        if keyname is not None:
            if c in prop[keyname]:
                msg += " (" + str(prop[keyname][c]) + ")"
        monitoring.to_log_and_console(msg, verbose)


def _get_nodes(prop, key):
    """

    Parameters
    ----------
    prop: property dictionary
    key: property for which the list of cell ids is required

    Returns
    -------
    The list of all cell ids for the given property, or None

    """
    proc = "_get_nodes"
    if key == 'cell_lineage' or key == 'cell_contact_surface':
        d = prop[key]
        nodes = list(set(d.keys()).union(set([v for values in list(d.values()) for v in values])))
        return nodes
    if key == 'cell_volume' or key == 'cell_surface' or key == 'cell_compactness' or key == 'cell_barycenter' or \
            key == 'cell_principal_values' or key == 'cell_name' or key == 'cell_principal_vectors':
        return list(prop[key].keys())
    if key == 'all_cells':
        return prop[key]
    monitoring.to_log_and_console(proc + ": property '" + str(key) + "' not handled yet")
    return None


def _find_node(prop, key, cell_id):

    keyname = None
    keynameset = set(prop.keys()).intersection(properties.keydictionary['name']['input_keys'])
    if len(keynameset) > 0:
        keyname = list(keynameset)[0]

    if key in properties.keydictionary['lineage']['input_keys']:
        pass
    elif key in properties.keydictionary['h_min']['input_keys']:
        pass
    elif key in properties.keydictionary['volume']['input_keys']:
        pass
    elif key in properties.keydictionary['surface']['input_keys']:
        pass
    elif key in properties.keydictionary['compactness']['input_keys']:
        pass
    elif key in properties.keydictionary['sigma']['input_keys']:
        pass
    elif key in properties.keydictionary['label_in_time']['input_keys']:
        pass
    elif key in properties.keydictionary['barycenter']['input_keys']:
        pass
    elif key in properties.keydictionary['fate']['input_keys']:
        pass
    elif key in properties.keydictionary['all-cells']['input_keys']:
        pass
    elif key in properties.keydictionary['principal-value']['input_keys']:
        pass
    elif key in properties.keydictionary['name']['input_keys']:
        pass
    elif key in properties.keydictionary['contact']['input_keys']:
        if cell_id in prop[key]:
            msg = "    - " + str(cell_id)
            if keyname is not None:
                if cell_id in prop[keyname]:
                    msg += " (" + str(prop[keyname][cell_id]) + ")"
            msg += " is indexed in '" + str(key) + "' dictionary"
            monitoring.to_log_and_console(msg)
        cells = []
        for c in prop[key]:
            if cell_id in prop[key][c]:
                cells += [c]
        if len(cells) > 0:
            msg = "    - " + str(cell_id)
            msg += " is a contact cell of"
            for c in cells:
                msg += " " + str(c) + ""
                if keyname is not None:
                    if c in prop[keyname]:
                        msg += " (" + str(prop[keyname][c]) + ")"
                if len(cells) == 2 and cells.index(c) == 0:
                    msg += " and"
                elif len(cells) > 2:
                    if cells.index(c) <= len(cells) - 2:
                        msg += ","
                    if cells.index(c) == len(cells) - 2:
                        msg += " and"

            monitoring.to_log_and_console(msg)

    elif key in properties.keydictionary['history']['input_keys']:
        pass
    elif key in properties.keydictionary['principal-vector']['input_keys']:
        pass
    else:
        monitoring.to_log_and_console("    unknown key '" + str(key) + "' for diagnosis", 1)


########################################################################################
#
#
#
########################################################################################

def _diagnosis_lineage(prop, description, diagnosis_parameters, time_digits_for_cell_id=4):

    keylineage = None
    keyset = set(prop.keys()).intersection(properties.keydictionary['lineage']['input_keys'])
    if len(keyset) > 0:
        keylineage = list(keyset)[0]

    keyname = None
    keyset = set(prop.keys()).intersection(properties.keydictionary['name']['input_keys'])
    if len(keyset) > 0:
        keyname = list(keyset)[0]

    keydiagnosis = 'selection_diagnosis_lineage'
    if keydiagnosis in prop:
        del prop[keydiagnosis]
    prop[keydiagnosis] = {}

    monitoring.to_log_and_console("  === " + str(description) + " diagnosis === ", 1)

    direct_lineage = prop[keylineage]
    first_time, last_time = _get_time_interval_from_lineage(direct_lineage,
                                                            time_digits_for_cell_id=time_digits_for_cell_id)
    monitoring.to_log_and_console("  - estimated time interval = [" + str(first_time) + ", " + str(last_time) + "]", 1)

    #
    # cell with more than 2 daughters
    #
    multiple_daughters = [cell for cell in direct_lineage if len(direct_lineage[cell]) > 2]
    divisions = [cell for cell in direct_lineage if len(direct_lineage[cell]) >= 2]

    #
    # get cells without daughters, remove cells from the last time point
    #
    direct_nodes = list(set(direct_lineage.keys()).union(set([v for values in list(direct_lineage.values())
                                                              for v in values])))
    leaves = set(direct_nodes) - set(direct_lineage.keys())
    early_leaves = [leave for leave in leaves if (leave/10**time_digits_for_cell_id) < last_time]

    #
    # count cells per time
    #
    div = 10 ** time_digits_for_cell_id
    cells_per_time = {}
    for c in direct_nodes:
        t = int(c) // div
        if t not in cells_per_time:
            cells_per_time[t] = [c]
        else:
            cells_per_time[t].append(c)

    monitoring.to_log_and_console("    at time " + str(first_time) + ", #cells = " +
                                  str(len(cells_per_time[first_time])), 1)
    monitoring.to_log_and_console("    at time " + str(last_time) + ", #cells = " +
                                  str(len(cells_per_time[last_time])), 1)

    #
    # build a reverse lineage
    #
    reverse_lineage = {}
    for k, values in direct_lineage.items():
        for v in values:
            if v not in reverse_lineage:
                reverse_lineage[v] = [k]
            else:
                reverse_lineage[v].append(k)

    #
    # branch lengths
    #
    length = {}
    daughters = [direct_lineage[c][0] for c in direct_lineage if len(direct_lineage[c]) == 2]
    daughters += [direct_lineage[c][1] for c in direct_lineage if len(direct_lineage[c]) == 2]
    daughters += [direct_lineage[c][2] for c in direct_lineage if len(direct_lineage[c]) == 3]
    for d in daughters:
        length[d] = 0
        c = d
        while c in direct_lineage and len(direct_lineage[c]) == 1:
            length[d] += 1
            c = direct_lineage[c][0]
        if int(c) // div == last_time:
            del length[d]

    short_length = [[c, length[c]] for c in length if length[c] <= diagnosis_parameters.minimal_length]
    if len(short_length) > 0:
        short_length = sorted(short_length, key=itemgetter(1))

    #
    # get cells with more than 1 mother
    #
    multiple_mothers = [cell for cell in reverse_lineage if len(reverse_lineage[cell]) > 1]

    #
    # get cells without mother, remove cells from the first time point
    #
    reverse_nodes = list(set(reverse_lineage.keys()).union(set([v for values in list(reverse_lineage.values())
                                                                for v in values])))
    orphans = set(reverse_nodes) - set(reverse_lineage.keys())
    late_orphans = [orphan for orphan in orphans if (orphan // 10 ** time_digits_for_cell_id) > first_time]

    monitoring.to_log_and_console("  - found " + str(len(direct_nodes)) + " cells", 1)

    nitems = 10
    if isinstance(diagnosis_parameters, DiagnosisParameters):
        nitems = int(diagnosis_parameters.items)

    if len(late_orphans) > 0:
        late_orphans.sort()
        cell_verbose = monitoring.verbose
        if len(late_orphans) <= nitems:
            cell_verbose = 1
        monitoring.to_log_and_console("  - " + str(len(late_orphans))
                                      + " lineage branches starting after the first time point", 1)
        _print_list(prop, late_orphans, time_digits_for_cell_id=time_digits_for_cell_id, verbose=cell_verbose)
        for c in late_orphans:
            prop[keydiagnosis][c] = 10

    if len(multiple_mothers) > 0:
        multiple_mothers.sort()
        cell_verbose = monitoring.verbose
        if len(multiple_mothers) <= nitems:
            cell_verbose = 1
        monitoring.to_log_and_console("  - " + str(len(multiple_mothers)) + " cells with multiple mother cells", 1)
        _print_list(prop, multiple_mothers, time_digits_for_cell_id=time_digits_for_cell_id, verbose=cell_verbose)
        for c in multiple_mothers:
            prop[keydiagnosis][c] = 20

    if len(leaves) > 0:
        monitoring.to_log_and_console("  - " + str(len(leaves))
                                      + " lineage terminal branches", 1)
    if len(early_leaves) > 0:
        early_leaves.sort()
        cell_verbose = monitoring.verbose
        if len(early_leaves) <= nitems:
            cell_verbose = 1
        monitoring.to_log_and_console("  - " + str(len(early_leaves))
                                      + " lineage terminal branches ending before the last time point", 1)
        _print_list(prop, early_leaves, time_digits_for_cell_id=time_digits_for_cell_id, verbose=cell_verbose)
        for c in early_leaves:
            prop[keydiagnosis][c] = 30

    if len(divisions) > 0:
        divisions.sort()
        monitoring.to_log_and_console("  - " + str(len(divisions)) + " cell divisions", 1)
    if len(multiple_daughters) > 0:
        multiple_daughters.sort()
        cell_verbose = monitoring.verbose
        if len(multiple_daughters) <= nitems:
            cell_verbose = 1
        monitoring.to_log_and_console("  - " + str(len(multiple_daughters))
                                      + " divisions yielding more than 2 branches", 1)
        _print_list(prop, multiple_daughters, time_digits_for_cell_id=time_digits_for_cell_id, verbose=cell_verbose)
        for c in multiple_daughters:
            prop[keydiagnosis][c] = 30

    if len(short_length) > 0:
        monitoring.to_log_and_console("  - " + str(len(short_length)) + " non-terminal branch length < " +
                                      str(diagnosis_parameters.minimal_length), 1)
        if len(short_length) <= diagnosis_parameters.items:
            nitems = len(short_length)
        else:
            nitems = int(diagnosis_parameters.items)
        for i in range(nitems):
            msg = "    - " + str(short_length[i][0])
            if keyname is not None:
                if short_length[i][0] in prop[keyname]:
                    msg += " (" + str(prop[keyname][short_length[i][0]]) + ")"
            msg += " has a branch/life length of " + str(short_length[i][1])
            monitoring.to_log_and_console(msg, 1)
        if nitems < len(short_length):
            msg = "    ... "
            monitoring.to_log_and_console(msg, 1)
        for c in short_length:
            prop[keydiagnosis][c[0]] = 40

    monitoring.to_log_and_console("")
    return prop


# def _diagnosis_h_min(d):
#     return


def _diagnosis_volume(prop, description, diagnosis_parameters, time_digits_for_cell_id=4):
    """

    :param prop:
    :param description:
    :param diagnosis_parameters:
    :return:
    """

    #     dictionary de int
    #     cell_volume.590002 = <type 'int'>
    #     590002: 236936

    keyname = None
    keyset = set(prop.keys()).intersection(properties.keydictionary['name']['input_keys'])
    if len(keyset) > 0:
        keyname = list(keyset)[0]

    keylineage = None
    keyset = set(prop.keys()).intersection(properties.keydictionary['lineage']['input_keys'])
    if len(keyset) > 0:
        keylineage = list(keyset)[0]

    keyvolume = None
    keyset = set(prop.keys()).intersection(properties.keydictionary['volume']['input_keys'])
    if len(keyset) > 0:
        keyvolume = list(keyset)[0]

    keydiagnosis = 'selection_diagnosis_volume'
    if keydiagnosis in prop:
        del prop[keydiagnosis]
    prop[keydiagnosis] = {}

    monitoring.to_log_and_console("  === " + str(description) + " diagnosis === ", 1)

    all_cell_with_volume = set(prop[keyvolume].keys())
    cell_with_volume = [c for c in all_cell_with_volume if _cell_id(c, time_digits_for_cell_id) != 1]
    background_with_volume = [c for c in all_cell_with_volume if _cell_id(c, time_digits_for_cell_id) == 1]

    volume = [[c, prop[keyvolume][c]] for c in cell_with_volume]
    volume = sorted(volume, key=itemgetter(1))

    monitoring.to_log_and_console("    found " + str(len(background_with_volume)) + " background cells with volume", 1)
    monitoring.to_log_and_console("    found " + str(len(cell_with_volume)) + " cells with volume", 1)

    monitoring.to_log_and_console("")
    monitoring.to_log_and_console("  - smallest volumes", 1)

    d = DiagnosisParameters()
    n = int(d.items)
    v = int(d.minimal_volume)

    if isinstance(diagnosis_parameters, DiagnosisParameters):
        n = int(diagnosis_parameters.items)
        v = int(diagnosis_parameters.minimal_volume)

    for i in range(len(list(prop[keyvolume].keys()))):
        if (n > 0 and i < n) or (v > 0 and int(volume[i][1]) <= v):
            msg = _decode_cell_id(volume[i][0])
            if keyname is not None:
                if volume[i][0] in prop[keyname]:
                    msg += " (" + str(prop[keyname][volume[i][0]]) + ")"
            msg += " has volume = " + str(volume[i][1])
            monitoring.to_log_and_console("    " + msg, 1)

    #
    # volume variation
    #
    lineage = prop[keylineage]
    div = 10 ** time_digits_for_cell_id
    first_time, last_time = _get_time_interval_from_lineage(lineage, time_digits_for_cell_id=time_digits_for_cell_id)
    first_cells = [lineage[c][0] for c in lineage if len(lineage[c]) == 2 and int(c) // div < last_time-1]
    first_cells += [lineage[c][1] for c in lineage if len(lineage[c]) == 2 and int(c) // div < last_time-1]

    volume_along_time = {}
    volume_derivative = {}

    for cell in first_cells:
        if cell not in prop[keyvolume] or cell not in lineage:
            continue
        volume_along_time[cell] = volume_along_time.get(cell, []) + [prop[keyvolume][cell]]
        volume_derivative[cell] = []
        pcell = cell
        while True:
            ncell = lineage[pcell][0]
            if ncell not in prop[keyvolume]:
                break
            volume_along_time[cell] += [prop[keyvolume][ncell]]
            volume_derivative[cell] += [100.0 * (prop[keyvolume][ncell] - prop[keyvolume][pcell]) /
                                        prop[keyvolume][pcell]]
            pcell = ncell
            if pcell not in lineage or len(lineage[pcell]) > 1:
                break

    volume_variation = {}
    volume_max_derivative = {}
    for cell in volume_along_time:
        if len(volume_along_time[cell]) <= 1:
            msg = "    * weird, " + _decode_cell_id(cell)
            if keyname is not None:
                if cell in prop[keyname]:
                    msg += " (" + str(prop[keyname][cell]) + ")"
            msg += " has a branch (life) length of " + str(len(volume_along_time[cell]))
            monitoring.to_log_and_console(msg, 1)
            continue
        volume_variation[cell] = 100.0 * (max(volume_along_time[cell]) - min(volume_along_time[cell])) / \
                                 statistics.median(volume_along_time[cell])
        if max(volume_derivative[cell]) > 0:
            volume_max_derivative[cell] = max(volume_derivative[cell])
        if min(volume_derivative[cell]) < 0:
            volume_max_derivative[cell] = max(-min(volume_derivative[cell]), volume_max_derivative.get(cell, 0.0))

    volume = [[c, volume_variation[c]] for c in volume_variation]
    volume = sorted(volume, key=itemgetter(1), reverse=True)

    monitoring.to_log_and_console("")
    monitoring.to_log_and_console("  - largest volume variation [100.0 * (max-min)/median]", 1)

    if isinstance(diagnosis_parameters, DiagnosisParameters):
        n = int(diagnosis_parameters.items)
        v = diagnosis_parameters.maximal_volume_variation

    for i in range(len(volume_variation.keys())):
        if (n > 0 and i < n) or (0 < v <= volume[i][1]):
            msg = _decode_cell_id(volume[i][0])
            if keyname is not None:
                if volume[i][0] in prop[keyname]:
                    msg += " (" + str(prop[keyname][volume[i][0]]) + ")"
            msg += " has volume variation = {:.2f}".format(volume[i][1])
            msg += " [branch length = " + str(len(volume_along_time[volume[i][0]])) + "]"
            monitoring.to_log_and_console("    " + msg, 1)

    volume = [[c, volume_max_derivative[c]] for c in volume_max_derivative]
    volume = sorted(volume, key=itemgetter(1), reverse=True)

    monitoring.to_log_and_console("")
    monitoring.to_log_and_console("  - largest volume derivative [100.0 * (v[t+1]-v[t])/v[t]]", 1)

    if isinstance(diagnosis_parameters, DiagnosisParameters):
        n = int(diagnosis_parameters.items)
        v = diagnosis_parameters.maximal_volume_derivative

    for i in range(len(volume_max_derivative.keys())):
        if (n > 0 and i < n) or (0 < v <= volume[i][1]):
            t = volume[i][0] // div
            msg = _decode_cell_id(volume[i][0])
            if keyname is not None:
                if volume[i][0] in prop[keyname]:
                    msg += " (" + str(prop[keyname][volume[i][0]]) + ")"
            msg += " has volume maximal derivative = {:.2f}".format(volume[i][1])
            msg += " [branch length = " + str(len(volume_along_time[volume[i][0]])) + "]"
            monitoring.to_log_and_console("    " + msg, 1)
            cell = volume[i][0]
            c = cell
            for j in range(len(volume_derivative[cell])):
                if (0 < v <= abs(volume_derivative[cell][j])) or abs(volume_derivative[cell][j]) > volume[i][1] - 0.1:
                    msg = "    - from time " + str(t+j) + " to " + str(t+j+1)
                    msg += ", volume evolves from " + str(volume_along_time[cell][j]) + " to "
                    msg += str(volume_along_time[cell][j+1])
                    msg += " (derivative = {:.2f}".format(volume_derivative[cell][j]) + ")"
                    monitoring.to_log_and_console(msg, 1)
                    if c not in prop[keydiagnosis]:
                        prop[keydiagnosis][c] = int(abs(volume_derivative[cell][j]))
                    else:
                        prop[keydiagnosis][c] = max(prop[keydiagnosis][c], int(abs(volume_derivative[cell][j])))
                    prop[keydiagnosis][lineage[c][0]] = int(abs(volume_derivative[cell][j]))
                    c = lineage[c][0]

    for c in prop[keydiagnosis]:
        if prop[keydiagnosis][c] > 255:
            prop[keydiagnosis][c] = 255

    monitoring.to_log_and_console("")
    return prop


# def _diagnosis_sigma(d):
#     return


# def _diagnosis_label_in_time(d):
#     return


# def _diagnosis_barycenter(d):
#     return


# def _diagnosis_fate(d):
#     return


# def _diagnosis_all_cells(d):
#     return


# def _diagnosis_principal_value(d):
#     return


def _diagnosis_name(prop, description, time_digits_for_cell_id=4):

    keyname = None
    keyset = set(prop.keys()).intersection(properties.keydictionary['name']['input_keys'])
    if len(keyset) > 0:
        keyname = list(keyset)[0]

    keylineage = None
    keyset = set(prop.keys()).intersection(properties.keydictionary['lineage']['input_keys'])
    if len(keyset) > 0:
        keylineage = list(keyset)[0]

    keydiagnosis = 'selection_diagnosis_name'
    if keydiagnosis in prop:
        del prop[keydiagnosis]
    prop[keydiagnosis] = {}


    monitoring.to_log_and_console("  === " + str(description) + " diagnosis === ", 1)

    name = prop[keyname]
    lineage = prop[keylineage]
    reverse_lineage = {v: k for k, values in lineage.items() for v in values}

    div = 10 ** time_digits_for_cell_id

    cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
    cells = sorted(cells)

    cells_per_time = {}
    names_per_time = {}
    missing_name = {}
    error_name = {}

    for c in cells:
        t = int(c) // div
        #
        # get cells and cell names at each time point
        #
        cells_per_time[t] = cells_per_time.get(t, []) + [c]
        if c in name:
            names_per_time[t] = names_per_time.get(t, []) + [name[c]]
        else:
            missing_name[t] = missing_name.get(t, []) + [c]

        if c not in name:
            if c in reverse_lineage:
                mother = reverse_lineage[c]
                if len(lineage[mother]) == 1 and mother in name:
                    prop[keydiagnosis][c] = 10
                    msg = "    cell " + str(c) + " has no name"
                    msg += ", but its mother cell " + str(mother) + " has a name " + str(name[mother])
                    monitoring.to_log_and_console(msg)
        elif c in reverse_lineage:
            mother = reverse_lineage[c]
            if mother not in name:
                prop[keydiagnosis][c] = 20
                msg = "    cell " + str(c) + " has a name '" + str(name[c]) + "'"
                msg += ", but its mother cell " + str(mother) + " has no name"
                monitoring.to_log_and_console(msg)
                error_name[t] = error_name.get(t, []) + [c]
            else:
                if len(lineage[mother]) == 1:
                    if name[mother] != name[c]:
                        prop[keydiagnosis][c] = 30
                        msg = "    cell " + str(c) + " has a name = " + str(name[c])
                        msg += " different than its mother cell " + str(mother) + " name = " + str(name[mother])
                        monitoring.to_log_and_console(msg)
                        error_name[t] = error_name.get(t, []) + [c]
                elif len(lineage[mother]) == 2:
                    daughter_names = uname.get_daughter_names(name[mother])
                    #
                    # check whether name[c] is in daughters' names
                    #
                    if name[c] not in daughter_names:
                        prop[keydiagnosis][c] = 40
                        msg = "    cell " + str(c) + " is named " + str(name[c])
                        msg += " but should be in " + str(daughter_names)
                        msg += " since its mother cell " + str(mother) + " is named " + str(name[mother])
                        monitoring.to_log_and_console(msg)
                        error_name[t] = error_name.get(t, []) + [c]
                    #
                    # check whether c's sister name is in daughters' names
                    #
                    else:
                        siblings = copy.deepcopy(lineage[mother])
                        siblings.remove(c)
                        daughter_names.remove(name[c])
                        if siblings[0] not in name:
                            prop[keydiagnosis][c] = 50
                            msg = "    cell " + str(siblings[0]) + " has no name "
                            msg += ", it should be " + str(daughter_names[0])
                            msg += " since its mother cell " + str(mother) + " is named " + str(name[mother])
                            msg += " and its sibling " + str(c) + " is named " + str(name[c])
                            monitoring.to_log_and_console(msg)
                            error_name[t] = error_name.get(t, []) + [siblings[0]]
                        elif name[siblings[0]] == name[c]:
                            prop[keydiagnosis][c] = 60
                            msg = "    cell " + str(siblings[0]) + " as named " + str(name[c])
                            msg += " as its sibling " + str(c)
                            msg += ", their mother cell " + str(mother) + " is named " + str(name[mother])
                            monitoring.to_log_and_console(msg)
                            error_name[t] = error_name.get(t, []) + [c]
                        elif name[siblings[0]] != daughter_names[0]:
                            prop[keydiagnosis][c] = 70
                            msg = "    cell " + str(siblings[0]) + " is named " + str(name[siblings[0]])
                            msg += " but should be " + str(daughter_names[0])
                            msg += " since its mother cell " + str(mother) + " is named " + str(name[mother])
                            msg += " and its sibling " + str(c) + " is named " + str(name[c])
                            monitoring.to_log_and_console(msg)
                            error_name[t] = error_name.get(t, []) + [siblings[0]]
                else:
                    prop[keydiagnosis][mother] = 80
                    msg = "    cell " + str(mother) + " has " + str(len(lineage[mother])) + " daughter cells"
                    msg += ": " + str(lineage[mother])
                    monitoring.to_log_and_console(msg)

    for t in names_per_time:
        repeats = {k: names_per_time[t].count(k) for k in set(names_per_time[t]) if names_per_time[t].count(k) > 1}
        if repeats != {}:
            msg = "  - there are " + str(len(repeats)) + " repeated names at time " + str(t)
            monitoring.to_log_and_console(msg)
            for n, p in repeats.items():
                cells = [c for c in name if name[c] == n]
                for c in cells:
                    prop[keydiagnosis][c] = 90
                msg = "    - " + str(n) + " is repeated " + str(p) + " times "
                monitoring.to_log_and_console(msg)

    monitoring.to_log_and_console("  - first time in lineage = " + str(min(cells_per_time.keys())))
    monitoring.to_log_and_console("  - last time in lineage = " + str(max(cells_per_time.keys())))
    monitoring.to_log_and_console("    - cells in lineage = " + str(len(cells)))
    monitoring.to_log_and_console("    - #cells at first time = " +
                                  str(len(cells_per_time[min(cells_per_time.keys())])))
    monitoring.to_log_and_console("    - #cells at last time = " +
                                  str(len(cells_per_time[max(cells_per_time.keys())])))
    monitoring.to_log_and_console("\t - names in lineage = " +
                                  str(len(list(collections.Counter(list(name.keys())).keys()))))
    monitoring.to_log_and_console("")

    return prop


# def _diagnosis_contact(d):
#     return

def _diagnosis_contact(prop, description, diagnosis_parameters, time_digits_for_cell_id=4):
    proc = "_diagnosis_contact"

    keycontact = None
    keyset = set(prop.keys()).intersection(properties.keydictionary['contact']['input_keys'])
    if len(keyset) > 0:
        keycontact = list(keyset)[0]

    keyname = None
    keyset = set(prop.keys()).intersection(properties.keydictionary['name']['input_keys'])
    if len(keyset) > 0:
        keyname = list(keyset)[0]

    keylineage = None
    keyset = set(prop.keys()).intersection(properties.keydictionary['lineage']['input_keys'])
    if len(keyset) > 0:
        keylineage = list(keyset)[0]

    keydiagnosis = 'selection_diagnosis_contact_surface'
    if keydiagnosis in prop:
        del prop[keydiagnosis]
    prop[keydiagnosis] = {}

    monitoring.to_log_and_console("  === " + str(description) + " diagnosis === ", 1)

    lineage = prop[keylineage]
    name = prop[keyname]
    contact = prop[keycontact]

    #
    # get the first time point
    #
    cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
    cells = sorted(cells)
    div = 10 ** time_digits_for_cell_id
    cells_per_time = {}
    for c in cells:
        t = int(c) // div
        #
        # get cells and cell names at each time point
        #
        cells_per_time[t] = cells_per_time.get(t, []) + [c]
    last_time = max(cells_per_time.keys())

    #
    # study single branches beginning right after a division, but not at the last time
    #
    first_cells = [lineage[c][0] for c in lineage if len(lineage[c]) == 2 and int(c) // div < last_time - 1]
    first_cells += [lineage[c][1] for c in lineage if len(lineage[c]) == 2 and int(c) // div < last_time - 1]

    reverse_lineage = {v: k for k, values in lineage.items() for v in values}

    score_along_time = {}
    cells_not_in_reversed_lineage = []
    cells_with_problem = []

    for cell in first_cells:
        if cell not in contact:
            msg = "    * weird, cell " + str(cell) + " is not in the 'contact surface' dictionary"
            monitoring.to_log_and_console(msg)
            continue
        first_time = int(cell) // div
        pcell = cell
        pneigh = copy.deepcopy(contact[cell])
        while True:
            ncell = lineage[pcell][0]
            t = int(ncell) // div
            nneigh = {}
            #
            # build a neighborhood with the same label than at the first time of the branch
            #
            if ncell not in contact:
                break
            for c in contact[ncell]:
                contrib = contact[ncell][c]
                if int(c) % div == 1 or int(c) % div == 0:
                    nneigh[1] = nneigh.get(1, 0.0) + contrib
                else:
                    for i in range(t - first_time):
                        if c not in reverse_lineage:
                            if c not in cells_not_in_reversed_lineage:
                                msg = "    * weird, cell " + str(c) + " is not in the reversed lineage"
                                monitoring.to_log_and_console(msg)
                                cells_not_in_reversed_lineage += [c]
                            if cell not in cells_with_problem:
                                msg = "    - diagnosis on "
                                msg += _decode_cell_id(cell, time_digits_for_cell_id=time_digits_for_cell_id)
                                msg += " should be handled with care"
                                monitoring.to_log_and_console(msg)
                                cells_with_problem += [cell]
                            break
                        c = reverse_lineage[c]
                    nneigh[c] = nneigh.get(c, 0.0) + contrib
            #
            #
            #
            score = ucontact.contact_distance(pneigh, nneigh, similarity=diagnosis_parameters.contact_similarity,
                                             change_contact_surfaces=False )
            score_along_time[cell] = score_along_time.get(cell, []) + [score]
            #
            #
            #
            pcell = ncell
            pneigh = copy.deepcopy(nneigh)
            if pcell not in lineage or len(lineage[pcell]) > 1:
                break

    if diagnosis_parameters.generate_figure:
        filename = 'figures_contact_score_evolution'
        file_suffix = ""
        if diagnosis_parameters.figurefile_suffix is not None and \
                isinstance(diagnosis_parameters.figurefile_suffix, str) and \
                len(diagnosis_parameters.figurefile_suffix) > 0:
            file_suffix = '_' + diagnosis_parameters.figurefile_suffix
        filename += file_suffix + '.py'

        f = open(filename, "w")

        f.write("import numpy as np\n")
        f.write("import matplotlib.pyplot as plt\n")
        f.write("import scipy.stats as stats\n")

        f.write("\n")
        f.write("savefig = True\n")

        f.write("\n")
        f.write("score_along_time = " + str(score_along_time) + "\n")
        f.write("lengths = [len(score_along_time[c]) for c in score_along_time]\n")

        f.write("\n")
        f.write("score_per_time = {}\n")
        f.write("for c in score_along_time:\n")
        f.write("    for i in range(len(score_along_time[c])):\n")
        f.write("        score_per_time[i] = score_per_time.get(i, []) + [score_along_time[c][i]]\n")
        f.write("scores_per_time = []\n")
        f.write("for i in range(max(lengths)):\n")
        f.write("    scores_per_time.append(score_per_time[i])\n")

        f.write("\n")
        f.write("fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(16, 6.5))\n")

        f.write("for c in score_along_time:\n")
        f.write("    ax1.plot(range(len(score_along_time[c])), score_along_time[c], 'k', markersize=1)\n")
        f.write("ax1.set_ylim(0, 0.7)\n")
        f.write("ax1.set_xlabel('time from division', fontsize=15)\n")
        f.write("ax1.set_ylabel('Score', fontsize=15)\n")
        f.write("title = 'Score/distance from t to t+1'\n")
        if diagnosis_parameters.figurefile_suffix is not None and \
                isinstance(diagnosis_parameters.figurefile_suffix, str) and \
                len(diagnosis_parameters.figurefile_suffix) > 0:
            f.write("title += ': ' + '" + str(diagnosis_parameters.figurefile_suffix) + "'\n")
        f.write("ax1.set_title(title, fontsize=15)\n")

        f.write("\n")
        f.write("ax2.boxplot(scores_per_time)\n")
        f.write("ax2.set_title(\"Box plot representation\", fontsize=15)\n")
        f.write("\n")
        f.write("for c in score_along_time:\n")
        f.write("    ax3.plot(range(len(score_along_time[c])), score_along_time[c], 'k', markersize=1)\n")
        f.write("ax3.set_xlabel('time from division', fontsize=15)\n")
        f.write("ax3.set_xlim(0, 5)\n")
        f.write("title = 'Score/distance from t to t+1'\n")
        f.write("ax3.set_title(title, fontsize=15)\n")
        f.write("\n")
        f.write("if savefig:\n")
        f.write("    plt.savefig('contact_score_evolution")
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '.png')\n")
        f.write("else:\n")
        f.write("    plt.show()\n")
        f.close()

    #
    # report
    #

    for cell in score_along_time:
        first_time = int(cell) // div
        msg = "    " + _decode_cell_id(cell, time_digits_for_cell_id=time_digits_for_cell_id)
        if cell in name:
            msg += " (" + name[cell] + ") "
        msg += " branch length of " + str(len(score_along_time[cell]))
        if cell in cells_with_problem:
            msg += " (warning: some cells were not in the reversed lineage)"
        #
        # branch of length 1
        #
        if len(score_along_time[cell]) <= 1:
            if first_time == last_time - 1:
                continue
            monitoring.to_log_and_console(msg)
            continue
        #
        # branch with all distances (between successive contact surface vectors) below the threshold
        #
        if max(score_along_time[cell][1:]) < diagnosis_parameters.maximal_contact_similarity:
            continue
        #
        # here, the branch has at least a couple of successive successive contact surface vectors has
        # a distance above the threshold
        #
        monitoring.to_log_and_console(msg)
        for i in range(1, len(score_along_time[cell])):
            if score_along_time[cell][i] < diagnosis_parameters.maximal_contact_similarity:
                continue
            c = cell
            for j in range(i):
                c = lineage[c][0]
            if c not in prop[keydiagnosis]:
                prop[keydiagnosis][c] = int(255 * score_along_time[cell][i])
            else:
                prop[keydiagnosis][c] = max(prop[keydiagnosis][c], int(255 * score_along_time[cell][i]))
            prop[keydiagnosis][lineage[c][0]] = int(255 * score_along_time[cell][i])
            msg = "    - distance of " + str(score_along_time[cell][i]) + " between "
            msg += "times " + str(first_time + i) + " and " + str(first_time + i + 1)
            msg += " (cells " + str(c) + " and " + str(lineage[c][0]) + ")"
            monitoring.to_log_and_console(msg)

    monitoring.to_log_and_console("")
    return prop

# def _diagnosis_history(d):
#     return


# def _diagnosis_principal_vector(d):
#     return

def _diagnosis_one_feature(prop, feature, diagnosis_parameters, time_digits_for_cell_id=4):
    proc = "_diagnosis_one_feature"
    
    #
    #
    #
    key = None
    for k in properties.keydictionary:
        if feature == k or feature in properties.keydictionary[k]['input_keys']:
            key = k
    if key is None:
        monitoring.to_log_and_console(proc + ": property '" + str(feature) + "' not found in 'keydictionary'", 1)
        return prop

    #
    #
    #
    propkeys = list(prop.keys())
    for k in propkeys:
        if k not in properties.keydictionary[key]['input_keys']:
            continue

        if k in properties.keydictionary['lineage']['input_keys']:
            prop = _diagnosis_lineage(prop, k, diagnosis_parameters=diagnosis_parameters,
                                      time_digits_for_cell_id=time_digits_for_cell_id)
        elif k in properties.keydictionary['h_min']['input_keys']:
            pass
        elif k in properties.keydictionary['volume']['input_keys']:
            prop = _diagnosis_volume(prop, k, diagnosis_parameters=diagnosis_parameters,
                                     time_digits_for_cell_id=time_digits_for_cell_id)
        elif k in properties.keydictionary['surface']['input_keys']:
            pass
        elif k in properties.keydictionary['compactness']['input_keys']:
            pass
        elif k in properties.keydictionary['sigma']['input_keys']:
            pass
        elif k in properties.keydictionary['label_in_time']['input_keys']:
            pass
        elif k in properties.keydictionary['barycenter']['input_keys']:
            pass
        elif k in properties.keydictionary['fate']['input_keys']:
            pass
        elif k in properties.keydictionary['all-cells']['input_keys']:
            pass
        elif k in properties.keydictionary['principal-value']['input_keys']:
            pass
        elif k in properties.keydictionary['name']['input_keys']:
            prop = _diagnosis_name(prop, k, time_digits_for_cell_id=time_digits_for_cell_id)
        elif k in properties.keydictionary['contact']['input_keys']:
            prop = _diagnosis_contact(prop, k, diagnosis_parameters=diagnosis_parameters,
                                      time_digits_for_cell_id=time_digits_for_cell_id)
        elif k in properties.keydictionary['history']['input_keys']:
            pass
        elif k in properties.keydictionary['principal-vector']['input_keys']:
            pass
        else:
            monitoring.to_log_and_console("    unknown key '" + str(k) + "' for diagnosis", 1)
    return prop


def diagnosis(prop, features=None, parameters=None, time_digits_for_cell_id=4):
    """

    :param prop: property dictionary
    :param features:
    :param parameters:
    :param time_digits_for_cell_id:
    :return:
    """

    if parameters is not None and isinstance(parameters, DiagnosisParameters):
        diagnosis_parameters = parameters
    else:
        diagnosis_parameters = DiagnosisParameters()

    keyname = None
    keynameset = set(prop.keys()).intersection(properties.keydictionary['name']['input_keys'])
    if len(keynameset) > 0:
        keyname = list(keynameset)[0]

    # monitoring.to_log_and_console("\n", 1)
    monitoring.to_log_and_console("... diagnosis", 1)

    monitoring.to_log_and_console("  === cell/key diagnosis === ", 1)
    div = 10 ** time_digits_for_cell_id
    # get nodes (ie cells) for each property
    # remove background
    nodes = {}
    for k in prop:
        nodes[k] = _get_nodes(prop, k)
        if nodes[k] is None:
            del nodes[k]
            continue

    #
    # background cells may remain, remove them again
    # check whether the label 0 exists
    #
    zeros = {}
    for k in nodes:
        nodes[k] = [c for c in nodes[k] if int(c) % div != 1]
        for c in nodes[k]:
            if int(c) % div == 0:
                zeros[k] = zeros.get(k, []) + [c]
                nodes[k].remove(c)
    if len(zeros) > 0:
        for k in zeros:
            msg = "  - property '" + str(k) + "': found " + str(len(zeros[k])) + " '0' labels"
            monitoring.to_log_and_console(msg, 1)
            monitoring.to_log_and_console("    " + str(zeros[k]), 1)
            for z in zeros[k]:
                _find_node(prop, k, z)
    del zeros

    #
    # check whether the same set of cell ids are indexed or used in each property
    #
    inodes = set()
    unodes = set()
    for k in nodes:
        unodes = unodes.union(set(nodes[k]))
    inodes = copy.deepcopy(unodes)
    for k in nodes:
        inodes = inodes.intersection(set(nodes[k]))
    if len(unodes) > len(inodes):
        diff = sorted(unodes - inodes)
        msg = "  - " + str(len(diff))
        msg += " cells (excluding 1s and 0s) are not in all features"
        monitoring.to_log_and_console(msg, 1)
        if len(diff) <= diagnosis_parameters.items:
            length = len(diff)
        else:
            length = diagnosis_parameters.items

        for i in range(length):
            msg = "    - " + str(diff[i])
            if keyname is not None:
                if diff[i] in prop[keyname]:
                    msg += " (" + str(prop[keyname][diff[i]]) + ")"
            insets = [k for k in nodes if diff[i] in nodes[k]]
            outsets = [k for k in nodes if diff[i] not in nodes[k]]
            if len(insets) > 0:
                msg += " is in " + str(insets)
            if len(outsets) > 0:
                if len(insets) > 0:
                    msg += " and"
                msg += " is not in " + str(outsets)
            monitoring.to_log_and_console(msg, 1)
        if length < len(diff):
            msg = "    ... "
            monitoring.to_log_and_console(msg, 1)

    monitoring.to_log_and_console("")

    #
    # diagnosis on all features of the property dictionary, or only on requested features
    #
    propkeys = list(prop.keys())
    if features is None or (isinstance(features, list) and len(features) == 0) or \
            (isinstance(features, str) and len(features) == 0):
        for k in propkeys:
            prop = _diagnosis_one_feature(prop, k, diagnosis_parameters,
                                          time_digits_for_cell_id=time_digits_for_cell_id)
    else:
        if isinstance(features, str):
            prop = _diagnosis_one_feature(prop, features, diagnosis_parameters,
                                   time_digits_for_cell_id=time_digits_for_cell_id)
        elif isinstance(features, list):
            for f in features:
                prop = _diagnosis_one_feature(prop, f, diagnosis_parameters,
                                              time_digits_for_cell_id=time_digits_for_cell_id)

    monitoring.to_log_and_console("  === end of diagnosis === ", 1)
    monitoring.to_log_and_console("")

    return prop





