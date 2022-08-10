
import copy
import os
import sys
import operator
import numpy as np
import scipy as sp
import scipy.cluster.hierarchy as sch

import astec.utils.common as common
import astec.utils.ascidian_name as uname
import astec.utils.neighborhood_distance as uneighborhood
import astec.utils.atlas_embryo as uatlase
import astec.utils.atlas_cell as uatlasc

monitoring = common.Monitoring()


###########################################################
#
#
#
############################################################

class DivisionParameters(uatlase.AtlasParameters):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix='atlas_'):

        if "doc" not in self.__dict__:
            self.doc = {}

        uatlase.AtlasParameters.__init__(self, prefix=prefix)

        #
        #
        #
        doc = "\t True or False. Exclude inner surfaces from the division-to-division distance calculation\n"
        self.doc['exclude_inner_surfaces'] = doc
        self.exclude_inner_surfaces = False

        #
        #
        #
        doc = "\t True or False. Performs some diagnosis after building the division atlas. \n"
        doc += "\t Incrementing the verboseness ('-v' in the command line) may give more details.\n"
        self.doc['division_diagnosis'] = doc
        self.division_diagnosis = False

        doc = "\t If True, will propose some daughters switches in the atlases. For a given division,\n"
        doc += "\t a global score is computed as the sum of all pairwise division similarity.\n"
        doc += "\t A switch is proposed for an atlas if it allows to decrease this global score.\n"
        self.doc['division_permutation_proposal'] = doc
        self.division_permutation_proposal = False

        #
        #
        #
        doc = "\t Cluster distance used to build dendrograms. Dendrograms are used either for\n"
        doc += "\t diagnosis purpose (if 'diagnosis_properties' is set to True) or to generate\n"
        doc += "\t figures (if 'generate_figure' is set to True)\n"
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

        doc = "\t Write out morphonet selection files."
        self.doc['write_selection'] = doc
        self.write_selection = False

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
        print('# DivisionParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        uatlase.AtlasParameters.print_parameters(self)

        self.varprint('exclude_inner_surfaces', self.exclude_inner_surfaces)

        self.varprint('division_diagnosis', self.division_diagnosis)
        self.varprint('division_permutation_proposal', self.division_permutation_proposal)

        self.varprint('dendrogram_cluster_distance', self.dendrogram_cluster_distance)
        self.varprint('write_selection', self.write_selection)

        self.varprint('cells_to_be_traced', self.cells_to_be_traced)
        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# DivisionParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        uatlase.AtlasParameters.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'exclude_inner_surfaces', self.exclude_inner_surfaces,
                      self.doc.get('exclude_inner_surfaces', None))

        self.varwrite(logfile, 'division_diagnosis', self.division_diagnosis, self.doc.get('division_diagnosis', None))
        self.varwrite(logfile, 'division_permutation_proposal', self.division_permutation_proposal,
                      self.doc.get('division_permutation_proposal', None))

        self.varwrite(logfile, 'dendrogram_cluster_distance', self.dendrogram_cluster_distance,
                      self.doc.get('dendrogram_cluster_distance', None))
        self.varwrite(logfile, 'write_selection', self.write_selection, self.doc.get('write_selection', None))

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

        uatlase.AtlasParameters.update_from_parameters(self, parameters)

        self.exclude_inner_surfaces = self.read_parameter(parameters, 'exclude_inner_surfaces',
                                                          self.exclude_inner_surfaces)

        self.division_diagnosis = self.read_parameter(parameters, 'division_diagnosis', self.division_diagnosis)
        self.division_diagnosis = self.read_parameter(parameters, 'diagnosis_properties', self.division_diagnosis)
        self.division_diagnosis = self.read_parameter(parameters, 'naming_diagnosis', self.division_diagnosis)
        self.division_diagnosis = self.read_parameter(parameters, 'diagnosis_naming', self.division_diagnosis)

        self.division_permutation_proposal = self.read_parameter(parameters, 'division_permutation_proposal',
                                                                 self.division_permutation_proposal)
        self.division_permutation_proposal = self.read_parameter(parameters, 'daughter_switch_proposal',
                                                                 self.division_permutation_proposal)

        self.dendrogram_cluster_distance = self.read_parameter(parameters, 'dendrogram_cluster_distance',
                                                               self.dendrogram_cluster_distance)
        self.write_selection = self.read_parameter(parameters, 'write_selection', self.write_selection)

        self.cells_to_be_traced = self.read_parameter(parameters, 'cells_to_be_traced', self.cells_to_be_traced)

    def update_from_parameter_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = common.load_source(parameter_file)
        self.update_from_parameters(parameters)


###########################################################
#
#
#
############################################################

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
        new_neighborhoods = uneighborhood.build_same_contact_surfaces(neighborhoods, [cell, sister])
        for n in new_neighborhoods:
            common_neighborhoods[n] = copy.deepcopy(new_neighborhoods[n])
    return common_neighborhoods


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
            if 2 * s > len(divisions[n]):
                majority[n][a] = s
            elif 2 * s == len(divisions[n]):
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
    neighborhoods = atlases.get_cell_neighborhood(delay_from_division=parameters.name_delay_from_division)

    summary = {}
    innersurfaces = []

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
        if parameters.exclude_inner_surfaces:
            innersurfaces = [d[0], d[1]]

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
                same_dist = division_distance(neighborhoods[d[0]][r1], neighborhoods[d[1]][r1],
                                              neighborhoods[d[0]][r2], neighborhoods[d[1]][r2],
                                              change_contact_surfaces=ccs, innersurfaces=innersurfaces)
                swit_dist = division_distance(neighborhoods[d[0]][r1], neighborhoods[d[1]][r1],
                                              switch_neighs[d[0]], switch_neighs[d[1]],
                                              change_contact_surfaces=ccs, innersurfaces=innersurfaces)
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


###########################################################
#
#
#
############################################################


def division_scipy_linkage(config, cluster_distance='single', change_contact_surfaces=True, innersurfaces=[],
                           distance='distance'):
    """

    Parameters
    ----------
    config: dictionary of dictionary of neighborhoods indexed by [reference] then by [0,1],
        where 0 stands for one daughter, and 1 for the other
    cluster_distance
    change_contact_surfaces
    innersurfaces

    Returns
    -------
    conddist: the squareform vector built from the distance matrice
       (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html)
    z: the hierarchical clustering encoded as a linkage matrix
       (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)
    labels: the list of atlas names

    """

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
            dist[i][j] = 0.0
            if distance == 'signature':
                dist[i][j] = 100.0 * division_signature(config[r][0], config[r][1], config[s][0], config[s][1],
                                                        change_contact_surfaces=change_contact_surfaces,
                                                        innersurfaces=innersurfaces)
            else:
                dist[i][j] = 100.0 * division_distance(config[r][0], config[r][1], config[s][0], config[s][1],
                                                       change_contact_surfaces=change_contact_surfaces,
                                                       innersurfaces=innersurfaces)
            dist[j][i] = dist[i][j]

    conddist = sp.spatial.distance.squareform(dist)
    z = sch.linkage(conddist, method=cluster_distance)

    return conddist, z, labels


def switched_division_neighborhoods(config, mother_name):
    #
    # copy neighborhoods from atlases and add switched neighborhoods
    #
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


def _diagnosis_linkage(atlases, parameters):
    proc = "_diagnosis_linkage"
    divisions = atlases.get_divisions()
    ccs = not parameters.use_common_neighborhood

    ref_atlases = atlases.get_atlases()
    output_selections = atlases.get_output_selections()

    merge_values = {}
    lastmerge_values = {}

    swmerge_values = {}
    swlastmerge_values = {}

    division_lastmerge_values = {}

    innersurfaces = []

    for n in divisions:
        stage = n.split('.')[0][1:]
        if len(divisions[n]) <= 2:
            continue

        d = uname.get_daughter_names(n)
        if parameters.exclude_inner_surfaces:
            innersurfaces = [d[0], d[1]]
        #
        # config is a dictionary indexed by [reference][0|1]
        # -> [reference][i] gives the neighborhood of daughter #i of division/mother n
        # swconfig contains the same neighborhoods than config plus the "switched" neighborhoods
        #
        config = atlases.extract_division_neighborhoods(n, delay_from_division=parameters.name_delay_from_division)
        swconfig = switched_division_neighborhoods(config, n)

        #
        # distance array for couples of atlases/references
        #
        conddist, z, labels = division_scipy_linkage(config, change_contact_surfaces=ccs, innersurfaces=innersurfaces)

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
            lineage = ref_atlases[r].cell_lineage
            name = ref_atlases[r].cell_name
            cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
            for c in cells:
                if c not in name:
                    continue
                if name[c] == n:
                    output_selections[keyselection][c] = round(lastmerge_value)

        #
        # distance array for couples of atlases/references plus the switched ones
        #
        swconddist, swz, swlabels = division_scipy_linkage(swconfig, change_contact_surfaces=ccs,
                                                          innersurfaces=innersurfaces)

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
            lineage = ref_atlases[r].cell_lineage
            name = ref_atlases[r].cell_name
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


def _dpp_global_generic_distance(mother, neighbors, references, parameters, debug=False):
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

    innersurfaces = []
    d = uname.get_daughter_names(mother)
    if parameters.exclude_inner_surfaces:
        innersurfaces = [d[0], d[1]]

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
            dist = division_distance(neighbors[0][r1], neighbors[1][r1], neighbors[0][r2], neighbors[1][r2],
                                     change_contact_surfaces=ccs, innersurfaces=innersurfaces)
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
    neighborhoods = atlases.get_cell_neighborhood(delay_from_division=parameters.name_delay_from_division)
    neighbors = {0: copy.deepcopy(neighborhoods[daughters[0]]), 1: copy.deepcopy(neighborhoods[daughters[1]])}

    # score before any changes
    debug = False
    if debug:
        print("")
        print("===== test division " + str(mother) + " : " + str(divisions[mother]))
    score = _dpp_global_generic_distance(mother, neighbors, divisions[mother], parameters, debug=debug)

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
            newscore[r] = _dpp_global_generic_distance(mother, tmp, divisions[mother], parameters, debug=debug)
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

            lineage = ref_atlases[ref].cell_lineage
            name = ref_atlases[ref].cell_name
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


def _division_distance(daughter00, daughter01, daughter10, daughter11, innersurfaces=[], debug=False):
    """
    Compute distance between two contact surface vectors. Do not compute 'common' neighborhood.
    Parameters
    ----------
    daughter00
    daughter01
    daughter10
    daughter11
    innersurfaces

    Returns
    -------

    """

    #
    # get the sum of difference, as well as the sums of contact surfaces
    # cells in 'innersurfaces' are excluded
    #
    nm0, n00, n10 = uneighborhood.cell_distance_elements(daughter00, daughter10, innersurfaces=innersurfaces)
    nm1, n01, n11 = uneighborhood.cell_distance_elements(daughter01, daughter11, innersurfaces=innersurfaces)
    score = (nm0 + nm1) / (n00 + n10 + n01 + n11)
    return score


def division_distance(daughter00, daughter01, daughter10, daughter11, change_contact_surfaces=True, innersurfaces=[],
                      debug=False):
    """

    Parameters
    ----------
    daughter00: dictionary depicting the neighborhood of daughter #0 of ref #0.
        Each key is a named neighbor, and the associated dictionary value give the contact surface.
    daughter01: dictionary depicting the neighborhood of daughter daughter #1 of ref #0
    daughter10: dictionary depicting the neighborhood of daughter daughter #0 of ref #1
    daughter11: dictionary depicting the neighborhood of daughter daughter #1 of ref #1
    change_contact_surfaces: True or False
    innersurfaces
    debug

    Returns
    -------

    """

    if change_contact_surfaces:
        tmp = {'foo': {0: daughter00, 1: daughter10}}
        v0 = uneighborhood.build_same_contact_surfaces(tmp, ['foo'], debug=debug)
        tmp = {'foo': {0: daughter01, 1: daughter11}}
        v1 = uneighborhood.build_same_contact_surfaces(tmp, ['foo'], debug=debug)
        score = _division_distance(v0['foo'][0], v1['foo'][0], v0['foo'][1], v1['foo'][1], innersurfaces=innersurfaces)
    else:
        score = _division_distance(daughter00, daughter01, daughter10, daughter11, innersurfaces=innersurfaces)

    return score


def _division_signature(daughter00, daughter01, daughter10, daughter11, innersurfaces=[]):
    neighbors = set(daughter00.keys()).union(set(daughter01.keys()), set(daughter10.keys()), set(daughter11.keys()))
    den = 0.0
    num = 0.0
    for k in neighbors:
        if k in innersurfaces:
            continue
        if k in daughter00 and k in daughter10:
            num += abs(daughter00[k] - daughter10[k])
        elif k in daughter00 and k not in daughter10:
            num += abs(daughter00[k])
        elif k not in daughter00 and k in daughter10:
            num += abs(daughter10[k])
        if k in daughter01 and k in daughter11:
            num += abs(daughter01[k] - daughter11[k])
        elif k in daughter01 and k not in daughter11:
            num += abs(daughter01[k])
        elif k not in daughter01 and k in daughter11:
            num += abs(daughter11[k])
        mnum0 = 0.0
        if k in daughter00:
            mnum0 += daughter00[k]
        if k in daughter01:
            mnum0 += daughter01[k]
        mnum1 = 0.0
        if k in daughter10:
            mnum1 += daughter10[k]
        if k in daughter11:
            mnum1 += daughter11[k]
        num -= abs(mnum0 - mnum1)
        den += mnum0 + mnum1
    score = num / den
    if score < 0.0:
        return 0.0
    return score


def division_signature(daughter00, daughter01, daughter10, daughter11, change_contact_surfaces=True, innersurfaces=[],
                       debug=False):
    """

    Parameters
    ----------
    daughter00: dictionary depicting the neighborhood of daughter #0 of ref #0.
        Each key is a named neighbor, and the associated dictionary value give the contact surface.
    daughter01: dictionary depicting the neighborhood of daughter daughter #1 of ref #0
    daughter10: dictionary depicting the neighborhood of daughter daughter #0 of ref #1
    daughter11: dictionary depicting the neighborhood of daughter daughter #1 of ref #1
    change_contact_surfaces: True or False
    innersurfaces
    debug

    Returns
    -------

    """

    if change_contact_surfaces:
        tmp = {'foo': {0: daughter00, 1: daughter10}}
        v0 = uneighborhood.build_same_contact_surfaces(tmp, ['foo'], debug=debug)
        tmp = {'foo': {0: daughter01, 1: daughter11}}
        v1 = uneighborhood.build_same_contact_surfaces(tmp, ['foo'], debug=debug)
        score = _division_signature(v0['foo'][0], v1['foo'][0], v0['foo'][1], v1['foo'][1], innersurfaces=innersurfaces)
    else:
        score = _division_signature(daughter00, daughter01, daughter10, daughter11, innersurfaces=innersurfaces)

    return score

###########################################################
#
#
#
############################################################

class DivisionAtlases(uatlasc.CellAtlases):
    def __init__(self, parameters=None):

        uatlasc.CellAtlases.__init__(self, parameters)

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

    ############################################################
    #
    # getters / setters
    #
    ############################################################

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
    #
    #
    ############################################################

    def _build_divisions(self, delay_from_division=None):
        """
        Build a dictionary index by mother cell name. Each entry contains reference names
        for which the both daughters exist
        Returns
        -------

        """
        proc = "_build_divisions"

        if self._divisions is not None:
            del self._divisions
            self._divisions = {}

        #
        # get all references per division/mother cells
        #
        neighborhoods = self.get_cell_neighborhood(delay_from_division=delay_from_division)
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

    #
    #
    #

    def _build_unequal_divisions(self, delay_from_division=None):
        """
        Build a dictionary index by mother cell name. Each entry contains reference names
        for which the both daughters exist
        Returns
        -------

        """

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

    #
    #
    #

    def build_division_atlases(self, parameters):

        #
        # cell based part
        # extract neighborhoods for required delays
        #
        delays = [parameters.name_delay_from_division]
        if parameters.confidence_delay_from_division is not None and \
                parameters.confidence_delay_from_division not in delays:
            delays += [parameters.confidence_delay_from_division]
        self.set_default_delay(parameters.name_delay_from_division)
        for d in delays:
            self.build_cell_atlases(parameters, delay_from_division=d)

        #
        # build common neighborhood reference if required
        #
        delays = self.get_cell_neighborhood_delays()
        if parameters.use_common_neighborhood:
            monitoring.to_log_and_console("... build common neighborhoods", 1)
            for d in delays:
                neighborhoods = self.get_cell_neighborhood(delay_from_division=d)
                self.set_cell_neighborhood(_build_common_neighborhoods(neighborhoods), delay_from_division=d)
            self._use_common_neighborhood = True
            monitoring.to_log_and_console("    done", 1)

        #
        # dictionary indexed by mother cell names,
        # give list of references for which both daughter cells exist
        #
        monitoring.to_log_and_console("... build division list", 1)
        self._build_divisions(delay_from_division=parameters.name_delay_from_division)
        monitoring.to_log_and_console("    done", 1)

        #
        # dictionary indexed by mother cell names,
        # give list of daughter cells, the largest one being the first one
        #
        monitoring.to_log_and_console("... build unequal division list", 1)
        for d in delays:
            self._build_unequal_divisions(delay_from_division=d)
        monitoring.to_log_and_console("    done", 1)

        if parameters.division_diagnosis:
            monitoring.to_log_and_console("")
            monitoring.to_log_and_console("============================================================")
            monitoring.to_log_and_console("===== diagnosis: atlases pairwise disagreements")
            _diagnosis_pairwise_switches(self, parameters)
            monitoring.to_log_and_console("===== diagnosis: dendrogram/linkage diagnosis")
            _diagnosis_linkage(self, parameters)
            monitoring.to_log_and_console("============================================================")
            monitoring.to_log_and_console("")

        #
        # look for daughter that may improve a global score
        # report it in the console/log file
        # as well as in morphonet selection file
        #
        if parameters.division_permutation_proposal:
            division_permutation_proposal(self, parameters)

    ############################################################
    #
    #
    #
    ############################################################

    def extract_division_neighborhoods(self, mother_name, delay_from_division=None):
        """

        Parameters
        ----------
        mother_name
        delay_from_division

        Returns
        -------
        A dictionary of cell neighborhoods indexed by [reference][0|1] where 0|1
        correspond to daughter names issued from uname.get_daughter_names()

        """

        delay = delay_from_division
        if delay is None:
            delay = self.get_default_delay()

        divisions = self.get_divisions()
        neighborhoods = self.get_cell_neighborhood(delay_from_division=delay)
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

    def generate_figure(self, parameters, time_digits_for_cell_id=4):

        uatlasc.CellAtlases.generate_figure(self, parameters, time_digits_for_cell_id=time_digits_for_cell_id)

        generate_figure = (isinstance(parameters.generate_figure, bool) and parameters.generate_figure) or \
                          (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'all') or \
                          (isinstance(parameters.generate_figure, list) and 'all' in parameters.generate_figure)

        #
        # draw histograms of both right pairing and wrong pairing
        # 2D histograms are at division level
        # 1D histograms are at cell (daughter) level
        #
        if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'distance-histograms') \
            or (isinstance(parameters.generate_figure, str)
                and parameters.generate_figure == 'cell-distance-histograms') \
            or (isinstance(parameters.generate_figure, list) and 'distance-histograms' in parameters.generate_figure) \
            or (isinstance(parameters.generate_figure, list)
                and 'cell-distance-histograms' in parameters.generate_figure) or generate_figure:
            monitoring.to_log_and_console("... generate cell distance histogram file", 1)
            _figures_cell_distance_histogram(self, parameters)
            monitoring.to_log_and_console("... done", 1)

        if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'distance-histograms') \
            or (isinstance(parameters.generate_figure, str)
                and parameters.generate_figure == 'division-distance-histograms') \
            or (isinstance(parameters.generate_figure, list) and 'distance-histograms' in parameters.generate_figure) \
            or (isinstance(parameters.generate_figure, list)
                and 'division-distance-histograms' in parameters.generate_figure) or generate_figure:
            monitoring.to_log_and_console("... generate division distance histogram file", 1)
            _figures_division_distance_histogram(self, parameters)
            monitoring.to_log_and_console("... done", 1)

        #
        # draw a dendrogram per division where atlases are grouped with similarity between division
        #
        if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'division-dendrograms') \
                or (isinstance(parameters.generate_figure, list)
                    and 'division-dendrograms' in parameters.generate_figure) \
                or generate_figure:
            monitoring.to_log_and_console("... generate division dendrogram figure file", 1)
            _figures_division_dendrogram(self, parameters)
            monitoring.to_log_and_console("...", 1)
            _figures_signature_dendrogram(self, parameters)
            monitoring.to_log_and_console("... done", 1)


################################################################################
#
# cell-to-cell and division-to-division distance histograms
#
################################################################################

def _write_raw_array(f, a, length=4):
    form = "{:1." + str(length) + "f}"
    last = len(a) - 1
    f.write("[")
    for i, v in enumerate(a):
        f.write(form.format(v))
        if i < last:
            f.write(", ")
    f.write("]")


def _write_array(f, name, a, length=4):
    f.write(str(name) + " = ")
    _write_raw_array(f, a, length=length)
    f.write("\n")


def _write_dict_of_arrays(f, name, a, length=4):
    last = len(a) - 1
    f.write(str(name) + " = {")
    for i, v in enumerate(a):
        f.write(str(v) + ": ")
        _write_raw_array(f, a[v], length=length)
        if i < last:
            f.write(", ")
    f.write("}\n")


def _figures_cell_distance_histogram(atlases, parameters):
    """
    Computes cell-to-cell and division-to-division distance histograms
    Parameters
    ----------
    atlases
    parameters

    Returns
    -------

    """
    proc = "_figures_cell_distance_histogram"

    filename = 'figures_cell_distance_histogram'
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

    compute_other_scores = True

    #
    # get the references per mother_name
    #
    divisions = atlases.get_divisions()
    ccs = not atlases.get_use_common_neighborhood()
    neighborhoods = atlases.get_cell_neighborhood()

    #
    # make generation-dependant calculation
    #
    division_per_generation = {}
    for n in divisions:
        generation = n.split('.')[0][1:]
        division_per_generation[generation] = division_per_generation.get(generation, []) + [n]
    for g in division_per_generation:
        print("    - generation " + str(g) + ": " + str(len(division_per_generation[g])) + " divisions")

    ndivision = 0
    for g in division_per_generation:
        # if int(g) > generationmax:
        #    continue
        ndivision += len(division_per_generation[g])

    right_cscores_per_generation = {}
    wrong_cscores_per_generation = {}

    other_scores_per_generation = {}

    #
    # compute cell-to-cell distances for
    # - similar cells (cell of same name across embryos)
    # - sister cells (only across embryos)
    # and division-to-division distances
    #
    for g in division_per_generation:
        for n in division_per_generation[g]:
            d = uname.get_daughter_names(n)
            #
            # if int(n.split('.')[0][1:]) > 6:
            #     continue
            for r1 in divisions[n]:
                for r2 in divisions[n]:
                    if r2 <= r1:
                        continue
                    d00 = uneighborhood.cell_distance(neighborhoods[d[0]][r1], neighborhoods[d[0]][r2],
                                                      change_contact_surfaces=ccs)
                    d11 = uneighborhood.cell_distance(neighborhoods[d[1]][r1], neighborhoods[d[1]][r2],
                                                      change_contact_surfaces=ccs)
                    d01 = uneighborhood.cell_distance(neighborhoods[d[0]][r1], neighborhoods[d[1]][r2],
                                                      change_contact_surfaces=True)
                    d10 = uneighborhood.cell_distance(neighborhoods[d[1]][r1], neighborhoods[d[0]][r2],
                                                      change_contact_surfaces=True)
                    right_cscores_per_generation[g] = right_cscores_per_generation.get(g, []) + [d00, d11]
                    wrong_cscores_per_generation[g] = wrong_cscores_per_generation.get(g, []) + [d01, d10]
    #
    # compute cell-to-cell distance for cells of the same generation
    # but not the same cell nor its sister
    #
    if compute_other_scores:
        i = 0
        for n in divisions:
            d = uname.get_daughter_names(n)
            generation = n.split('.')[0][1:]
            # if int(generation) > generationmax:
            #     continue
            if i % 10 == 0:
                print("      " + str(i) + "/" + str(ndivision))
            other_scores = []
            for m in division_per_generation[generation]:
                if m <= n:
                    continue
                f = uname.get_daughter_names(m)
                for r1 in divisions[n]:
                    for r2 in divisions[m]:
                        if r2 == r1:
                            continue
                        d00 = uneighborhood.cell_distance(neighborhoods[d[0]][r1], neighborhoods[f[0]][r2],
                                                          change_contact_surfaces=True)
                        d11 = uneighborhood.cell_distance(neighborhoods[d[1]][r1], neighborhoods[f[1]][r2],
                                                          change_contact_surfaces=True)
                        d01 = uneighborhood.cell_distance(neighborhoods[d[0]][r1], neighborhoods[f[1]][r2],
                                                          change_contact_surfaces=True)
                        d10 = uneighborhood.cell_distance(neighborhoods[d[1]][r1], neighborhoods[f[0]][r2],
                                                          change_contact_surfaces=True)
                        other_scores += [d00, d11, d01, d10]
            i += 1
            other_scores_per_generation[generation] = other_scores_per_generation.get(generation, []) + other_scores

    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")
    f.write("import scipy.stats as stats\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    _write_dict_of_arrays(f, "right_cscores_per_generation", right_cscores_per_generation, length=4)
    f.write("right_cscores = []\n")
    f.write("for g in right_cscores_per_generation:\n")
    f.write("    right_cscores += right_cscores_per_generation[g]\n")

    f.write("\n")
    _write_dict_of_arrays(f, "wrong_cscores_per_generation", wrong_cscores_per_generation, length=4)
    f.write("wrong_cscores = []\n")
    f.write("for g in wrong_cscores_per_generation:\n")
    f.write("    wrong_cscores += wrong_cscores_per_generation[g]\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(7.5, 7.5))\n")
    f.write("labels = ['same cell', 'sister cell']\n")
    f.write("ax.hist([right_cscores, wrong_cscores], 100, histtype='bar', label=labels)\n")
    f.write("ax.legend(prop={'size': 10})\n")
    f.write("ax.set_title('atlas-to-atlas cell distances', fontsize=12)\n")
    f.write("ax.tick_params(labelsize=10)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('histogram1D_cell")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    if compute_other_scores:
        f.write("\n")
        _write_dict_of_arrays(f, "other_scores_per_generation", other_scores_per_generation, length=4)
        f.write("other_scores = []\n")
        f.write("for g in other_scores_per_generation:\n")
        f.write("    other_scores += other_scores_per_generation[g]\n")

        f.write("\n")
        f.write("fig, ax = plt.subplots(figsize=(7.5, 7.5))\n")
        f.write("labels = ['same cell', 'sister cell', 'other cell']\n")
        f.write("ax.hist([right_cscores, wrong_cscores, other_scores], 100, histtype='bar', label=labels, density=True)\n")
        f.write("ax.legend(prop={'size': 10})\n")
        f.write("ax.set_title('atlas-to-atlas cell distances', fontsize=12)\n")
        f.write("ax.tick_params(labelsize=10)\n")

        f.write("\n")
        f.write("if savefig:\n")
        f.write("    plt.savefig('density1D_cell")
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '.png')\n")
        f.write("else:\n")
        f.write("    plt.show()\n")
        f.write("    plt.close()\n")

    f.write("\n")

    f.close()


def _figures_division_distance_histogram(atlases, parameters):
    """
    Computes cell-to-cell and division-to-division distance histograms
    Parameters
    ----------
    atlases
    parameters

    Returns
    -------

    """
    proc = "_figures_division_distance_histogram"

    filename = 'figures_division_distance_histogram'
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
    # get the references per mother_name
    #
    divisions = atlases.get_divisions()
    ccs = not atlases.get_use_common_neighborhood()
    neighborhoods = atlases.get_cell_neighborhood()

    #
    # make generation-dependant calculation
    #
    division_per_generation = {}
    for n in divisions:
        generation = n.split('.')[0][1:]
        division_per_generation[generation] = division_per_generation.get(generation, []) + [n]
    for g in division_per_generation:
        division_per_generation[g] = sorted(division_per_generation[g])
        print("    - generation " + str(g) + ": " + str(len(division_per_generation[g])) + " divisions")

    ndivision = 0
    for g in division_per_generation:
        ndivision += len(division_per_generation[g])

    right_dscores_per_generation = {}
    wrong_dscores_per_generation = {}
    diff_dscores_per_generation = {}

    #
    # compute cell-to-cell distances for
    # - similar cells (cell of same name across embryos)
    # - sister cells (only across embryos)
    # and division-to-division distances
    #
    innersurfaces = []
    for g in division_per_generation:
        for n in division_per_generation[g]:
            d = uname.get_daughter_names(n)
            if parameters.exclude_inner_surfaces:
                innersurfaces = [d[0], d[1]]
            #
            # if int(n.split('.')[0][1:]) > 6:
            #     continue
            for r1 in divisions[n]:
                for r2 in divisions[n]:
                    if r2 <= r1:
                        continue

                    div00 = division_distance(neighborhoods[d[0]][r1], neighborhoods[d[1]][r1], neighborhoods[d[0]][r2],
                                              neighborhoods[d[1]][r2], change_contact_surfaces=ccs,
                                              innersurfaces=innersurfaces)

                    # daughter neighborhood have the same neighbors if atlases.get_use_common_neighborhood() is True
                    # there is then no need to change the contact surfaces
                    div01 = division_distance(neighborhoods[d[0]][r1], neighborhoods[d[1]][r1], neighborhoods[d[1]][r2],
                                              neighborhoods[d[0]][r2], change_contact_surfaces=ccs,
                                              innersurfaces=innersurfaces)

                    # trace = n == "b7.0003_" or n == "b7.0008*"
                    # if trace:
                    #     print("division distance of " + n + " between " + r1 + " and " + r2 + ":")
                    #     print("\t right pairing = " + str(div00))
                    #     print("\t wrong pairing = " + str(div01))
                    right_dscores_per_generation[g] = right_dscores_per_generation.get(g, []) + [div00]
                    wrong_dscores_per_generation[g] = wrong_dscores_per_generation.get(g, []) + [div01]
                    diff_dscores_per_generation[g] = diff_dscores_per_generation.get(g, []) + [div01-div00]

    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")
    f.write("import scipy.stats as stats\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    _write_dict_of_arrays(f, "right_dscores_per_generation", right_dscores_per_generation, length=4)
    f.write("right_dscores = []\n")
    f.write("for g in right_dscores_per_generation:\n")
    f.write("    right_dscores += right_dscores_per_generation[g]\n")

    f.write("\n")
    _write_dict_of_arrays(f, "wrong_dscores_per_generation", wrong_dscores_per_generation, length=4)
    f.write("wrong_dscores = []\n")
    f.write("for g in wrong_dscores_per_generation:\n")
    f.write("    wrong_dscores += wrong_dscores_per_generation[g]\n")

    f.write("\n")
    f.write("gens = set(right_dscores_per_generation.keys())\n")
    f.write("gens.intersection(set(wrong_dscores_per_generation.keys()))\n")
    f.write("gens = sorted(list(gens))\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(7.5, 7.5))\n")
    f.write("labels = ['right pairing', 'wrong pairing']\n")
    f.write("ax.hist([right_dscores, wrong_dscores], 100, histtype='bar', label=labels)\n")
    f.write("ax.legend(prop={'size': 10})\n")
    f.write("ax.set_title('atlas-to-atlas division distances', fontsize=12)\n")
    f.write("ax.tick_params(labelsize=10)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('histogram1D_division")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    f.write("fig, axs = plt.subplots(2, 2, figsize=(7.5, 7.5), sharex=True, sharey=True)\n")
    f.write("labels = ['right pairing', 'wrong pairing']\n")
    f.write("if len(gens) >= 1:\n")
    f.write("    axs[0, 0].hist([right_dscores_per_generation[gens[0]], wrong_dscores_per_generation[gens[0]]]")
    f.write(", 100, histtype='bar', label=labels)\n")
    f.write("    axs[0, 0].legend(prop={'size': 10})\n")
    f.write("    title = 'generation = ' + str(gens[0])\n")
    f.write("    axs[0, 0].set_title(title, fontsize=10)\n")
    f.write("    axs[0, 0].tick_params(labelsize=10)\n")
    f.write("if len(gens) >= 2:\n")
    f.write("    axs[0, 1].hist([right_dscores_per_generation[gens[1]], wrong_dscores_per_generation[gens[1]]]")
    f.write(", 100, histtype='bar', label=labels)\n")
    f.write("    axs[0, 1].legend(prop={'size': 10})\n")
    f.write("    title = 'generation = ' + str(gens[1])\n")
    f.write("    axs[0, 1].set_title(title, fontsize=10)\n")
    f.write("    axs[0, 1].tick_params(labelsize=10)\n")
    f.write("if len(gens) >= 3:\n")
    f.write("    axs[1, 0].hist([right_dscores_per_generation[gens[2]], wrong_dscores_per_generation[gens[2]]]")
    f.write(", 100, histtype='bar', label=labels)\n")
    f.write("    axs[1, 0].legend(prop={'size': 10})\n")
    f.write("    title = 'generation = ' + str(gens[2])\n")
    f.write("    axs[1, 0].set_title(title, fontsize=10)\n")
    f.write("    axs[1, 0].tick_params(labelsize=10)\n")
    f.write("if len(gens) >= 4:\n")
    f.write("    axs[1, 1].hist([right_dscores_per_generation[gens[3]], wrong_dscores_per_generation[gens[3]]]")
    f.write(", 100, histtype='bar', label=labels)\n")
    f.write("    axs[1, 1].legend(prop={'size': 10})\n")
    f.write("    title = 'generation = ' + str(gens[3])\n")
    f.write("    axs[1, 1].set_title(title, fontsize=10)\n")
    f.write("    axs[1, 1].tick_params(labelsize=10)\n")
    f.write("fig.suptitle('atlas-to-atlas division distances', fontsize=12)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('histogram1D_division_generation")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    _write_dict_of_arrays(f, "diff_dscores_per_generation", diff_dscores_per_generation, length=4)
    f.write("diff_dscores = []\n")
    f.write("for g in diff_dscores_per_generation:\n")
    f.write("    diff_dscores += diff_dscores_per_generation[g]\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(7.5, 7.5))\n")
    f.write("ax.hist(diff_dscores, 100, histtype='bar')\n")
    f.write("ax.set_title('atlas-to-atlas division distance difference', fontsize=12)\n")
    f.write("ax.tick_params(labelsize=10)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('histogram1D_division_difference")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")

    f.close()


################################################################################
#
# cell neighbors number wrt total number of cells in the embryo
#
################################################################################

def _linkage_balance(z, nlabels):
    """

    Parameters
    ----------
    z: the hierarchical clustering encoded as a linkage matrix
       (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)
    nlabels

    Returns
    -------
    balance: min( #elements of the two last clusters) / max( #elements of the two last clusters)
        1.0 means the two last clusters are balanced
    """
    #
    # z[-1, 0] and z[-1, 1] are the labels of the two last clusters to be fused
    # if label < nlabels, it means the label is the one of an single element (one original observation)
    # else the cluster was made/described at line 'round(z[-1, 0]) - nlabels' of z, and the 3rd
    # element gives the number of original observations of this cluster
    #
    if z[-1, 0] < nlabels:
        ncluster0 = 1
    else:
        ncluster0 = z[round(z[-1, 0]) - nlabels, 3]
    if z[-1, 1] < nlabels:
        ncluster1 = 1
    else:
        ncluster1 = z[round(z[-1, 1]) - nlabels, 3]
    balance = min(ncluster0, ncluster1) / max(ncluster0, ncluster1)
    return balance


def _figures_division_dendrogram(atlases, parameters):
    """
    Parameters
    ----------
    atlases: nested dictionary of neighborhood, where the keys are ['cell name']['reference name']
        where 'reference name' is the name of the reference lineage, and neighborhood a dictionary of contact surfaces
        indexed by cell names (only for the first time point after the division)
    parameters

    Returns
    -------

    """
    proc = "figures_division_dendrogram"

    filename = 'figures_division_dendrogram'
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
    # get the references per mother_name
    #
    divisions = atlases.get_divisions()
    ccs = not parameters.use_common_neighborhood
    cluster_distance = parameters.dendrogram_cluster_distance

    #
    #
    #

    dendro_values = {}
    merge_values = {}

    swmerge_values = {}
    swlastmerge_values = {}

    f = open(filename, "w")

    f.write("import sys\n")
    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")
    f.write("import scipy.cluster.hierarchy as sch\n")

    f.write("\n")
    f.write("cluster_distance = '" + str(cluster_distance) + "'\n")
    f.write("\n")

    cellidentifierlist = []
    innersurfaces = []

    for n in divisions:
        stage = n.split('.')[0][1:]
        if len(divisions[n]) <= 2:
            continue

        #
        # config is a dictionary indexed par [reference][0 or 1]
        # config[r][0] = neighborhoods[daughters[0]][r]
        # config[r][1] = neighborhoods[daughters[1]][r]
        #
        config = atlases.extract_division_neighborhoods(n)
        #
        # swconfig = config + switched daughters
        #
        swconfig = switched_division_neighborhoods(config, n)

        #
        # distance array for couples of atlases/references
        #
        d = uname.get_daughter_names(n)
        if parameters.exclude_inner_surfaces:
            innersurfaces = [d[0], d[1]]
        conddist, z, labels = division_scipy_linkage(config, cluster_distance=cluster_distance,
                                                    change_contact_surfaces=ccs, innersurfaces=innersurfaces)

        merge_values[stage] = merge_values.get(stage, []) + list(z[:, 2])
        # lastmerge_value is the distance between the two last clusters
        lastmerge_value = z[:, 2][-1]
        # balance in [0, 1] reflects the balance between the two last cluster. 1.0 means the two last clusters
        # have the same number of original observations
        balance = _linkage_balance(z, len(labels))
        dendro_values[stage] = dendro_values.get(stage, []) + [(lastmerge_value, balance, len(labels))]

        #
        # distance array for couples of atlases/references plus the switched ones
        # daughter neighborhoods may be in a common reference, if so there is no need to change
        # the contact surfaces
        #
        swconddist, swz, swlabels = division_scipy_linkage(swconfig, cluster_distance=cluster_distance,
                                                          change_contact_surfaces=ccs, innersurfaces=innersurfaces)

        swmerge_values[stage] = swmerge_values.get(stage, []) + list(swz[:, 2])
        swlastmerge_value = swz[:, 2][-1]
        swlastmerge_values[stage] = swlastmerge_values.get(stage, []) + [swlastmerge_value]

        #
        # identifier for mother cell
        #
        cellname = n.split('.')[0] + "_" + n.split('.')[1][0:4]
        if n.split('.')[1][4] == '_':
            cellname += 'U'
        elif n.split('.')[1][4] == '*':
            cellname += 'S'
        else:
            cellname += 'S'
        fileidentifier = 'HC{:03d}_'.format(int(lastmerge_value)) + cellname
        cellidentifier = cellname + '_HC{:03d}'.format(int(lastmerge_value))
        cellidentifier += '_BAL{:03d}'.format(round(100.0 * balance))
        cellidentifierlist.append((cellidentifier, n))

        f.write("\n")
        f.write("savefig = True\n")

        f.write("\n")
        f.write("\n")
        f.write("def draw_" + cellidentifier + "(savefig=False):\n")
        f.write("    " + "\n")
        f.write("    " + "cdist = " + str(list(conddist)) + "\n")
        f.write("    " + "labels = " + str(labels) + "\n")
        f.write("    " + "\n")
        f.write("    " + "cswdist = " + str(list(swconddist)) + "\n")
        f.write("    " + "swlabels = " + str(swlabels) + "\n")
        f.write("\n")
        f.write("    " + "title = '" + str(n) + " (linkage=' + cluster_distance + '), ")
        f.write("delay={:d}'\n".format(parameters.name_delay_from_division))
        f.write("\n")
        f.write("    " + "Z = sch.linkage(cdist, method=cluster_distance)\n")
        f.write("    " + "fig = plt.figure(figsize=(16, 8))\n")
        f.write("    " + "dn = sch.dendrogram(Z, labels=labels, orientation='right')\n")
        f.write("    " + "plt.title(title, fontsize=24)\n")
        f.write("    " + "plt.xticks(fontsize=16)\n")
        f.write("    " + "plt.yticks(fontsize=14)\n")
        f.write("    " + "plt.xlim([0, 100])\n")

        f.write("    if savefig:\n")
        f.write("        plt.savefig('" + str(fileidentifier) + "_' + cluster_distance +'")
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '.png')\n")
        f.write("        plt.savefig('" + str(cellidentifier) + "_' + cluster_distance +'")
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '.png')\n")
        f.write("    else:\n")
        f.write("        plt.show()\n")
        f.write("    plt.close()\n")
        f.write("\n")

        f.write("    " + "Z = sch.linkage(cswdist, method=cluster_distance)\n")
        f.write("    " + "fig = plt.figure(figsize=(18, 8))\n")
        f.write("    " + "dn = sch.dendrogram(Z, labels=swlabels, orientation='right')\n")
        f.write("    " + "plt.title(title, fontsize=24)\n")
        f.write("    " + "plt.xticks(fontsize=16)\n")
        f.write("    " + "plt.yticks(fontsize=12)\n")
        f.write("    " + "plt.xlim([0, 100])\n")

        f.write("    if savefig:\n")
        f.write("        plt.savefig('" + str(fileidentifier) + "_' + cluster_distance +'")
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '_SW.png')\n")
        f.write("        plt.savefig('" + str(cellidentifier) + "_' + cluster_distance +'")
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '_SW.png')\n")
        f.write("    else:\n")
        f.write("        plt.show()\n")
        f.write("    plt.close()\n")

        f.write("\n")

    f.write("\n")
    f.write("\n")
    f.write("def draw_all(celllist=[], savefig=True):\n")
    for cellid in cellidentifierlist:
        f.write("    if celllist == [] or '" + cellid[1] + "' in celllist:\n")
        f.write("        draw_" + cellid[0] + "(savefig=savefig)\n")

    f.write("\n")
    f.write("\n")

    f.write("celllist = []\n")
    f.write("if len(sys.argv) > 1 and sys.argv[1][:2] == '-h':\n")
    f.write("    print('Usage: ' + str(sys.argv[0]) + ' list of cell names (between quotes)')\n")
    f.write("if len(sys.argv) > 1:\n")
    f.write("    for i in range(1, len(sys.argv)):\n")
    f.write("        celllist += [sys.argv[i]]\n")
    f.write("\n")
    f.write("# print(\"celllist=\"+str(celllist))\n")
    f.write("draw_all(celllist=celllist, savefig=savefig)\n")
    f.write("\n")

    #
    # Statistics from the dendrograms
    #

    generations = list(merge_values.keys())
    generations = sorted(generations)
    f.write("generations = " + str(generations) + "\n")
    f.write("merge_values = [")
    for i, g in enumerate(generations):
        f.write(str(list(merge_values[g])))
        if i < len(generations) - 1:
            f.write(", ")
    f.write("]\n")
    f.write("\n")

    f.write("dendro_values = [")
    for i, g in enumerate(generations):
        f.write(str(list(dendro_values[g])))
        if i < len(generations) - 1:
            f.write(", ")
    f.write("]\n")
    f.write("\n")

    f.write("merge_labels = [")
    for i, g in enumerate(generations):
        f.write("'generation " + str(g) + "'")
        if i < len(generations) - 1:
            f.write(", ")
    f.write("]\n")

    f.write("\n")
    f.write("lastmerge_values = [[v[0] for v in dv] for dv in dendro_values]\n")
    f.write("balance_values = [[v[1] for v in dv] for dv in dendro_values]\n")

    f.write("\n")
    f.write("fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(16, 6.5))\n")

    f.write("ax1.hist(merge_values, bins=list(range(0, 101, 2)), histtype='bar', stacked=True, label=merge_labels)\n")
    f.write("ax1.set_title('dendrogram merge values', fontsize=12)\n")
    f.write("ax1.legend(prop={'size': 10})\n")
    f.write("ax1.tick_params(labelsize=10)\n")

    f.write("ax2.hist(lastmerge_values, bins=list(range(0, 101, 2)), histtype='bar', stacked=True, label=merge_labels)\n")
    f.write("ax2.set_title('dendrogram last merge values', fontsize=12)\n")
    f.write("ax2.legend(prop={'size': 10})\n")
    f.write("ax2.tick_params(labelsize=10)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('dendrogram_merge_histogram_' + cluster_distance +'")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(8, 8))\n")

    f.write("ax.hist(balance_values, bins=50, range=(0, 1), histtype='bar', stacked=True, label=merge_labels)\n")
    f.write("ax.set_title('dendrogram balance values', fontsize=12)\n")
    f.write("ax.legend(prop={'size': 10})\n")
    f.write("ax.tick_params(labelsize=10)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('dendrogram_balance_histogram_' + cluster_distance +'")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    f.write("x = [lv for lvg in lastmerge_values for lv in lvg]\n")
    f.write("y = [bv for bvg in balance_values for bv in bvg]\n")
    f.write("c_values = [[int(generations[i])] * len(lvg) for i, lvg in enumerate(lastmerge_values)]\n")
    f.write("c = [cv for cvg in c_values for cv in cvg]\n")
    f.write("s_values = [[v[2] for v in dv] for dv in dendro_values]\n")
    f.write("s = [10 * sv for svg in s_values for sv in svg]\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(8, 8))\n")
    f.write("scatter = ax.scatter(x, y, c=c, s=s, alpha=0.5)\n")
    f.write("ax.set_xlabel('Last merge value (cluster distance * 100)')\n")
    f.write("ax.set_ylabel('Linkage balance')\n")
    f.write("legend1 = ax.legend(*scatter.legend_elements(alpha=0.8), loc='upper left', title='Generations')\n")
    f.write("ax.add_artist(legend1)\n")
    f.write("handles, labels = scatter.legend_elements(prop='sizes', alpha=0.6)\n")
    f.write("labels = [l.replace('0}', '}') for l in labels]\n")
    f.write("legend2 = ax.legend(handles, labels, loc='upper right', title='Sizes')\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('dendrogram_balance_merge_scatter_' + cluster_distance +'")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")

    f.close()


def _figures_signature_dendrogram(atlases, parameters):
    """
    Parameters
    ----------
    atlases: nested dictionary of neighborhood, where the keys are ['cell name']['reference name']
        where 'reference name' is the name of the reference lineage, and neighborhood a dictionary of contact surfaces
        indexed by cell names (only for the first time point after the division)
    parameters

    Returns
    -------

    """
    proc = "figures_signature_dendrogram"

    filename = 'figures_signature_dendrogram'
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
    # get the references per mother_name
    #
    divisions = atlases.get_divisions()
    ccs = not parameters.use_common_neighborhood
    cluster_distance = parameters.dendrogram_cluster_distance

    #
    #
    #

    dendro_values = {}
    merge_values = {}

    f = open(filename, "w")

    f.write("import sys\n")
    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")
    f.write("import scipy.cluster.hierarchy as sch\n")

    f.write("\n")
    f.write("cluster_distance = '" + str(cluster_distance) + "'\n")
    f.write("\n")

    cellidentifierlist = []
    innersurfaces = []

    for n in divisions:
        stage = n.split('.')[0][1:]
        if len(divisions[n]) <= 2:
            continue

        #
        # config is a dictionary indexed par [reference][0 or 1]
        # config[r][0] = neighborhoods[daughters[0]][r]
        # config[r][1] = neighborhoods[daughters[1]][r]
        #
        config = atlases.extract_division_neighborhoods(n)

        #
        # distance array for couples of atlases/references
        #
        d = uname.get_daughter_names(n)
        if parameters.exclude_inner_surfaces:
            innersurfaces = [d[0], d[1]]
        conddist, z, labels = division_scipy_linkage(config, cluster_distance=cluster_distance,
                                                     change_contact_surfaces=ccs, innersurfaces=innersurfaces,
                                                     distance='signature')

        merge_values[stage] = merge_values.get(stage, []) + list(z[:, 2])
        # lastmerge_value is the distance between the two last clusters
        lastmerge_value = z[:, 2][-1]
        # balance in [0, 1] reflects the balance between the two last cluster. 1.0 means the two last clusters
        # have the same number of original observations
        balance = _linkage_balance(z, len(labels))
        dendro_values[stage] = dendro_values.get(stage, []) + [(lastmerge_value, balance, len(labels))]

        #
        # identifier for mother cell
        #
        cellname = n.split('.')[0] + "_" + n.split('.')[1][0:4]
        if n.split('.')[1][4] == '_':
            cellname += 'U'
        elif n.split('.')[1][4] == '*':
            cellname += 'S'
        else:
            cellname += 'S'
        fileidentifier = 'HC{:03d}_'.format(int(lastmerge_value)) + cellname
        cellidentifier = cellname + '_HC{:03d}'.format(int(lastmerge_value))
        cellidentifier += '_BAL{:03d}'.format(round(100.0 * balance))
        cellidentifierlist.append((cellidentifier, n))

        f.write("\n")
        f.write("savefig = True\n")

        f.write("\n")
        f.write("\n")
        f.write("def draw_" + cellidentifier + "(savefig=False):\n")
        f.write("    " + "\n")
        f.write("    " + "cdist = " + str(list(conddist)) + "\n")
        f.write("    " + "labels = " + str(labels) + "\n")
        f.write("\n")
        f.write("    " + "title = '" + str(n) + " (linkage=' + cluster_distance + '), ")
        f.write("delay={:d}'\n".format(parameters.name_delay_from_division))
        f.write("\n")
        f.write("    " + "Z = sch.linkage(cdist, method=cluster_distance)\n")
        f.write("    " + "fig = plt.figure(figsize=(16, 8))\n")
        f.write("    " + "dn = sch.dendrogram(Z, labels=labels, orientation='right')\n")
        f.write("    " + "plt.title(title, fontsize=24)\n")
        f.write("    " + "plt.xticks(fontsize=16)\n")
        f.write("    " + "plt.yticks(fontsize=14)\n")
        f.write("    " + "plt.xlim([0, 100])\n")

        f.write("    if savefig:\n")
        f.write("        plt.savefig('SIG_" + str(fileidentifier) + "_' + cluster_distance +'")
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '.png')\n")
        f.write("        plt.savefig('SIG_" + str(cellidentifier) + "_' + cluster_distance +'")
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '.png')\n")
        f.write("    else:\n")
        f.write("        plt.show()\n")
        f.write("    plt.close()\n")
        f.write("\n")

    f.write("\n")
    f.write("\n")
    f.write("def draw_all(celllist=[], savefig=True):\n")
    for cellid in cellidentifierlist:
        f.write("    if celllist == [] or '" + cellid[1] + "' in celllist:\n")
        f.write("        draw_" + cellid[0] + "(savefig=savefig)\n")

    f.write("\n")
    f.write("\n")

    f.write("celllist = []\n")
    f.write("if len(sys.argv) > 1 and sys.argv[1][:2] == '-h':\n")
    f.write("    print('Usage: ' + str(sys.argv[0]) + ' list of cell names (between quotes)')\n")
    f.write("if len(sys.argv) > 1:\n")
    f.write("    for i in range(1, len(sys.argv)):\n")
    f.write("        celllist += [sys.argv[i]]\n")
    f.write("\n")
    f.write("# print(\"celllist=\"+str(celllist))\n")
    f.write("draw_all(celllist=celllist, savefig=savefig)\n")
    f.write("\n")

    #
    # Statistics from the dendrograms
    #

    generations = list(merge_values.keys())
    generations = sorted(generations)
    f.write("generations = " + str(generations) + "\n")
    f.write("merge_values = [")
    for i, g in enumerate(generations):
        f.write(str(list(merge_values[g])))
        if i < len(generations) - 1:
            f.write(", ")
    f.write("]\n")
    f.write("\n")

    f.write("dendro_values = [")
    for i, g in enumerate(generations):
        f.write(str(list(dendro_values[g])))
        if i < len(generations) - 1:
            f.write(", ")
    f.write("]\n")
    f.write("\n")

    f.write("merge_labels = [")
    for i, g in enumerate(generations):
        f.write("'generation " + str(g) + "'")
        if i < len(generations) - 1:
            f.write(", ")
    f.write("]\n")

    f.write("\n")
    f.write("lastmerge_values = [[v[0] for v in dv] for dv in dendro_values]\n")
    f.write("balance_values = [[v[1] for v in dv] for dv in dendro_values]\n")

    f.write("\n")
    f.write("fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(16, 6.5))\n")

    f.write("ax1.hist(merge_values, bins=list(range(0, 101, 2)), histtype='bar', stacked=True, label=merge_labels)\n")
    f.write("ax1.set_title('dendrogram merge values', fontsize=12)\n")
    f.write("ax1.legend(prop={'size': 10})\n")
    f.write("ax1.tick_params(labelsize=10)\n")

    f.write("ax2.hist(lastmerge_values, bins=list(range(0, 101, 2)), histtype='bar', stacked=True, label=merge_labels)\n")
    f.write("ax2.set_title('dendrogram last merge values', fontsize=12)\n")
    f.write("ax2.legend(prop={'size': 10})\n")
    f.write("ax2.tick_params(labelsize=10)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('SIG_dendrogram_merge_histogram_' + cluster_distance +'")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(8, 8))\n")

    f.write("ax.hist(balance_values, bins=50, range=(0, 1), histtype='bar', stacked=True, label=merge_labels)\n")
    f.write("ax.set_title('dendrogram balance values', fontsize=12)\n")
    f.write("ax.legend(prop={'size': 10})\n")
    f.write("ax.tick_params(labelsize=10)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('SIG_dendrogram_balance_histogram_' + cluster_distance +'")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    f.write("x = [lv for lvg in lastmerge_values for lv in lvg]\n")
    f.write("y = [bv for bvg in balance_values for bv in bvg]\n")
    f.write("c_values = [[int(generations[i])] * len(lvg) for i, lvg in enumerate(lastmerge_values)]\n")
    f.write("c = [cv for cvg in c_values for cv in cvg]\n")
    f.write("s_values = [[v[2] for v in dv] for dv in dendro_values]\n")
    f.write("s = [10 * sv for svg in s_values for sv in svg]\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(8, 8))\n")
    f.write("scatter = ax.scatter(x, y, c=c, s=s, alpha=0.5)\n")
    f.write("ax.set_xlabel('Last merge value (cluster distance * 100)')\n")
    f.write("ax.set_ylabel('Linkage balance')\n")
    f.write("legend1 = ax.legend(*scatter.legend_elements(alpha=0.8), loc='upper left', title='Generations')\n")
    f.write("ax.add_artist(legend1)\n")
    f.write("handles, labels = scatter.legend_elements(prop='sizes', alpha=0.6)\n")
    f.write("labels = [l.replace('0}', '}') for l in labels]\n")
    f.write("legend2 = ax.legend(handles, labels, loc='upper right', title='Sizes')\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('SIG_dendrogram_balance_merge_scatter_' + cluster_distance +'")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")

    f.close()
