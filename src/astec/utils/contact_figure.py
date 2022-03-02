
import sys
import copy
import os
from collections import Counter

import astec.utils.common as common
import astec.utils.ascidian_name as uname
import astec.utils.contact as ucontact
import astec.utils.contact_atlas as ucontacta
import astec.utils.ioproperties as ioproperties

monitoring = common.Monitoring()


def _write_array(f, name, a, length=4):
    form = "{:1." + str(length) + "f}"
    last = len(a) - 1
    f.write(str(name) + " = [")
    for i, v in enumerate(a):
        f.write(form.format(v))
        if i < last:
            f.write(", ")
    f.write("]\n")


################################################################################
#
# cell-to-cell distance between successive cells in a branch with respect to distance from first cell
# (from first time point/division to last time point/division)
#
################################################################################

def figures_distance_along_branch(atlases, parameters, time_digits_for_cell_id=4):
    """

    Parameters
    ----------
    atlases
    parameters
    time_digits_for_cell_id

    Returns
    -------

    """
    proc = "figures_distance_along_branch"

    if not isinstance(parameters, ucontacta.AtlasParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: " +
                                      str(type(parameters)))
        sys.exit(1)

    filename = 'figures_distance_along_branch'
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

    ref_atlases = atlases.get_atlases()
    div = 10 ** time_digits_for_cell_id

    contact_distance_along_time = {}

    for ref in ref_atlases:

        contact_distance_along_time[ref] = {}
        lineage = ref_atlases[ref]['cell_lineage']
        contact = ref_atlases[ref]['cell_contact_surface']
        name = ref_atlases[ref]['cell_name']
        reverse_lineage = {v: k for k, values in lineage.items() for v in values}

        #
        # get the first and last time points
        #
        cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
        cells = sorted(cells)
        cells_per_time = {}
        for c in cells:
            t = int(c) // div
            #
            # get cells and cell ids at each time point
            #
            cells_per_time[t] = cells_per_time.get(t, []) + [c]
        last_time = max(cells_per_time.keys())

        #
        # study single branches beginning right after a division, but not at the last time
        #
        first_cells = [lineage[c][0] for c in lineage if len(lineage[c]) == 2 and int(c) // div < last_time - 1]
        first_cells += [lineage[c][1] for c in lineage if len(lineage[c]) == 2 and int(c) // div < last_time - 1]

        for cell in first_cells:
            if cell not in contact:
                msg = "    * weird, cell " + str(cell) + " is not in the 'contact surface' dictionary"
                msg += " of atlas '" + str(ref) + "'"
                monitoring.to_log_and_console(msg, 3)
                continue

            if cell not in name:
                msg = "    * weird, cell " + str(cell) + " is not in the 'name' dictionary"
                msg += " of atlas '" + str(ref) + "'"
                monitoring.to_log_and_console(msg, 3)
                keyd = cell
            else:
                keyd = name[cell]

            first_time = int(cell) // div
            pcell = cell
            pneigh = copy.deepcopy(contact[cell])
            #
            # extract next neighborhood, change neighbors wrt first cell and compute distance
            #
            emergency_break = False
            while True:
                if pcell not in lineage or len(lineage[pcell]) > 1:
                    break
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
                    # background
                    if int(c) % div == 1 or int(c) % div == 0:
                        nneigh[1] = nneigh.get(1, 0.0) + contrib
                    else:
                        for i in range(t - first_time):
                            if c not in reverse_lineage:
                                msg = "    * weird, cell " + str(c) + " is not in the reversed lineage"
                                msg += " of atlas '" + str(ref) + "'"
                                monitoring.to_log_and_console(msg)
                                emergency_break = True
                                break
                            c = reverse_lineage[c]
                        if emergency_break:
                            break
                        nneigh[c] = nneigh.get(c, 0.0) + contrib
                if emergency_break:
                    break
                #
                #
                #
                d = ucontact.cell_contact_distance(pneigh, nneigh, distance=parameters.cell_contact_distance,
                                                   change_contact_surfaces=False)
                contact_distance_along_time[ref][keyd] = contact_distance_along_time[ref].get(keyd, []) + [d]
                pcell = ncell
                pneigh = copy.deepcopy(nneigh)
    #
    #
    #
    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")
    f.write("import scipy.stats as stats\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    f.write("contact_distance = " + str(contact_distance_along_time) + "\n")
    f.write("lengths = [len(contact_distance[r][c]) for r in contact_distance for c in contact_distance[r]]\n")

    f.write("\n")
    f.write("dict_dist_per_time = {}\n")
    f.write("for r in contact_distance:\n")
    f.write("    for c in contact_distance[r]:\n")
    f.write("        for i, v in enumerate(contact_distance[r][c]):\n")
    f.write("            dict_dist_per_time[i] = dict_dist_per_time.get(i, []) + [contact_distance[r][c][i]]\n")
    f.write("dist_per_time = []\n")
    f.write("for i in range(max(lengths)):\n")
    f.write("    dist_per_time.append(dict_dist_per_time[i])\n")

    f.write("\n")
    f.write("ticks = [x*10 for x in set([x//10 for x in list(range(max(lengths)))])]\n")

    f.write("\n")
    f.write("fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True, figsize=(16, 6.5))\n")

    f.write("\n")
    f.write("ax1.set_ylim(0, 0.7)\n")
    f.write("ax1.set_xlabel('time from division', fontsize=15)\n")
    f.write("ax1.set_ylabel('distance', fontsize=15)\n")
    f.write("ax1.boxplot(dist_per_time)\n")
    f.write("ax1.set_title(\"[t, t+1] distances\", fontsize=15)\n")
    f.write("ax1.set_xticks(ticks)\n")
    f.write("ax1.set_xticklabels(ticks)\n")

    f.write("\n")
    f.write("ax2.set_ylim(0, 0.7)\n")
    f.write("ax2.set_xlim(0, 20.5)\n")
    f.write("ax2.set_xlabel('time from division', fontsize=15)\n")
    f.write("ax2.set_ylabel('distance', fontsize=15)\n")
    f.write("ax2.boxplot(dist_per_time)\n")
    f.write("ax2.set_title(\"[t, t+1] distances (close-up)\", fontsize=15)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('distance_along_branch")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")
    f.close()


################################################################################
#
# cell number wrt time without and with temporal registration
#
################################################################################

def figures_temporal_registration(atlases, parameters, time_digits_for_cell_id=4):
    """
    Plot cell number curves without and with linear temporal alignment
    Parameters
    ----------
    atlases
    parameters
    time_digits_for_cell_id

    Returns
    -------

    """
    proc = "figures_temporal_registration"

    filename = 'figures_temporal_registration'
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

    ref_atlases = atlases.get_atlases()
    atlas_names = list(ref_atlases.keys())

    cells_per_time = {}
    temporal_coefficients = {}
    for n in atlas_names:
        lineage = ref_atlases[n]['cell_lineage']
        cells = list(set(lineage.keys()).union(set([v for values in list(lineage.values()) for v in values])))
        cells = sorted(cells)
        div = 10 ** time_digits_for_cell_id
        cells_per_time[n] = {}
        for c in cells:
            t = int(c) // div
            cells_per_time[n][t] = cells_per_time[n].get(t, 0) + 1
        temporal_coefficients[n] = ref_atlases[n]['temporal_alignment']

    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    f.write("cells_per_time = " + str(cells_per_time) + "\n")
    f.write("temporal_coefficients = " + str(temporal_coefficients) + "\n")

    f.write("\n")
    f.write("fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True, figsize=(16, 8))\n")
    f.write("labels = []\n")
    f.write("for n in cells_per_time:\n")
    f.write("    labels += [n]\n")
    f.write("    x = list(cells_per_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    y = [cells_per_time[n][i] for i in x]\n")
    f.write("    ax1.plot(x, y)\n")
    f.write("ax1.set_title(\"cell number (without alignment)\", fontsize=15)\n")
    f.write("ax1.legend(labels, prop={'size': 10})\n")

    f.write("for n in cells_per_time:\n")
    f.write("    x = list(cells_per_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    t = [temporal_coefficients[n][0] * i + temporal_coefficients[n][1] for i in x]\n")
    f.write("    y = [cells_per_time[n][i] for i in x]\n")
    f.write("    ax2.plot(t, y)\n")
    f.write("ax2.set_title(\"cell number (with alignment)\", fontsize=15)\n")
    f.write("ax2.legend(labels, prop={'size': 10})\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('temporal_alignment_cell_number")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")
    f.close()


################################################################################
#
# cell neighbors number wrt total number of cells in the embryo
#
################################################################################

def _neighbor_histogram(neighbors, prop, filename, threshold=0.05, time_digits_for_cell_id=4):
    proc = "_neighbor_histogram"

    if 'cell_contact_surface' not in prop:
        monitoring.to_log_and_console(str(proc) + ": 'cell_contact_surface' is not in dictionary: ")
        return neighbors
    contact = prop['cell_contact_surface']

    cell_number_threshold = 10
    #
    # nodespertime is a dictionary
    # nodespertime[t] = #nodes at time t
    #
    nodes = list(contact.keys())
    div = 10 ** time_digits_for_cell_id
    times = [n // div for n in nodes]
    nodespertime = Counter(times)

    for n in contact:
        surface = 0
        for k in contact[n]:
            surface += contact[n][k]
        frac = []
        for k in contact[n]:
            if k % div == 0 or k % div == 1:
                continue
            frac += [contact[n][k]/surface]
        #
        # frac: array of surface fraction (without the background)
        #
        l = len([x for x in frac if x >= threshold])
        #
        # print a message for large neighbor number
        #
        if l > cell_number_threshold:
            msg = "cell " + str(n)
            if 'cell_name' in prop:
                if n in prop['cell_name']:
                    msg += " (" + str(prop['cell_name'][n]) + ")"
            msg += " of properties '" + str(filename) + "' has " + str(l) + " neighbors"
            msg += " (above " + str(threshold*100) + "%)"
            monitoring.to_log_and_console("\t " + msg)
        t = n // div
        neighbors[nodespertime[t]] = neighbors.get(nodespertime[t], []) + [l]
    return neighbors


def _figures_neighbor_histogram_period(f):

    f.write("\n")
    f.write("axtwin = ax.secondary_xaxis('top')\n")
    f.write("axtwin.tick_params('x', direction='inout', length=10, width=2)\n")
    f.write("axtwin.set_xticks([76, 110, 250, 332, 624])\n")

    f.write("\n")
    f.write("x = [2, 76, 76, 2]\n")
    f.write("y = [0, 0, 14, 14]\n")
    f.write("plt.fill(x, y, color='grey', alpha=0.4)\n")
    f.write("plt.text(15, 13.5, 'cleavage', ha='left', va='center', fontsize=12, wrap=True)\n")

    f.write("\n")
    f.write("x = [110, 250, 250, 110]\n")
    f.write("y = [0, 0, 14, 14]\n")
    f.write("plt.fill(x, y, color='grey', alpha=0.4)\n")
    f.write("plt.text(180, 13.5, 'gastrula', ha='center', va='center', fontsize=12, wrap=True)\n")

    f.write("\n")
    f.write("x = [332, 624, 624, 332]\n")
    f.write("y = [0, 0, 14, 14]\n")
    f.write("plt.fill(x, y, color='grey', alpha=0.4)\n")
    f.write("plt.text(480, 13.5, 'neurula', ha='center', va='center', fontsize=12, wrap=True)\n")


def figures_neighbor_histogram(atlasfiles, parameters, time_digits_for_cell_id=4):
    """
    Plot cell neighbor number with respect to total cell number
    Parameters
    ----------
    atlasfiles
    parameters
    time_digits_for_cell_id

    Returns
    -------

    """
    proc = "figures_neighbor_histogram"

    neighbors = {}
    if isinstance(atlasfiles, str):
        prop = ioproperties.read_dictionary(atlasfiles, inputpropertiesdict={})
        neighbors = _neighbor_histogram(neighbors, prop, atlasfiles, time_digits_for_cell_id=time_digits_for_cell_id)
        del prop
    elif isinstance(atlasfiles, list):
        for f in atlasfiles:
            prop = ioproperties.read_dictionary(f, inputpropertiesdict={})
            neighbors = _neighbor_histogram(neighbors, prop, f, time_digits_for_cell_id=time_digits_for_cell_id)
            del prop

    #
    # neighbors is a dictionary indexed by the number of cells of an embryo
    # neighbors[64] is an array containing the number of neighbors
    #
    filename = 'figures_neighbor_histogram'
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
    #
    #
    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")
    f.write("import scipy.stats as stats\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    f.write("neighbors = " + str(neighbors) + "\n")
    f.write("ncells = sorted(list(neighbors.keys()))\n")

    f.write("\n")
    f.write("neighbors_per_ncell = []\n")
    f.write("ticks = []\n")
    f.write("xp = []\n")
    f.write("mp = []\n")
    f.write("sp = []\n")
    f.write("for i in ncells:\n")
    f.write("    xp += [i]\n")
    f.write("    mp += [np.mean(neighbors[i])]\n")
    f.write("    sp += [np.std(neighbors[i])]\n")
    f.write("    neighbors_per_ncell.append(neighbors[i])\n")
    f.write("    if i % 20 == 0:\n")
    f.write("        ticks += [str(i)]\n")
    f.write("    else:\n")
    f.write("        ticks += ['']\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(16, 6.5))\n")

    f.write("\n")
    f.write("ax.set_xlabel('number of cells', fontsize=15)\n")
    f.write("ax.set_ylabel('number of neighbors', fontsize=15)\n")
    f.write("ax.boxplot(neighbors_per_ncell, positions=ncells)\n")
    f.write("ax.set_title(\"cell neighbors\", fontsize=15)\n")
    f.write("ax.set_xticklabels(ticks)\n")
    f.write("ax.set_xlim(left=0)\n")

    _figures_neighbor_histogram_period(f)

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('neighbor_histogram_boxplot")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")

    f.write("\n")
    f.write("fig, ax = plt.subplots(figsize=(16, 6.5))\n")

    f.write("\n")
    f.write("ax.set_xlabel('number of cells', fontsize=15)\n")
    f.write("ax.set_ylabel('average number of neighbors (+/- std dev)', fontsize=15)\n")
    f.write("ax.errorbar(xp, mp, yerr=sp, color='red', fmt='-', ecolor='blue')\n")
    f.write("ax.set_title(\"cell neighbors\", fontsize=15)\n")
    f.write("ax.set_xlim(left=0)\n")

    _figures_neighbor_histogram_period(f)

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('neighbor_histogram_average")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")

    f.close()


################################################################################
#
# cell-to-cell and division-to-division distance histograms
#
################################################################################

def figures_distance_histogram(atlases, parameters):
    """
    Computes cell-to-cell and division-to-division distance histograms
    Parameters
    ----------
    atlases
    parameters

    Returns
    -------

    """
    proc = "figures_distance_histogram"

    filename = 'figures_distance_histogram'
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
    neighborhoods = atlases.get_neighborhoods()
    similarity = atlases.get_division_contact_similarity()

    right_cscores = []
    wrong_cscores = []

    other_scores = []

    right_dscores = []
    wrong_dscores = []
    diff_dscores = []

    #
    # compute cell-to-cell distances for
    # - similar cells (cell of same name across embryos)
    # - sister cells (only across embryos)
    # and division-to-division distances
    #
    for n in divisions:
        d = uname.get_daughter_names(n)
        #
        # if int(n.split('.')[0][1:]) > 6:
        #     continue
        for r1 in divisions[n]:
            for r2 in divisions[n]:
                if r2 <= r1:
                    continue
                d00 = atlases.get_cell_distance(d[0], r1, d[0], r2, change_contact_surfaces=ccs)
                d11 = atlases.get_cell_distance(d[1], r1, d[1], r2, change_contact_surfaces=ccs)
                d01 = atlases.get_cell_distance(d[0], r1, d[1], r2, change_contact_surfaces=True)
                d10 = atlases.get_cell_distance(d[1], r1, d[0], r2, change_contact_surfaces=True)
                right_cscores += [d00, d11]
                wrong_cscores += [d01, d10]

                div00 = ucontacta.division_contact_generic_distance(atlases, neighborhoods[d[0]][r1],
                                                                    neighborhoods[d[1]][r1], neighborhoods[d[0]][r2],
                                                                    neighborhoods[d[1]][r2],
                                                                    similarity=similarity, change_contact_surfaces=ccs)
                # daughter neighborhood have the same neighbors if atlases.get_use_common_neighborhood() is True
                # there is then no need to change the contact surfaces
                div01 = ucontacta.division_contact_generic_distance(atlases, neighborhoods[d[0]][r1],
                                                                    neighborhoods[d[1]][r1], neighborhoods[d[1]][r2],
                                                                    neighborhoods[d[0]][r2],
                                                                    similarity=similarity, change_contact_surfaces=ccs)
                # trace = n == "b7.0003_" or n == "b7.0008*"
                # if trace:
                #     print("division distance of " + n + " between " + r1 + " and " + r2 + ":")
                #     print("\t right pairing = " + str(div00))
                #     print("\t wrong pairing = " + str(div01))
                right_dscores += [div00]
                wrong_dscores += [div01]
                diff_dscores += [div01-div00]

    #
    # compute cell-to-cell distance for cells of the same generation
    # but not the same cell nor its sister
    #

    if compute_other_scores:

        # generationmax = 7

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

        i = 0
        for n in divisions:
            d = uname.get_daughter_names(n)
            generation = n.split('.')[0][1:]
            # if int(generation) > generationmax:
            #     continue
            if i % 10 == 0:
                print("      " + str(i) + "/" + str(ndivision))
            for m in division_per_generation[generation]:
                if m <= n:
                    continue
                f = uname.get_daughter_names(m)
                for r1 in divisions[n]:
                    for r2 in divisions[m]:
                        if r2 == r1:
                            continue
                        d00 = atlases.get_cell_distance(d[0], r1, f[0], r2, change_contact_surfaces=True)
                        d11 = atlases.get_cell_distance(d[1], r1, f[1], r2, change_contact_surfaces=True)
                        d01 = atlases.get_cell_distance(d[0], r1, f[1], r2, change_contact_surfaces=True)
                        d10 = atlases.get_cell_distance(d[1], r1, f[0], r2, change_contact_surfaces=True)
                        other_scores += [d00, d11, d01, d10]
            i += 1

    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")
    f.write("import scipy.stats as stats\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    _write_array(f, "right_cscores", right_cscores, length=4)

    f.write("\n")
    _write_array(f, "wrong_cscores", wrong_cscores, length=4)

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
        _write_array(f, "other_scores", other_scores, length=4)
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
    _write_array(f, "right_dscores", right_dscores, length=4)
    _write_array(f, "wrong_dscores", wrong_dscores, length=4)
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
    _write_array(f, "diff_dscores", diff_dscores, length=4)
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


def figures_division_dendrogram(atlases, parameters):
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

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")
    f.write("import scipy.cluster.hierarchy as sch\n")

    f.write("\n")
    f.write("cluster_distance = '" + str(cluster_distance) + "'\n")
    f.write("\n")

    cellidentifierlist = []

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
        swconfig = ucontacta.switched_division_neighborhoods(config, n)

        #
        # distance array for couples of atlases/references
        #
        conddist, z, labels = ucontacta.call_to_scipy_linkage(atlases, config, cluster_distance=cluster_distance,
                                                              change_contact_surfaces=ccs)

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
        swconddist, swz, swlabels = ucontacta.call_to_scipy_linkage(atlases, swconfig,
                                                                    cluster_distance=cluster_distance,
                                                                    change_contact_surfaces=ccs)

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
        cellidentifierlist.append(cellidentifier)

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
    f.write("def draw_all(savefig=False):\n")
    for cellidentifier in cellidentifierlist:
        f.write("    draw_" + cellidentifier + "(savefig=savefig)\n")

    f.write("\n")
    f.write("\n")
    for cellidentifier in cellidentifierlist:
        f.write("if False:\n")
        f.write("    draw_" + cellidentifier + "(savefig=False)\n")
    f.write("\n")
    f.write("if True:\n")
    f.write("    draw_all(savefig=savefig)\n")
    f.write("\n")

    generations = list(merge_values.keys())
    generations = sorted(generations)
    f.write("generations = " + str(generations) + "\n")
    f.write("merge_values = [")
    for i, g in enumerate(generations):
        f.write(str(list(merge_values[g])))
        if i < len(generations) - 1:
            f.write(", ")
    f.write("]\n")
    dendro_values[stage] = dendro_values.get(stage, []) + [(lastmerge_value, balance, len(labels))]
    f.write("dendro_values = [")
    for i, g in enumerate(generations):
        f.write(str(list(dendro_values[g])))
        if i < len(generations) - 1:
            f.write(", ")
    f.write("]\n")
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


################################################################################
#
# cell-to-cell distance between successive cells in a branch with respect to distance from first cell
# (from first time point/division to last time point/division)
#
################################################################################

def figures_embryo_volume(atlases, parameters, time_digits_for_cell_id=4):
    """

    Parameters
    ----------
    atlases
    parameters
    time_digits_for_cell_id

    Returns
    -------

    """
    proc = "figures_embryo_volume"

    if not isinstance(parameters, ucontacta.AtlasParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: " +
                                      str(type(parameters)))
        sys.exit(1)

    filename = 'figures_embryo_volume'
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

    ref_atlases = atlases.get_atlases()
    atlas_names = list(ref_atlases.keys())
    div = 10 ** time_digits_for_cell_id

    volume_along_time = {}
    temporal_coefficients = {}
    for n in atlas_names:
        volume_along_time[n] = {}
        volume = ref_atlases[n]['cell_volume']
        for c in volume:
            t = int(c) // div
            volume_along_time[n][t] = volume_along_time[n].get(t, 0) + volume[c]
        temporal_coefficients[n] = ref_atlases[n]['temporal_alignment']

    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    f.write("volume_along_time = " + str(volume_along_time) + "\n")
    f.write("temporal_coefficients = " + str(temporal_coefficients) + "\n")

    f.write("\n")
    f.write("fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True, figsize=(16, 8))\n")
    f.write("labels = []\n")
    f.write("for n in volume_along_time:\n")
    f.write("    labels += [n]\n")
    f.write("    x = list(volume_along_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    y = [volume_along_time[n][i] for i in x]\n")
    f.write("    ax1.plot(x, y)\n")
    f.write("ax1.set_title(\"embryo volume (without alignment)\", fontsize=15)\n")
    f.write("ax1.legend(labels, prop={'size': 10})\n")

    f.write("for n in volume_along_time:\n")
    f.write("    x = list(volume_along_time[n].keys())\n")
    f.write("    x = sorted(x)\n")
    f.write("    t = [temporal_coefficients[n][0] * i + temporal_coefficients[n][1] for i in x]\n")
    f.write("    y = [volume_along_time[n][i] for i in x]\n")
    f.write("    ax2.plot(t, y)\n")
    f.write("ax2.set_title(\"embryo volume (with alignment)\", fontsize=15)\n")
    f.write("ax2.legend(labels, prop={'size': 10})\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('temporal_alignment_embryo_volume")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")
    f.close()
