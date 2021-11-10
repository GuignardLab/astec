
import sys
import copy
import numpy as np
import scipy as sp
import scipy.cluster.hierarchy as sch

import astec.utils.common as common
import astec.utils.ascidian_name as uname
import astec.utils.contact_atlas as ucontacta

monitoring = common.Monitoring()


def figures_division_dendrogram(atlases, parameters, linkage_method='single'):
    """
    Parameters
    ----------
    atlases: nested dictionary of neighborhood, where the keys are ['cell name']['reference name']
        where 'reference name' is the name of the reference lineage, and neighborhood a dictionary of contact surfaces
        indexed by cell names (only for the first time point after the division)
    parameters
    linkage_method: 'single', ’complete’, ’average’, ’weighted’, ’centroid’,  ’median’, ’ward’

    Returns
    -------

    """

    filename = 'figures_division_dendrogram'
    file_suffix = None
    if parameters.figurefile_suffix is not None and isinstance(parameters.figurefile_suffix, str) and \
            len(parameters.figurefile_suffix) > 0:
        file_suffix = '_' + parameters.figurefile_suffix
    if file_suffix is not None:
        filename += file_suffix
    filename += '.py'

    #
    # get the references per mother_name
    #
    divisions = atlases.get_divisions()
    ccs = not parameters.use_common_neighborhood
    neighborhoods = atlases.get_neighborhoods()
    similarity = atlases.get_division_contact_similarity()

    #
    #
    #

    merge_values = {}
    lastmerge_values = {}

    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")
    f.write("import scipy.cluster.hierarchy as sch\n")

    f.write("\n")
    f.write("linkage_method = '" + str(linkage_method) + "'\n")
    f.write("\n")

    cellidentifierlist = []

    for n in divisions:
        stage = n.split('.')[0][1:]
        if len(divisions[n]) <= 2:
            continue
        daughters = uname.get_daughter_names(n)

        if False and (n != 'a7.0002_' and n != 'a8.0003_' and n != 'a8.0004_'):
            continue

        #
        #
        #
        config = {}
        swconfig = {}
        for r in divisions[n]:
            config[r] = {}
            config[r][0] = copy.deepcopy(neighborhoods[daughters[0]][r])
            config[r][1] = copy.deepcopy(neighborhoods[daughters[1]][r])
            swconfig[r] = {}
            swconfig[r][0] = copy.deepcopy(config[r][0])
            swconfig[r][1] = copy.deepcopy(config[r][1])
            sr = 'switched-' + str(r)
            swconfig[sr] = {}
            swconfig[sr][0] = copy.deepcopy(neighborhoods[daughters[1]][r])
            swconfig[sr][1] = copy.deepcopy(neighborhoods[daughters[0]][r])
            if daughters[1] in swconfig[sr][0]:
                msg = "  weird, " + str(daughters[1]) + " was found in its neighborhood for reference " + str(r)
                print("      " + msg)
            if daughters[0] in swconfig[sr][0]:
                swconfig[sr][0][daughters[1]] = swconfig[sr][0][daughters[0]]
                del swconfig[sr][0][daughters[0]]
            if daughters[0] in swconfig[sr][1]:
                msg = "  weird, " + str(daughters[0]) + " was found in its neighborhood for reference " + str(r)
                print("      " + msg)
            if daughters[1] in swconfig[sr][1]:
                swconfig[sr][1][daughters[0]] = swconfig[sr][1][daughters[1]]
                del swconfig[sr][1][daughters[1]]

        #
        # distance array for couples of atlases/references
        #
        labels = []
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
                dist[i][j] = 100.0 * ucontacta.division_contact_generic_distance(atlases, config[r][0], config[r][1],
                                                                                 config[s][0], config[s][1],
                                                                                 similarity=similarity,
                                                                                 change_contact_surfaces=ccs)
                dist[j][i] = dist[i][j]

        conddist = sp.spatial.distance.squareform(dist)
        z = sch.linkage(conddist, method=linkage_method)

        merge_values[stage] = merge_values.get(stage, []) + list(z[:, 2])
        lastmerge_value = z[:, 2][-1]
        lastmerge_values[stage] = lastmerge_values.get(stage, []) + [lastmerge_value]

        swdist = np.zeros((len(swconfig), len(swconfig)))
        swlabels = []
        for i, r in enumerate(swconfig):
            swlabels += [r]
            for j, s in enumerate(swconfig):
                if r == s:
                    swdist[i][i] = 0.0
                    continue
                if r > s:
                    continue
                # if r == 'switched-' + str(s) or s == 'switched-' + str(r):
                #    continue
                swdist[i][j] = 100.0 * ucontacta.division_contact_generic_distance(atlases, swconfig[r][0],
                                                                                   swconfig[r][1], swconfig[s][0],
                                                                                   swconfig[s][1],
                                                                                   similarity=similarity,
                                                                                   change_contact_surfaces=ccs)
                swdist[j][i] = swdist[i][j]
        swconddist = sp.spatial.distance.squareform(swdist)

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
        # method='single' (default)
        # method='complete'
        # method='average'
        # method='weighted'
        # method='centroid'
        # method='median'
        # method='ward'
        f.write("\n")
        f.write("    " + "title = '" + str(n) + " (linkage=' + linkage_method + '), ")
        f.write("delay={:d}'\n".format(parameters.delay_from_division))
        f.write("\n")
        f.write("    " + "Z = sch.linkage(cdist, method=linkage_method)\n")
        f.write("    " + "fig = plt.figure(figsize=(16, 8))\n")
        f.write("    " + "dn = sch.dendrogram(Z, labels=labels, orientation='right')\n")
        f.write("    " + "plt.title(title, fontsize=24)\n")
        f.write("    " + "plt.xticks(fontsize=16)\n")
        f.write("    " + "plt.yticks(fontsize=14)\n")
        f.write("    " + "plt.xlim([0, 100])\n")

        f.write("    if savefig:\n")
        f.write("        plt.savefig('" + str(fileidentifier) + "_' + linkage_method +'")
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '.png')\n")
        f.write("        plt.savefig('" + str(cellidentifier) + "_' + linkage_method +'")
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '.png')\n")
        f.write("    else:\n")
        f.write("        plt.show()\n")
        f.write("    plt.close()\n")
        f.write("\n")

        f.write("    " + "Z = sch.linkage(cswdist, method=linkage_method)\n")
        f.write("    " + "fig = plt.figure(figsize=(18, 8))\n")
        f.write("    " + "dn = sch.dendrogram(Z, labels=swlabels, orientation='right')\n")
        f.write("    " + "plt.title(title, fontsize=24)\n")
        f.write("    " + "plt.xticks(fontsize=16)\n")
        f.write("    " + "plt.yticks(fontsize=12)\n")
        f.write("    " + "plt.xlim([0, 100])\n")

        f.write("    if savefig:\n")
        f.write("        plt.savefig('" + str(fileidentifier) + "_' + linkage_method +'")
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '_SW.png')\n")
        f.write("        plt.savefig('" + str(cellidentifier) + "_' + linkage_method +'")
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
    f.write("merge_values = [")
    for i, g in enumerate(generations):
        f.write(str(list(merge_values[g])))
        if i < len(generations) - 1:
            f.write(", ")
    f.write("]\n")
    f.write("lastmerge_values = [")
    for i, g in enumerate(generations):
        f.write(str(list(lastmerge_values[g])))
        if i < len(generations) - 1:
            f.write(", ")
    f.write("]\n")
    f.write("merge_labels = [")
    for i, g in enumerate(generations):
        f.write("'" + str(g) + "'")
        if i < len(generations) - 1:
            f.write(", ")
    f.write("]\n")

    f.write("\n")
    f.write("fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(16, 6.5))\n")

    f.write("ax1.hist(merge_values, bins=list(range(101)), histtype='bar', stacked=True, label=merge_labels)\n")
    f.write("ax1.set_title('dendrogram merge values', fontsize=12)\n")
    f.write("ax1.legend(prop={'size': 10})\n")
    f.write("ax1.tick_params(labelsize=10)\n")

    f.write("ax2.hist(lastmerge_values, bins=list(range(101)), histtype='bar', stacked=True, label=merge_labels)\n")
    f.write("ax2.set_title('dendrogram last merge values', fontsize=12)\n")
    f.write("ax2.legend(prop={'size': 10})\n")
    f.write("ax2.tick_params(labelsize=10)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('dendrogram_merge_histogram_' + linkage_method +'")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")

    f.close()


def figures_division_graph(atlases, parameters):
    """
    Build a file to draw graphs for each division where a node is a reference and
    Parameters
    ----------
    atlases: nested dictionary of neighborhood, where the keys are ['cell name']['reference name']
        where 'reference name' is the name of the reference lineage, and neighborhood a dictionary of contact surfaces
        indexed by cell names (only for the first time point after the division)
    parameters

    Returns
    -------

    """

    proc = "figures_division_graph"

    filename = 'figures_division_graph'
    file_suffix = None

    if not isinstance(parameters, ucontacta.AtlasParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: " +
                                      str(type(parameters)))
        sys.exit(1)

    if parameters.figurefile_suffix is not None and isinstance(parameters.figurefile_suffix, str) and \
            len(parameters.figurefile_suffix) > 0:
        file_suffix = '_' + parameters.figurefile_suffix
    if file_suffix is not None:
        filename += file_suffix
    filename += '.py'

    #
    # get the references per mother_name
    #
    divisions = atlases.get_divisions()
    ccs = not parameters.use_common_neighborhood

    #
    # edgeValues[division][ref1][ref2]: set an edge value in [0, 1]
    # 1 means perfect equality, while 0 means total difference
    # values[division]: all edge values for division
    #
    mother_names = sorted(divisions.keys())
    edgeValues = {}
    values = {}
    for n in mother_names:
        edgeValues[n] = {}
        for r1 in divisions[n]:
            edgeValues[n][r1] = {}
            for r2 in divisions[n]:
                if r2 <= r1:
                    continue
                edgeValues[n][r1][r2] = atlases.get_division_similarity(n, r1, r2, change_contact_surfaces=ccs)
            if edgeValues[n][r1] != {}:
                values[n] = values.get(n, []) + list(edgeValues[n][r1].values())

    #
    # values[n] is a list of all the score values for any couple of atlases exhibiting the division of 'n'
    # max_values[n] is the worst value for a comparison of two atlases exhibiting the division of 'n'
    #
    min_values = {}
    for n in values.keys():
        min_values[n] = min(values[n])

    f = open(filename, "w")

    f.write("import networkx as nx\n")
    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")

    sublist = ['a7.0002_', 'a8.0003_', 'a8.0004_']

    cellidentifierlist = []
    for mother_name in min_values:

        if False and mother_name not in sublist:
            continue

        min_value = 100.0 * min_values[mother_name]

        # identifier for mother cell
        cellname = mother_name.split('.')[0] + "_" + mother_name.split('.')[1][0:4]
        if mother_name.split('.')[1][4] == '_':
            cellname += 'U'
        elif mother_name.split('.')[1][4] == '*':
            cellname += 'S'
        else:
            cellname += 'S'
        fileidentifier = 'GRAPH{:03d}_'.format(int(min_value)) + cellname
        cellidentifier = cellname + '_GRAPH{:03d}'.format(int(min_value))
        cellidentifierlist.append(cellidentifier)

        #
        # edges[n1][n2] in [1,5] is an integer giving the weight of an edge
        # the smallest the weight, the weakest the edge
        #
        edges = {}
        for n1 in divisions[mother_name]:
            if n1 not in edgeValues[mother_name]:
                continue
            edges[n1] = {}
            for n2 in divisions[mother_name]:
                if n2 not in edgeValues[mother_name][n1]:
                    continue
                if n2 <= n1:
                    continue
                edges[n1][n2] = int(1.0 + 4.0 * edgeValues[mother_name][n1][n2])

        f.write("\n")
        f.write("\n")
        f.write("def draw_" + cellidentifier + "(savefig=False):\n")

        nodes = list(divisions[mother_name])
        nodes.sort()

        f.write("\n")
        f.write("    G_" + cellidentifier + " = nx.Graph()\n")
        f.write("    G_" + cellidentifier + ".add_nodes_from(" + str(nodes) + ")\n")
        f.write("\n")

        #
        # for each edge weight, get the edge weight, and use the right edge value for colorization
        #
        edgelist = {}
        edgecolor = {}
        for n1 in edges:
            for n2 in edges[n1]:
                if n2 <= n1:
                    continue
                edgelist[edges[n1][n2]] = edgelist.get(edges[n1][n2], []) + [(n1, n2)]
                edgecolor[edges[n1][n2]] = edgecolor.get(edges[n1][n2], []) + [edgeValues[mother_name][n1][n2]]

        for i in edgelist:
            f.write("    G_" + cellidentifier + ".add_edges_from(" + str(edgelist[i]) + ", weight=" + str(i) + ")\n")

        f.write("\n")
        f.write("    node_labels = {}\n")
        for n in nodes:
            f.write("    node_labels['" + str(n) + "'] = '" + str(n) + "'\n")
        f.write("\n")

        f.write("    fig, ax = plt.subplots(figsize=(7.5, 7.5))\n")
        f.write("    pos = nx.circular_layout(G_" + cellidentifier + ", scale=1.8)\n")
        f.write("    nx.draw_networkx_nodes(G_" + cellidentifier + ", pos=pos, ax=ax, node_size=1000)\n")
        f.write("    nx.draw_networkx_labels(G_" + cellidentifier + ", ax=ax, pos=pos, labels=node_labels," +
                " font_weight='bold')\n")
        f.write("    ax.set_xlim([-2.2, 2.2])\n")
        f.write("    ax.set_ylim([-2.1, 2.1])\n")
        for i in edgelist:
            f.write("    nx.draw_networkx_edges(G_" + cellidentifier + ", pos=pos, ax=ax, edgelist=" + str(edgelist[i]))
            f.write(", edge_color=" + str(edgecolor[i]))
            f.write(", edge_vmin=0.0, edge_vmax=1.0")
            # '_r' at the end of the colormap name means reversed colormap
            f.write(", edge_cmap = plt.cm.brg_r")
            f.write(", width=" + str(i))
            f.write(")\n")
        title = "division of " + str(mother_name)
        title += ", edge value in [{:.2f}, {:.2f}]".format(min_values[mother_name], max(values[mother_name]))

        f.write("    ax.set_title(\"" + str(title) + "\")\n")

        f.write("    if savefig:\n")
        f.write("        plt.savefig('" + str(fileidentifier))
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '.png')\n")
        f.write("        plt.savefig('" + str(cellidentifier))
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '.png')\n")
        f.write("    else:\n")
        f.write("        plt.show()\n")
        f.write("    plt.close()\n")

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
    f.write("    draw_all(savefig=True)\n")
    f.write("\n")

    f.close()


def figures_distance_histogram(atlases, parameters):

    filename = 'figures_distance_histogram'
    file_suffix = None
    if parameters.figurefile_suffix is not None and isinstance(parameters.figurefile_suffix, str) and \
            len(parameters.figurefile_suffix) > 0:
        file_suffix = '_' + parameters.figurefile_suffix
    if file_suffix is not None:
        filename += file_suffix
    filename += '.py'

    #
    # get the references per mother_name
    #
    divisions = atlases.get_divisions()
    ccs = not atlases.get_use_common_neighborhood()
    neighborhoods = atlases.get_neighborhoods()
    similarity = atlases.get_division_contact_similarity()

    right_cscores = []
    right_sscores = []
    wrong_cscores = []
    wrong_sscores = []

    compute_other_scores = True
    other_scores = []

    right_dscores = []
    wrong_dscores = []

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
                if d01 > 2.0 or d10 > 2.0:
                    print("mother = " + str(n))
                    print("  d[0,1] = " + str(d[0]) + ", " + str(d[1]))
                    print("  r1, r2 = " + str(r1) + ", " + str(r2))
                    print("  d01, d10 = " + str(d01) + ", " + str(d10))
                right_cscores += [d00, d11]
                right_sscores += [d11, d00]
                wrong_cscores += [d01, d10]
                wrong_sscores += [d10, d01]

                div00 = ucontacta.division_contact_generic_distance(atlases, neighborhoods[d[0]][r1],
                                                                    neighborhoods[d[1]][r1], neighborhoods[d[0]][r2],
                                                                    neighborhoods[d[1]][r2],
                                                                    similarity=similarity, change_contact_surfaces=ccs)
                div01 = ucontacta.division_contact_generic_distance(atlases, neighborhoods[d[0]][r1],
                                                                    neighborhoods[d[1]][r1], neighborhoods[d[1]][r2],
                                                                    neighborhoods[d[0]][r2],
                                                                    similarity=similarity, change_contact_surfaces=ccs)
                right_dscores += [div00]
                wrong_dscores += [div01]

    if compute_other_scores:
        for n in divisions:
            d = uname.get_daughter_names(n)
            generation = n.split('.')[0][1:]
            if int(generation) > 7:
                continue
            for m in divisions:
                if m <= n:
                    continue
                if m.split('.')[0][1:] != generation:
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

    step = atlases.get_probability_step()

    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")
    f.write("import scipy.stats as stats\n")

    f.write("\n")
    f.write("savefig = True\n")

    f.write("\n")
    f.write("right_cscores = " + str(right_cscores) + "\n")
    f.write("right_sscores = " + str(right_sscores) + "\n")
    f.write("\n")
    f.write("values = np.array([right_cscores, right_sscores])\n")
    f.write("kernel = stats.gaussian_kde(values)\n")
    f.write("\n")
    f.write("step = " + str(step) + "\n")
    f.write("X, Y = np.mgrid[0:1:step, 0:1:step]\n")
    f.write("positions = np.vstack([X.ravel(), Y.ravel()])\n")
    f.write("Z = np.reshape(kernel(positions).T, X.shape)\n")
    f.write("\n")
    f.write("scale = 100.0 / Z.sum()\n")
    f.write("N = np.reshape([[(Z[Z <= Z[j][i]].sum()) * scale for i in range(Z.shape[1])] for j in range(Z.shape[0])], Z.shape)\n")

    f.write("\n")
    f.write("fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(16, 6.5))\n")
    f.write("ax1.plot(right_cscores, right_sscores, 'k.', markersize=2)\n")
    f.write("ax1.set_box_aspect(1)\n")
    f.write("ax1.set_xlim(0, 1)\n")
    f.write("ax1.set_xlabel('Daughter #1/daughter #1 distance', fontsize=15)\n")
    f.write("ax1.set_ylabel('Daughter #2/daughter #2 distance', fontsize=15)\n")
    f.write("title = 'distance couples'\n")
    f.write("ax1.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("ax2.imshow(np.rot90(Z), cmap=plt.cm.rainbow_r, extent=[0, 1, 0, 1], vmin=Z.min(), vmax=Z.max())\n")
    f.write("ax2.set_box_aspect(1)\n")
    f.write("ax2.set_xlabel('Daughter #1/daughter #1 distance', fontsize=15)\n")
    f.write("title = 'kernel estimation in [{:.2f}, {:.2f}]'.format(Z.min(), Z.max())\n")
    f.write("ax2.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("ax3.imshow(np.rot90(N), cmap=plt.cm.rainbow_r, extent=[0, 1, 0, 1])\n")
    f.write("ax3.set_box_aspect(1)\n")
    f.write("CS = ax3.contour(X, Y, N, levels=np.arange(10, 100, 10))\n")
    f.write("ax3.clabel(CS, inline=True, fontsize=10)\n")
    f.write("ax3.set_xlabel('Daughter #1/daughter #1 distance', fontsize=15)\n")
    f.write("title = 'distance couple probability'\n")
    f.write("ax3.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('histogram2D_rightpairing")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

    f.write("\n")
    f.write("wrong_cscores = " + str(wrong_cscores) + "\n")
    f.write("wrong_sscores = " + str(wrong_sscores) + "\n")
    f.write("\n")
    f.write("values = np.array([wrong_cscores, wrong_sscores])\n")
    f.write("kernel = stats.gaussian_kde(values)\n")
    f.write("\n")
    f.write("step = " + str(step) + "\n")
    f.write("X, Y = np.mgrid[0:1:step, 0:1:step]\n")
    f.write("positions = np.vstack([X.ravel(), Y.ravel()])\n")
    f.write("Z = np.reshape(kernel(positions).T, X.shape)\n")
    f.write("\n")
    f.write("scale = 100.0 / Z.sum()\n")
    f.write("N = [[(Z[Z <= Z[j][i]].sum()) * scale for i in range(Z.shape[1])] for j in range(Z.shape[0])]\n")

    f.write("\n")
    f.write("fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(16, 6.5))\n")
    f.write("ax1.plot(wrong_cscores, wrong_sscores, 'k.', markersize=2)\n")
    f.write("ax1.set_box_aspect(1)\n")
    f.write("ax1.set_xlim(0, 1)\n")
    f.write("ax1.set_xlabel('Daughter #1/daughter #2 distance', fontsize=15)\n")
    f.write("ax1.set_ylabel('Daughter #2/daughter #1 distance', fontsize=15)\n")
    f.write("title = 'score couples'\n")
    f.write("ax1.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("ax2.imshow(np.rot90(Z), cmap=plt.cm.rainbow_r, extent=[0, 1, 0, 1], vmin=Z.min(), vmax=Z.max())\n")
    f.write("ax2.set_box_aspect(1)\n")
    f.write("ax2.set_xlabel('Daughter #1/daughter #2 distance', fontsize=15)\n")
    f.write("title = 'kernel estimation in [{:.2f}, {:.2f}]'.format(Z.min(), Z.max())\n")
    f.write("ax2.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("ax3.imshow(np.rot90(N), cmap=plt.cm.rainbow_r, extent=[0, 1, 0, 1])\n")
    f.write("ax3.set_box_aspect(1)\n")
    f.write("CS = ax3.contour(X, Y, N, levels=np.arange(10, 100, 10))\n")
    f.write("ax3.clabel(CS, inline=True, fontsize=10)\n")
    f.write("ax3.set_xlabel('Daughter #1/daughter #2 distance', fontsize=15)\n")
    f.write("title = 'score probability'\n")
    f.write("ax3.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('histogram2D_wrongpairing")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("    plt.close()\n")

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
        f.write("other_scores = " + str(other_scores) + "\n")
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
    f.write("right_dscores = " + str(right_dscores) + "\n")
    f.write("wrong_dscores = " + str(wrong_dscores) + "\n")
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

    f.close()
