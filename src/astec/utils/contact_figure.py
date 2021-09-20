
import astec.utils.ascidian_name as uname
import astec.utils.contact_atlas as ucontacta
import astec.utils.contact as ucontact


def figures_division_score(neighborhoods, parameters):
    """
    Build a file to draw figures of 'consistency' between reference embryos for mother cells
    Parameters
    ----------
    neighborhoods: nested dictionary of neighborhood, where the keys are ['cell name']['reference name']
        where 'reference name' is the name of the reference lineage, and neighborhood a dictionary of contact surfaces
        indexed by cell names (only for the first time point after the division)
    parameters

    Returns
    -------

    """

    filename = 'figures_division_score'
    file_suffix = None
    if parameters.figurefile_suffix is not None and isinstance(parameters.figurefile_suffix, str) and \
            len(parameters.figurefile_suffix) > 0:
        file_suffix = '_' + parameters.figurefile_suffix
    filename += file_suffix + '.py'

    #
    # get the references per mother_name
    #
    cell_names = sorted(list(neighborhoods.keys()))
    references = {}
    for cell_name in cell_names:
        mother_name = uname.get_mother_name(cell_name)
        references[mother_name] = references.get(mother_name, set()).union(set(neighborhoods[cell_name].keys()))

    #
    # get an edge value in [0, 1]
    # 1 means perfect equality, while 0 means total difference
    #
    mother_names = sorted(references.keys())
    edgeValues = {}
    values = {}
    for n in mother_names:
        d = uname.get_daughter_names(n)
        edgeValues[n] = {}
        for n1 in references[n]:
            edgeValues[n][n1] = {}
            for n2 in references[n]:
                if n2 <= n1:
                    continue
                if d[0] in neighborhoods and (n1 in neighborhoods[d[0]] and n2 in neighborhoods[d[0]]):
                    s0 = ucontact.contact_distance(neighborhoods[d[0]][n1], neighborhoods[d[0]][n2],
                                                   similarity=parameters.contact_similarity)
                else:
                    s0 = None
                if d[1] in neighborhoods and (n1 in neighborhoods[d[1]] and n2 in neighborhoods[d[1]]):
                    s1 = ucontact.contact_distance(neighborhoods[d[1]][n1], neighborhoods[d[1]][n2],
                                                   similarity=parameters.contact_similarity)
                else:
                    s1 = None
                #
                # get_score() is in [0, 1]:
                # 0: perfect agreement
                # 1: full disagreement
                #
                # edgeValues will be in [0, 1]
                # 0: weak edge
                # 1: string edge
                #
                if s0 is not None and s1 is not None:
                    edgeValues[n][n1][n2] = 1.0 - 0.5 * (s0 + s1)
                # elif s0 is not None:
                #     edges[n][n1][n2] = s0
                # elif s1 is not None:
                #     edges[n][n1][n2] = s1
            if edgeValues[n][n1] != {}:
                values[n] = values.get(n, []) + list(edgeValues[n][n1].values())
    #
    # values[n] is a list of all the score values for any couple of atlases exhibiting the division of 'n'
    # max_values[n] is the worst value for a comparison of two atlases exhibiting the division of 'n'
    #
    max_values = {}
    for n in values.keys():
        max_values[n] = 1.0 - min(values[n])

    f = open(filename, "w")

    f.write("import networkx as nx\n")
    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")

    cellidentifierlist = []
    for mother_name in max_values:
        max_value = 100.0 * max_values[mother_name]

        # identifier for mother cell
        cellname = mother_name.split('.')[0] + "_" + mother_name.split('.')[1][0:4]
        if mother_name.split('.')[1][4] == '_':
            cellname += 'U'
        elif mother_name.split('.')[1][4] == '*':
            cellname += 'S'
        else:
            cellname += 'S'
        fileidentifier = 'S{:03d}_'.format(int(max_value)) + cellname
        cellidentifier = cellname + '_S{:03d}'.format(int(max_value))
        cellidentifierlist.append(cellidentifier)

        #
        # edges[n1][n2] in [1,5] is an integer giving the weight of an edge
        # the smallest the weight, the weakest the edge
        #
        edges = {}
        for n1 in references[mother_name]:
            if n1 not in edgeValues[mother_name]:
                continue
            edges[n1] = {}
            for n2 in references[mother_name]:
                if n2 not in edgeValues[mother_name][n1]:
                    continue
                if n2 <= n1:
                    continue
                edges[n1][n2] = int(1.0 + 4.0 * edgeValues[mother_name][n1][n2])

        f.write("\n")
        f.write("\n")
        f.write("def draw_" + cellidentifier + "(savefig=False):\n")

        nodes = list(references[mother_name])
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
        title += ", max score value = " + '{:.2f}'.format(max_values[mother_name])

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

    f.write("\n")
    f.write("\n")
    f.write("def draw_all(savefig=False):\n")
    for cellidentifier in cellidentifierlist:
        f.write("    draw_" + cellidentifier + "(savefig=savefig)\n")
        f.write("    plt.close()\n")

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


def figures_division_probability(atlases, parameters):
    """
    Build a file to draw figures of 'probability' between reference embryos for mother cells
    Parameters
    ----------
    atlases: nested dictionary of neighborhood, where the keys are ['cell name']['reference name']
        where 'reference name' is the name of the reference lineage, and neighborhood a dictionary of contact surfaces
        indexed by cell names (only for the first time point after the division)
    parameters

    Returns
    -------

    """

    filename = 'figures_division_probability'
    file_suffix = None
    if parameters.figurefile_suffix is not None and isinstance(parameters.figurefile_suffix, str) and \
            len(parameters.figurefile_suffix) > 0:
        file_suffix = '_' + parameters.figurefile_suffix
    filename += file_suffix + '.py'

    #
    # get the references per mother_name
    #
    neighborhoods = atlases.get_neighborhoods()
    cell_names = sorted(list(neighborhoods.keys()))
    references = {}
    for cell_name in cell_names:
        mother_name = uname.get_mother_name(cell_name)
        references[mother_name] = references.get(mother_name, set()).union(set(neighborhoods[cell_name].keys()))

    #
    # get an edge value in [0, 1]
    # 1 means perfect equality, while 0 means total difference
    #
    mother_names = sorted(references.keys())
    edgeValues = {}
    values = {}
    for n in mother_names:
        d = uname.get_daughter_names(n)
        #
        # check whether each reference has the two daughters
        #
        refs = list(references[n])
        for r in refs:
            if r in neighborhoods[d[0]] and r in neighborhoods[d[1]]:
                continue
            references[n].remove(r)
        if len(references[n]) <= 1:
            continue

        edgeValues[n] = {}
        for n1 in references[n]:
            edgeValues[n][n1] = {}
            for n2 in references[n]:
                if n2 <= n1:
                    continue
                s0 = ucontact.contact_distance(neighborhoods[d[0]][n1], neighborhoods[d[0]][n2],
                                               similarity=parameters.contact_similarity)
                s1 = ucontact.contact_distance(neighborhoods[d[1]][n1], neighborhoods[d[1]][n2],
                                               similarity=parameters.contact_similarity)
                #
                # probability is in [0, 100]
                # edgeValues will be in [0, 1]
                # 0: weak edge
                # 1: strong edge
                #
                edgeValues[n][n1][n2] = atlases.get_probability(s0, s1) / 100.0
            if edgeValues[n][n1] != {}:
                values[n] = values.get(n, []) + list(edgeValues[n][n1].values())
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

    cellidentifierlist = []
    for mother_name in min_values:
        min_value = 100.0 * min_values[mother_name]

        # identifier for mother cell
        cellname = mother_name.split('.')[0] + "_" + mother_name.split('.')[1][0:4]
        if mother_name.split('.')[1][4] == '_':
            cellname += 'U'
        elif mother_name.split('.')[1][4] == '*':
            cellname += 'S'
        else:
            cellname += 'S'
        fileidentifier = 'PROB{:03d}_'.format(int(min_value)) + cellname
        cellidentifier = cellname + '_PROB{:03d}'.format(int(min_value))
        cellidentifierlist.append(cellidentifier)

        #
        # edges[n1][n2] in [1,5] is an integer giving the weight of an edge
        # the smallest the weight, the weakest the edge
        #
        edges = {}
        for n1 in references[mother_name]:
            if n1 not in edgeValues[mother_name]:
                continue
            edges[n1] = {}
            for n2 in references[mother_name]:
                if n2 not in edgeValues[mother_name][n1]:
                    continue
                if n2 <= n1:
                    continue
                edges[n1][n2] = int(1.0 + 4.0 * edgeValues[mother_name][n1][n2])

        f.write("\n")
        f.write("\n")
        f.write("def draw_" + cellidentifier + "(savefig=False):\n")

        nodes = list(references[mother_name])
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
        title += ", min probability value = " + '{:.2f}'.format(min_values[mother_name])

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

    f.write("\n")
    f.write("\n")
    f.write("def draw_all(savefig=False):\n")
    for cellidentifier in cellidentifierlist:
        f.write("    draw_" + cellidentifier + "(savefig=savefig)\n")
        f.write("    plt.close()\n")

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


def figures_cell_neighborhood_pca(neighborhoods, parameters, min_samples=4):

    cells = sorted(neighborhoods, key=lambda key: len(neighborhoods[key]), reverse=True)

    filename = 'figures_cell_neighborhood_pcafigures_cell_neighborhood_pca'
    file_suffix = None
    if parameters.figurefile_suffix is not None and isinstance(parameters.figurefile_suffix, str) and \
            len(parameters.figurefile_suffix) > 0:
        file_suffix = '_' + parameters.figurefile_suffix
    filename += file_suffix + '.py'

    f = open(filename, "w")

    f.write("import matplotlib.pyplot as plt\n")
    f.write("import matplotlib.colors as clr\n")
    f.write("from sklearn.decomposition import PCA\n")

    cellidentifierlist = []
    for cell in cells:
        if len(neighborhoods[cell]) <= min_samples:
            continue
        cellname = cell.split('.')[0] + "_" + cell.split('.')[1][0:4]
        if cell.split('.')[1][4] == '_':
            cellname += 'U'
        elif cell.split('.')[1][4] == '*':
            cellname += 'S'
        else:
            cellname += 'S'
        cellidentifier = cellname + '_NA{:03d}'.format(int(len(neighborhoods[cell])))
        cellidentifierlist.append(cellidentifier)

        atlases = list(neighborhoods[cell].keys())
        neighbors = neighborhoods[cell][atlases[0]]
        arr = []
        for a in atlases:
            vec = []
            norm = 0.0
            for n in neighbors:
                vec += [neighborhoods[cell][a][n]]
                norm += neighborhoods[cell][a][n]
            vec = [i/norm for i in vec]
            arr += [vec]

        f.write("\n")
        f.write("\n")
        f.write("def draw_" + cellidentifier + "(savefig=False, cmap=plt.cm.gist_rainbow):\n")

        f.write("\n")
        f.write("    L_" + cellname + " = " + str(list(neighborhoods[cell].keys())) + "\n")
        f.write("    N_" + cellname + " = " + str(arr) + "\n")

        f.write("\n")
        f.write("    pca = PCA(n_components=3)\n")
        f.write("    pca.fit(N_" + cellname + ")\n")
        f.write("    explained_variance_2D = pca.explained_variance_ratio_[0] + pca.explained_variance_ratio_[1]\n")
        f.write("    explained_variance_3D = sum(pca.explained_variance_ratio_)\n")
        f.write("    components = pca.fit_transform(N_" + cellname + ")\n")
        f.write("    sortedcomp = sorted(list(zip(L_" + cellname + ", components)), key=lambda x: x[1][0])\n")
        f.write("    fig = plt.figure(figsize=(15, 6.5))\n")
        f.write("    ax = fig.add_subplot(1, 2, 1)\n")
        f.write("    ax.set_xlabel('Principal Component 1', fontsize=15)\n")
        f.write("    ax.set_ylabel('Principal Component 2', fontsize=15)\n")
        f.write("    title = '" + str(cell) + "' + ', explained variance = {:.2f}'.format(explained_variance_2D)\n")
        f.write("    ax.set_title(title, fontsize=20)\n")
        f.write("    for i in range(len(N_" + cellname + ")):\n")
        f.write("        ax.scatter(x=[sortedcomp[i][1][0]], y=[sortedcomp[i][1][1]], c=[i], vmin=0, \\\n")
        f.write("            vmax=len(L_" + cellname + "), label=sortedcomp[i][0], cmap=cmap)\n")
        f.write("    ax.grid()\n")
        f.write("    ax.legend()\n")
        f.write("    ymin, ymax = ax.get_ylim()\n")
        f.write("\n")

        f.write("    ax1 = fig.add_subplot(1, 2, 2, projection='3d')\n")
        f.write("    ax1.set_xlabel('Principal Component 2', fontsize=15)\n")
        f.write("    ax1.set_ylabel('Principal Component 1', fontsize=15)\n")
        f.write("    ax1.set_zlabel('Principal Component 3', fontsize=15)\n")
        f.write("    ax1.set_xlim(ymax, ymin)\n")
        f.write("    title = '" + str(cell) + "' + ', explained variance = {:.2f}'.format(explained_variance_3D)\n")
        f.write("    ax1.set_title(title, fontsize=20)\n")
        f.write("    for i in range(len(N_" + cellname + ")):\n")
        f.write("        ax1.scatter([sortedcomp[i][1][1]], [sortedcomp[i][1][0]], [sortedcomp[i][1][2]], c=[i], \\\n")
        f.write("            vmin=0, vmax=len(L_" + cellname + "), label=sortedcomp[i][0], cmap=cmap)\n")
        f.write("    norm = clr.Normalize(vmin=0, vmax=len(L_" + cellname + "))\n")
        f.write("    zmin, zmax = ax1.get_zlim()\n")
        f.write("    for i in range(len(N_" + cellname + ")):\n")
        f.write("        color = cmap(norm(i))\n")
        f.write("        ax1.plot([sortedcomp[i][1][1], sortedcomp[i][1][1]], \\\n")
        f.write("            [sortedcomp[i][1][0], sortedcomp[i][1][0]], [zmin, sortedcomp[i][1][2]], color=color)\n")
        f.write("\n")

        f.write("    if savefig:\n")
        f.write("        plt.savefig('" + str(cellidentifier))
        if file_suffix is not None:
            f.write(file_suffix)
        f.write("'" + " + '.png')\n")
        f.write("    else:\n")
        f.write("        plt.show()\n")

    f.write("\n")
    f.write("\n")
    f.write("def draw_all(savefig=False, cmap=plt.cm.gist_rainbow):\n")
    for cellidentifier in cellidentifierlist:
        f.write("    draw_" + cellidentifier + "(savefig=savefig, cmap=cmap)\n")
        f.write("    plt.close()\n")

    f.write("\n")
    f.write("\n")
    f.write("cmap = plt.cm.gist_rainbow\n")
    for cellidentifier in cellidentifierlist:
        f.write("if False:\n")
        f.write("    draw_" + cellidentifier + "(savefig=False, cmap=cmap)\n")
    f.write("\n")
    f.write("if True:\n")
    f.write("    draw_all(savefig=True, cmap=cmap)\n")
    f.write("\n")

    f.close()


def figures_histogram_scores(neighborhoods, parameters):

    filename = 'figures_histogram_scores'
    file_suffix = None
    if parameters.figurefile_suffix is not None and isinstance(parameters.figurefile_suffix, str) and \
            len(parameters.figurefile_suffix) > 0:
        file_suffix = '_' + parameters.figurefile_suffix
    filename += file_suffix + '.py'

    cell_score_by_stage = {}
    sister_score_by_stage = {}

    for cell in neighborhoods:
        stage = int(cell.split('.')[0][1:])
        sister = uname.get_sister_name(cell)
        for ref in neighborhoods[cell]:
            #
            # cell score
            #
            for r in neighborhoods[cell]:
                score = ucontact.contact_distance(neighborhoods[cell][ref], neighborhoods[cell][r],
                                                  similarity=parameters.contact_similarity)
                cell_score_by_stage[stage] = cell_score_by_stage.get(stage, []) + [score]
            if sister not in neighborhoods:
                continue
            if sister < cell:
                continue
            for r in neighborhoods[sister]:
                score = ucontact.contact_distance(neighborhoods[cell][ref], neighborhoods[sister][r],
                                                  similarity=parameters.contact_similarity)
                sister_score_by_stage[stage] = sister_score_by_stage.get(stage, []) + [score]

    cell_scores = []
    sister_scores = []
    for s in cell_score_by_stage:
        cell_scores += cell_score_by_stage[s]
    for s in sister_score_by_stage:
        sister_scores += sister_score_by_stage[s]

    f = open(filename, "w")

    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")
    f.write("\n")
    f.write("savefig = True\n")
    f.write("\n")
    f.write("cell_score_by_stage = " + str(cell_score_by_stage))
    f.write("\n")
    f.write("sister_score_by_stage = " + str(sister_score_by_stage))
    f.write("\n")
    f.write("cell_scores = []\n")
    f.write("sister_scores = []\n")
    f.write("for s in cell_score_by_stage:\n")
    f.write("    cell_scores += cell_score_by_stage[s]\n")
    f.write("for s in sister_score_by_stage:\n")
    f.write("    sister_scores += sister_score_by_stage[s]\n")
    f.write("\n")

    f.write("fig, ax = plt.subplots(figsize=(7.5, 7.5))\n")
    f.write("labels = ['same cell', 'sister cell']\n")
    f.write("ax.hist([cell_scores, sister_scores], 100, histtype='bar', label=labels)\n")
    f.write("ax.legend(prop={'size': 10})\n")
    f.write("ax.set_title('all scores', fontsize=12)\n")
    f.write("ax.tick_params(labelsize=10)\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('histogram_scores_all")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.write("\n")

    f.write("for s in cell_score_by_stage:\n")
    f.write("    fig, ax = plt.subplots(figsize=(7.5, 7.5))\n")
    f.write("    labels = ['same cell', 'sister cell']\n")
    f.write("    ax.hist([cell_score_by_stage[s], sister_score_by_stage[s]], 100, histtype='bar', label=labels)\n")
    f.write("    ax.legend(prop={'size': 10})\n")
    f.write("    title = \"scores for stage #{:02d}\".format(s)\n")
    f.write("    ax.set_title(title, fontsize=15)\n")
    f.write("    ax.tick_params(labelsize=15)\n")
    f.write("    if savefig:\n")
    f.write("        plt.savefig('histogram_scores_by_stage_S{:02d}")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'.format(s)" + " + '.png')\n")
    f.write("    else:\n")
    f.write("        plt.show()\n")
    f.write("\n")

    f.close()


def figures_histogram2D_scores(atlases, parameters):

    filename = 'figures_histogram2D_scores'
    file_suffix = None
    if parameters.figurefile_suffix is not None and isinstance(parameters.figurefile_suffix, str) and \
            len(parameters.figurefile_suffix) > 0:
        file_suffix = '_' + parameters.figurefile_suffix
    filename += file_suffix + '.py'

    neighborhoods = atlases.get_neighborhoods()

    right_cscores = []
    right_sscores = []

    wrong_cscores = []
    wrong_sscores = []

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

                right_cscores.append(ucontact.contact_distance(neighborhoods[cell][ref], neighborhoods[cell][r],
                                                               similarity=parameters.contact_similarity))
                right_sscores.append(ucontact.contact_distance(neighborhoods[sister][ref], neighborhoods[cell][r],
                                                               similarity=parameters.contact_similarity))
                wrong_cscores.append(ucontact.contact_distance(neighborhoods[cell][ref], neighborhoods[cell][r],
                                                               similarity=parameters.contact_similarity))
                wrong_sscores.append(ucontact.contact_distance(neighborhoods[sister][ref], neighborhoods[cell][r],
                                                               similarity=parameters.contact_similarity))

    step = 0.01

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
    f.write("N = [[(Z[Z <= Z[j][i]].sum()) * scale for i in range(Z.shape[1])] for j in range(Z.shape[0])]\n")

    f.write("\n")
    f.write("fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(16, 6.5))\n")
    f.write("ax1.plot(right_cscores, right_sscores, 'k.', markersize=2)\n")
    f.write("ax1.set_box_aspect(1)\n")
    f.write("ax1.set_xlim(0, 1)\n")
    f.write("ax1.set_xlabel('Daughter #1/daughter #1 score', fontsize=15)\n")
    f.write("ax1.set_ylabel('Daughter #2/daughter #2 score', fontsize=15)\n")
    f.write("title = 'score couples'\n")
    f.write("ax1.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("ax2.imshow(np.rot90(Z), cmap=plt.cm.rainbow_r, extent=[0, 1, 0, 1], vmin=Z.min(), vmax=Z.max())\n")
    f.write("ax2.set_box_aspect(1)\n")
    f.write("ax2.set_xlabel('Daughter #1/daughter #1 score', fontsize=15)\n")
    f.write("title = 'kernel estimation in [{:.2f}, {:.2f}]'.format(Z.min(), Z.max())\n")
    f.write("ax2.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("ax3.imshow(np.rot90(N), cmap=plt.cm.rainbow_r, extent=[0, 1, 0, 1])\n")
    f.write("ax3.set_box_aspect(1)\n")
    f.write("CS = ax3.contour(X, Y, N, levels=np.arange(10, 100, 10))\n")
    f.write("ax3.clabel(CS, inline=True, fontsize=10)\n")
    f.write("ax3.set_xlabel('Daughter #1/daughter #1 score', fontsize=15)\n")
    f.write("title = 'score probability'\n")
    f.write("ax3.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('kde_estimation_rightpairing")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")

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
    f.write("ax1.set_xlabel('Daughter #1/daughter #2 score', fontsize=15)\n")
    f.write("ax1.set_ylabel('Daughter #2/daughter #1 score', fontsize=15)\n")
    f.write("title = 'score couples'\n")
    f.write("ax1.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("ax2.imshow(np.rot90(Z), cmap=plt.cm.rainbow_r, extent=[0, 1, 0, 1], vmin=Z.min(), vmax=Z.max())\n")
    f.write("ax2.set_box_aspect(1)\n")
    f.write("ax2.set_xlabel('Daughter #1/daughter #2 score', fontsize=15)\n")
    f.write("title = 'kernel estimation in [{:.2f}, {:.2f}]'.format(Z.min(), Z.max())\n")
    f.write("ax2.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("ax3.imshow(np.rot90(N), cmap=plt.cm.rainbow_r, extent=[0, 1, 0, 1])\n")
    f.write("ax3.set_box_aspect(1)\n")
    f.write("CS = ax3.contour(X, Y, N, levels=np.arange(10, 100, 10))\n")
    f.write("ax3.clabel(CS, inline=True, fontsize=10)\n")
    f.write("ax3.set_xlabel('Daughter #1/daughter #2 score', fontsize=15)\n")
    f.write("title = 'score probability'\n")
    f.write("ax3.set_title(title, fontsize=15)\n")

    f.write("\n")
    f.write("if savefig:\n")
    f.write("    plt.savefig('kde_estimation_wrongpairing")
    if file_suffix is not None:
        f.write(file_suffix)
    f.write("'" + " + '.png')\n")
    f.write("else:\n")
    f.write("    plt.show()\n")
    f.close()
