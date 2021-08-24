
import astec.utils.neighborhood as uneighborhood

#
# figures
#

def figures_neighborhood_consistency(neighborhoods, parameters, filename='figures_consistency.py'):
    """
    Build a file to draw figures of 'consistency' between reference embryos for mother cells
    Parameters
    ----------
    neighborhoods: nested dictionary of neighborhood, where the keys are ['cell name']['reference name']
        where 'reference name' is the name of the reference lineage, and neighborhood a dictionary of contact surfaces
        indexed by cell names (only for the first time point after the division)
    filename
    min_discrepancy

    Returns
    -------

    """
    references, tested_couples, discrepancy = uneighborhood.get_neighborhood_consistency(neighborhoods, parameters)
    if len(discrepancy) == 0:
        return

    first_discrepancy = {}
    first_percents = []
    second_discrepancy = {}
    second_percents = []

    mother_names = sorted(discrepancy.keys())
    for n in mother_names:
        if uneighborhood.get_mother_name(n) in mother_names or uneighborhood.get_mother_name(n) in first_discrepancy:
            second_discrepancy[n] = discrepancy[n]
            second_percents.append(100.0 * float(len(discrepancy[n])) / float(tested_couples[n]))
        else:
            first_discrepancy[n] = discrepancy[n]
            first_percents.append(100.0 * float(len(discrepancy[n])) / float(tested_couples[n]))

    percents = []
    for mother_name in mother_names:
        percents.append(100.0 * float(len(discrepancy[mother_name])) / float(tested_couples[mother_name]))
    [sorted_percents, sorted_mothers] = list(zip(*sorted(zip(percents, mother_names), reverse=True)))

    f = open(filename, "w")

    f.write("import networkx as nx\n")
    f.write("import numpy as np\n")
    f.write("import matplotlib.pyplot as plt\n")

    cellidentifierlist = []
    for mother_name in sorted_mothers:
        percent = 100.0 * float(len(discrepancy[mother_name])) / float(tested_couples[mother_name])

        # identifier for mother cell
        cellidentifier = mother_name[0:2]+mother_name[3:7]
        if mother_name[7] == '_':
            cellidentifier += 'U'
        elif mother_name[7] == '*':
            cellidentifier += 'S'
        else:
            cellidentifier += 'S'
        if mother_name in first_discrepancy:
            cellidentifier += "_1"
        elif mother_name in second_discrepancy:
            cellidentifier += "_2"
        fileidentifier = 'D{:03d}_'.format(int(percent)) + cellidentifier
        cellidentifier += '_D{:03d}'.format(int(percent))

        cellidentifierlist.append(cellidentifier)

        edges = {}
        for n1 in references[mother_name]:
            edges[n1] = {}
            for n2 in references[mother_name]:
                if n2 <= n1:
                    continue
                edges[n1][n2] = 4
        for d in discrepancy[mother_name]:
            if d[0] < d[1]:
                edges[d[0]][d[1]] -= 1
            elif d[1] < d[0]:
                edges[d[1]][d[0]] -= 1

        f.write("\n")
        f.write("\n")
        f.write("def draw_" + cellidentifier + "(savefig=False, cutnodes=True):\n")

        nodes = list(references[mother_name])
        nodes.sort()

        f.write("\n")
        f.write("    G_" + cellidentifier + " = nx.Graph()\n")
        f.write("    G_" + cellidentifier + ".add_nodes_from(" + str(nodes) + ")\n")
        f.write("\n")

        edgelist = {}
        for n1 in nodes:
            for n2 in nodes:
                if n2 <= n1:
                    continue
                edgelist[edges[n1][n2]] = edgelist.get(edges[n1][n2], []) + [(n1, n2)]

        for i in range(1, 5):
            if i not in edgelist:
                continue
            if i == 2:
                f.write("    G_" + cellidentifier + ".add_edges_from(" + str(edgelist[i]) + ", weight=0)\n")
            else:
                f.write("    G_" + cellidentifier + ".add_edges_from(" + str(edgelist[i]) + ", weight=" + str(i) +
                        ")\n")
        f.write("\n")

        f.write("    snodes = nx.spectral_ordering(G_" + cellidentifier + ")\n")
        f.write("    scutmin, pnodes = nx.minimum_cut(G_" + cellidentifier + ", snodes[0], snodes[-1]," +
                " capacity='weight')\n")

        f.write("\n")
        f.write("    node_labels = {}\n")
        for n in nodes:
            f.write("    node_labels['" + str(n) + "'] = '" + str(n) + "'\n")
        f.write("\n")

        f.write("    fig, ax = plt.subplots(figsize=(7.5, 7.5))\n")
        f.write("    pos = nx.circular_layout(G_" + cellidentifier + ", scale=1.8)\n")
        f.write("    if cutnodes:\n")
        f.write("        nx.draw_networkx_nodes(G_" + cellidentifier + ", nodelist=pnodes[0] " +
                ", node_color='#b41f78', pos=pos, ax=ax, node_size=1000)\n")
        f.write("        nx.draw_networkx_nodes(G_" + cellidentifier + ", nodelist=pnodes[1] " +
                ", node_color='#1f78b4', pos=pos, ax=ax, node_size=1000)\n")
        f.write("    else:\n")
        f.write("        nx.draw_networkx_nodes(G_" + cellidentifier + ", pos=pos, ax=ax, node_size=1000)\n")
        f.write("    nx.draw_networkx_labels(G_" + cellidentifier + ", ax=ax, pos=pos, labels=node_labels," +
                " font_weight='bold')\n")
        f.write("    ax.set_xlim([-2.2, 2.2])\n")
        f.write("    ax.set_ylim([-2.1, 2.1])\n")
        for i in range(1, 5):
            if i not in edgelist:
                continue
            f.write("    nx.draw_networkx_edges(G_" + cellidentifier + ", pos=pos, ax=ax, edgelist=" + str(edgelist[i]))
            f.write(", edge_color='#1f78b4'")
            f.write(", width=" + str(i))
            if i == 1:
                f.write(", style='dotted'")
            elif i == 2:
                f.write(", style='dashed'")
            f.write(")\n")
        title = "division of " + str(mother_name)
        title += ", discrepancy = " + '{:.2f}'.format(percent) + "%"
        if mother_name in first_discrepancy:
            title += ", 1st order"
        elif mother_name in second_discrepancy:
            title += ", 2nd order"
        title += ", cut = {:2d}"
        f.write("    cstr = \"C{:02d}\".format(scutmin)\n")

        f.write("    ax.set_title(\"" + str(title) + "\".format(scutmin))\n")
        f.write("    if savefig:\n")
        f.write("        plt.savefig('" + str(fileidentifier) + "' + '_'" + " + str(cstr)" + " + '.png')\n")
        f.write("        plt.savefig('" + str(cellidentifier) + "' + '_'" + " + str(cstr)" + " + '.png')\n")
        f.write("        plt.savefig(str(cstr)" + " + '_" + str(cellidentifier) + ".png')\n")
        f.write("    else:\n")
        f.write("        plt.show()\n")

    f.write("\n")
    f.write("\n")
    f.write("def draw_all(savefig=False, cutnodes=True):\n")
    for cellidentifier in cellidentifierlist:
        f.write("    draw_" + cellidentifier + "(savefig=savefig, cutnodes=cutnodes)\n")
        f.write("    plt.close()\n")

    f.write("\n")
    f.write("\n")
    for cellidentifier in cellidentifierlist:
        f.write("if False:\n")
        f.write("    draw_" + cellidentifier + "(savefig=False, cutnodes=True)\n")
    f.write("\n")
    f.write("if True:\n")
    f.write("    draw_all(savefig=True, cutnodes=True)\n")
    f.write("\n")

    f.write("if True:\n")
    f.write("    first_percents = " + str(first_percents) + "\n")
    f.write("    second_percents = " + str(second_percents) + "\n")
    f.write("    fig, ax = plt.subplots(figsize=(5, 5))\n")
    f.write("    bins = range(0,100,5)\n")
    f.write("    label = ['1st order', '2nd order']\n")
    f.write("    ax.hist([first_percents, second_percents], bins, histtype='bar', label=label)\n")
    f.write("    ax.legend(prop={'size': 10})\n")
    f.write("    ax.set_title('neighborhood discrepancies')\n")
    f.write("    fig.tight_layout()\n")
    f.write("    if True:\n")
    f.write("        plt.savefig('neighborhood_consistency_histogram.png')\n")
    f.write("    else:\n")
    f.write("        plt.show()\n")
    f.write("\n")
    f.close()


def figures_neighborhood_score(neighborhoods, parameters):
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

    filename = 'figures_score'
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
        mother_name = uneighborhood.get_mother_name(cell_name)
        references[mother_name] = references.get(mother_name, set()).union(set(neighborhoods[cell_name].keys()))

    #
    # get an edge value in [0, 1]
    # 1 means perfect equality, while 0 means total difference
    #
    mother_names = sorted(references.keys())
    edgeValues = {}
    values = {}
    for n in mother_names:
        d = uneighborhood.get_daughter_names(n)
        edgeValues[n] = {}
        for n1 in references[mother_name]:
            edgeValues[n][n1] = {}
            for n2 in references[mother_name]:
                if n2 <= n1:
                    continue
                if d[0] in neighborhoods and (n1 in neighborhoods[d[0]] and n2 in neighborhoods[d[0]]):
                    s0 = uneighborhood.get_score(neighborhoods[d[0]][n1], neighborhoods[d[0]][n2],
                                                 parameters.neighborhood_comparison)
                else:
                    s0 = None
                if d[1] in neighborhoods and (n1 in neighborhoods[d[1]] and n2 in neighborhoods[d[1]]):
                    s1 = uneighborhood.get_score(neighborhoods[d[1]][n1], neighborhoods[d[1]][n2],
                                                 parameters.neighborhood_comparison)
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
        cellidentifier = mother_name[0:2]+mother_name[3:7]
        if mother_name[7] == '_':
            cellidentifier += 'U'
        elif mother_name[7] == '*':
            cellidentifier += 'S'
        else:
            cellidentifier += 'S'
        fileidentifier = 'S{:03d}_'.format(int(max_value)) + cellidentifier
        cellidentifier += '_S{:03d}'.format(int(max_value))

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
