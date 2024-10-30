import copy
import os
import sys

import astec.utils.common as common
import astec.utils.ascidian_name as uname
import astec.utils.atlas_embryo as uatlase
import astec.utils.neighborhood_distance as uneighborhood

monitoring = common.Monitoring()


###########################################################
#
#
#
############################################################


def _get_branch_length(cell, lineage):
    length = 0
    c = cell
    while c in lineage and len(lineage[c]) == 1:
        length += 1
        c = lineage[c][0]
    return length


def _get_symmetric_neighborhood(neighborhood):
    """
    Changes the names of the cells in the neighborhood to get the symmetric neighborhood
    :param neighborhood:
    :return:
    """
    symneighborhood = {}
    for n in neighborhood:
        if n == "background":
            sn = "background"
        elif n == "other-half":
            sn = "other-half"
        else:
            sn = uname.get_symmetric_name(n)
        symneighborhood[sn] = neighborhood[n]
    return symneighborhood


def neighborhood_normalization(neighborhood, atlas, timepoint, cell_normalization=None):
    if cell_normalization == "local":
        total = 0
        for n in neighborhood:
            total += neighborhood[n]
        for n in neighborhood:
            neighborhood[n] /= total
    elif cell_normalization == "global":
        vs = atlas.get_voxelsize_correction(timepoint)
        for n in neighborhood:
            neighborhood[n] *= vs * vs
    return neighborhood


###########################################################
#
# closed set: defined by a set of interior cells
#
############################################################


class CellAtlases(uatlase.Atlases):
    def __init__(self, parameters=None):
        uatlase.Atlases.__init__(self, parameters)

        self._default_delay = None

        # nested dictionary of neighborhoods, where the keys are
        # [delay_from_division]['cell name(s)']['reference name']
        # - delay_from_division: allows to have different neighborhoods for the same reference
        # - 'cell name(s)': a tuple of the cell names (Conklin) of the interior of the structure
        #    one name correspond to a single cell
        #    names of two sister cells correspond to a division
        # - 'reference name' is the file name,
        # a neighborhood is a dictionary of contact surfaces indexed by cell names
        self._cell_neighborhood = {}

        # nested dictionary of neighborhoods, where the keys are
        # [delay_from_division]['cell name(s)']['reference name']
        # give list of the cell name(s) volumes
        self._volumes = {}

    ############################################################
    #
    # getters / setters
    #
    ############################################################

    def get_default_delay(self):
        return self._default_delay

    def set_default_delay(self, delay=0):
        self._default_delay = delay

    def get_cell_neighborhood_delays(self):
        if self._cell_neighborhood == {}:
            return None
        delays = sorted(list(self._cell_neighborhood.keys()))
        return delays

    def get_cell_neighborhood(self, delay_from_division=None):
        if delay_from_division is None:
            delay = self.get_default_delay()
        else:
            delay = delay_from_division
        if delay not in self._cell_neighborhood:
            self._cell_neighborhood[delay] = {}
        return self._cell_neighborhood[delay]

    def set_cell_neighborhood(self, neighborhoods, delay_from_division=0):
        self._cell_neighborhood[delay_from_division] = neighborhoods

    def get_volumes(self, delay_from_division=None):
        if delay_from_division is None:
            delay = self.get_default_delay()
        else:
            delay = delay_from_division
        if delay not in self._volumes:
            self._volumes[delay] = {}
        return self._volumes[delay]

    def set_volumes(self, volumes, delay_from_division=0):
        self._volumes[delay_from_division] = volumes

    ############################################################
    #
    # Build single cell atlas
    #
    ############################################################

    def _add_cell_atlas(
        self,
        atlas,
        atlas_name,
        parameters,
        delay_from_division=0,
        time_digits_for_cell_id=4,
    ):
        """

        Parameters
        ----------
        atlas: embryo properties (atlas to be added)
        atlas_name: file or atlas name (will indexed the added neighborhoods for each cell name)
        parameters:
        delay_from_division:
        time_digits_for_cell_id:

        Returns
        -------

        """
        proc = "add_cell_atlas"

        if not isinstance(parameters, uatlase.AtlasParameters):
            monitoring.to_log_and_console(
                str(proc)
                + ": unexpected type for 'parameters' variable: "
                + str(type(parameters))
            )
            sys.exit(1)

        #
        # build a nested dictionary of neighborhood, where the keys are
        # ['cell name']['reference name']
        # where 'reference name' is the name
        # of the reference lineage, and neighborhood a dictionary of contact surfaces indexed by cell names
        # only consider the first time point after the division
        #

        neighborhoods = self.get_cell_neighborhood(
            delay_from_division=delay_from_division
        )
        volumes = self.get_volumes(delay_from_division=delay_from_division)

        #
        # clean cell names:
        # - remove empty names
        # - leading or trailing spaces
        # was necessary for original property files (science ones stored on fileshare)
        #
        cell_name = atlas.cell_name
        cells = list(cell_name.keys())
        for c in cells:
            if cell_name[c] == "":
                del cell_name[c]
                continue
            cell_name[c] = cell_name[c].strip()

        #
        cell_volume = atlas.cell_volume
        cell_lineage = atlas.cell_lineage
        cell_contact = atlas.cell_contact_surface

        div = 10**time_digits_for_cell_id

        #
        # get the daughter cells just after division
        #
        reverse_lineage = {v: k for k, values in cell_lineage.items() for v in values}
        daughters = [
            cell_lineage[c][0] for c in cell_lineage if len(cell_lineage[c]) == 2
        ]
        daughters += [
            cell_lineage[c][1] for c in cell_lineage if len(cell_lineage[c]) == 2
        ]

        ancestor_name = []
        missing_name = []
        missing_contact = []
        missing_neighbors = []

        #
        # parse the cells
        #
        for daugh in daughters:
            #
            # mother cell should be named
            #
            if reverse_lineage[daugh] not in cell_name:
                continue
            #
            # check whether the cell is in dictionaries
            #
            if daugh not in cell_name:
                if daugh not in missing_name:
                    missing_name.append(daugh)
                    monitoring.to_log_and_console(
                        "\t"
                        + str(proc)
                        + ": daughter cell #"
                        + str(daugh)
                        + " was not found in 'cell_name' dictionary. Skip it",
                        6,
                    )
                continue

            if daugh not in cell_contact:
                if daugh not in missing_contact:
                    missing_contact.append(daugh)
                    monitoring.to_log_and_console(
                        "\t"
                        + str(proc)
                        + ": daughter cell #"
                        + str(daugh)
                        + " was not found in 'cell_contact_surface' dictionary. Skip it"
                    )
                continue

            #
            # check whether the mother name is the right one
            #
            if cell_name[reverse_lineage[daugh]] != uname.get_mother_name(
                cell_name[daugh]
            ):
                msg = (
                    "weird, name of daughter cell #"
                    + str(daugh)
                    + " is "
                    + str(cell_name[daugh])
                )
                msg += (
                    " while its mother #" + str(reverse_lineage[daugh]) + " is named "
                )
                msg += str(cell_name[reverse_lineage[daugh]]) + ". Skip it"
                monitoring.to_log_and_console("\t" + str(proc) + ": " + msg)
                continue

            #
            # get the daughter cell after some additional delay
            # positive delay: count from the division
            # negative delay: count from the end of the shortest branch for the two sisters
            #   to get the same delay for the two sisters
            #   kept for historical reason (and also for comparison purpose)
            #
            d = daugh
            local_delay_from_division = 0
            if delay_from_division >= 0:
                local_delay_from_division = delay_from_division
            elif delay_from_division < 0:
                length0 = _get_branch_length(d, cell_lineage)
                sisters = copy.deepcopy(cell_lineage[reverse_lineage[d]])
                sisters.remove(d)
                length1 = _get_branch_length(sisters[0], cell_lineage)
                local_delay_from_division = min(length0, length1) + delay_from_division
                if local_delay_from_division < 0:
                    local_delay_from_division = 0

            for i in range(local_delay_from_division):
                if d not in cell_lineage:
                    break
                if len(cell_lineage[d]) > 1:
                    break
                nextd = cell_lineage[d][0]
                if nextd not in cell_name or nextd not in cell_contact:
                    break
                d = nextd

            #
            # build the neighborhood of cell of id 'd'
            #
            neighbor = {}
            neighbor_is_complete = True
            # half_id = '*' or '_'
            half_id = cell_name[d][-1]
            for c in cell_contact[d]:
                n = int(c) % div
                if n == 1 or n == 0:
                    neighbor["background"] = (
                        neighbor.get("background", 0) + cell_contact[d][c]
                    )
                else:
                    cname = c
                    #
                    # if c is not named, get a named ancestor
                    #
                    if cname not in cell_name:
                        while cname in reverse_lineage and cname not in cell_name:
                            cname = reverse_lineage[cname]
                    if cname not in cell_name:
                        neighbor_is_complete = False
                        if c not in missing_neighbors:
                            missing_neighbors.append(c)
                            msg = (
                                "cell #"
                                + str(c)
                                + " (nor any ancestor) was not found in 'cell_name' dictionary."
                            )
                            monitoring.to_log_and_console("\t" + proc + ": " + msg)
                        continue
                    #
                    # cname is either the cell id or the id of a named ancestor
                    #
                    this_cell_name = cell_name[cname]
                    if cname != c:
                        #
                        # ancestor_name is just used to keep trace of 'ancestor' replacement
                        #
                        if c not in ancestor_name:
                            ancestor_name += [c]
                            msg = (
                                "use name '"
                                + str(this_cell_name)
                                + "' of cell #"
                                + str(cname)
                            )
                            msg += " for cell #" + str(c)
                            monitoring.to_log_and_console("\t" + proc + ": " + msg)
                    if parameters.differentiate_other_half:
                        neighbor[this_cell_name] = (
                            neighbor.get(this_cell_name, 0) + cell_contact[d][c]
                        )
                    else:
                        if this_cell_name[-1] == half_id:
                            neighbor[this_cell_name] = (
                                neighbor.get(this_cell_name, 0) + cell_contact[d][c]
                            )
                        else:
                            neighbor["other-half"] = (
                                neighbor.get("other-half", 0) + cell_contact[d][c]
                            )

            #
            # check whether all neighboring cells have a name
            #
            if not neighbor_is_complete:
                msg = (
                    ": neighborhood of "
                    + str(cell_name[d])
                    + " is not complete. Skip it"
                )
                monitoring.to_log_and_console("\t" + str(proc) + msg)
                continue

            #
            # surface normalization
            #
            neighbor = neighborhood_normalization(
                neighbor,
                atlas,
                int(d) // div,
                cell_normalization=parameters.cell_normalization,
            )

            #
            # add cell neighborhood and volume
            #
            key_cell_name = cell_name[d]

            if key_cell_name not in neighborhoods:
                neighborhoods[key_cell_name] = {}
            if atlas_name in neighborhoods[key_cell_name]:
                msg = (
                    "weird, "
                    + str(atlas_name)
                    + " was already indexed for neighbors of cell "
                    + str(cell_name[d])
                )
                monitoring.to_log_and_console("\t" + str(proc) + ": " + msg)
            neighborhoods[key_cell_name][atlas_name] = neighbor

            if key_cell_name not in volumes:
                volumes[key_cell_name] = {}
            if atlas_name in volumes[key_cell_name]:
                msg = (
                    "weird, "
                    + str(atlas_name)
                    + " was already indexed for volume of cell "
                    + str(cell_name[d])
                )
                monitoring.to_log_and_console("\t" + str(proc) + ": " + msg)
            if d not in cell_volume:
                msg = (
                    "\t cell #" + str(d) + " was not found in 'cell_volume' dictionary."
                )
                monitoring.to_log_and_console("\t" + str(proc) + ": " + msg)
            else:
                volumes[key_cell_name][atlas_name] = cell_volume[d]
            #
            # add symmetric neighborhood if asked
            #
            if parameters.add_symmetric_neighborhood:
                sname = uname.get_symmetric_name(cell_name[d])
                sreference = "sym-" + atlas_name
                sneighbor = _get_symmetric_neighborhood(neighbor)
                key_sym_name = sname

                if key_sym_name not in neighborhoods:
                    neighborhoods[key_sym_name] = {}
                if sreference in neighborhoods[key_sym_name]:
                    msg = (
                        "weird, "
                        + str(sreference)
                        + " was already indexed for cell "
                        + str(sname)
                    )
                    monitoring.to_log_and_console("\t" + str(proc) + ": " + msg)
                neighborhoods[key_sym_name][sreference] = sneighbor

                if key_sym_name not in volumes:
                    volumes[key_sym_name] = {}
                if sreference in volumes[key_sym_name]:
                    msg = (
                        "weird, "
                        + str(sreference)
                        + " was already indexed for volume of cell "
                        + str(sname)
                    )
                    monitoring.to_log_and_console("\t" + str(proc) + ": " + msg)
                if d not in cell_volume:
                    msg = (
                        "\t cell #"
                        + str(d)
                        + " was not found in 'cell_volume' dictionary."
                    )
                    monitoring.to_log_and_console("\t" + str(proc) + ": " + msg)
                else:
                    volumes[key_sym_name][sreference] = cell_volume[d]

        self.set_cell_neighborhood(
            neighborhoods, delay_from_division=delay_from_division
        )
        self.set_volumes(volumes, delay_from_division=delay_from_division)

        if len(missing_name) > 0:
            msg = (
                ": daughter cells without names = "
                + str(len(missing_name))
                + "/"
                + str(len(daughters))
            )
            msg += " in '" + str(atlas_name) + "'"
            monitoring.to_log_and_console(str(proc) + msg)

        if len(missing_contact) > 0:
            msg = ": daughter cells without contact surfaces  = " + str(
                len(missing_contact)
            )
            msg += "/" + str(len(daughters)) + " in '" + str(atlas_name) + "'"
            monitoring.to_log_and_console(str(proc) + msg)

        if len(missing_neighbors) > 0:
            msg = ": neighboring cells without names  = " + str(len(missing_neighbors))
            msg += " in '" + str(atlas_name) + "'"
            monitoring.to_log_and_console(str(proc) + msg)

        # monitoring.to_log_and_console("")

        return

    def build_cell_atlases(
        self, parameters, delay_from_division=0, time_digits_for_cell_id=4
    ):
        ref_atlases = self.get_atlases()
        for a in ref_atlases:
            self._add_cell_atlas(
                ref_atlases[a],
                a,
                parameters,
                delay_from_division=delay_from_division,
                time_digits_for_cell_id=time_digits_for_cell_id,
            )

    ############################################################
    #
    #
    #
    ############################################################

    def generate_figure(self, parameters, time_digits_for_cell_id=4):
        uatlase.Atlases.generate_figure(self, parameters)

        do_generate_figure = (
            (
                isinstance(parameters.generate_figure, bool)
                and parameters.generate_figure
            )
            or (
                isinstance(parameters.generate_figure, str)
                and parameters.generate_figure == "all"
            )
            or (
                isinstance(parameters.generate_figure, list)
                and "all" in parameters.generate_figure
            )
        )

        #
        # cell-to-cell distance between successive cells in a branch with respect to distance from first cell
        # (from first time point/division to last time point/division)
        #
        if (
            (
                isinstance(parameters.generate_figure, str)
                and parameters.generate_figure == "cell-distance-along-branch"
            )
            or (
                isinstance(parameters.generate_figure, list)
                and "cell-distance-along-branch" in parameters.generate_figure
            )
            or do_generate_figure
        ):
            monitoring.to_log_and_console(
                "... generate cell distance along branch file", 1
            )
            _figures_distance_along_branch(
                self, parameters, time_digits_for_cell_id=time_digits_for_cell_id
            )
            monitoring.to_log_and_console("... done", 1)


################################################################################
#
# cell-to-cell distance between successive cells in a branch with
# respect to distance from first cell
# (from first time point/division to last time point/division)
#
################################################################################


def _figures_distance_along_branch(atlases, parameters, time_digits_for_cell_id=4):
    """

    Parameters
    ----------
    atlases
    parameters
    time_digits_for_cell_id

    Returns
    -------

    """
    proc = "_figures_distance_along_branch"

    filename = "figures_distance_along_branch"
    file_suffix = None
    if (
        parameters.figurefile_suffix is not None
        and isinstance(parameters.figurefile_suffix, str)
        and len(parameters.figurefile_suffix) > 0
    ):
        file_suffix = "_" + parameters.figurefile_suffix
    if file_suffix is not None:
        filename += file_suffix
    filename += ".py"

    if parameters.outputDir is not None and isinstance(parameters.outputDir, str):
        if not os.path.isdir(parameters.outputDir):
            if not os.path.exists(parameters.outputDir):
                os.makedirs(parameters.outputDir)
            else:
                monitoring.to_log_and_console(
                    proc + ": '" + str(parameters.outputDir) + "' is not a directory ?!"
                )
        if os.path.isdir(parameters.outputDir):
            filename = os.path.join(parameters.outputDir, filename)

    ref_atlases = atlases.get_atlases()
    div = 10**time_digits_for_cell_id

    contact_distance_along_time = {}

    for ref in ref_atlases:
        contact_distance_along_time[ref] = {}
        lineage = ref_atlases[ref].cell_lineage
        contact = ref_atlases[ref].cell_contact_surface
        name = ref_atlases[ref].cell_name
        reverse_lineage = {v: k for k, values in lineage.items() for v in values}

        #
        # get the first and last time points
        #
        cells = list(
            set(lineage.keys()).union(
                set([v for values in list(lineage.values()) for v in values])
            )
        )
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
        first_cells = [
            lineage[c][0]
            for c in lineage
            if len(lineage[c]) == 2 and int(c) // div < last_time - 1
        ]
        first_cells += [
            lineage[c][1]
            for c in lineage
            if len(lineage[c]) == 2 and int(c) // div < last_time - 1
        ]

        for cell in first_cells:
            if cell not in contact:
                msg = (
                    "    * weird, cell "
                    + str(cell)
                    + " is not in the 'contact surface' dictionary"
                )
                msg += " of atlas '" + str(ref) + "'"
                monitoring.to_log_and_console(msg, 3)
                continue

            if cell not in name:
                msg = (
                    "    * weird, cell "
                    + str(cell)
                    + " is not in the 'name' dictionary"
                )
                msg += " of atlas '" + str(ref) + "'"
                monitoring.to_log_and_console(msg, 3)
                keyd = cell
            else:
                keyd = name[cell]

            first_time = int(cell) // div
            pcell = cell
            pneigh = copy.deepcopy(contact[cell])
            pneigh = neighborhood_normalization(
                pneigh,
                ref_atlases[ref],
                first_time,
                cell_normalization=parameters.cell_normalization,
            )
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
                                msg = (
                                    "    * weird, cell "
                                    + str(c)
                                    + " is not in the reversed lineage"
                                )
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
                nneigh = neighborhood_normalization(
                    nneigh,
                    ref_atlases[ref],
                    t,
                    cell_normalization=parameters.cell_normalization,
                )
                d = uneighborhood.cell_distance(
                    pneigh, nneigh, change_contact_surfaces=False
                )
                contact_distance_along_time[ref][keyd] = contact_distance_along_time[
                    ref
                ].get(keyd, []) + [d]
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
    f.write(
        "lengths = [len(contact_distance[r][c]) for r in contact_distance for c in contact_distance[r]]\n"
    )

    f.write("\n")
    f.write("dict_dist_per_time = {}\n")
    f.write("for r in contact_distance:\n")
    f.write("    for c in contact_distance[r]:\n")
    f.write("        for i, v in enumerate(contact_distance[r][c]):\n")
    f.write(
        "            dict_dist_per_time[i] = dict_dist_per_time.get(i, []) + [contact_distance[r][c][i]]\n"
    )
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
    f.write('ax1.set_title("[t, t+1] distances", fontsize=15)\n')
    f.write("ax1.set_xticks(ticks)\n")
    f.write("ax1.set_xticklabels(ticks)\n")

    f.write("\n")
    f.write("ax2.set_ylim(0, 0.7)\n")
    f.write("ax2.set_xlim(0, 20.5)\n")
    f.write("ax2.set_xlabel('time from division', fontsize=15)\n")
    f.write("ax2.set_ylabel('distance', fontsize=15)\n")
    f.write("ax2.boxplot(dist_per_time)\n")
    f.write('ax2.set_title("[t, t+1] distances (close-up)", fontsize=15)\n')

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
