import copy
import numpy as np

import astec.utils.common as common
import astec.utils.ascidian_name as uname

monitoring = common.Monitoring()


###########################################################
#
#
#
############################################################


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


###########################################################
#
#
#
############################################################


def _is_ancestor_in_stages(neigh, neighbors_by_stage):
    """

    Parameters
    ----------
    neigh: the cell name to be consider
    neighbors_by_stage: cell names to be compared to

    Returns
    -------

    """

    if (
        neigh == "background"
        or neigh == "other-half"
        or isinstance(neigh, int)
        or isinstance(neigh, np.int64)
    ):
        return False

    stage = int(neigh.split(".")[0][1:])
    stages = list(neighbors_by_stage.keys())
    stages = sorted(stages, key=lambda x: int(x), reverse=True)

    for s in stages:
        if int(s) >= int(stage) or int(s) == 0:
            continue
        diff = int(stage) - int(s)
        mother = neigh
        for i in range(diff):
            mother = uname.get_mother_name(mother)
        if mother in neighbors_by_stage[s]:
            return True
    return False


def build_same_contact_surfaces(
    neighborhoods, celllist, celltobeexcluded=None, maximal_generation=None, debug=False
):
    """

    Parameters
    ----------
    neighborhoods: dictionary of dictionary of contact surface vectors
    celllist: build same set of contact surfaces for all neighborhoods of cell from the list
    celltobeexcluded: cells that will not forced to be replaced by an ancestor
        typically, the two daughter cells when the outer cells are homogenized at their mother generation
    maximal_generation: if not None, it is the maximal generation that has to be found in the neighborhood,
        except for the cells listed in celltobeexcluded
    debug

    Returns
    -------
    dictionary of contact surface vectors, with the same contact cells
    """
    common_neighborhoods = {}
    for c in celllist:
        if c not in neighborhoods:
            continue
        common_neighborhoods[c] = copy.deepcopy(neighborhoods[c])

    #
    # get the neighbors for each atlas
    #
    neighbors_by_stage = {}
    for c in celllist:
        if c not in neighborhoods:
            continue
        for r in neighborhoods[c]:
            for n in neighborhoods[c][r]:
                #
                # 'background' and 'other-half' are names used for named neighborhood
                # partially named neighborhood will have cell ids as keys
                #
                if (
                    n == "background"
                    or n == "other-half"
                    or isinstance(n, int)
                    or isinstance(n, np.int64)
                ):
                    neighbors_by_stage[0] = neighbors_by_stage.get(0, []) + [n]
                    continue
                stage = int(n.split(".")[0][1:])
                neighbors_by_stage[stage] = neighbors_by_stage.get(stage, []) + [n]
    # suppress duplicate
    for s in neighbors_by_stage:
        if s == 0:
            continue
        neighbors_by_stage[s] = list(sorted(set(neighbors_by_stage[s])))

    #
    # get the stages in decreasing order
    #
    stages = list(neighbors_by_stage.keys())
    stages = sorted(stages, key=lambda x: int(x), reverse=True)
    neighbor_already_processed = []

    #
    # process stages, the latest first
    #
    for s in stages:
        if s == 0:
            continue
        if debug:
            print("stage " + str(s))
        for neigh in neighbors_by_stage[s]:
            if debug:
                print("   neighbor " + str(neigh))
            if neigh in neighbor_already_processed:
                if debug:
                    print("   --- already processed ")
                continue
            #
            # should the neighbor be forced to be replaced by its mother?
            # - not if the cells are to be excluded
            # - not in case of generic name (eg 'background')
            #
            max_test = (
                (isinstance(celltobeexcluded, list) and neigh not in celltobeexcluded)
                and s > 0
                and isinstance(maximal_generation, int)
                and int(neigh.split(".")[0][1:]) > maximal_generation
            )
            if _is_ancestor_in_stages(neigh, neighbors_by_stage) or max_test:
                #
                # replace the neighbor by its mother
                #
                # 1. add mother to previous stage neighbors (if required)
                #
                mother = uname.get_mother_name(neigh)
                if debug:
                    print("      replace by mother " + str(mother))
                previous_stage = int(s) - 1
                if previous_stage not in neighbors_by_stage:
                    neighbors_by_stage[previous_stage] = []
                if mother not in neighbors_by_stage[previous_stage]:
                    neighbors_by_stage[previous_stage] += [mother]

                d = uname.get_daughter_names(mother)
                #
                # 2. add mother in the contact surfaces
                #    add daughter contact surfaces to the mother one, and remove the daughter contact surfaces
                #
                for c in celllist:
                    if c not in common_neighborhoods:
                        continue
                    for r in common_neighborhoods[c]:
                        if mother not in common_neighborhoods[c][r]:
                            common_neighborhoods[c][r][mother] = 0.0
                        for i in range(2):
                            if d[i] in common_neighborhoods[c][r]:
                                common_neighborhoods[c][r][
                                    mother
                                ] += common_neighborhoods[c][r][d[i]]
                                del common_neighborhoods[c][r][d[i]]

                #
                # 3. both daughters have been processed
                #
                neighbor_already_processed += d

            else:
                #
                # keep the neighbor
                # add it when it is missing
                #
                if debug:
                    print("      keep neighbor " + str(neigh))

                for c in celllist:
                    if c not in common_neighborhoods:
                        continue
                    for r in common_neighborhoods[c]:
                        if neigh not in common_neighborhoods[c][r]:
                            common_neighborhoods[c][r][neigh] = 0.0

        #
        # end of loop "for neigh in neighbors_by_stage[s]:
        #

    return common_neighborhoods


###########################################################
#
#
#
############################################################


def cell_distance_elements(vect0, vect1, innersurfaces=[]):
    n0 = 0.0
    n1 = 0.0
    nm = 0.0
    # sum of eligible surfaces for vect0
    # remove innersurfaces if any
    for k in vect0:
        if k not in innersurfaces:
            n0 += vect0[k]
    # sum of eligible surfaces for vectq
    # remove innersurfaces if any
    for k in vect1:
        if k not in innersurfaces:
            n1 += vect1[k]
    # sum of differences
    # 1. common neighbors
    # 2. neighbors only in vect0
    # 3. neighbors only in vect1
    for k in set(vect0.keys()) & set(vect1.keys()):
        if k not in innersurfaces:
            nm += abs(vect0[k] - vect1[k])
    for k in set(vect0.keys()) - set(vect1.keys()):
        if k not in innersurfaces:
            nm += abs(vect0[k])
    for k in set(vect1.keys()) - set(vect0.keys()):
        if k not in innersurfaces:
            nm += abs(vect1[k])
    return nm, n0, n1


def _l1_normalized_modulus(vect0, vect1, innersurfaces=[]):
    nm, n0, n1 = cell_distance_elements(vect0, vect1, innersurfaces=innersurfaces)
    return nm / (n0 + n1)


def cell_distance(neigh0, neigh1, change_contact_surfaces=True, title=None):
    """
    Compute the similarity score of two neighborhoods, as the normalized scalar
    product of the vectors of contact surfaces.
    Parameters
    ----------
    neigh0: dictionary depicting a cell neighborhood. Each key is a neighbor, and the
            associated dictionary value give the contact surface
    neigh1: dictionary depicting a cell neighborhood
    change_contact_surfaces: if True, build vectors with neighboring cell of same name
        (eg may fuse daughters into a fake mother cell, etc)
    title: if not None, print out the neighborhoods as well as the score. For debug purpose.

    Returns
    -------
    the similarity score.

    """

    if change_contact_surfaces:
        #
        # build a dictionary of one element
        # build_same_contact_surfaces() build a common reference for a list of keys of the dictionary
        #
        tmp = {"foo": {0: neigh0, 1: neigh1}}
        vect = build_same_contact_surfaces(tmp, ["foo"])
        score = _l1_normalized_modulus(vect["foo"][0], vect["foo"][1])
        if title is not None:
            _print_common_neighborhoods(vect["foo"][0], vect["foo"][1], title=title)
            monitoring.to_log_and_console("\t score = " + str(score) + "\n")
    else:
        score = _l1_normalized_modulus(neigh0, neigh1)

    return score
