
import copy
import functools
import math
import os
import sys

import astec.utils.common as common
import astec.utils.ascidian_name as uname

monitoring = common.Monitoring()

class ContactSurfaceParameters(common.PrefixedParameter):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, prefix=None):

        common.PrefixedParameter.__init__(self, prefix=prefix)

        if "doc" not in self.__dict__:
            self.doc = {}

        doc = "\t Defines the contact surface similarity. Contact surface vectors are normalized before"
        doc += "comparison. Possible values are:\n"
        doc += "\t - 'inner_product': normalization by the inner product.\n"
        doc += "\t - 'l1_distance': normalization by the l1-norm. \n"
        doc += "\t - 'l2_distance': normalization by the l1-norm.\n"
        doc += "\t Default is 'l2_distance'"
        doc += "\t This measure is normalized into [0, 1]: 0 means perfect equality, 1 means total dissimilarity"
        self.doc['contact_similarity'] = doc
        self.contact_similarity = 'l2_distance'

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('#')
        print('# ContactSurfaceParameters')
        print('#')
        print("")

        common.PrefixedParameter.print_parameters(self)

        self.varprint('contact_similarity', self.contact_similarity)

        print("")

    def write_parameters_in_file(self, logfile):
        logfile.write("\n")
        logfile.write("# \n")
        logfile.write("# ContactSurfaceParameters\n")
        logfile.write("# \n")
        logfile.write("\n")

        common.PrefixedParameter.write_parameters_in_file(self, logfile)

        self.varwrite(logfile, 'contact_similarity', self.contact_similarity)

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
        self.contact_similarity = self.read_parameter(parameters, 'contact_similarity', self.contact_similarity)
        self.contact_similarity = self.read_parameter(parameters, 'neighborhood_comparison', self.contact_similarity)

    def update_from_parameter_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = common.load_source(parameter_file)
        self.update_from_parameters(parameters)


#######################################################################################
#
# for debug purposes
#
#######################################################################################

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


#######################################################################################
#
#
#
#######################################################################################

def _is_ancestor_in_stages(neigh, neighbors_by_stage):
    """

    Parameters
    ----------
    neigh: the cell name to be consider
    neighbors_by_stage: cell names to be compared to

    Returns
    -------

    """

    if neigh == 'background' or neigh == 'other-half':
        return False

    stage = int(neigh.split('.')[0][1:])
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


def _compare_cell(a, b):
    #
    # sort
    # - stages in decreasing values
    # - [a, b] in increasing value
    # - p in in increasing value
    # - [*, _] in increasing value
    #
    if a == 'background':
        return -1
    if b == 'background':
        return 1
    if a == 'other-half':
        return -1
    if b == 'other-half':
        return 1
    astage = int(a.split('.')[0][1:])
    bstage = int(b.split('.')[0][1:])
    if astage < bstage:
        return 1
    elif astage > bstage:
        return -1
    aabvalue = a.split('.')[0][0]
    babvalue = b.split('.')[0][0]
    if aabvalue < babvalue:
        return -1
    elif aabvalue > babvalue:
        return 1
    ap = int(a.split('.')[1][0:4])
    bp = int(b.split('.')[1][0:4])
    if ap < bp:
        return -1
    elif ap > bp:
        return 1
    alrvalue = a.split('.')[1][4]
    blrvalue = b.split('.')[1][4]
    if alrvalue < blrvalue:
        return -1
    elif alrvalue > blrvalue:
        return 1
    return 0


def build_same_contact_surfaces(neighborhoods, debug=False):
    """

    Parameters
    ----------
    neighborhoods: dictionary of contact surface vectors
    debug

    Returns
    -------
    dictionary of contact surface vectors, with the same contact cells
    """

    common_neighborhoods = copy.deepcopy(neighborhoods)

    #
    # get the neighbors for each atlas
    #
    neighbors_by_stage = {}
    for r in neighborhoods:
        for n in neighborhoods[r]:
            if n == 'background' or n == 'other-half':
                neighbors_by_stage[0] = neighbors_by_stage.get(0, []) + [n]
                continue
            stage = int(n.split('.')[0][1:])
            neighbors_by_stage[stage] = neighbors_by_stage.get(stage, []) + [n]
    # suppress duplicate
    for s in neighbors_by_stage:
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
            # should the neighbor be replaced by its mother?
            #
            if _is_ancestor_in_stages(neigh, neighbors_by_stage):
                #
                # replace the neighbor by its mother
                #
                # 1. add mother to previous stage neighbors (if required)
                #
                mother = uname.get_mother_name(neigh)
                if debug:
                    print("      replace by mother " + str(mother))
                previous_stage = int(s)-1
                if previous_stage not in neighbors_by_stage:
                    neighbors_by_stage[previous_stage] = []
                if mother not in neighbors_by_stage[previous_stage]:
                    neighbors_by_stage[previous_stage] += [mother]

                d = uname.get_daughter_names(mother)
                #
                # 2. add mother in the contact surfaces
                #    add daughter contact surfaces to the mother one, and remove the daughter contact surfaces
                #
                for r in common_neighborhoods:
                    if mother not in common_neighborhoods[r]:
                        common_neighborhoods[r][mother] = 0.0
                    for i in range(2):
                        if d[i] in common_neighborhoods[r]:
                            common_neighborhoods[r][mother] += common_neighborhoods[r][d[i]]
                            del common_neighborhoods[r][d[i]]

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

                for r in common_neighborhoods:
                    if neigh not in common_neighborhoods[r]:
                        common_neighborhoods[r][neigh] = 0.0

        #
        # end of loop "for neigh in neighbors_by_stage[s]:
        #

    return common_neighborhoods

#######################################################################################
#
# compute distance between two vectors of contact surfaces
#
########################################################################################

def _scalar_product(vect0, vect1):
    n0 = 0.0
    n1 = 0.0
    ps = 0.0
    for k in vect0:
        n0 += vect0[k] * vect0[k]
    for k in vect1:
        n1 += vect1[k] * vect1[k]
    # & = set intersection
    for k in set(vect0.keys()) & set(vect1.keys()):
        ps += vect0[k] * vect1[k]
    return ps / (math.sqrt(n0 * n1))


def _l1_normalized_modulus(vect0, vect1):
    n0 = 0.0
    n1 = 0.0
    nm = 0.0
    for k in vect0:
        n0 += vect0[k]
    for k in vect1:
        n1 += vect1[k]
    for k in set(vect0.keys()) & set(vect1.keys()):
        nm += abs(vect0[k]/n0 - vect1[k]/n1)
    for k in set(vect0.keys()) - set(vect1.keys()):
        nm += abs(vect0[k] / n0)
    for k in set(vect1.keys()) - set(vect0.keys()):
        nm += abs(vect1[k] / n0)
    return nm


def _l2_normalized_modulus(vect0, vect1):
    n0 = 0.0
    n1 = 0.0
    nm = 0.0
    for k in vect0:
        n0 += vect0[k]
    for k in vect1:
        n1 += vect1[k]
    for k in set(vect0.keys()) & set(vect1.keys()):
        nm += (vect0[k]/n0 - vect1[k]/n1) * (vect0[k]/n0 - vect1[k]/n1)
    for k in set(vect0.keys()) - set(vect1.keys()):
        nm += vect0[k]/n0 * vect0[k]/n0
    for k in set(vect1.keys()) - set(vect0.keys()):
        nm += vect1[k]/n1 * vect1[k]/n1
    return math.sqrt(nm)


def _contact_distance(neigh0, neigh1, similarity='l2_distance'):
    """
    Compute distance between two contact surface vectors. Do not compute 'common' neighborhood.
    Parameters
    ----------
    neigh0
    neigh1
    similarity

    Returns
    -------

    """
    proc = "_contact_distance"
    #
    # scalar product in [0, 1] (0: full disagreement, 1: perfect agreement)
    # _l1_normalized_modulus in [0, 2] (0: perfect agreement, 2:full disagreement)
    # _l2_normalized_modulus in [0, sqrt(2)] (0: perfect agreement, sqrt(2):full disagreement)
    #
    if similarity.lower() == 'l1_distance' or similarity.lower() == 'l1-distance' or similarity.lower() == 'l1_norm' \
            or similarity.lower() == 'l1-norm':
        score = _l1_normalized_modulus(neigh0, neigh1) / 2.0
    elif similarity.lower() == 'l2_distance' or similarity.lower() == 'l2-distance' or similarity.lower() == 'l2_norm' \
            or similarity.lower() == 'l2-norm':
        score = _l2_normalized_modulus(neigh0, neigh1) / 1.4142136
    elif similarity.lower() == 'scalar_product' or similarity.lower() == 'scalar-product' \
            or similarity.lower() == 'inner_product' or similarity.lower() == 'inner-product':
        score = 1.0 - _scalar_product(neigh0, neigh1)
    else:
        monitoring.to_log_and_console(proc + ": unhandled score computation '" + str(similarity) + "'\n")
        sys.exit(1)
    return score


def contact_distance(neigh0, neigh1, similarity='l2_distance', change_contact_surfaces=True, title=None, debug=False):
    """
    Compute the similarity score of two neighborhoods, as the normalized scalar
    product of the vectors of contact surfaces.
    Parameters
    ----------
    neigh0: dictionary depicting a cell neighborhood. Each key is a neighbor, and the
            associated dictionary value give the contact surface
    neigh1: dictionary depicting a cell neighborhood
    similarity:
    change_contact_surfaces: if True, build vectors with neighboring cell of same name
        (eg may fuse daughters into a fake mother cell, etc)
    title: if not None, print out the neighborhoods as well as the score. For debug purpose.
    debug

    Returns
    -------
    the similarity score.

    """

    if change_contact_surfaces:
        tmp = {}
        tmp[0] = neigh0
        tmp[1] = neigh1
        vect = build_same_contact_surfaces(tmp, debug=debug)
        score = _contact_distance(vect[0], vect[1], similarity)
        if title is not None:
            _print_common_neighborhoods(vect[0], vect[1], title=title)
            monitoring.to_log_and_console("\t score = " + str(score) + "\n")
    else:
        score = _contact_distance(neigh0, neigh1, similarity)

    return score
