import os

import numpy as np
import pickle as pkl
import xml.etree.ElementTree as ElementTree


class ToyMonitoring(object):
    def __init__(self):
        self.verbose = -1
        self.debug = 0
        self.log_filename = None

    @staticmethod
    def to_console(text):
        print(text)

    def to_log(self, text):
        if self.log_filename is not None:
            with open(self.log_filename, "a") as logfile:
                logfile.write(text + "\n")

    def to_log_and_console(self, text, verboseness=0):
        if self.verbose >= verboseness or self.debug > 0:
            self.to_console(text)
        self.to_log(text)


import pkg_resources

try:
    dist = pkg_resources.get_distribution("astec")
    from astec.utils import common

    monitoring = common.Monitoring()
except pkg_resources.DistributionNotFound:
    monitoring = ToyMonitoring()

########################################################################################
#
# key correspondences
#
# Examples
# - from 'full_properties_Samson_MN20.pkl', keys are
#   ['volumes information',
#    'Cells labels in time',
#    'Barycenters',
#    'Fate',
#    'All Cells',
#    'Principal values',
#    'Names',
#    'cell_cell_contact_information',
#    'Lineage tree',
#    'Cells history',
#    'Principal vectors']
# - from 'new_lineage_tree_MN20.pkl', keys are cell labels
# - from '140317-Patrick-St8_seg_lineage.pkl', keys are
#   ['h_mins_information', 'lin_tree', 'volumes_information', 'sigmas_information']
#
#
########################################################################################

keydictionary = {
    "lineage": {
        "output_key": "cell_lineage",
        "input_keys": ["lineage_tree", "lin_tree", "Lineage tree", "cell_lineage"],
    },
    "h_min": {
        "output_key": "cell_h_min",
        "input_keys": ["cell_h_min", "h_mins_information"],
    },
    "volume": {
        "output_key": "cell_volume",
        "input_keys": [
            "cell_volume",
            "volumes_information",
            "volumes information",
            "vol",
        ],
    },
    "surface": {
        "output_key": "cell_surface",
        "input_keys": ["cell_surface", "cell surface"],
    },
    "compactness": {
        "output_key": "cell_compactness",
        "input_keys": [
            "cell_compactness",
            "Cell Compactness",
            "compacity",
            "cell_sphericity",
        ],
    },
    "sigma": {
        "output_key": "cell_sigma",
        "input_keys": ["cell_sigma", "sigmas_information", "sigmas"],
    },
    "label_in_time": {
        "output_key": "cell_labels_in_time",
        "input_keys": ["cell_labels_in_time", "Cells labels in time", "time_labels"],
    },
    "barycenter": {
        "output_key": "cell_barycenter",
        "input_keys": ["cell_barycenter", "Barycenters", "barycenters"],
    },
    "fate": {"output_key": "cell_fate", "input_keys": ["cell_fate", "Fate"]},
    "fate2": {"output_key": "cell_fate_2", "input_keys": ["cell_fate_2", "Fate2"]},
    "fate3": {"output_key": "cell_fate_3", "input_keys": ["cell_fate_3", "Fate3"]},
    "fate4": {"output_key": "cell_fate_4", "input_keys": ["cell_fate_4", "Fate4"]},
    "all-cells": {
        "output_key": "all_cells",
        "input_keys": ["all_cells", "All Cells", "All_Cells", "all cells", "tot_cells"],
    },
    "principal-value": {
        "output_key": "cell_principal_values",
        "input_keys": ["cell_principal_values", "Principal values"],
    },
    "apicobasal-length": {
        "output_key": "cell_apicobasal_length",
        "input_keys": ["cell_apicobasal_length"],
    },
    "name": {
        "output_key": "cell_name",
        "input_keys": ["cell_name", "Names", "names", "cell_names"],
    },
    "contact": {
        "output_key": "cell_contact_surface",
        "input_keys": ["cell_contact_surface", "cell_cell_contact_information"],
    },
    "contact-edge": {
        "output_key": "cell_contact_edge",
        "input_keys": ["cell_contact_edge"],
    },
    "contact-edge-length": {
        "output_key": "cell_contact_edge_length",
        "input_keys": ["cell_contact_edge_length"],
    },
    "contact-edge-segment": {
        "output_key": "cell_contact_edge_segment",
        "input_keys": ["cell_contact_edge_segment"],
    },
    "history": {
        "output_key": "cell_history",
        "input_keys": ["cell_history", "Cells history", "cell_life", "life"],
    },
    "principal-vector": {
        "output_key": "cell_principal_vectors",
        "input_keys": ["cell_principal_vectors", "Principal vectors"],
    },
    "name-score": {
        "output_key": "cell_naming_score",
        "input_keys": ["cell_naming_score", "Scores", "scores"],
    },
    "problems": {
        "output_key": "problematic_cells",
        "input_keys": ["problematic_cells"],
    },
    "urchin_apicobasal_length": {
        "output_key": "urchin_cell_apicobasal_length",
        "input_keys": ["urchin_cell_apicobasal_length"],
    },
    "urchin_apicobasal_segment": {
        "output_key": "urchin_cell_apicobasal_segment",
        "input_keys": ["urchin_cell_apicobasal_segment"],
    },
    "urchin_adjacency": {
        "output_key": "urchin_cell_adjacency",
        "input_keys": ["urchin_cell_adjacency"],
    },
    "urchin_apical_surface": {
        "output_key": "urchin_apical_surface",
        "input_keys": ["urchin_apical_surface"],
    },
    "urchin_apical_surface_barycenter": {
        "output_key": "urchin_apical_surface_barycenter",
        "input_keys": ["urchin_apical_surface_barycenter"],
    },
    "urchin_basal_surface": {
        "output_key": "urchin_basal_surface",
        "input_keys": ["urchin_basal_surface"],
    },
    "urchin_basal_surface_barycenter": {
        "output_key": "urchin_basal_surface_barycenter",
        "input_keys": ["urchin_basal_surface_barycenter"],
    },
    "urchin_apical_contact_edge_length": {
        "output_key": "urchin_apical_contact_edge_length",
        "input_keys": ["urchin_apical_contact_edge_length"],
    },
    "urchin_apical_contact_edge_segment": {
        "output_key": "urchin_apical_contact_edge_segment",
        "input_keys": ["urchin_apical_contact_edge_segment"],
    },
    "urchin_basal_contact_edge_length": {
        "output_key": "urchin_basal_contact_edge_length",
        "input_keys": ["urchin_basal_contact_edge_length"],
    },
    "urchin_basal_contact_edge_segment": {
        "output_key": "urchin_basal_contact_edge_segment",
        "input_keys": ["urchin_basal_contact_edge_segment"],
    },
    "urchin_vegetal_distance": {
        "output_key": "urchin_vegetal_distance",
        "input_keys": ["urchin_vegetal_distance"],
    },
    "unknown": {"output_key": "unknown_key", "input_keys": ["unknown_key"]},
}


def _normalize_dictionary_keys(inputdict):
    """

    :param inputdict:
    :return:
    """

    if inputdict == {}:
        return {}

    outputdict = {}

    for inputkey in inputdict:
        foundkey = False
        for k in keydictionary:
            # print "       compare '" + str(tmpkey) + "' with '" + str(k) + "'"
            if inputkey in keydictionary[k]["input_keys"]:
                outputkey = keydictionary[k]["output_key"]
                # monitoring.to_log_and_console("   ... recognized key '" + str(outputkey) + "'", 4)
                #
                # update if key already exists, else just create the dictionary entry
                #
                outputdict[outputkey] = inputdict[inputkey]
                foundkey = True
                break

        if foundkey is False:
            outputdict[inputkey] = inputdict[inputkey]

    return outputdict


def get_dictionary_entry(inputdict, keystring):
    proc = "get_dictionary_entry"
    if keystring not in keydictionary:
        monitoring.to_log_and_console(
            str(proc) + ": keystring must be in " + str(list(keydictionary.keys())), 1
        )
        return {}
    for k in keydictionary[keystring]["input_keys"]:
        if k in inputdict:
            return inputdict[k]
    else:
        monitoring.to_log_and_console(
            str(proc) + ": '" + str(keystring) + "' was not found in input dictionary",
            1,
        )
        monitoring.to_log_and_console(
            "    keys were: " + str(list(inputdict.keys())), 1
        )
        return {}


########################################################################################
#
# to translate a dictionary into XML
#
########################################################################################

#
# from stackoverflow.com
# questions/3095434/inserting-newlines-in-xml-file-generated-via-xml-etree-elementtree-in-python
#
#
# used for pretty printing
#


def _indent(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            _indent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


########################################################################################
#
# types
# 'lineage':  'lineage_tree' :
#     dictionary de liste de int
#     lineage_tree.590002 = <type 'list'>
# 'h_min' : 'cell_h_min' :
# 'volume' : 'cell_volume' :
#     dictionary de int
#     cell_volume.590002 = <type 'int'>
#     590002: 236936
# 'sigma': 'cell_sigma':
# 'label_in_time': 'cell_labels_in_time'
#     dictionary de liste de numpy.int64
#     cell_labels_in_time.1 = <type 'list'>
#     1: [10002, 10003, 10004, ..., 10082]
# 'barycenter': 'cell_barycenter'
#     dictionary de numpy.ndarray de numpy.float64
#     cell_barycenter.590002 = <type 'numpy.ndarray'>
#     590002: array([ 258.41037242,  226.74975943,  303.67167927])
# 'fate': 'cell_fate'
#     dictionary de str
#     cell_fate.590002 = <type 'str'>
#     590002: 'Mesoderm Notochord 1'
# 'all-cells': 'all_cells'  # liste de toutes les cellules ?
#     liste de numpy.int64
#     all_cells = <type 'list'>
# 'principal-value': 'cell_principal_values'
#     dictionary de liste de numpy.float64
#     cell_principal_values.590002 = <type 'list'>
#     590002: [1526.0489371146978, 230.60881177650205, 91.063513300019849]
# 'name': 'cell_name'
#     dictionary de str
#     cell_name.590002 = <type 'str'>
#     590002: 'a9.0012_'
# 'contact': 'cell_contact_surface',
#     dictionary de dictionary de int
#     cell_contact_surface.590002.590019 = <type 'int'>
#     590002: {590001: 1808, 590003: 1436, 590004: 5012, ..., 590225: 2579}
# 'history': 'cell_history'
#     dictionary de numpy.ndarray de numpy.int64
#     cell_history.590002 = <type 'numpy.ndarray'>
#     590002: array([510002, 520002, 530002, 540002, 550002, 560002, 570002, 580002,
#         590002, 600002, 610002, 620002, 630002, 640002, 650002, 660002,
#         670002, 680002, 690002, 700002, 710002, 720002, 730002, 740002,
#         750002, 760002, 770002, 780002, 790002, 800002, 810002, 820002,
#         830002, 840002, 850002, 860002, 870002, 880002])
# 'principal-vector': 'cell_principal_vectors'    # liste de numpy.ndarray
#     dictionary de liste de numpy.ndarray de numpy.float64
#     cell_principal_vectors.590002 = <type 'list'>
#     590002: [array([ 0.17420991, -0.74923203,  0.63898534]),
#         array([-0.24877611,  0.59437038,  0.7647446 ]),
#         array([ 0.95276511,  0.29219037,  0.08284582])]
#
########################################################################################


def _set_xml_element_text(element, value):
    """

    :param element:
    :param value:
    :return:
    """
    proc = "_set_xml_element_text"

    #
    # dictionary : recursive call
    #   dictionary element may be list, int, numpy.ndarray, str
    # list : may be int, numpy.int64, numpy.float64, numpy.ndarray
    #

    if type(value) == dict:
        # print proc + ": type is dict"
        keylist = list(value.keys())
        keylist.sort()
        for k in keylist:
            _dict2xml(element, k, value[k])

    elif type(value) == list:
        #
        # empty list
        #

        if len(value) == 0:
            element.text = repr(value)
        #
        # 'lineage', 'label_in_time', 'all-cells', 'principal-value'
        #

        elif type(value[0]) in (int, float, np.int64, np.float64):
            # element.text = str(value)
            value.sort()
            element.text = repr(value)

        elif isinstance(value[0], str):
            value.sort()
            text = "["
            for i in range(len(value)):
                text += "'" + value[i] + "'"
                if i < len(value) - 1:
                    text += ", "
                    if i > 0 and i % 5 == 0:
                        text += "\n  "
            text += "]"
            element.text = text
            del text
        #
        # 'principal-vector' case
        #  liste de numpy.ndarray de numpy.float64
        #
        elif type(value[0]) == np.ndarray:
            text = "["
            for i in range(len(value)):
                # text += str(list(value[i]))
                text += repr(list(value[i]))
                if i < len(value) - 1:
                    text += ", "
                    if i > 0 and i % 10 == 0:
                        text += "\n  "
            text += "]"
            element.text = text
            del text

        else:
            monitoring.to_log_and_console(
                proc
                + ": error, element list type ('"
                + str(type(value[0]))
                + "') not handled yet"
            )

    #
    # 'barycenter', 'cell_history'
    #
    elif type(value) == np.ndarray:
        # element.text = str(list(value))
        element.text = repr(list(value))

    #
    # 'volume', 'contact'
    #
    elif type(value) in (int, float, np.int64, np.float64):
        # element.text = str(value)
        element.text = repr(value)

    #
    # 'fate', 'name'
    #
    elif type(value) == str:
        element.text = repr(value)

    else:
        monitoring.to_log_and_console(
            proc
            + ": element type '"
            + str(type(value))
            + "' not handled yet, uncomplete translation"
        )


#
#
#


def _dict2xml(parent, tag, value):
    """

    :param parent:
    :param tag:
    :param value:
    :return:
    """

    #
    # integers can not be XML tags
    #
    if type(tag) in (int, np.int64):
        child = ElementTree.Element("cell", attrib={"cell-id": str(tag)})
    else:
        mn_type, name = _get_morphonet_type_name(str(tag))
        if mn_type is not None:
            child = ElementTree.Element(name, attrib={"mn_type": str(mn_type)})
        else:
            child = ElementTree.Element(str(tag))

    _set_xml_element_text(child, value)

    parent.append(child)
    return parent


#
# procedure d'appel et creation de la racine
#


def dict2xml(dictionary, defaultroottag="data"):
    """

    :param dictionary:
    :param defaultroottag:
    :return:
    """

    proc = "dict2xml"

    if type(dictionary) is not dict:
        monitoring.to_log_and_console(
            proc + ": error, input is of type '" + str(type(dictionary)) + "'"
        )
        return None

    #
    # s'il n'y a qu'un seul element dans le dictionnaire, on appelle la racine
    # d'apres celui-ci (eg lineage_tree), sinon on cree une racine qui contient
    # tous les elements
    #

    if len(dictionary) == 1:
        roottag = list(dictionary.keys())[0]
        root = ElementTree.Element(roottag)
        _set_xml_element_text(root, dictionary[roottag])

    elif len(dictionary) > 1:
        root = ElementTree.Element(defaultroottag)
        for k, v in dictionary.items():
            _dict2xml(root, k, v)

    else:
        monitoring.to_log_and_console(proc + ": error, empty dictionary ?!")
        return None

    _indent(root)
    tree = ElementTree.ElementTree(root)

    return tree


########################################################################################
#
# to translate a XML tree into dictionary
#
########################################################################################


def _set_dictionary_value(root):
    """

    :param root:
    :return:
    """

    if len(root) == 0:
        #
        # pas de branche, on renvoie la valeur
        #

        # return ast.literal_eval(root.text)
        if root.text is None:
            return None
        else:
            return eval(root.text)

    else:
        dictionary = {}

        for child in root:
            # print("child.tag=" + str(child.tag))
            # print "len(child)=" + str(len(child))
            # print "child.text=" + str(child.text)

            key = child.tag
            if child.tag == "cell":
                key = np.int64(child.attrib["cell-id"])
            elif "mn_type" in child.attrib:
                key = "morphonet_" + child.attrib["mn_type"] + "_" + str(child.tag)
            dictionary[key] = _set_dictionary_value(child)

    return dictionary


def xml2dict(tree):
    """

    :param tree:
    :return:
    """

    proc = "xml2dict"

    root = tree.getroot()

    dictionary = {}

    for k, v in keydictionary.items():
        if root.tag == v["output_key"]:
            monitoring.to_log_and_console(
                "   ... " + proc + ": process root.tag = '" + str(root.tag) + "'", 3
            )
            dictionary[str(root.tag)] = _set_dictionary_value(root)
            break
    else:
        for child in root:
            monitoring.to_log_and_console(
                "   ... " + proc + ": process child.tag = '" + str(child.tag) + "'", 3
            )
            value = _set_dictionary_value(child)
            if value is None:
                monitoring.to_log_and_console(
                    "       "
                    + proc
                    + ": empty property '"
                    + str(child.tag)
                    + "' ?! "
                    + " ... skip it",
                    1,
                )
            else:
                key = child.tag
                if child.tag == "cell":
                    key = np.int64(child.attrib["cell-id"])
                elif "mn_type" in child.attrib:
                    key = "morphonet_" + child.attrib["mn_type"] + "_" + str(child.tag)
                dictionary[key] = value

    return dictionary


########################################################################################
#
# to read a set of files into a dictionary
#
########################################################################################


#
# update dictionary from what has been read
#


def _update_read_dictionary(propertiesdict, tmpdict, filename):
    """

    :param propertiesdict:
    :param tmpdict:
    :return:
    """
    proc = "_update_read_dictionary"
    unknownkeys = []

    for tmpkey in tmpdict:
        foundkey = False

        for k in keydictionary:
            # print "       compare '" + str(tmpkey) + "' with '" + str(k) + "'"
            if tmpkey in keydictionary[k]["input_keys"]:
                outputkey = keydictionary[k]["output_key"]
                monitoring.to_log_and_console(
                    "   ... recognized key '" + str(outputkey) + "'", 2
                )
                #
                # update if key already exists, else just create the dictionary entry
                #
                if outputkey in propertiesdict:
                    if (
                        type(propertiesdict[outputkey]) is dict
                        and type(tmpdict[tmpkey]) is dict
                    ):
                        propertiesdict[outputkey].update(tmpdict[tmpkey])
                    elif (
                        type(propertiesdict[outputkey]) is list
                        and type(tmpdict[tmpkey]) is list
                    ):
                        propertiesdict[outputkey] += tmpdict[tmpkey]
                    else:
                        monitoring.to_log_and_console(
                            proc
                            + ": error, can not update property '"
                            + str(outputkey)
                            + "'"
                        )
                else:
                    propertiesdict[outputkey] = tmpdict[tmpkey]
                foundkey = True
                break

        if foundkey is False:
            unknownkeys.append(tmpkey)

    if len(unknownkeys) > 0 and len(unknownkeys) == len(list(tmpdict.keys())):
        #
        # no key was found
        # it is assumed it's a lineage tree: add some test here ?
        #
        monitoring.to_log_and_console(
            "   ... assume '" + str(filename) + "' is a lineage", 1
        )
        outputkey = keydictionary["lineage"]["output_key"]
        if outputkey in propertiesdict:
            if type(propertiesdict[outputkey]) is dict and type(tmpdict) is dict:
                propertiesdict[outputkey].update(tmpdict)
            else:
                monitoring.to_log_and_console(
                    proc + ": error, can not update property '" + str(outputkey) + "'"
                )
        else:
            propertiesdict[outputkey] = tmpdict

    elif len(unknownkeys) > 0:
        #
        # some unknown keys were found
        #
        monitoring.to_log_and_console(
            "   ... unrecognized key(s) are '" + str(unknownkeys) + "'", 1
        )

        # previous behavior: use keydictionary['unknown']['output_key'] as key
        # for *one* unknown property
        #
        #  outputkey = keydictionary['unknown']['output_key']
        # if len(unknownkeys) == 1:
        #     tmpkey = unknownkeys[0]
        #     if outputkey in propertiesdict:
        #         if type(propertiesdict[outputkey]) is dict and type(tmpdict[tmpkey]) is dict:
        #             propertiesdict[outputkey].update(tmpdict[tmpkey])
        #         elif type(propertiesdict[outputkey]) is list and type(tmpdict[tmpkey]) is list:
        #             propertiesdict[outputkey] += tmpdict[tmpkey]
        #         else:
        #             monitoring.to_log_and_console(proc + ": error, can not update property '" + str(outputkey)
        #                                           + "'")
        #     else:
        #         propertiesdict[outputkey] = tmpdict[tmpkey]
        # else:
        #     monitoring.to_log_and_console(proc + ": error, can not update many unknown properties")

        #
        # use unknown keys as such
        #
        for k in unknownkeys:
            propertiesdict[k] = tmpdict[k]

    return propertiesdict


#
# types issued from the reading of xml files may be erroneous
# fix it
#


def _set_types_from_xml(propertiesdict):
    """

    :param propertiesdict:
    :return:
    """

    if propertiesdict == {}:
        return {}

    if "cell_barycenter" in propertiesdict:
        monitoring.to_log_and_console("   ... translate types of 'cell_barycenter'", 3)
        for c in propertiesdict["cell_barycenter"]:
            propertiesdict["cell_barycenter"][c] = np.array(
                propertiesdict["cell_barycenter"][c]
            )

    if "cell_history" in propertiesdict:
        monitoring.to_log_and_console("   ... translate types of 'cell_history'", 3)
        for c in propertiesdict["cell_history"]:
            propertiesdict["cell_history"][c] = np.array(
                propertiesdict["cell_history"][c]
            )

    if "cell_principal_vectors" in propertiesdict:
        monitoring.to_log_and_console(
            "   ... translate types of 'cell_principal_vectors'", 3
        )
        for c in propertiesdict["cell_principal_vectors"]:
            for v in range(len(propertiesdict["cell_principal_vectors"][c])):
                propertiesdict["cell_principal_vectors"][c][v] = np.array(
                    propertiesdict["cell_principal_vectors"][c][v]
                )

    return propertiesdict


#
#
#


def _read_xml_file(filename, propertiesdict):
    monitoring.to_log_and_console("... reading '" + str(filename) + "'", 1)
    inputxmltree = ElementTree.parse(filename)
    tmpdict = xml2dict(inputxmltree)
    propertiesdict = _update_read_dictionary(propertiesdict, tmpdict, filename)
    del tmpdict
    return propertiesdict


#
#
#


def _read_pkl_file(filename, propertiesdict):
    monitoring.to_log_and_console("... reading '" + str(filename) + "'", 1)
    inputfile = open(filename, "rb")
    tmpdict = pkl.load(inputfile)
    inputfile.close()
    propertiesdict = _update_read_dictionary(propertiesdict, tmpdict, filename)
    del tmpdict
    return propertiesdict


#
#
#


def read_dictionary(inputfilenames, inputpropertiesdict={}):
    """

    :param inputfilenames:
    :param inputpropertiesdict:
    :return:
    """
    proc = "read_dictionary"

    if inputfilenames is None:
        monitoring.to_log_and_console(proc + ": error, no input files")
        return {}

    propertiesdict = inputpropertiesdict

    #
    #
    #

    if type(inputfilenames) == str:
        if not os.path.isfile(inputfilenames):
            monitoring.to_log_and_console(
                proc + ": error, file '" + str(inputfilenames) + "' does not exist"
            )
            return propertiesdict

        if inputfilenames.endswith("xml") is True:
            propertiesdict = _read_xml_file(inputfilenames, propertiesdict)
            propertiesdict = _set_types_from_xml(propertiesdict)
        elif inputfilenames.endswith("pkl") is True:
            propertiesdict = _read_pkl_file(inputfilenames, propertiesdict)
        else:
            monitoring.to_log_and_console(
                proc
                + ": error: extension not recognized for '"
                + str(inputfilenames)
                + "'"
            )

        propertiesdict = _normalize_dictionary_keys(propertiesdict)
        return propertiesdict

    #
    # here, we assume type(inputfilenames) == list
    #

    #
    # read xml files
    #

    for filename in inputfilenames:
        if not os.path.isfile(filename):
            monitoring.to_log_and_console(
                proc + ": error, file '" + str(filename) + "' does not exist"
            )
            continue

        if filename.endswith("xml") is True:
            propertiesdict = _read_xml_file(filename, propertiesdict)

    #
    # translation of xml may take place here
    #

    propertiesdict = _set_types_from_xml(propertiesdict)

    #
    # read pkl files
    #

    for filename in inputfilenames:
        if not os.path.isfile(filename):
            monitoring.to_log_and_console(
                proc + ": error, file '" + str(filename) + "' does not exist"
            )
            continue

        if filename.endswith("pkl") is True:
            propertiesdict = _read_pkl_file(filename, propertiesdict)

    #
    #
    #

    for filename in inputfilenames:
        if filename[len(filename) - 3 : len(filename)] == "xml":
            continue
        elif filename[len(filename) - 3 : len(filename)] == "pkl":
            continue
        else:
            monitoring.to_log_and_console(
                proc + ": error: extension not recognized for '" + str(filename) + "'"
            )

    propertiesdict = _normalize_dictionary_keys(propertiesdict)
    return propertiesdict


def write_dictionary(inputfilename, inputpropertiesdict):
    """

    :param inputfilename:
    :param inputpropertiesdict:
    :return:
    """
    proc = "write_dictionary"

    if inputfilename.endswith("pkl") is True:
        lineagefile = open(inputfilename, "wb")
        pkl.dump(inputpropertiesdict, lineagefile)
        lineagefile.close()
    elif inputfilename.endswith("xml") is True:
        xmltree = dict2xml(inputpropertiesdict)
        xmltree.write(inputfilename)
        del xmltree
    elif inputfilename.endswith("tlp") is True:
        write_tlp_file(inputfilename, inputpropertiesdict)
    else:
        monitoring.to_log_and_console(
            str(proc)
            + ": error when writing lineage file. Extension not recognized for '"
            + os.path.basename(inputfilename)
            + "'",
            1,
        )
    return


########################################################################################
#
#
#
########################################################################################


def _get_morphonet_type_name(key):
    mn_type = None
    if key[:10] == "morphonet_" or key[:10] == "selection_":
        if key[10:20] == "selection_":
            mn_type = "selection"
            name = key[20:]
        elif key[10:16] == "float_":
            mn_type = "float"
            name = key[16:]
        else:
            name = key[10:]
    else:
        name = key
    return mn_type, name


def write_morphonet_selection(d, time_digits_for_cell_id=4, directory=None):
    proc = "write_morphonet_selection"
    div = int(10**time_digits_for_cell_id)
    for key in d:
        if not isinstance(key, str):
            # print("skip key '" + str(key) + "', not a string")
            continue
        if len(key) < 10:
            # print("skip key '" + str(key) + "', too short")
            continue
        if key[:9] != "selection" and key[:9] != "morphonet":
            # print("skip key '" + str(key) + "', not a selection")
            continue

        mn_type, name = _get_morphonet_type_name(key)

        # print("write key '" + str(key) + "'")
        filename = name + ".txt"
        if directory is not None and isinstance(directory, str):
            if not os.path.isdir(directory):
                if not os.path.exists(directory):
                    os.makedirs(directory)
                else:
                    monitoring.to_log_and_console(
                        proc + ": '" + str(directory) + "' is not a directory ?!"
                    )
            if os.path.isdir(directory):
                filename = os.path.join(directory, filename)
        f = open(filename, "w")
        f.write("# " + str(name) + "\n")
        if mn_type == "float":
            f.write("type:float\n")
        elif mn_type == "selection":
            f.write("type:selection\n")
        else:
            f.write("type:float\n")
        for c in d[key]:
            #
            # object tuple (OTP): t, id, ch,
            # specifying the time point (t),
            # the visualisation channel (ch) [can be skipped]
            # and the id of the specific object given in the obj file (id).
            #
            otp = "{:d}".format(int(c) // div) + ", {:d}".format(int(c) % div) + ": "
            if isinstance(d[key][c], (int, float, np.int64, np.float64)):
                f.write(otp + str(d[key][c]) + "\n")
            elif isinstance(d[key][c], (list, np.ndarray)):
                if isinstance(d[key][c][0], (int, float, np.int64, np.float64)):
                    for v in d[key][c]:
                        f.write(otp + str(v) + "\n")
                else:
                    msg = ": list element type'" + str(type(d[key][c][0])) + "'"
                    msg += "(key = '" + str(key) + "')"
                    msg += " not handled yet. Skip it."
                    monitoring.to_log_and_console(proc + msg)
            else:
                msg = ": element type'" + str(type(d[key][c])) + "'"
                msg += "(key = '" + str(key) + "')"
                msg += " not handled yet. Skip it."
                monitoring.to_log_and_console(proc + msg)
        f.close()


########################################################################################
#
# write tlp file
# this was inspired from pkl2tlp from L. Guignard
#
########################################################################################


def write_tlp_file(tlpfilename, dictionary):
    """

    :param tlpfilename:
    :param dictionary:
    :return:
    """

    proc = "write_tlp_file"

    #
    # is there a lineage
    #
    if keydictionary["lineage"]["output_key"] in dictionary:
        lineage = dictionary[keydictionary["lineage"]["output_key"]]
    else:
        monitoring.to_log_and_console(proc + ": no lineage was found.")
        return

    #
    # open file
    #
    f = open(tlpfilename, "w")
    f.write('(tlp "2.0"\n')

    #
    # write nodes = lineage.keys() + lineage.values()
    #
    nodes = set(lineage.keys()).union(
        set([v for values in list(lineage.values()) for v in values])
    )
    f.write("(nodes ")
    for n in nodes:
        f.write(str(n) + " ")
    f.write(")\n")

    #
    # write edges
    #
    count_edges = 0
    for m, ds in lineage.items():
        count_edges += 1
        for d in ds:
            f.write("(edge " + str(count_edges) + " " + str(m) + " " + str(d) + ")\n")

    #
    # write node ids
    #
    f.write('(property 0 int "id"\n')
    f.write('\t(default "0" "0")\n')
    for node in nodes:
        f.write("\t(node " + str(node) + str(' "') + str(node) + '")\n')
    f.write(")\n")

    #
    #
    #
    for p in dictionary:
        if p == keydictionary["lineage"]["output_key"]:
            pass
        elif p == keydictionary["all-cells"]["output_key"]:
            pass
        #
        # property as single double
        #
        elif (
            p == keydictionary["volume"]["output_key"]
            or p == keydictionary["surface"]["output_key"]
            or p == keydictionary["apicobasal-length"]["output_key"]
            or p == keydictionary["compactness"]["output_key"]
        ):
            prop = dictionary[p]
            default = np.median(list(prop.values()))
            f.write('(property 0 double "' + str(p) + '"\n')
            f.write('\t(default "' + str(default) + '" "0")\n')
            for node in nodes:
                f.write(
                    "\t(node "
                    + str(node)
                    + str(' "')
                    + str(prop.get(node, default))
                    + '")\n'
                )
            f.write(")\n")
        #
        # property as string
        #
        elif (
            p == keydictionary["fate"]["output_key"]
            or p == keydictionary["fate2"]["output_key"]
            or p == keydictionary["fate3"]["output_key"]
            or p == keydictionary["fate4"]["output_key"]
            or p == keydictionary["name"]["output_key"]
        ):
            prop = dictionary[p]
            f.write('(property 0 string "' + str(p) + '"\n')
            f.write('\t(default "' + "no string" + '" "0")\n')
            for node in nodes:
                f.write(
                    "\t(node "
                    + str(node)
                    + str(' "')
                    + str(prop.get(node, "no string"))
                    + '")\n'
                )
            f.write(")\n")
        #
        #
        #
        elif (
            p == keydictionary["h_min"]["output_key"]
            or p == keydictionary["sigma"]["output_key"]
            or p == keydictionary["label_in_time"]["output_key"]
            or p == keydictionary["barycenter"]["output_key"]
            or p == keydictionary["principal-value"]["output_key"]
            or p == keydictionary["contact"]["output_key"]
            or p == keydictionary["history"]["output_key"]
            or p == keydictionary["principal-vector"]["output_key"]
            or p == keydictionary["name-score"]["output_key"]
        ):
            pass
        else:
            monitoring.to_log_and_console(
                proc + ": property '" + str(p) + "' not handled yet for writing."
            )

    #
    # close file
    #
    f.write(")")
    f.write("(nodes ")

    f.close()
