
import numpy as np
from scipy.ndimage import binary_closing as nd_binary_closing
from scipy.ndimage import find_objects
from skimage.morphology import binary_dilation
import networkx as nx
import raster_geometry as rg
import os
import pandas as pd
import pickle



def calculate_incremented_bounding_boxes(watershed_labels, increment = 1):

    """ calculate incremented bounding boxes
    Args:
        watershed_labels(array): labelled image of cell positions
        increment (int): factor by how much bounding boxes should be incremented

    Returns: 
        incremented_b_boxes (): A list of tuples, with 3 slices each. Bounding boxes are incremented by one voxel in each dimension. Labels correspond to the indexes of the bounding box.
    """

    #calculate bounding boxes
    b_boxes = find_objects(watershed_labels)
    z_max, x_max, y_max = watershed_labels.shape

    #increment bounding boxes by one voxel
    incremented_b_boxes = {}
    for i, b_box in enumerate(b_boxes):
        if b_box is not None:
            new_entry = (slice(max(b_box[0].start - increment, 0), min(b_box[0].stop + increment, z_max-1), b_box[0].step),
                          slice(max(b_box[1].start - increment, 0), min(b_box[1].stop + increment, x_max-1), b_box[1].step),
                          slice(max(b_box[2].start - increment, 0), min(b_box[2].stop + increment, y_max-1), b_box[2].step))

            #adapt index to match the label of the object
            incremented_b_boxes[i+1] = new_entry

    return incremented_b_boxes

def extract_touching_surfaces(watershed_labels):
    
    """ find membranes between all predicted cells (based on watershed labels)
    Args:
        watershed_labels (3d array): watershed labels as 3D matrix of interegers, where 0 is background

    Returns:
        membrane_labels (3D array): 3D matrix of integers, showing all overlapping membranes
        mapp_mem (dict): dictionary that maps membrane id to tuples of neighboring cells, as labelled in watershed image

    """
    
    #create list of all labels and remove the background label (= 1 in astec output)
    unique_labels = sorted(list(np.unique(watershed_labels)))
    
    #find whether 0 or 1 is largest label and remove the bigger one (assuming it is the background)
    zero_count = np.count_nonzero(watershed_labels == 0)
    one_count = np.count_nonzero(watershed_labels == 1)
     
    if zero_count > one_count:
        if 0 in unique_labels:
            unique_labels.remove(0)
    elif zero_count < one_count:        
        if 1 in unique_labels:
            unique_labels.remove(1)
    else:
        print("neither 0 nor 1 is a label in the segmentation image")
        
    #calculate bounding boxes for all cells, dilate bounding box by 1 voxel
    incremented_b_boxes = calculate_incremented_bounding_boxes(watershed_labels)
    
    #initiate output_image and mapper dictionary
    membrane_labels = np.zeros_like(watershed_labels)
    mapp_mem = {}
    
    #iterate through labels and find neighbors, calculate overlap
    membrane_id = 1
    for label in unique_labels:
        b_box = incremented_b_boxes[label]
        sub_image = watershed_labels[b_box]
        
        #find neighbours in bounding box
        mask_dilation = binary_dilation(sub_image == label)
        neighbors = np.unique(sub_image[mask_dilation])
        
        #put all overlaps in one image as labels, map labels to cell combinations in mapper dictionary   
        for neighbor in neighbors:
            #because we sorted unique_labels, we know that if neighbor is smaller than label, 
            #we have already processed this combination
            if neighbor > label:
                n_mask = mask_dilation & (sub_image == neighbor)
                if (n_mask.sum() > 0) and (tuple(sorted([label, neighbor])) not in mapp_mem.keys()):
                    membrane_labels[b_box][n_mask] = membrane_id
                    mapp_mem[tuple(sorted([label, neighbor]))] = int(membrane_id)
                    membrane_id += 1
   
    #some very small membrane fragment might have gotten covered by bigger ones, so we make sure only existing membranes exist in the mapper
    final_membranes = np.unique(membrane_labels) 
    lost_membranes = set(mapp_mem.values()).difference(set(final_membranes))
    if len(lost_membranes) > 0:
        mapp_mem = {pair: mem_id for pair, mem_id in mapp_mem.items() if mem_id in set(final_membranes)}
                 
    return membrane_labels, mapp_mem



def volume_ratio_after_closing(interface_image, mapper, voxel_size = (1, 0.173, 0.173), iterations = 1):
    
    '''
    calculate ratio of the volumes of each interface before and after using scipy.ndimage.binary_closing 
    --> volume increase points to abrupt shape changes in membrane = not an actual membrane
    Args:
        interface_image
        mapper (dict): dictionary mapping each membrane to the pair of cells it seperates: {(label cell1, label cell2) : membrane_id}
        appr_int_scale(tuple): voxel_size in order z, y, x
        iterations (int): how many times should the closing algorithm be applied
    '''
    
    volume_ratios = {}
    volumes = {}

    # test if all values are equal
    if len(set(voxel_size)) > 1:
           
        # create 3D kernel for dilation which takes the anisotropy of the image into account
        z_dim, y_dim, x_dim = voxel_size
           
        # use scaling factor to find the smallest possible ellipsoid
        semiaxes = np.array([1/z_dim, 1/y_dim, 1/x_dim])*np.max(voxel_size)
        shape = [int(c) for c in np.ceil(semiaxes*2+1)]
        structure = rg.ellipsoid(shape, semiaxes)
        increment = np.ceil(iterations * np.max(semiaxes)).astype(int)
    else:
        structure = None
    # calculate bounding boxes for all cells, dilate bounding box by amount needied for the binary closing
    incremented_b_boxes = calculate_incremented_bounding_boxes(interface_image, increment = increment)
    
    for label in mapper.values():
        # subset image for region around membrane of interest
        b_box = incremented_b_boxes[label]
        sub_image = interface_image[b_box]
        interface = sub_image == label
        if not interface.any():
            continue
            
        interface_closed = nd_binary_closing(interface, structure = structure, iterations = iterations)
        
        vol_before = interface.sum()
        vol_after = interface_closed.sum()
        volume_ratios[label] = float(vol_after/vol_before)
        volumes[label] = int(vol_before)

    return volume_ratios, volumes



def find_connected_components(false_pairs_list): 
    
    """find connected components among pairs of cells to be merged, create sets of labels that are to be merged into one cell
    Args:
        false_labels_list (list): list with all pairs of cells that need to be merged
    Returns:
        cc_list (list): list of sets of cells that are to be merged
    """

    G = nx.Graph()
    for edge in false_pairs_list:
        G.add_edge(*edge)
        
    cc_list = list(nx.connected_components(G))
    
    return cc_list


def merge_labels_with_false_membranes(false_pairs_list, original_watershed_labels):
    
    """combine labels that are seperated by interfaces that are detected as false membranes
    Args:
        false_labels_list (list): list with all pairs of cells that need to be merged
        original_watershed_labels(array): array of integers, labelled image
        output_name(str): name for saving the new watershed image with merged cells
        
    Returns:
        merged_watershed (3D array): new image with merged cells, according to labels in merging_dictionary
    """
    
    #find connected components among cells that need to be merged (group all labels for merging into one big cell)
    cc_list = find_connected_components(false_pairs_list)
    #create image with merged cells based on sets of cells and original watershed image
    merged_watershed = original_watershed_labels.copy()
    merging_dict = {}
    for cell_set in cc_list:
        smallest_id = min(cell_set)
        merging_dict[smallest_id] = cell_set
        for label in cell_set:
            merged_watershed[merged_watershed == label] = smallest_id
    
    return merged_watershed, merging_dict

def translate_cell_pair_to_previous (cell_pair, reversed_correspondences):
    return tuple(sorted([reversed_correspondences[cell_pair[0]], reversed_correspondences[cell_pair[1]]]))