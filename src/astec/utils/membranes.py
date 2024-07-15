
import numpy as np
from scipy.ndimage import binary_closing as nd_binary_closing
from scipy.ndimage import find_objects
from skimage.morphology import binary_dilation
import networkx as nx
import raster_geometry as rg


def calculate_incremented_bounding_boxes(watershed_labels, increment = 1):

    """ calculate incremented bounding boxes
    Args:
        watershed_labels(array): labelled image of cell positions
        increment (int): factor by how much bounding boxes should be incremented

    Returns: 
        incremented_b_boxes (): A list of tuples, with 3 slices each. Bounding boxes are incremented by one voxel in each dimension. Labels correspond to the indexes of the bounding box.
    """
    # in case 1 is the background, set it to 0
    zero_count = np.count_nonzero(watershed_labels == 0)
    one_count = np.count_nonzero(watershed_labels == 1)
     
    if one_count > zero_count:
        watershed_labels[watershed_labels == 1] = 0
        
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



def volume_ratio_after_closing(interface_image, mapper, voxelsize, iterations = 1):
    
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
    if len(set(voxelsize)) > 1:
        print(f"{voxelsize=}")
        # create 3D kernel for dilation which takes the anisotropy of the image into account
        x_dim, y_dim, z_dim = voxelsize
           
        # use scaling factor to find the smallest possible ellipsoid
        semiaxes = np.array([1/z_dim, 1/y_dim, 1/x_dim])*np.max(voxelsize)
        shape = [int(c) for c in np.ceil(semiaxes*2+1)]
        structure = rg.ellipsoid(shape, semiaxes)
        increment = np.ceil(iterations * np.max(semiaxes)).astype(int)
    else:
        structure = None
        increment = 1
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


def merge_labels_with_false_membranes(false_pairs_list, original_watershed_labels, selected_seeds, correspondences):
    
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
    changed_cells = {}
    #find mother_cell id for each cell, if there are several use the smallest id for the merged cell
    for cell_set in cc_list:
        list_of_mothers = [key for key, value in correspondences.items() if len(set(value).intersection(cell_set)) > 0]
        #TODO use reversed corr here and further down as well
        smallest_id = min(list_of_mothers)
        #TODO should this be one the labels of the daughters instead? this way we can make sure to not have this value assigned to another cell already
        for label in cell_set:
            mother = [key for key, value in correspondences.items() if label in value][0]
            #TODO only change if label != smallest_id
            #TODO use bounding boxes
            merged_watershed[merged_watershed == label] = smallest_id
            selected_seeds[selected_seeds == label] = smallest_id
            changed_cells[mother] = smallest_id
            # update correspondences, while not just replacing the old with new, but keeping sibling relationships
            old_daugthers = correspondences[mother]
            new_daughters = [smallest_id if x == label else x for x in old_daugthers]
            # in case we merged both of the daughters (pretty likely) we reduce the list using a set transformation
            correspondences[mother] = list(set(new_daughters))
    return merged_watershed, selected_seeds, changed_cells

def translate_cell_pair_to_previous (cell_pair, reversed_correspondences):
    return tuple(sorted([reversed_correspondences[cell_pair[0]], reversed_correspondences[cell_pair[1]]]))



def update_correspondences_dictionary(correspondences, changed_cells):
    ''' 
    update correspondences dictionary to remove mothers that have the same daughters (important for downstream compatibility)

    ''' 

    new_correspondences = {}
    processed_daughters = set()
    for mother_cell in correspondences:
        if mother_cell in changed_cells:
            new_daughter = changed_cells[mother_cell]
            if new_daughter not in processed_daughters:
                # filter correspondences for all entries containing the new daughter cell
                cut_correspondences = {key: value for key, value in correspondences.items() if new_daughter in value}
                if len(cut_correspondences) == 1: # because if there is just one, the necessary change was already doine during the merging
                    new_correspondences[mother_cell] = correspondences[mother_cell]
                else:
                    print(f"{len(cut_correspondences)=}")
                    len_daughters = [len(daughter_list) for daughter_list in cut_correspondences.values()]
                    # check whether any of those mothers have divided and not been merged
                    if max(len_daughters) == 1 and len(len_daughters) > 1:
                        #pick mother with the smallest label (that should be the same label as the daughter), discard the others
                        only_mother = min(cut_correspondences.keys())
                        new_correspondences[only_mother] = [new_daughter]
                    elif max(len_daughters) == 2: # that means there must be cells that have divided
                        # is there more than one daughter pair? 
                        daughter_pairs = [x for x in cut_correspondences.values() if len(x) == 2]
                        if len(daughter_pairs) > 1: # if there is just one, the necessary changes have already been done in the merging
                            for i, daughter_pair_list in enumerate(daughter_pairs):
                                mother = [key for key, value in correspondences.items() if value == daughter_pair_list][0]
                                # TODO find the right mother to keep - currently it is arbitrarily chosen as the first in line
                                if i == 0:
                                    new_daughters = daughter_pair_list
                                else:
                                    # remove entry from other daughter lists
                                    new_daughters = [x for x in daughter_pair_list if x != new_daughter]
                                new_correspondences[mother] = new_daughters
            processed_daughters.add(new_daughter)
        else:
            new_correspondences[mother_cell] = correspondences[mother_cell]
    return new_correspondences