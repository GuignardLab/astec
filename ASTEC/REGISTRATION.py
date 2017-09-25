import numpy as np
import sys
import os
from scipy import ndimage as nd
import pickle as pkl
sys.path.append(os.path.join(os.path.dirname(__file__),"CommunFunctions"))
from ImageHandling import imread, imsave, SpatialImage
from cpp_wrapping import *

#def associate_labels(path_ref, path_flo, path_pairs_ref_flo, path_ref_out, path_flo_out, path_labels_out, zeros=False, ref=False, force=False, verbose=False):

def readLUT(file):
	'''
	Return a dictionnary of integer key-to-key correspondances 
	'''
	lut={}
	if os.path.exists(file):
		f=open(file)
		for line in f:
			li=line.strip()
			if not li.startswith('#'):
				info=li.split()
				if len(info) == 2:
					if not lut.has_key(int(info[0])):
						lut[int(info[0])]=None
					lut[int(info[0])]=int(info[1])
	return lut

def symmetry_plane(image_bin, image_seg=None, image_direction_histogram="tmp_sphere.inr", 
	equation_output=None, trsf_output=None, plane_output=None, 
	frac=None, maxima=None, 
	d=10, dmin=2, delta=2, p=0.1, sigma=None, iterations=None, realSize=True, keep_all=False, verbose=False):
	"""
	maxima: pour specifier le ou les maxima (0=plus gros, 1=deuxieme plus gros, etc.) que l'on souhaite utiliser pour l'extraction du plan de symetrie.
			Si non specifie, on prend par defaut la liste des maxima de hauteur >= hauteur max * frac 
	frac: pour specifier la tolerance de hauteur que l'on se donne pour l'extraction des maxima en fonction de la hauteur maximale de l'histogramme des directions.
			Si non specifie, le programme par defaut prend la valeur 0.5 pour ce parametre.
	"""
    ###### Exemple de chaine d'execution de commandes shell #######
    # symmetryPlane ../BIN/bin_t0${i}_on_t099.inr -sphere ../HISTO/bin_t0${i}_on_t099_R15A32.inr -sigma 0.1 -weq ../SYM/planeEq_t0${i}_on_t099_max_0.txt -d 0 -dmin 1 
    #
    # diceMaximisation ../WAT/OUT_t0${i}_on_t099-wat.inr ../SYM/planeEq_t0${i}_on_t099_max_0_dmax.txt -n `cat ../SYM/planeEq_t0${i}_on_t099_max_0.txt | grep -v "#" | awk -F':' '{print $2}'` -delta 10 
    #
    # symmetryPlane ../BIN/bin_t0${i}_on_t099.inr -n `cat ../SYM/planeEq_t0${i}_on_t099_max_0_dmax.txt | grep new | awk -F':' '{print $2}'` -weq ../SYM/planeEq_t0${i}_on_t099_max_0_dmax_a.txt -d 10 -dmin 10 -p 0.1 

	#assert(os.path.exists(flo_file_bin) and os.path.exists(flo_file_hist))
	#symmetryPlane(flo_file_bin, flo_file_hist, equation_output=flo_file_sym_eq, trsf_output=flo_file_sym_alignment_trsf, plane_output=flo_file_sym_plane, 
	#			  maximum=maximum_ref, d=d_ref, dmin=dmin_ref, delta=delta_ref, p=p_ref, sigma=sigma_ref, realSize=True, verbose=verbose)

	# Histogramme directionnel
	if not os.path.exists(image_direction_histogram):
		assert(os.path.exists(image_bin))
		if not os.path.dirname(image_direction_histogram):
			try:
				os.makedirs(os.path.dirname(image_direction_histogram))
			except:
				print "Unable to create experiences path. Check the complete path..."
				raise
		directionHistogram(image_bin, image_direction_histogram, verbose=verbose)

	# Extraction de la liste des directions candidates

	
	if image_seg:
		# Histogramme maxima extraction
		candidates=directionHistogramMaxima(image_direction_histogram, file_out=None, maxima=maxima, frac=frac, verbose=False)

		n=candidates.shape[0]
		means_max=[]
		i_max=0
		tmp_eq_files=[]
		tmp_plane_files=[]
		tmp_distribution_files=[]
		tmp_trsf_files=[]
		# boucle sur les candidats
		for i in range(n):
			tmp_eq_files[i]='tmp_symmetryPlane_'+str(i)+'.eq'
			tmp_plane_files[i]=None
			tmp_distribution_files[i]=None
			tmp_trsf_files[i]=None
			if plane_output:
				tmp_plane_files[i]='tmp_symmetryPlane_'+str(i)+'.inr'
			if distribution_output:
				tmp_distribution_files[i]='tmp_symmetryPlane_'+str(i)+'.distribution'
			if trsf_output:
				tmp_trsf_files[i]='tmp_symmetryPlane_'+str(i)+'.trsf'
			symmetryPlane(image_bin, path_input_sphere=None, normal=candidates[i], equation_output=tmp_eq_files[i], plane_output=tmp_plane_files[i], distribution_output=tmp_distribution_files[i], trsf_output=tmp_trsf_files[i], 
                  d=d, dmin=dmin, delta=delta, p=p, sigma=sigma, realSize=realSize, iterations=iterations, verbose=verbose, lazy=True)
			d,m,s=dice(image_seg, symmetry=tmp_eq_files[i])
			if means_max==[] or m>max(means_max):
				i_max=i
			means_max[i]=m
		# Copie des fichiers correspondants au plan de symetrie du candidat optimal
		cmd="cp " + tmp_eq_files[i_max] + " " + 'tmp_symmetryPlane.eq'
		if verbose:
			print cmd
		os.system(cmd)
		if equation_output:
			cmd="cp " + tmp_eq_files[i_max] + " " + equation_output
			if verbose:
				print cmd
			os.system(cmd)
		if plane_output:
			copy(tmp_plane_files[i_max],plane_output, verbose=verbose)
		if distribution_output:
			cmd="cp " + tmp_distribution_files[i_max] + " " + distribution_output
			if verbose:
				print cmd
			os.system(cmd)
		if trsf_output:
			cmd="cp " + tmp_trsf_files[i_max] + " " + trsf_output
			if verbose:
				print cmd
			os.system(cmd)
		# Effacement des donnees temporaires
		if not keep_all:
			for i in range(n):
				cmd="rm " + tmp_eq_files[i]
				if tmp_plane_files[i]:
					cmd += " " + tmp_plane_files[i]
				if tmp_distribution_files[i]:
					cmd += " " + tmp_distribution_files[i]
				if tmp_trsf_files[i]:
					cmd += " " + tmp_trsf_files[i]
			if verbose:
				print cmd
			os.system(cmd)
	else:
		maximum=None
		if len(maxima)==1 and type(maxima[0])==int:
			maximum=maxima[0]
		symmetryPlane(image_bin, path_input_sphere=image_direction_histogram, normal=None, equation_output='tmp_symmetryPlane.eq', plane_output=plane_output, distribution_output=distribution_output, trsf_output=trsf_output, 
                  maximum=maximum, d=d, dmin=dmin, delta=delta, p=p, sigma=sigma, realSize=realSize, iterations=iterations, verbose=verbose, lazy=True)
		# Copie du fichier correspondants a l'equation du plan de symetrie
		if equation_output:
			cmd="cp tmp_symmetryPlane.eq " + equation_output
			if verbose:
				print cmd
			os.system(cmd)
	# Lecture de l'equation du plan de symetrie
	f=open("tmp_symmetryPlane.eq")
	n=[]
	for line in f:
		li=line.strip()
		#print li
		if not li.startswith('#'):
			#print li
			l=li.split()
			for value in l[1:5]:
				n.append(float(value))
	# Effacement du fichier d'equation temporaire
	cmd="rm tmp_symmetryPlane.eq" 
	if verbose:
		print cmd
	os.system(cmd)
	return n



def spatial_registration(ref_seg_post_file, flo_seg_post_file, ref_fused_file=None, flo_fused_file=None,  # Input images
	path_trsf_flo_ref=None, path_pairs_ref_flo=None, path_dices_ref_flo=None, path_residuals_ref_flo=None, # Outputs registration
	embryoKey_ref=None, embryoKey_flo=None, folder_tmp='WORKSPACE/', # Temporary files and paths
	init_ref=0.9, init_flo=0.9, realScale_ref=True, realScale_flo=True, # arguments for membrane_renforcement ref and flo
	sensitivity_ref=0.99, sensitivity_flo=0.99, # arguments for anisotropicHist ref and flo
	maximum_ref=None, d_ref=None, dmin_ref=None, delta_ref=None, p_ref=None, sigma_ref=None, # arguments for symmetryPlane ref
	maximum_flo=None, d_flo=None, dmin_flo=None, delta_flo=None, p_flo=None, sigma_flo=None, # arguments for symmetryPlane flo
	trsf_type='rigid', estimator='lts', lts_fraction=0.9, # arguments for planeRegistration
	background_ref = 1, background_flo = 1, # Background labels for planeRegistration
	keep_mem=False,
	keep_bin=False,
	keep_hist=False,
	keep_sym=False,
	keep_inter=False,
	verbose=False):
	'''
	# Temporary folders:

	folder_mem = folder_tmp + 'mem/'
	folder_bin = folder_tmp + 'bin/'
	folder_hist = folder_tmp + 'hist/'
	folder_sym = folder_tmp + 'sym/'
	folder_inter = folder_tmp + 'inter/'

	# Keep temporary folders ?

	keep_mem = False,
	keep_bin = False,
	keep_hist = False,
	keep_sym = False,
	keep_inter = False,

	# Parameters:

	ref_seg_post_file, flo_seg_post_file, ref_fused_file=None, flo_fused_file=None,  				# Input images
	path_trsf_flo_ref=None, path_pairs_ref_flo=None, path_dices_ref_flo=None, path_residuals_ref_flo=None, 		# Outputs registration
	embryoKey_ref=None, embryoKey_flo=None, folder_tmp='WORKSPACE/', 			# Temporary files and paths
	init_ref=0.9, init_flo=0.9, realScale_ref=True, realScale_flo=True, 		# arguments for membrane_renforcement ref and flo
	sensitivity_ref=0.99, sensitivity_flo=0.99, 		# arguments for anisotropicHist ref and flo
	maximum_ref=None, d_ref=None, dmin_ref=None, delta_ref=None, p_ref=None, sigma_ref=None, 		# arguments for symmetryPlane ref
	maximum_flo=None, d_flo=None, dmin_flo=None, delta_flo=None, p_flo=None, sigma_flo=None, 		# arguments for symmetryPlane flo
	trsf_type='rigid', estimator='lts', lts_fraction=0.9,		 				# arguments for planeRegistration
	background_ref = 1, background_flo = 1, 			# Background labels for planeRegistration

	# Verbosity:

	verbose=False
	'''

	delete_workspace = not (keep_mem or keep_bin or keep_hist or keep_sym or keep_inter)
	
	### Ref ###

	if not embryoKey_ref:
		embryoKey_ref=ref_seg_post_file.split(os.path.sep)[-1]
		embryoKey_ref=embryoKey_ref.split('-')[1]+'_'+embryoKey_ref.split('_')[-1].split('.')[0]

	### Flo ###

	if not embryoKey_flo:
		embryoKey_flo=flo_seg_post_file.split(os.path.sep)[-1]
		embryoKey_flo=embryoKey_flo.split('-')[1]+'_'+embryoKey_flo.split('_')[-1].split('.')[0]
	
	### TEMPORARY PATH VARIABLES ###

	folder_mem = folder_tmp + 'mem/'
	folder_bin = folder_tmp + 'bin/'
	folder_hist = folder_tmp + 'hist/'
	folder_sym = folder_tmp + 'sym/'
	folder_inter = folder_tmp + 'inter/'

	### Flo ###

	flo_file_mem_prefix = folder_mem + embryoKey_flo + "_mem"
	flo_file_bin_prefix = folder_bin + embryoKey_flo + "_bin"
	flo_file_bin = flo_file_bin_prefix+'.inr'
	flo_file_theta = flo_file_bin_prefix+'.theta.inr'
	flo_file_phi = flo_file_bin_prefix+'.phi.inr'
	flo_file_hist = folder_hist + embryoKey_flo + "_directionHistogram.inr"
	flo_file_sym_prefix = folder_sym + embryoKey_flo + "_sym"
	flo_file_sym_eq = flo_file_sym_prefix + '.eq'
	flo_file_sym_plane = flo_file_sym_prefix + '.plane.inr'
	flo_file_sym_alignment_trsf = flo_file_sym_prefix + '.trsf'

	### Ref ###

	ref_file_mem_prefix = folder_mem + embryoKey_ref + "_mem"
	ref_file_bin_prefix = folder_bin + embryoKey_ref + "_bin"
	ref_file_bin = ref_file_bin_prefix+'.inr'
	ref_file_theta = ref_file_bin_prefix+'.theta.inr'
	ref_file_phi = ref_file_bin_prefix+'.phi.inr'
	ref_file_hist = folder_hist + embryoKey_ref + "_directionHistogram.inr"
	ref_file_sym_prefix = folder_sym + embryoKey_ref + "_sym"
	ref_file_sym_eq = ref_file_sym_prefix + '.eq'
	ref_file_sym_plane = ref_file_sym_prefix + '.plane.inr'
	ref_file_sym_alignment_trsf = ref_file_sym_prefix + '.trsf'

	### Registration ###

	ref_flo_file_trsf = folder_inter + embryoKey_flo + '_' + embryoKey_ref + '.trsf'
	ref_flo_file_pairs = folder_inter + embryoKey_ref + '_' + embryoKey_flo + '.pairs'
	ref_flo_file_res = folder_inter + embryoKey_ref + '_' + embryoKey_flo + '_residuals.m'
	ref_flo_file_dices = folder_inter + embryoKey_ref + '_' + embryoKey_flo + '.dices'


	### STUFF ###

	if not os.path.isdir(folder_tmp):
		try:
			os.makedirs(folder_tmp)
		except:
			print "Unable to create experiences path. Check the complete path..."
			raise

	# Fonction de rehaussement de membranes :
	if not (os.path.exists(ref_file_bin) or os.path.exists(ref_file_mem_prefix+'.ext.inr')):
		assert(os.path.exists(ref_fused_file))
		if not os.path.isdir(folder_mem):
			try:
				os.makedirs(folder_mem)
			except:
				print "Unable to create experiences path. Check the complete path..."
				raise
		membrane_renforcement(ref_fused_file, prefix_output=ref_file_mem_prefix, init=init_ref, realScale=realScale_ref, verbose=verbose)

	if not (os.path.exists(flo_file_bin) or os.path.exists(flo_file_mem_prefix+'.ext.inr')):
		assert(os.path.exists(flo_fused_file))
		if not os.path.isdir(folder_mem):
			try:
				os.makedirs(folder_mem)
			except:
				print "Unable to create experiences path. Check the complete path..."
				raise
		membrane_renforcement(flo_fused_file, prefix_output=flo_file_mem_prefix, init=init_flo, realScale=realScale_flo, verbose=verbose)

	# Binarisation de membranes :
	if not os.path.exists(ref_file_bin):
		assert(os.path.exists(ref_file_mem_prefix+'.ext.inr'))
		if not os.path.isdir(folder_bin):
			try:
				os.makedirs(folder_bin)
			except:
				print "Unable to create experiences path. Check the complete path..."
				raise
		anisotropicHist(ref_file_mem_prefix+'.ext.inr', path_output=ref_file_bin, sensitivity=sensitivity_ref, verbose=verbose)

	if not os.path.exists(flo_file_bin):
		assert(os.path.exists(flo_file_mem_prefix+'.ext.inr'))
		if not os.path.isdir(folder_bin):
			try:
				os.makedirs(folder_bin)
			except:
				print "Unable to create experiences path. Check the complete path..."
				raise
		anisotropicHist(flo_file_mem_prefix+'.ext.inr', path_output=flo_file_bin, sensitivity=sensitivity_flo, verbose=verbose)

	if not keep_mem:
		cmd='rm -rf ' + folder_mem
		if verbose:
			print cmd
		os.system(cmd)

	# Calcul de l'histogramme directionnel
	if not os.path.exists(ref_file_hist):
		assert(os.path.exists(ref_file_bin))
		if not os.path.isdir(folder_hist):
			try:
				os.makedirs(folder_hist)
			except:
				print "Unable to create experiences path. Check the complete path..."
				raise
		directionHistogram(ref_file_bin, ref_file_hist, verbose=verbose)
	if not os.path.exists(flo_file_hist):
		assert(os.path.exists(flo_file_bin))
		if not os.path.isdir(folder_hist):
			try:
				os.makedirs(folder_hist)
			except:
				print "Unable to create experiences path. Check the complete path..."
				raise
		directionHistogram(flo_file_bin, flo_file_hist, verbose=verbose)

	# Calcul de l'equation du plan de symetrie
	assert(os.path.exists(ref_file_bin))
	assert(os.path.exists(flo_file_bin))
	assert(os.path.exists(ref_file_hist))
	assert(os.path.exists(flo_file_hist))
	if not os.path.isdir(folder_sym):
		try:
			os.makedirs(folder_sym)
		except:
			print "Unable to create experiences path. Check the complete path..."
			raise
	if not os.path.exists(flo_file_sym_eq):
		assert(os.path.exists(flo_file_bin) and os.path.exists(flo_file_hist))
		symmetryPlane(flo_file_bin, flo_file_hist, equation_output=flo_file_sym_eq, trsf_output=flo_file_sym_alignment_trsf, plane_output=flo_file_sym_plane, 
			maximum=maximum_ref, d=d_ref, dmin=dmin_ref, delta=delta_ref, p=p_ref, sigma=sigma_ref, realSize=True, verbose=verbose)
	if not os.path.exists(ref_file_sym_eq):
		assert(os.path.exists(ref_file_bin) and os.path.exists(ref_file_hist))
		symmetryPlane(ref_file_bin, ref_file_hist, equation_output=ref_file_sym_eq, trsf_output=ref_file_sym_alignment_trsf, plane_output=ref_file_sym_plane, 
			maximum=maximum_flo, d=d_flo, dmin=dmin_flo, delta=delta_flo, p=p_flo, sigma=sigma_flo, realSize=True, verbose=verbose)

	if not keep_bin:
		cmd='rm -rf ' + folder_bin
		if verbose:
			print cmd
		os.system(cmd)

	if not keep_hist:
		cmd='rm -rf ' + folder_hist
		if verbose:
			print cmd
		os.system(cmd)

	# Registration
	assert(os.path.exists(ref_seg_post_file))
	assert(os.path.exists(flo_seg_post_file))
	assert(os.path.exists(ref_file_sym_eq))
	assert(os.path.exists(flo_file_sym_eq))
	if not os.path.isdir(folder_inter):
		try:
			os.makedirs(folder_inter)
		except:
			print "Unable to create experiences path. Check the complete path..."
			raise
	planeRegistration(ref_seg_post_file , flo_seg_post_file, ref_file_sym_eq, flo_file_sym_eq,
						  path_trsf_ref_to_flo=ref_flo_file_trsf, path_residuals=ref_flo_file_res, 
						  path_dices=ref_flo_file_dices, path_pairs_ref_flo=ref_flo_file_pairs, 
	                      background_ref = background_ref, background_flo = background_flo,
	                      trsf_type=trsf_type, estimator=estimator, lts_fraction=lts_fraction, verbose=verbose, lazy=True)

	if not keep_sym:
		cmd='rm -rf ' + folder_sym
		if verbose:
			print cmd
		os.system(cmd)


	if not path_trsf_flo_ref:
		path_trsf_flo_ref = os.path.basename(ref_flo_file_trsf)
	if path_trsf_flo_ref:
		cmd="cp " + ref_flo_file_trsf + ' ' + path_trsf_flo_ref
		if verbose:
			print cmd
		os.system(cmd)
	
	if not path_pairs_ref_flo:
		path_pairs_ref_flo = os.path.basename(ref_flo_file_pairs)
	if path_pairs_ref_flo:
		cmd="cp " + ref_flo_file_pairs + ' ' + path_pairs_ref_flo
		if verbose:
			print cmd
		os.system(cmd)

	if not path_dices_ref_flo:
		path_dices_ref_flo = os.path.basename(ref_flo_file_dices)
	if path_dices_ref_flo:
		cmd="cp " + ref_flo_file_dices + ' ' + path_dices_ref_flo
		if verbose:
			print cmd
		os.system(cmd)

	if not path_residuals_ref_flo:
		path_residuals_ref_flo = os.path.basename(ref_flo_file_res)
	if path_residuals_ref_flo:
		cmd="cp " + ref_flo_file_res + ' ' + path_residuals_ref_flo
		if verbose:
			print cmd
		os.system(cmd)




	if not keep_inter:
		cmd='rm -rf ' + folder_inter
		if verbose:
			print cmd
		os.system(cmd)

	if delete_workspace and folder_tmp:
		cmd='rm -rf ' + folder_tmp
		if verbose:
			print cmd
		os.system(cmd)

	return path_trsf_flo_ref, path_pairs_ref_flo, path_dices_ref_flo, path_residuals_ref_flo


def to_256_bits_key(n, background=[0,1], key_no_correspondences='no_correspondences', lut=None):
	'''
	Returns :
		0 if n in background
		255 if n is equal to key_no_correspondences or not lut.has_key(n) or not n in lut
		((n - d )% 255 ) + 1 if key_no_correspondences is None and lut is None
		((n - d )% 254 ) + 1 otherwise (default)
		where n <- lut[n] if lut is of type dict (may be of type list instead if one does not expect a lut but a set of cells with existing correspondences)
		(default = None), otherwise original parameter value is kept.
	This method is appropriated for Fiji LUT "glasbeyBGD"
	'''
	modulo = 254
	if key_no_correspondences == None and not lut:
		modulo = 255
	else:
		if type(n) == type(key_no_correspondences) and n == key_no_correspondences:
			return 255

	assert(type(n)==int)
	if n in background:
		return 0
	if lut:
		if type(lut)==dict:
			if not lut.has_key(n):
				return 255
			n = lut[n]
		else:
			if not n in lut:
				return 255

	d = len([i for i in background if i < n])# dividende
	l = ((n - d )% modulo ) + 1
	return l


def segmentationRelabellingAtTime(path_segmentation, path_output, lineage, time_point, reverse_lut=None, lut=None, trsf=None, template=None, voxelsize=None, iso=None, dimensions=None, backgroundLabels=[0,1], visu=True, verbose=False):

	'''
	Builds the segmented image path_output encoded in 8 bits from input segmentation path_segmentation, with a relabellization and, if visu=True, a thin erosion of the cells for vizualization 
	(the user may choose the 'glasbeyBGD' LUT in Fiji for path_output vizualisation). The relabelling is processed with respect to the lineage relabelling, to lut (or reverse_lut), and to time_point.

	Parameters:

		lut: if specified, do not specify reverse_lut. 
		     The lut option can be of two types:
		     	- list (or set or tuple): in that case, associates value 255 if the relabelled value l of path_segmentation is not in lut. Otherwise, associate the "to_256_bits_key" conversion of l.
		     	- dict: in that case, associates value 255 if the relabelled value l of path_segmentation is not in lut.keys(). Otherwise, associate the "to_256_bits_key" conversion of lut[l].
		     This parameter may be useful when the user wants to synchronize a reference image considering a floating sequence registration onto the reference and when the reference-to-floating correspondences dictionary is known.
		     In that case of use, the lut parameter should be set with the list given by correspondences.keys(), so that cells with no correspondences will be automatically set to output value 255.
		
		reverse_lut: if specified, do not specify lut. Must be a dict type and is equivalent than to set lut option with the value {v:k for k,v in reverse_luts.items()} (ie the reverse version of reverse_lut).
					 May be a useful parameter when the user wants to synchronise a floating image onto a reference sequence and when the reference-to-floating correspondences dictionary is known. 
					 In that case of use, the reverse_lut parameter should be set with the correspondences dictionary.

	Optional parameters:

		trsf (default=None):  if None or False, path_output will have the same field of view as path_segmentation.
							  If trsf is a path for a transformation or a transformation matrix 4*4 (np.array or list) that transforms spatially a target field into path_segmentation (with convention trsf = T_segmentation<-target) 
							  then the function applies this transformation to path_segmentation (with optional dimensions, voxelsize, iso and template parameters).

		voxelsize (default=None): forces path_output voxelsizes to the given values (must be a list or tuple of length 3) (may be used independently to trsf parameter).

		iso (default=None): forces path_output to have isotropic voxel size of given value (float type) (may be used independently to trsf parameter).

		dimensions (default=None): forces the path_output dimensions (must be an integer list or tuple of length 3) (may be used independently to trsf parameter).

		template (default=None): path to a template image that forces the path_output dimensions and voxelsize to correspond to the tmeplate image (must be a str type). If used, do not use iso, voxelsize neither dimensions parameters.

		backgroundLabels (default=[0,1]): labels of the input image to be set to 0 in the output image

		visu (default=True): option for vizualisation. If True, then a thin erosion of original cells borders is processed to see well cell borders in the output image.

		verbose (defaut=False): option for verbosity of the function.
	'''
	assert os.path.exists(path_segmentation)
	if template:
		assert os.path.exists(template)

	if not os.path.dirname(path_output):
		try:
			os.makedirs(os.path.dirname(path_output))
		except:
			print "Unable to create output path. Check the complete path..."
			raise

	if type(trsf) != np.ndarray and trsf == None:
		trsf=False

	out_original_labels={}
	rev_lut=lineage.reverse_luts[time_point]

	if reverse_lut:
		assert not lut 
		lut={v:k for k,v in reverse_lut.items()}

	if not lut:
		for k, v in rev_lut.items():
			k_256_bits=to_256_bits_key(k, backgroundLabels)
			if not out_original_labels.has_key(k_256_bits):
				out_original_labels[k_256_bits]=set()
			out_original_labels[k_256_bits].add(v)
	else:
		for k, v in rev_lut.items():
			k_256_bits=to_256_bits_key(k, backgroundLabels, lut=lut)
			if not out_original_labels.has_key(k_256_bits):
				out_original_labels[k_256_bits]=set()
			out_original_labels[k_256_bits].add(v)

	# We add the background correspondence...
	if not out_original_labels.has_key(0):
		out_original_labels[0]=set()
	for v in backgroundLabels:
		out_original_labels[0].add(v)
		out_original_labels[0].add(v)

	txt=""
	for k, v in out_original_labels.items():
		for val in v:
			txt += "%d %d\n"%(val,k)

	path_tmp_lut = "tmp_segmentationRelabellingAtTime_t%03d.lut"%time_point
	if verbose:
		print "Writing temporary file %s"%path_tmp_lut
	f=open(path_tmp_lut, "w")
	f.write(txt)
	f.close()




	tmp=path_output + ".tmp_segmentationRelabellingAtTime_t%03d.inr"%time_point

	trsf_tmp_file='tmp_segmentationRelabellingAtTime_t%03d.trsf'%time_point
	# Check if need to compute the trsf (->set cpt_trsf) or if it is a filename (-> copy to temporary file)
	if type(trsf)!=bool:
		if type(trsf)==str:
			assert os.path.exists(trsf), "Error: unexpected value '%s' for trsf parameter (file not found). See help."%trsf
			cmd="cp %s %s"%(trsf, trsf_tmp_file)
			if verbose:
				print cmd
			os.system(cmd)
	# If trsf is a tabular, then write the temporary trsf file
	if type(trsf)!=bool:
		if type(trsf)==list:
			trsf=np.array(trsf)
		if type(trsf)==np.ndarray:
			f=open(trsf_tmp_file, "w")
			f.write(str(trsf).replace('[','').replace(']',''))
			f.close()


	if template:
		assert (not dimensions and not voxelsize and not iso) 
	assert type(trsf) != bool or not trsf


	if voxelsize or dimensions or iso or template or type(trsf) != bool :
		if type(trsf) != bool: 
			assert os.path.exists(trsf_tmp_file), "Unexpected problem while setting the transformation T_flo<-ref."
			apply_trsf(path_segmentation, path_trsf=trsf_tmp_file, path_output=path_output, nearest=True, template=template, voxelsize=voxelsize, dimensions=dimensions, iso=iso, verbose=verbose)
		else:
			apply_trsf(path_segmentation, path_output=path_output, nearest=True, template=template, voxelsize=voxelsize, dimensions=dimensions, iso=iso, verbose=verbose)

	if visu:
		if voxelsize or dimensions or iso or type(trsf) != bool:
			labelBorders(path_output, tmp, verbose)
		else:
			labelBorders(path_segmentation, tmp, verbose)
		Logic(tmp, tmp, Mode='inv', verbose=verbose)


	if voxelsize or dimensions or iso or type(trsf) != bool:
		fuseLabelsWithLUT(path_output, path_output, path_tmp_lut, U8=True, verbose=verbose)
	else:
		fuseLabelsWithLUT(path_segmentation, path_output, path_tmp_lut, U8=True, verbose=verbose)

	if os.path.exists(path_tmp_lut):
		cmd='rm '+path_tmp_lut
		if verbose:
			print cmd
		os.system(cmd)

	if visu:
		Logic(tmp, path_output, path_output, Mode='mask', verbose=verbose)
		cmd="rm " + tmp
		if verbose:
			print cmd
		os.system(cmd)

	# RM TEMPORARY FILE
	if type(trsf) != bool:
		if os.path.exists(trsf_tmp_file):
			cmd="rm %s"%trsf_tmp_file
			if verbose:
				print cmd
			os.system(cmd)


def associateSegmentationsAtTimesFromLineageCorrespondences(path_ref, path_flo, path_ref_out, path_flo_out, 
															lineage_ref, lineage_flo, correspondences, time_ref, time_flo, 
															trsf=False, estimator='lts', lts_fraction=1.0,
															#template=None,
															voxelsize=None, iso=None, dimensions=None,
															backgroundLabels=[0,1], visu=True, verbose=False):
	'''
	Builds the segmented images path_ref_out and path_flo_out from input segmentations path_ref and path_flo, with a relabellization and thin erosion of the cells for vizualization 
	(the user may choose the 'glasbeyBGD' LUT in Fiji). The relabelling is processed with respect to the lineage_ref, lineage_flo, correspondences, time_ref and time_flo parameters.
	Optional parameters:
		trsf (default=False): if False, path_flo_out will have the same dimensions as path_flo.
							  If True or 'rigid', the function will compute automatically a rigid transformation that registers lineage_ref with lineage_flo at given time_ref and time_flo,
							  under the constraint of given correspondences, and apply the transformation to path_flo_out so that dimension(path_flo_out)=dimension(path_ref_out).
							  If 'affine', the same as if True, but with an affine transformation instead of a rigid.
							  In these two cases, a transformation will be automatically computed with the parameters 'estimator' and 'lts_fraction' that are related to the pointCloudRegistration function.
							  If trsf is a path for a transformation or a transformation matrix 4*4 (np.array or list) that registers spatially path_flo into path_ref (with convention trsf = T_flo<-ref) 
							  then the function applies this transformation to path_flo so that dimension(path_flo_out)=dimension(path_ref_out).
		voxelsize (default=None): forces path_ref_out and path_flo_out voxelsizes to the given values (must be a list or tuple of length 3).
		iso (default=None): forces path_ref_out and path_flo_out to have isotropic voxel size of given value (float type).
		dimensions (default=None): forces the path_ref_out and path_flo_out output dimensions (must be an integer list or tuple of length 3).
		backgroundLabels (default=[0,1]): labels of the input images to be set to 0 in the output images
		visu (default=True): option for vizualisation. If True, then a thin erosion of original cells borders is processed to see well cell borders in the output images.
		verbose (defaut=False): option for verbosity of the function.
	'''
	

	assert(os.path.exists(path_ref))
	assert(os.path.exists(path_flo))
	key_no_correspondences='no_correspondences'
	
	if type(trsf) != np.ndarray and trsf == None:
		trsf=False

	def get_correspondence(lineage_ref, lineage_flo, time_ref, time_flo, correspondences, key_no_correspondences='no_correspondences'):
		'''
		Returns a dict structure with each item as follows:
		keys, values: 
			- key_no_correspondences, (set_ref, set_flo) :
				key for labels that belong to lineages at given times and that do not exist in correspondences dictionary
				set_ref for relabellings that belong to lineages_ref at time_ref and not belong to correspondences (ie correspondences.has_key(l) returns False for each element l of set_ref)
				set_flo for relabellings that belong to lineages_flo at time_flo and not belong to correspondences (ie l in correspondences.values() returns False for each element l of set_flo)
			- label : (set_ref, set_flo) :
				set_ref for set of relabellings from lineage_ref that is or are associated to the set of relabellings set_flo from lineage_flo (all elements of set_ref are existing cells of lineage_ref at time_ref)
				set_flo for set of relabellings from lineage_flo that is or are associated to the set of relabellings set_ref from lineage_ref (all elements of set_flo are existing cells of lineage_flo at time_flo)
				A property that must be respected is that at least one of the two sets is a singleton (ie has one and only one element). The other set can have 0, 1 or more elements.
				The key label is a label which value is the "greatest common ancestor" of elements of set_ref so that all the elements of set_flo are associated to this ancestor or its descendents.
					
		'''
		ref_correspondences=correspondences.keys()
		flo_correspondences=correspondences.values()


		relabels_ref = lineage_ref.reverse_luts[time_ref].keys()
		relabels_flo = lineage_flo.reverse_luts[time_flo].keys()

		not_corresponding_ref = set(relabels_ref).difference(ref_correspondences).intersection(relabels_ref)
		not_corresponding_flo = set(relabels_flo).difference(flo_correspondences).intersection(relabels_flo)


		out={key_no_correspondences:(not_corresponding_ref, not_corresponding_flo)} # ref and flo sets of cells that have to be visualized as cells with no correspondences


		relabels_ref = list(set(relabels_ref).difference(not_corresponding_ref)) # we remove from relabels_ref the cells with no correspondences
		#relabels_ref.sort()
		relabels_flo = list(set(relabels_flo).difference(not_corresponding_flo)) # we remove from relabels_flo the cells with no correspondences
		#relabels_flo.sort()

		relabels_flo_as_ancestors=set()
		scanned_flo_set=set()

		for relabel_ref in relabels_ref:
			assert relabel_ref in ref_correspondences
			# Cas ou relabel_ref existe dans le dict 'correspondences'
			corres_flo=correspondences[relabel_ref]
			if corres_flo in relabels_flo:
				# Cas simple : correspondance 1 a 1 trouvee
				assert not out.has_key(relabel_ref)
				out[relabel_ref]= ({relabel_ref},{corres_flo})
				scanned_flo_set.add(corres_flo)
			else:
				# Cas possibles :
				# 1. relabels_flo contient des descendants de corres_flo -> on associe relabel_ref a tous les descendants de corres_flo qui appartiennent au dictionnaire 'correspondences'
				# 2. relabels_flo contient un ancetre de corres_flo -> on traite le cas dans le sens "inverse" en cherchant tous les correspondants de l'ancetre dans relabels_ref qui appartiennent au dictionnaire 'correspondences'
				#	 -> on procede a la recherche "inverse" afin de determiner les correspondances de l'ancetre dans relabels_ref (boucle finale sur relabels_flo_loop)
				# 3. relabels_flo ne contient aucune cellule de la lignee de corres_flo -> out[relabel_ref]=(relabel_ref, set())
				
				corres_flo_descendents = lineage_flo.relabelled_descendents(corres_flo)
				corres_flo_ancestors = lineage_flo.relabelled_ancestors(corres_flo)

				# Test cas 1:
				intersection=set(relabels_flo).intersection(corres_flo_descendents)
				if intersection:
					# Cas 1:
					assert(not set(relabels_flo).intersection(corres_flo_ancestors)) # cohabitation logiquement impossible entre descendants et ancetres d'une meme cellule a un instant donne du developpement d'un embryon

					corresponding=intersection.intersection(correspondences.values())# extraction des descendants du label en correspondance avec relabel_ref qui appartiennent au dictionnaire 'correspondences' -> out key=relabel_ref
					not_corresponding=intersection.difference(correspondences.values()).intersection(intersection) # extraction des descendants du label en correspondance avec relabel_ref qui n'appartiennent pas au dictionnaire 'correspondences' -> out key=key_no_correspondences
					assert not not_corresponding
					assert not out.has_key(relabel_ref)
					out[relabel_ref] = ({relabel_ref},corresponding) 
					scanned_flo_set=scanned_flo_set.union(corresponding)
				else:
					# Cas 2 ou 3:
					# Test cas 2:
					intersection=set(relabels_flo).intersection(corres_flo_ancestors)
					assert(len(intersection)<=1) # on ne devrait pas avoir une cohabitation de plusieurs ancetres d'un unique label a un pas de temps donne
					if intersection:
						# Cas 2: 
						intersection_value=intersection.pop()
						assert intersection_value in flo_correspondences
						relabels_flo_as_ancestors.add(intersection_value)
					else:
						# Cas 3:
						assert not out.has_key(relabel_ref)
						out[relabel_ref]=({relabel_ref}, set())

		rev_correspondences={v:k for k,v in correspondences.items()}
		for relabel_flo in relabels_flo_as_ancestors:
			# Cette boucle concerne tous les labels de relabels_flo qui s'associent a un ancetre de cellules existentes dans relabels_ref...
			assert rev_correspondences.has_key(relabel_flo)
			corres_ref=rev_correspondences[relabel_flo]
			corres_ref_descendents = lineage_ref.relabelled_descendents(corres_ref)
			intersection=set(relabels_ref).intersection(corres_ref_descendents)
			assert len(intersection)>0
			corresponding=intersection.intersection(rev_correspondences.values())# extraction des descendants du label en correspondance avec relabel_flo qui appartiennent au dictionnaire 'rev_correspondences' -> out key=corres_ref
			not_corresponding=intersection.difference(rev_correspondences.values()).intersection(intersection) # extraction des descendants du label en correspondance avec relabel_flo qui n'appartiennent pas au dictionnaire 'rev_correspondences' -> out key=key_no_correspondences (theoriquement, deja ajoute a ce stade)
			assert len(not_corresponding) == len(not_corresponding.intersection(out[key_no_correspondences][0]))
			assert not out.has_key(corres_ref)
			out[corres_ref]=(corresponding, {relabel_flo})
			scanned_flo_set.add(relabel_flo)

		# Etape finale : aucun element de relabels_ref ni de relabels_flo ne doit avoir ete omis
		last_loop_flo = scanned_flo_set.difference(relabels_flo) 
		for relabel_flo in last_loop_flo: # labels de relabels_flo qui n'ont pas encore ete scannes, i.e. qui n'ont soit pas de correspondants dans rev_correspondences (->key_no_correspondences), soit pas de correspondants dont la lignee appartient a relabels_ref (->key=rev_correspondences[relabel_flo])
			assert relabel_flo in flo_correspondences
			assert not out.has_key(rev_correspondences[relabel_flo])
			out[rev_correspondences[relabel_flo]]=(set(), {relabel_flo})

		return out

	# extraction des correspondences entre cellules de reference et flottantes avec nouvel etiquetage
	out_relabelled=get_correspondence(lineage_ref, lineage_flo, time_ref, time_flo, correspondences, key_no_correspondences)

	out_original_labels={}
	rev_lut_ref=lineage_ref.reverse_luts[time_ref]
	rev_lut_flo=lineage_flo.reverse_luts[time_flo]

	for k, (v_ref, v_flo) in out_relabelled.items():
		k_256_bits=to_256_bits_key(k, backgroundLabels, key_no_correspondences)
		if not out_original_labels.has_key(k_256_bits):
			out_original_labels[k_256_bits]=(set(), set())
		for relabel_ref in v_ref:
			assert(rev_lut_ref.has_key(relabel_ref))
			out_original_labels[k_256_bits][0].add(rev_lut_ref[relabel_ref])
		for relabel_flo in v_flo:
			assert(rev_lut_flo.has_key(relabel_flo))
			out_original_labels[k_256_bits][1].add(rev_lut_flo[relabel_flo])

	# We add the background correspondence...
	if not out_original_labels.has_key(0):
		out_original_labels[0]=(set(), set())
	for v in backgroundLabels:
		out_original_labels[0][0].add(v)
		out_original_labels[0][1].add(v)

	txt_ref=""
	txt_flo=""
	for k, (v_ref, v_flo) in out_original_labels.items():
		for val in v_ref:
			txt_ref += "%d %d\n"%(val,k)
		for val in v_flo:
			txt_flo += "%d %d\n"%(val,k)

	path_tmp_lut_ref = "tmp_associateSegmentationsAtTimesFromLineageCorrespondences_t%03d_t%03d_ref.lut"%(time_ref,time_flo)
	path_tmp_lut_flo = 'tmp_associateSegmentationsAtTimesFromLineageCorrespondences_t%03d_t%03d_flo.lut'%(time_ref,time_flo)
	if verbose:
		print "Wrinting temporary file %s"%path_tmp_lut_ref
	f=open(path_tmp_lut_ref, "w")
	f.write(txt_ref)
	f.close()
	if verbose:
		print "Wrinting temporary file %s"%path_tmp_lut_flo
	g=open(path_tmp_lut_flo, "w")
	g.write(txt_flo)
	g.close()


	tmp_ref=path_ref_out + ".tmp_t%03d_t%03d.inr"%(time_ref,time_flo)
	tmp_flo=path_flo_out + ".tmp_t%03d_t%03d.inr"%(time_ref,time_flo)


	cpt_trsf=None
	trsf_tmp_file='tmp_associateSegmentationsAtTimesFromLineageCorrespondences_flo_ref.trsf'
	# Check if need to compute the trsf (->set cpt_trsf) or if it is a filename (-> copy to temporary file)
	if type(trsf)==bool and trsf:
		cpt_trsf='rigid' 
	if type(trsf)!=bool:
		if type(trsf)==str:
			if trsf=='rigid':
				cpt_trsf='rigid'
			if trsf=='affine':
				cpt_trsf='affine'
			else:
				assert os.path.exists(trsf), "Error: unexpected value '%s' for trsf parameter (file not found). See help."%trsf
				cmd="cp %s %s"%(trsf, trsf_tmp_file)
				if verbose:
					print cmd
				os.system(cmd)
	# If need to compute trsf (->set trsf as a np.array)
	if cpt_trsf:
		point_cloud_ref=lineage_ref.relabelled_barycenters_at_time(time_ref)
		point_cloud_flo=lineage_flo.relabelled_barycenters_at_time(time_flo)
		if type(backgroundLabels)==list or type(backgroundLabels)==tuple or type(backgroundLabels)==set:
			for v in backgroundLabels:
				point_cloud_ref.pop(v, None)
		if type(backgroundLabels)==int:
			point_cloud_ref.pop(backgroundLabels, None)

		out=pointCloudRegistration( point_cloud_ref, point_cloud_flo, correspondences, skip_not_found=True, 
						trsf_type=cpt_trsf, estimator=estimator, lts_fraction=lts_fraction, 
						lazy=False, verbose=verbose)
		trsf=np.array(out['trsf']) 
	# If trsf is a tabular, then write the temporary trsf file
	if type(trsf)!=bool:
		if type(trsf)==list:
			trsf=np.array(trsf)
		if type(trsf)==np.ndarray:
			f=open(trsf_tmp_file, "w")
			f.write(str(trsf).replace('[','').replace(']',''))
			f.close()


	if voxelsize or dimensions or iso:
		assert (not dimensions and not voxelsize) or (not dimensions and not iso) or (not voxelsize and not iso) # not the two options at the same time
		apply_trsf(path_ref, path_output=path_ref_out, nearest=True, voxelsize=voxelsize, dimensions=dimensions, iso=iso, verbose=verbose)
		#apply_trsf(path_flo, path_output=path_flo_out, nearest=True, voxelsize=voxelsize, dimensions=dimensions, iso=iso, verbose=verbose)

	#Checks the existence of trsf_tmp_file if trsf was mentionned
	assert type(trsf) != bool or not trsf
	if type(trsf) != bool:
		assert os.path.exists(trsf_tmp_file), "Unexpected problem while setting the transformation T_flo<-ref."
		if voxelsize or dimensions or iso:
			apply_trsf(path_flo, path_trsf=trsf_tmp_file, path_output=path_flo_out, template=path_ref_out, nearest=True, verbose=verbose)
		else:
			apply_trsf(path_flo, path_trsf=trsf_tmp_file, path_output=path_flo_out, template=path_ref, nearest=True, verbose=verbose)
	else:
		if voxelsize or dimensions or iso:
			apply_trsf(path_flo, path_output=path_flo_out, nearest=True, voxelsize=voxelsize, dimensions=dimensions, iso=iso, verbose=verbose)

	if visu:
		if voxelsize or dimensions or iso:
			labelBorders(path_ref_out, tmp_ref, verbose)
			#cmd=path_labelborders + " " + path_ref_out + ' ' + tmp_ref 
		else:
			labelBorders(path_ref, tmp_ref, verbose)
		if type(trsf) != bool or voxelsize or dimensions or iso:
			labelBorders(path_flo_out, tmp_flo, verbose)
		else:
			labelBorders(path_flo, tmp_flo, verbose)
		Logic(tmp_ref, tmp_ref, Mode='inv', verbose=verbose)
		Logic(tmp_flo, tmp_flo, Mode='inv', verbose=verbose)


	if voxelsize or dimensions or iso:
		fuseLabelsWithLUT(path_ref_out, path_ref_out, path_tmp_lut_ref, U8=True, verbose=verbose)
	else:
		fuseLabelsWithLUT(path_ref, path_ref_out, path_tmp_lut_ref, U8=True, verbose=verbose)
	if type(trsf) != bool:
		fuseLabelsWithLUT(path_flo_out, path_flo_out, path_tmp_lut_flo, U8=True, verbose=verbose)
	else:
		fuseLabelsWithLUT(path_flo, path_flo_out, path_tmp_lut_flo, U8=True, verbose=verbose)

	if os.path.exists(path_tmp_lut_ref):
		cmd='rm '+path_tmp_lut_ref
		if verbose:
			print cmd
		os.system(cmd)
	if os.path.exists(path_tmp_lut_flo):
		cmd='rm '+path_tmp_lut_flo
		if verbose:
			print cmd
		os.system(cmd)

	if visu:
		Logic(tmp_ref, path_ref_out, path_ref_out, Mode='mask', verbose=verbose)
		Logic(tmp_flo, path_flo_out, path_flo_out, Mode='mask', verbose=verbose)
		cmd="rm " + tmp_ref + ' ' + tmp_flo
		if verbose:
			print cmd
		os.system(cmd)

	# RM TEMPORARY FILE
	if type(trsf) != bool:
		if os.path.exists(trsf_tmp_file):
			cmd="rm %s"%trsf_tmp_file
			if verbose:
				print cmd
			os.system(cmd)


def associateRefFloLabels(path_ref, path_flo, path_pairs_ref_flo, path_ref_out, path_flo_out, path_labels_out=None, path_trsf_flo_ref=None, zeros=True, backgroundLabels=[0,1], ref=False, force=False, visu=False, verbose=False):    
  '''
   path_ref : correspond au parametre path_label_ref de la methode planeRegistration
   path_flo : correspond au parametre path_label_flo de la methode planeRegistration
   path_pairs_ref_flo : fichier de correspondances label-ref / label-flo (tel qu'en sortie de planeRegistration / pointCloudRegistration path_pairs_ref_flo %s) (accepte aussi les dictionnaires avec les items de correspondences sous forme de key=ref_label: value=flo_label)
   path_trsf_flo_ref : fichier de transformation flo<-ref (tel qu'en sortie de planeRegistration path_trsf_flo_ref %s) permettant de calculer un repositionnement de l'image flottante dans le referentiel de l'image de reference (accepte aussi les transformations sous forme de numpy.ndarray ou de "list of lists")
   zeros : les zeros en entree sont conserves en sortie si True
   backgroundLabels : si zeros a True, met a zero les labels stipules dans la liste passee en argument
   force : force la valeur de sortie des labels (tableau d'entree n*3) si True
   ref : utilise les labels de l'image in comme labels de reference si True
   visu : (defaut: False) mise a zero des voxels a la frontiere entre deux labels differents, utile pour la visualisation des images
   verbose : option de verbosite
  '''
  flag_remove_path_correspondences=False
  flag_remove_path_trsf=False
  assert(os.path.exists(path_ref))
  assert(os.path.exists(path_flo))
  if type(path_pairs_ref_flo)==dict:
    from morpheme_lineage import write_correspondences
    correspondences=path_pairs_ref_flo
    path_pairs_ref_flo="tmp_ref_flo_associateRefFloLabels.pairs"
    write_correspondences(correspondences, path_pairs_ref_flo)
    flag_remove_path_correspondences=True
  assert( type(path_pairs_ref_flo)==str and os.path.exists(path_pairs_ref_flo))

  tmp_ref=path_ref_out + ".tmp.inr"
  tmp_flo=path_flo_out + ".tmp.inr"

  if (type(path_trsf_flo_ref)==list or type(path_trsf_flo_ref)==np.ndarray):
    trsf=np.array(path_trsf_flo_ref)
    path_trsf_flo_ref="tmp_ref_flo_associateRefFloLabels.trsf"
    f=open(path_trsf_flo_ref, "w")
    f.write(str(trsf).replace('[','').replace(']',''))
    f.close()
    flag_remove_path_trsf=True


  if path_trsf_flo_ref:
    assert(os.path.exists(path_trsf_flo_ref))
    apply_trsf(path_flo, path_trsf=path_trsf_flo_ref, path_output=path_flo_out, template=path_ref, nearest=True, verbose=verbose)

  if zeros and len(backgroundLabels):
    resetLabels(path_ref, path_ref_out, backgroundLabels, verbose=verbose)
    if path_trsf_flo_ref:
      resetLabels(path_flo_out, path_flo_out, backgroundLabels, verbose=verbose)
    else:
      resetLabels(path_flo, path_flo_out, backgroundLabels, verbose=verbose)

  if visu:
    if zeros and len(backgroundLabels):
      labelBorders(path_ref_out, tmp_ref, verbose)
      #cmd=path_labelborders + " " + path_ref_out + ' ' + tmp_ref 
    else:
      labelBorders(path_ref_out, tmp_ref, verbose)
    if path_trsf_flo_ref or (zeros and len(backgroundLabels)):
      labelBorders(path_flo_out, tmp_flo, verbose)
    else:
      labelBorders(path_flo, tmp_flo, verbose)
    Logic(tmp_ref, tmp_ref, Mode='inv', verbose=verbose)
    Logic(tmp_flo, tmp_flo, Mode='inv', verbose=verbose)


  if zeros and len(backgroundLabels):
    associateLabels(path_ref_out, path_flo_out, path_pairs_ref_flo, path_ref_out, path_flo_out, path_labels_out=path_labels_out, zeros=zeros, ref=ref, force=force, verbose=verbose)
  else:
    if path_trsf_flo_ref:
      associateLabels(path_ref, path_flo_out, path_pairs_ref_flo, path_ref_out, path_flo_out, path_labels_out=path_labels_out, zeros=zeros, ref=ref, force=force, verbose=verbose)
    else:
      associateLabels(path_ref, path_flo, path_pairs_ref_flo, path_ref_out, path_flo_out, path_labels_out=path_labels_out, zeros=zeros, ref=ref, force=force, verbose=verbose)
  
  if flag_remove_path_correspondences:
    cmd='rm '+path_pairs_ref_flo
    if verbose:
      print cmd
    os.system(cmd)

  if flag_remove_path_trsf:
    cmd='rm '+path_trsf_flo_ref
    if verbose:
      print cmd
    os.system(cmd)

  if visu:
    Logic(tmp_ref, path_ref_out, path_ref_out, Mode='mask', verbose=verbose)
    Logic(tmp_flo, path_flo_out, path_flo_out, Mode='mask', verbose=verbose)
    cmd="rm " + tmp_ref + ' ' + tmp_flo
    if verbose:
      print cmd
    os.system(cmd)


def sisters_association(point_cloud_ref, point_cloud_flo, correspondences, sisters_ref, sisters_flo, threshold_angle=90.0,
				background_ref = 1, background_flo = 1, trsf_type='affine', estimator='lts', lts_fraction=1.0, bash_options=None,
				verbose=False):
	'''
	Function that determines correspondences between reference and floating pairs of labels, given a reference 
	and a floating point-clouds (pairs of labels to be associated must belong to these point-clouds) and given a set of correspondences.
	Inputs:
		point_cloud_ref : dict
		point_cloud_flo : dict
		correspondences : dict
		sisters_ref     : list or tuple of length 2
		sisters_flo     : list or tuple of length 2
	Optional association parameter:
		threshold_angle : positive angle in degrees (default value is 90). 
						  Sisters association is kept for the result iif sisters axes form an absolute angle that is smaller than the 
						  given threshold

	Optional parameters (for pointCloudRegistration cpp_wrapping imported function):
		background_ref  : label de fond a ignorer dans l'image ref correspondante
		background_flo  : label de fond a ignorer dans l'image flo correspondante
		trsf_type       : computes 'affine' (set as default) or 'rigid' transformation
		estimator       : transformation estimator:
   							'wlts': weighted least trimmed squares
   							'lts': least trimmed squares (default)
   							'wls': weighted least squares
							'ls': least squares
		lts_fraction    : for trimmed estimations, fraction of pairs that are kept (default: 1.0)
		bash_options    : parameter of type str (None by default) to set iif the user wants to add specific bash options for the 
						  execution of the program pointCloudRegistration from library vt (morpheme_privat)
	Outputs:
		new_correspondences : dict that is a copy of the input correspondences with the new keys corresponding to sisters_ref labels
							  and values corresponding to associated sisters_flo labels so that it minimizes the squares distances
							  between associated labels with respect to the T_flo<-ref transformation computed in function of the input 
							  correspondences
		(psn, volume_ratio_ref, volume_ratio_flo) : tuple of length 3:
							psn: "produit scalaire normalise", indeed normalized dot product value between division axes of sisters_ref 
								 and sisters_flo
							volume_ratio_ref: volume ratio between reference sisters under constraint that ratio <=1 (if the volumes 
											  are not known, returns 1) 
							volume_ratio_flo: volume ratio between floating sisters following the association determined by this function
											  with the reference sisters (means that it may be > 1) (if the volumes are not known, returns 1)
	'''
	from math import pi
	if len(sisters_ref) != 2 or len(sisters_flo) != 2:
		if verbose:
			print "sisters_ref or sisters_flo does not have the expected number of elements (2)."
		return {}, (None, None, None)

	out=pointCloudRegistration( point_cloud_ref, point_cloud_flo, correspondences, skip_not_found=True, 
								background_ref=background_ref, background_flo=background_flo,
								trsf_type=trsf_type, estimator=estimator, lts_fraction=lts_fraction, 
								lazy=False, verbose=verbose)
	trsf=np.array(out['trsf']) 
	daughter_ref_0=np.array(point_cloud_ref[sisters_ref[0]][0:3]+(1,))
	daughter_ref_1=np.array(point_cloud_ref[sisters_ref[1]][0:3]+(1,))
	daughter_flo_0=np.array(point_cloud_flo[sisters_flo[0]][0:3]+(1,))
	daughter_flo_1=np.array(point_cloud_flo[sisters_flo[1]][0:3]+(1,))
	trsf_daughter_ref_0=np.dot(trsf, daughter_ref_0)
	trsf_daughter_ref_1=np.dot(trsf, daughter_ref_1)
	trsf_vector_ref=trsf_daughter_ref_1-trsf_daughter_ref_0
	vector_flo=daughter_flo_1-daughter_flo_0
	ps=np.vdot(trsf_vector_ref, vector_flo)
	psn=ps/(np.linalg.norm(trsf_vector_ref)*np.linalg.norm(vector_flo)) # psn est le produit scalaire normalise entre les deux vecteurs definis par les centres de gravite des cellules filles juste apres la division
	angle_degree= np.arccos(abs(psn))*180/pi
	# Pour l'appariement de cellules, la minimisation au sens des moindres carres de distance equivaut a associer daughter_ref_0 a daughter_flo_0 et daughter_ref_1 a daughter_ref_1 ssi le psn > 0, et sinon on inverse les associations
	volume_ratio_ref=1
	volume_ratio_flo=1
	new_correspondences = {}
	#new_correspondences = correspondences.copy()
	if psn > 0:
		if angle_degree < threshold_angle:
			new_correspondences[sisters_ref[0]]=sisters_flo[0]
			new_correspondences[sisters_ref[1]]=sisters_flo[1]
		if len(point_cloud_ref[sisters_ref[0]])==len(point_cloud_ref[sisters_ref[1]])==len(point_cloud_flo[sisters_flo[0]])==len(point_cloud_flo[sisters_flo[1]])==4: 
			volume_ratio_ref = point_cloud_ref[sisters_ref[1]][3] / point_cloud_ref[sisters_ref[0]][3]
			volume_ratio_flo = point_cloud_flo[sisters_flo[1]][3] / point_cloud_flo[sisters_flo[0]][3]
	else:
		if angle_degree < threshold_angle:
			new_correspondences[sisters_ref[0]]=sisters_flo[1]
			new_correspondences[sisters_ref[1]]=sisters_flo[0]
		if len(point_cloud_ref[sisters_ref[0]])==len(point_cloud_ref[sisters_ref[1]])==len(point_cloud_flo[sisters_flo[0]])==len(point_cloud_flo[sisters_flo[1]])==4: 
			volume_ratio_ref = point_cloud_ref[sisters_ref[1]][3] / point_cloud_ref[sisters_ref[0]][3]
			volume_ratio_flo = point_cloud_flo[sisters_flo[0]][3] / point_cloud_flo[sisters_flo[1]][3]
	if volume_ratio_ref > 1:
		volume_ratio_ref = 1/volume_ratio_ref
		volume_ratio_flo = 1/volume_ratio_flo
	return new_correspondences, (angle_degree, volume_ratio_ref, volume_ratio_flo)



def propagate_all_correspondences(lineage_ref, lineage_flo, correspondences, starting_time_ref=None, threshold_angle=90.0, stop_when_blocked=False, verbose=False ):
	"""
	Function for propagation of all the cell-to-cell correspondences between given lineage_ref and lineage_flo (as Morpheme_lineage instances). 
	Propagation will start at specified starting_time_ref (by default, it starts at the first time-point of lineage_ref).
	The parameter 'correspondences' must be a dictionary with keys (resp. values) as relabellings of reference (resp. floating) embryo cells.
	This parameter can also be a filename at readable format for the "readLUT" function that contains the relabelled cell-to-cell 
	correspondences between ref and flo embryos.
	Optional parameter for sisters_association:
			threshold_angle : positive angle in degrees (default value is 90). 
						      Sisters association is kept for the result iif sisters axes form an absolute angle that is smaller than the 
						      given threshold


	Returns:
		new_correspondences    : dict structure with all the cell-to-cell correspondences built by this function.
		unpropagated_cells_ref : list of relabelled cells of lineage_ref that have not been propagated due to unrespected threshold condition.
		scalars                : dict structure with keys as relabelled cells that have been propagated (or tried to be propagated) with 
								 their associated scalar measures computed by the function "sisters_association".
	"""
	import morpheme_lineage as ml
	if type(correspondences)==str:
		correspondences=readLUT(correspondences)
	assert type(correspondences)==dict
	new_correspondences = correspondences.copy()
	if starting_time_ref==None:
		starting_time_ref=lineage_ref.timepoints()[0]

	scalars={}
	unpropagated_cells_ref=[]
	next_time_ref=starting_time_ref
	while next_time_ref < lineage_ref.timepoints()[-1]:
		if verbose:
			print "Time-point ref t%03d" % next_time_ref
		next_labels_ref, next_time_ref=lineage_ref.relabelled_next_death(next_time_ref)
		# Propagation of correspondences for all cells of lineage_ref that die at time next_time_ref
		for next_label_ref in next_labels_ref:
			if next_label_ref in new_correspondences.keys():
				if verbose:
					print "relabel ref %04d" % next_label_ref
				#correspondences = new_correspondences
				point_cloud_ref, point_cloud_flo, daughters_ref, daughters_flo=ml.lineages_relabelled_correspondences_propagation(lineage_ref, lineage_flo, 
																new_correspondences, 
																next_label_ref)
				if daughters_ref and daughters_flo:
					if len(daughters_ref) == 2 and (new_correspondences.has_key(daughters_ref[0]) or new_correspondences.has_key(daughters_ref[1])):
						if verbose:
							print "-> daughters already in the correspondence dictionary"
					else:
						daughters_correspondences, scalars[next_label_ref] = sisters_association(point_cloud_ref, point_cloud_flo, new_correspondences, daughters_ref, daughters_flo, threshold_angle=threshold_angle)
						new_correspondences.update(daughters_correspondences)
						if not daughters_correspondences:
							unpropagated_cells_ref.append(next_label_ref)
							if verbose:
								print "-> daughters division axes form an angle %.1f which is higher than the tolerance" % scalars[next_label_ref][0]
							if stop_when_blocked:
								return new_correspondences, unpropagated_cells_ref, scalars

				else:
					if verbose:
						print "-> end of lineage correspondences"
			else:
				if verbose:
					print "-> skipping relabel ref %04d (not in the correspondence dictionary)" % next_label_ref
		next_time_ref=next_time_ref+1

	return new_correspondences, unpropagated_cells_ref, scalars

def temporal_affine_registration(lineage_ref, lineage_flo, correspondences, include_end_of_time=False):
	"""
	Affine temporal registration of a floating lineage lineage_flo onto a reference lineage lineage_ref 
	with respect to a cell-to-cell correspondence dictionary correspondences.
	The registration is processed with a basic linear regression of data given by corresponding cell death time-points.
	Inputs:
		lineage_ref
		lineage_flo
		correspondences
		linclude_end_of_time (default=False) : excludes from the regression the data from cells whose death time-points 
											   coincide with the last time-point of the reference or floating lineage.
	Outputs:
		a, b : providing the relation t_flo = a * t_ref + b 
	"""

	t_ref_max=lineage_ref.timepoints()[-1]
	t_flo_max=lineage_flo.timepoints()[-1]

	times_death_ref=[]
	times_death_flo=[]

	for k, v in correspondences.items():
		if not lineage_ref.exists(k):
			print "Warning : relabel %d not found in lineage_ref."%k
		else:
			if not lineage_flo.exists(v):
				print "Warning : relabel %d not found in lineage_flo."%v
			else:
				t_ref=lineage_ref.relabelled_death(k)
				t_flo=lineage_flo.relabelled_death(v)
				if include_end_of_time or (t_ref < t_ref_max and t_flo < t_flo_max):
					times_death_ref.append(t_ref)
					times_death_flo.append(t_flo)

	X=times_death_ref
	Y=times_death_flo

	# Linear regression on (X,Y)
	# -> find (a,b) that minimize square error given by sum(|y-(ax+b)|^2) for 

	c=np.cov(X,Y, bias=False)
	a=c[0,1]/c[0,0];
	b=np.mean(Y)-a*np.mean(X)

	return a,b



