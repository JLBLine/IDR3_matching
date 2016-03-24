#!/usr/bin/python
from sys import path
path.append('/home/jline/Documents/cataloguing/PUMA/scripts')
import make_table_lib as mkl
import atpy
import numpy as np
from itertools import combinations
import optparse
from astropy.table import Table, Column, MaskedColumn
from make_table_function import *

dr= np.pi/180.0

parser = optparse.OptionParser()

parser.add_option('-f', '--cat_freqs', 
	help='Enter the frequencies of each catalogue')

parser.add_option('-c', '--choices', 
	help='Enter the text file of source choices')

parser.add_option('-p', '--pref_cats',
	help='Enter names of matched cataloges, in order of preferable position')

parser.add_option('-m', '--matched_cats', 
	help='Enter the names of each catalogue, base catalogue first')

parser.add_option('-i', '--input_bayes', 
	help='Enter name of eyeball bayes file')

parser.add_option('-o', '--output_name', 
	help='Enter name for output catalogue')

options, args = parser.parse_args()

cat_fs = options.cat_freqs.split(',')
cat_freqs= []
for fs in cat_fs:
	split = fs.split('~')
	split = map(float,split)
	cat_freqs.append(split)

matched_cats = options.matched_cats.split(',')

pref_cats = []
for pref in options.pref_cats.split(','): pref_cats.append(pref)

num_freqs = []
for freq in cat_fs: num_freqs.append(len(freq.split('~')))

mkl.closeness = 1.15/60.0
mkl.high_prob = 0.95
mkl.low_prob = 0.8
mkl.chi_thresh = 10.0
mkl.jstat_thresh = 0.1
mkl.num_freqs = num_freqs
mkl.split = 0

def get_allinfo(all_info,good_names):
	'''Takes a list of strings. Each string is a line containing the information for a single source
	in a matched group from an output file of calculate_bayes.py. Gets all of the information
	from each entry and returns them in a source_group() class'''
	src_all = source_group()
	for entry in all_info:
		info = entry.split()
		if info[1] in good_names or all_info.index(entry) == 0:
			src_all.cats.append(info[0])
			src_all.names.append(info[1])
			src_all.ras.append(float(info[2]))
			src_all.rerrs.append(float(info[3]))
			src_all.decs.append(float(info[4]))
			src_all.derrs.append(float(info[5]))
			src_all.majors.append(float(info[-5]))
			src_all.minors.append(float(info[-4]))
			src_all.PAs.append(float(info[-3]))
			src_all.flags.append(info[-2])
			src_all.IDs.append(info[-1])
			
			##If the source only has one frequency. Append as an array
			##so that all entries to src_all.freqs etc are of the same
			##type. This deals with cats with multiple freqs, otherwise
			##the position, name etc will have to be repeated for each
			##frequency
			if len(info)==14:
				src_all.freqs.append(np.array([float(info[6])]))          
				src_all.fluxs.append(np.array([float(info[7])]))
				src_all.ferrs.append(np.array([float(info[8])]))
				#src_all.freqs.append(float(info[6]))  ##Left here in case ever want to append just the freq, not
				#src_all.fluxs.append(float(info[7]))  ##an array
				#src_all.ferrs.append(float(info[8]))
				
			##If not, work out how many freqs there are and append to lists
			else:
				extra = (len(info)-14) / 3
				freqs = []
				fluxs = []
				ferrs = []
				for i in xrange(extra+1):
					##Test to see if flux is a nan or -100000.0; make sure all flux/freq info is -100000.0 if so
					if np.isnan(float(info[7+(3*i)])) == False or float(info[7+(3*i)]) != -100000.0:
						freqs.append(float(info[6+(3*i)]))
						fluxs.append(float(info[7+(3*i)]))
						ferrs.append(float(info[8+(3*i)]))
					else:
						freqs.append(-100000.0)
						fluxs.append(-100000.0)
						ferrs.append(-100000.0)
				src_all.freqs.append(np.array(freqs))
				src_all.fluxs.append(np.array(fluxs))
				src_all.ferrs.append(np.array(ferrs))
					#src_all.freqs.append(float(info[6+(3*i)]))
					#src_all.fluxs.append(float(info[7+(3*i)]))
					#src_all.ferrs.append(float(info[8+(3*i)]))
	return src_all
	
def get_srcg(src_all,info,good_names):
	'''Takes a string which contains the information for all sources in particular combination
	Uses num_freqs to work out where each piece of information is, and then return the relevant
	information in a source_group() class. Any catalogues with no matches are entered as -10000.0'''
	src_g = source_group()
	##Work out where in the string each different catalogue will start, using num_freqs
	##to work out how many entries each catalogue will have
	indexes = [(14+((i-1)*3)) for i in num_freqs]
	starts = [0]
	for i in xrange(len(indexes)-1): starts.append(sum(indexes[:i+1]))
	#all_freqs = []
	#all_fluxs = []
	#all_ferrs = []
	for j in xrange(len(starts)): 
		num_freq = num_freqs[j]
		ind = starts[j]
		#cat = info[ind]
		
		#print info[ind+1],good_names
		
		cat_slot = matched_cats[j]
		
		if cat_slot in src_all.cats:
			srcall_ind = src_all.cats.index(cat_slot)
			freqss = []
			fluxss = []
			ferrss = []
			for k in xrange(num_freq):
				##Test to see if flux is a nan or -100000.0; make sure all flux/freq info is -100000.0 if so
				if np.isnan(src_all.fluxs[srcall_ind][k]) == False or float(src_all.fluxs[srcall_ind][k]) != -100000.0:
					freqss.append(float(src_all.freqs[srcall_ind][k]))
					fluxss.append(float(src_all.fluxs[srcall_ind][k]))
					ferrss.append(float(src_all.ferrs[srcall_ind][k]))
				else:
					freqss.append(-100000.0)
					fluxss.append(-100000.0)
					ferrss.append(-100000.0)
				#if float(info[6+ind+(3*k)])!=-100000.0: all_freqs.append(float(info[6+ind+(3*k)]))
				#if float(info[7+ind+(3*k)])!=-100000.0: all_fluxs.append(float(info[7+ind+(3*k)]))
				#if float(info[8+ind+(3*k)])!=-100000.0: all_ferrs.append(float(info[8+ind+(3*k)]))
			src_g.freqs.append(freqss)
			src_g.fluxs.append(fluxss)
			src_g.ferrs.append(ferrss)
			src_g.cats.append(src_all.cats[srcall_ind])
			src_g.names.append(src_all.names[srcall_ind])
			src_g.ras.append(src_all.ras[srcall_ind])
			src_g.rerrs.append(src_all.rerrs[srcall_ind])
			src_g.decs.append(src_all.decs[srcall_ind])
			src_g.derrs.append(src_all.derrs[srcall_ind])
			src_g.majors.append(src_all.majors[srcall_ind])
			src_g.minors.append(src_all.minors[srcall_ind])
			src_g.PAs.append(src_all.PAs[srcall_ind])
			src_g.flags.append(src_all.flags[srcall_ind])
			src_g.IDs.append(src_all.IDs[srcall_ind])
			src_g.prob = float(info[-1])
			#freqss = []
			#fluxss = []
			#ferrss = []
			#for k in xrange(num_freq):
				###Test to see if flux is a nan or -100000.0; make sure all flux/freq info is -100000.0 if so
				#if np.isnan(float(info[7+ind+(3*k)])) == False or float(info[7+(3*i)]) != -100000.0:
					#freqss.append(float(info[6+ind+(3*k)]))
					#fluxss.append(float(info[7+ind+(3*k)]))
					#ferrss.append(float(info[8+ind+(3*k)]))
				#else:
					#freqss.append(-100000.0)
					#fluxss.append(-100000.0)
					#ferrss.append(-100000.0)
				##if float(info[6+ind+(3*k)])!=-100000.0: all_freqs.append(float(info[6+ind+(3*k)]))
				##if float(info[7+ind+(3*k)])!=-100000.0: all_fluxs.append(float(info[7+ind+(3*k)]))
				##if float(info[8+ind+(3*k)])!=-100000.0: all_ferrs.append(float(info[8+ind+(3*k)]))
			#src_g.freqs.append(freqss)
			#src_g.fluxs.append(fluxss)
			#src_g.ferrs.append(ferrss)
			#src_g.cats.append(info[ind])
			#src_g.names.append(info[ind+1])
			#src_g.ras.append(float(info[ind+2]))
			#src_g.rerrs.append(float(info[ind+3]))
			#src_g.decs.append(float(info[ind+4]))
			#src_g.derrs.append(float(info[ind+5]))
			#src_g.majors.append(float(info[ind+9+((num_freq-1)*3)]))
			#src_g.minors.append(float(info[ind+10+((num_freq-1)*3)]))
			#src_g.PAs.append(float(info[ind+11+((num_freq-1)*3)]))
			#src_g.flags.append(info[ind+12+((num_freq-1)*3)])
			#src_g.IDs.append(info[ind+13+((num_freq-1)*3)])
			#src_g.prob = float(info[-1])
		else:
			freqss = []
			fluxss = []
			ferrss = []
			for k in xrange(num_freq):
				freqss.append(-100000.0)
				fluxss.append(-100000.0)
				ferrss.append(-100000.0)
				#if float(info[6+ind+(3*k)])!=-100000.0: all_freqs.append(float(info[6+ind+(3*k)]))
				#if float(info[7+ind+(3*k)])!=-100000.0: all_fluxs.append(float(info[7+ind+(3*k)]))
				#if float(info[8+ind+(3*k)])!=-100000.0: all_ferrs.append(float(info[8+ind+(3*k)]))
			src_g.freqs.append(freqss)
			src_g.fluxs.append(fluxss)
			src_g.ferrs.append(ferrss)
			src_g.cats.append('-100000.0')
			src_g.names.append('-100000.0')
			src_g.ras.append(-100000.0)
			src_g.rerrs.append(-100000.0)
			src_g.decs.append(-100000.0)
			src_g.derrs.append(-100000.0)
			src_g.majors.append(-100000.0)
			src_g.minors.append(-100000.0)
			src_g.PAs.append(-100000.0)
			src_g.flags.append('-100000.0')
			src_g.IDs.append('-100000.0')
			src_g.prob = float(info[-1])
	return src_g

def combine_flux(src_all,src_g):
	'''Takes a src_group() class that contains all group info. Indentifies which catalogue is repeated, combines
	the fluxes and fits a new line. Returns the reduced frequency list, the combined flux and flux error 
	arrays, as well as the line_fit object and residuals'''
	
	accepted_inds = xrange(len(src_all.names))
	
	##Find repeated catalogues
	repeated_cats = set([src_all.cats[ind] for ind in accepted_inds if src_all.cats.count(src_all.cats[ind])>1])
	##This won't neccesarily be in the that the cats appear in src_all.cats so reorder
	repeat_indexs = [src_all.cats.index(cat) for cat in repeated_cats]
	repeated_cats = [cat for ind,cat in sorted(zip(repeat_indexs,repeated_cats),key=lambda pair: pair[0])]
	
	##These are used to test the combined spectrum
	temp_freqs = [src_all.freqs[i] for i in xrange(len(src_all.freqs)) if src_all.cats[i] not in repeated_cats]
	temp_fluxs = [src_all.fluxs[i] for i in xrange(len(src_all.fluxs)) if src_all.cats[i] not in repeated_cats]
	temp_ferrs = [src_all.ferrs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	
	base_freqs = [src_all.freqs[i] for i in xrange(len(src_all.freqs)) if src_all.cats[i] not in repeated_cats]
	base_fluxs = [src_all.fluxs[i] for i in xrange(len(src_all.fluxs)) if src_all.cats[i] not in repeated_cats]
	base_ferrs = [src_all.ferrs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_cats = [src_all.cats[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_ras = [src_all.ras[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_rerrs = [src_all.rerrs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_decs = [src_all.decs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_derrs = [src_all.derrs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_majors = [src_all.majors[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_minors = [src_all.minors[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_PAs = [src_all.PAs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_names = [src_all.names[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	
	##Will need these for fitting/passing on to plotting function
	comb_freqs = []
	comb_fluxs = []
	comb_ferrs = []
	ra_ws = []
	dec_ws = []
	rerr_ws = []
	derr_ws = []
	
	##Need these in case of doing the split test
	resolved_diff_inds = []
	unrepeat_dists = []
	num_of_repeats = []
	
	for repeat_cat in repeated_cats:
		num_of_repeat = [i for i in xrange(len(src_all.names)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		num_of_repeats.append(num_of_repeat)
		
	for cat,name in zip(src_g.cats[1:],src_g.names[1:]):
		if cat not in repeated_cats and cat != '100000.0':
			if cat == 'vlssr':
				src_g.vlssr = name
			elif cat == 'mrc':
				src_g.mrc = name
			elif cat == 'sumss':
				src_g.sumss = name
			elif cat == 'nvss':
				src_g.nvss = name
			elif cat == 'atg20':
				src_g.atg20 = name
			elif cat == 'TGSS':
				src_g.TGSS = name
			else:
				pass
	
	##For each repeated catalogue:
	for repeat_cat in repeated_cats:
		##Find the frequency/ies of repeated cat
		comb_freq = src_all.freqs[src_all.cats.index(repeat_cat)]
		##Find the flux/es of the repeat_cat sources that were accepted by retained_sources()
		flux_to_comb = [src_all.fluxs[i] for i in xrange(len(src_all.fluxs)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		##Find the flux error/s of the repeat_cat sources that were accepted by retained_sources()
		ferr_to_comb = [src_all.ferrs[i] for i in xrange(len(src_all.ferrs)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		
		##ALso find all of the positional information to combine
		ras_to_comb = [src_all.ras[i] for i in xrange(len(src_all.ras)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		rerrs_to_comb = [src_all.rerrs[i] for i in xrange(len(src_all.rerrs)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		decs_to_comb = [src_all.decs[i] for i in xrange(len(src_all.decs)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		derrs_to_comb = [src_all.derrs[i] for i in xrange(len(src_all.derrs)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		
		##This is to write down all the names of the sources combined
		names_to_comb = [src_all.names[i] for i in xrange(len(src_all.names)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		name_string = ''
		for name in names_to_comb: name_string += ','+name
		
		if repeat_cat == 'vlssr':
			src_g.vlssr = name_string[1:]
		elif repeat_cat == 'mrc':
			src_g.mrc = name_string[1:]
		elif repeat_cat == 'sumss':
			src_g.sumss = name_string[1:]
		elif repeat_cat == 'nvss':
			src_g.nvss = name_string[1:]
		elif repeat_cat == 'atg20':
			src_g.atg20 = name_string[1:]
		elif repeat_cat == 'TGSS':
			src_g.TGSS = name_string[1:]
		else:
			pass
		
		##TEST TO SEE IF THE REPEATED SOURCES ARE RESOLVED (BY A GIVEN RESOLUTION THRESHOLD)
		##Need to do the split test here even if not propagating to the final catalogue
		##-------------------------------------------------------------------------------------------------------
		
		comb_flux = sum(flux_to_comb) #Sum the fluxes
		comb_ferr = np.zeros(1)
		for ferr in ferr_to_comb: comb_ferr+= ferr**2 #Add the errors in quadrature
		comb_ferr = comb_ferr**0.5

		##Need these later to fit the combined flux and to populate a new src_g if
		##the fit passes
		temp_freqs.append(comb_freq)
		temp_fluxs.append([comb_flux])
		temp_ferrs.append([comb_ferr[0]])
		comb_freqs.append(comb_freq)
		comb_fluxs.append(comb_flux)
		comb_ferrs.append(comb_ferr)
		
		##Weight by the first flux in the flux list (in case catalogue has multiple frequencies)
		flux_s = [flux[0] for flux in flux_to_comb]

		##A a little code in case sources are very close to RA=0, and some are reporting 359.9,
		##and others 0.001 - if you don't account for this, the weigthed position goes mental
		wrap = 'no'
		for combo in combinations(ras_to_comb,2):
			diff = combo[0]-combo[1]
			if abs(diff)>180.0:
				wrap='yes'
		if wrap=='yes':
			for i in xrange(len(ras_to_comb)):
				if ras_to_comb[i]<180.0:
					ras_to_comb[i]+=360.0
		##Weight the sources by their flux
		weights = [flux/sum(flux_s) for flux in flux_s]

		##Do the weighting, if weighted RA is above 360 deg, rescale, 
		##add the errors as shown in the write up
		ra_w = np.dot(ras_to_comb,weights)
		if ra_w>360.0: ra_w-=360.0
		dec_w = np.dot(decs_to_comb,weights)
		rerr_w = (np.dot(rerrs_to_comb,weights)**2)**0.5
		derr_w = (np.dot(derrs_to_comb,weights)**2)**0.5
		
		ra_ws.append(ra_w)
		dec_ws.append(dec_w)
		rerr_ws.append(rerr_w)
		derr_ws.append(derr_w)
	
	log_temp_freqs = []
	log_temp_fluxs = []
	log_temp_ferrs = []
	
	##Get the sources out of the array in list format (which is used later when making the sources
	##to add to the final table)
	#for freqs in temp_freqs:
		#for freq in freqs: log_temp_freqs.append(np.log(freq))
	#for fluxs in temp_fluxs:
		#for flux in fluxs: log_temp_fluxs.append(np.log(flux))
	#for i in xrange(len(temp_ferrs)):
		#for j in xrange(len(temp_ferrs[i])): log_temp_ferrs.append(temp_ferrs[i][j]/temp_fluxs[i][j])
		
	#Get the sources out of the array in list format (which is used later when making the sources
	##to add to the final table)
	for i in xrange(len(temp_freqs)):
		for j in xrange(len(temp_freqs[i])):
			if temp_fluxs[i][j] == -100000.0 or np.isnan(temp_fluxs[i][j])==True:
				pass
			else:
				log_temp_freqs.append(np.log(temp_freqs[i][j]))
	for i in xrange(len(temp_freqs)):
		for j in xrange(len(temp_freqs[i])):
			if temp_fluxs[i][j] == -100000.0 or np.isnan(temp_fluxs[i][j])==True:
				pass
			else:
				log_temp_fluxs.append(np.log(temp_fluxs[i][j]))
	for i in xrange(len(temp_freqs)):
		for j in xrange(len(temp_freqs[i])):
			if temp_fluxs[i][j] == -100000.0 or np.isnan(temp_fluxs[i][j])==True:
				pass
			else:
				log_temp_ferrs.append(temp_ferrs[i][j]/temp_fluxs[i][j])
		
	print temp_freqs,log_temp_freqs
	print temp_fluxs, log_temp_fluxs
	print temp_ferrs, log_temp_ferrs

	##Fit and find residuals to the combined spectrum
	comb_fit,comb_jstat,comb_bse,comb_chi_red = mkl.fit_line(np.array(log_temp_freqs),np.array(log_temp_fluxs),np.array(log_temp_ferrs))
	print comb_chi_red
	print "============================================="
	
	
	
	
	##Find out where in srg_g the repeated cats appear
	repeat_cat_inds = [src_g.cats.index(cat) for cat in repeated_cats]
	
	split_flag=''
	##Make labels for when we're plotting a put in combined_names
	combined_names = []
	#if comb_jstat<=jstat_thresh or comb_chi_red<=chi_thresh:
		##Create the combined source no matter what for plotting purposes
		##Loop over all the combined sources and repopulate the entries of a src_g
		##at the point where the repeated catalogues appear
	for i in xrange(len(comb_fluxs)):
		srcg_ind = repeat_cat_inds[i]
		src_g.ras[srcg_ind] = ra_ws[i]
		src_g.rerrs[srcg_ind] = rerr_ws[i]
		src_g.decs[srcg_ind] = dec_ws[i]
		src_g.derrs[srcg_ind] = derr_ws[i]
		src_g.PAs[srcg_ind] = -100000.0
		src_g.majors[srcg_ind] = -100000.0
		src_g.minors[srcg_ind] = -100000.0
		src_g.names[srcg_ind] = "Combined-%s" %src_g.cats[srcg_ind]
		combined_names.append("Combined-%s" %src_g.cats[srcg_ind])
		src_g.fluxs[srcg_ind] = [comb_fluxs[i]]
		src_g.ferrs[srcg_ind] = comb_ferrs[i]
		
	##Need to make sure the single sources in src_all are used rather than
	##whatever was in src_g
	for i in xrange(len(base_fluxs)):
		srcg_ind = src_g.cats.index(base_cats[i])
		src_g.ras[srcg_ind] = base_ras[i]
		src_g.rerrs[srcg_ind] = base_rerrs[i]
		src_g.decs[srcg_ind] = base_decs[i]
		src_g.derrs[srcg_ind] = base_derrs[i]
		src_g.PAs[srcg_ind] = base_PAs[i]
		src_g.majors[srcg_ind] = base_majors[i]
		src_g.minors[srcg_ind] = base_minors[i]
		src_g.names[srcg_ind] = base_names[i]
		src_g.fluxs[srcg_ind] = base_fluxs[i]
		src_g.ferrs[srcg_ind] = base_ferrs[i]
		
	#srg_g.freqs = temp_freqs
	src_g.SI = float(comb_fit.params[0])
	src_g.intercept = comb_fit.params[1]
	src_g.SI_err = comb_bse[0]
	src_g.intercept_err = comb_bse[1]
	src_g.chi_resid = comb_chi_red
	src_g.epsilon_red = comb_jstat

	##If good fit, report that in the final stats object
	if comb_chi_red<=2:
		src_g.low_resids = 0
	else:
		src_g.low_resids = 1
			
	return src_g

sources = []
sources_stats = []

##Input eyeball file
bayes_comp = open(options.input_bayes).read().split('END_GROUP')
del bayes_comp[-1]

choices = open(options.choices,'r').read().split('\n')
if choices[-1] == '': del choices[-1]

choice_dict = {}
choice_names = []

for choice in choices:
	name,inc_names = choice.split()[0],choice.split()[1]
	inc_names = inc_names.split(',')
	choice_dict[name] = inc_names
	choice_names.append(name)
	
for comp in bayes_comp:

	##Get the information into nice usable forms, and get rid of empty/pointless
	##entries
	chunks = comp.split('START_COMP')
	all_info = chunks[0].split('\n')
	
	if 'START' in all_info[1]:
		src_name = all_info[2].split()[1]
	else:
		src_name = all_info[1].split()[1]
	
	for entry in all_info:
		info = entry.split()
	
	##FOR SOME REASON CAN'T DO BOTH OF THESE LINES IN THE SAME FOR LOOP?!?!?!
	for entry in all_info:   
		if entry=='': del all_info[all_info.index(entry)]
	for entry in all_info:
		if 'START' in entry: del all_info[all_info.index(entry)]

	matches = chunks[1].split('\n')
	
	del matches[0],matches[-1]
	stats = matches[-1]
	del matches[-2:]
	
	match1 = matches[0].split()
	
	##There may be sources in the mod pile that don't actually need editing
	##Only use what is in the list
	if src_name in choice_names:
		#print src_name
		##Need to only include sources from src_all and src_g we actually want
		if choice_dict[src_name] == ['all']:
			src_all = mkl.get_allinfo(all_info)
			src_g = get_srcg(src_all,match1,src_all.names)
			#name_choice = [name for name in src]
			#src_g = mkl.get_srcg(match1)
		else:
			src_all = get_allinfo(all_info,choice_dict[src_name])
			src_g = get_srcg(src_all,match1,choice_dict[src_name])
		
		##Get some info and find which catalogues are present 
		
		present_cats = [cat for cat in src_all.cats if cat!='-100000.0']
		
		flags = src_all.flags[0].split('+')
				
		meh,num_matches,accept_matches,accepted_inds,accept_type,stage = stats.split()
		g_stats = mkl.group_stats()
		g_stats.num_matches = int(num_matches)
		g_stats.retained_matches = int(accept_matches)
		
		src_comb = combine_flux(src_all,src_g)
				
		num_matches = len(matches)
		src_comb.inspected = 2
		src_comb.accept = accept_type
		src_comb.accept_type = stage
		sources.append(src_comb)
		#g_stats.accept_type = 'multiple'
		
		g_stats.outcome = accept_type
		if 'combine' in stage: 
			stage = 'multiple'
		elif 'position' in stage:
			stage = 'isolated'
		g_stats.accept_type = stage
		
		##See how many unique catalogues there are
		num_cats = len(set([cat for cat in src_all.cats if cat!='-100000.0']))
		g_stats.num_cats = num_cats
		
		sources_stats.append(g_stats)
	
t=Table(masked=True)

out_t = make_the_table(t,sources,sources_stats)

out_t.write("%s.vot" %options.output_name,format='votable',overwrite=True)
#out_t.write("%s.fits" %options.output_name,format='fits',overwrite=True)
	
##inpsec = 3	
	
