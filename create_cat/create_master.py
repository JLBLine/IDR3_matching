from numpy import *
from astropy.table import Table, vstack, Column
import optparse
import subprocess

parser = optparse.OptionParser()

parser.add_option('-t', '--tags',
	help='Enter tags of the checked lists')
	
parser.add_option('-v', '--v5_tag',
	help='Enter names v5 match tag')
	
parser.add_option('-e', '--extras', 
	help='Enter name of any extra tables to append')
	
parser.add_option('-b', '--base_table', 
	help='Enter name of PUMA matched base table')

parser.add_option('-d', '--delete_names', 
	help='Enter names of sources to be deleted')

parser.add_option('-o', '--output_name', 
	help='Enter name of output vot')

options, args = parser.parse_args()

tags = options.tags.split(',')

main_table = Table.read(options.base_table)
accept_names = list(main_table['KGS_ID'])

#print len(accept_names),accept_names.index('351.0')

cluster_rads = array(['h' for i in xrange(len(accept_names))])

clus_col = Column(name='R_cluster',data=cluster_rads,description='The clustering radius use in the KGS process; q for one-quarter beam radius, h for one-half beam radius',dtype=str)

#main_table.add_columns([clus_col])

mask_list = []

def delete_rows(table_name,main_table):
	skip = 'no'
	try:
		tag_table = Table.read(table_name)
		tag_names = list(tag_table['KGS_ID'])
		for name in tag_names:
			accept_names = list(main_table['KGS_ID'])
			ind = accept_names.index(name)
			main_table.remove_row(ind)
		
	except IndexError:
		pass
	
def append_table(table_name,main_table):
	try:
		tag_table = Table.read(table_name)
		main_table = vstack([main_table,tag_table])
	except IndexError:
		pass

##For each of the modified combine catagories, find and replace the update
##entry in the main table
for tag in tags:
	print tag
	delete_rows('%s_choices.vot' %tag,main_table)
	#delete_rows('%s_agrees.vot' %tag,main_table)
	delete_rows('%s_edited.vot' %tag,main_table)
	
main_table.write('KATALOGSS_PUMA_final_1.vot',format='votable',overwrite=True)


main_table = Table.read('KATALOGSS_PUMA_final_1.vot')

##I don't know why but cannot do this in a function like the one above
for tag in tags:
	for label in ['choices','edited']:
		try:
			tag_table = Table.read('%s_%s.vot' %(tag,label))
			print 'Adding table %s_%s.vot' %(tag,label)
			main_table = vstack([main_table,tag_table])
		except IndexError:
			pass

if options.v5_tag:
	for label in ['choices','edited']:
		try:
			tag_table = Table.read('%s_%s.vot' %(options.v5_tag,label))
			print 'Adding table %s_%s.vot' %(options.v5_tag,label)
			main_table = vstack([main_table,tag_table])
		except IndexError:
			pass
		
if options.extras:
	tables = options.extras.split(',')
	for table in tables:
		try:
			print 'Adding table '+table
			tag_table = Table.read(table)
			main_table = vstack([main_table,tag_table])
		except IndexError:
			pass
		
if options.delete_names:
	delete_names = open(options.delete_names,'r').read().split('\n')
	if delete_names[-1]=='': del delete_names[-1]
			
	for name in delete_names:
		if '#' in name:
			pass
		else:
			accept_names = list(main_table['KGS_ID'])
			ind = accept_names.index(name)
			main_table.remove_row(ind)
		

##Old way of doing it
####Read in Pattis variables
#import pickle

###titles ['sig_flux', 'sig_dec', 'R_class', 'snr_mean', 'ra_corr', 'N_det', 'EB_factor', 'flux', 'sig_ra', 'snr', 'ra', 'dec', 'id', 'dec_corr']

#v11_details = pickle.load( open( "source_cand_v11.p", "rb" ) )

#p_names = list(v11_details['id'])
#rclass = v11_details['R_class']
#ndet = v11_details['N_det']
#eb_factor =  v11_details['EB_factor']
##nexp = v11_details['N_exp']
##meanbeam = v11_details['Mean_Beam']
##locdens = v11_details['Local_Density']


###Find names of all sources in final table
#names = list(main_table['KGS_ID'])

###Create arrays to hold data
#t_rclass = ones(len(names))*-1.0
#t_eb = ones(len(names))*-1.0
#t_ndet = ones(len(names))*-1.0
#t_meanbeam = ones(len(names))*-1.0
##t_nexp = ones(len(names))*-1.0
##t_locdens = ones(len(names))*-1.0

#p_names = [str(float(name)) for name in p_names]

#for name in names:
	#try:
		#ind = names.index(name)
		#if 'v10' in name:
			#name = str(float(name.split('_')[0]))
		#p_ind = p_names.index(name)
		#t_rclass[ind] = rclass[p_ind]
		#t_ndet[ind] = ndet[p_ind]
		#t_eb[ind] = eb_factor[p_ind]
		##t_meanbeam[ind] = meanbeam[p_ind]
		##t_locdens[ind] = locdens[p_ind]
		##t_nexp[ind] = nexp[p_ind]
	#except ValueError:
		#pass
	
	
###Now add in the v5 shizz
#v5_details = pickle.load( open( "sources_v5.p" , "rb" ) )
##['sig_flux', 'rms_ra', 'sig_dec', 'rms_dec', 'id', 'flux0_mean_sig', 'flux0', 'flux_mean_sig', 'ra', 'bm', 'nobs2', 'nobs', 'fd', 'pdet', 'multi', 'rdet', 'sig_flux0', 'flux', 'ext', 'sig_ra', 'dec', 'ndet']
	
#p_names = list(v5_details['id'])
##rclass = v5_details['rdet']
#ndet = v5_details['ndet']
##eb_factor =  v11_details['EB_factor']
#meanbeam = v5_details['bm']
	
#p_names = [str(float(name)) for name in p_names]

#for name in names:
	#try:
		#ind = names.index(name)
		#if 'v5' in name:
			#name = str(float(name.split('_')[0]))
			#p_ind = p_names.index(name)
			##t_rclass[ind] = rclass[p_ind]
			#t_ndet[ind] = ndet[p_ind]
			##t_eb[ind] = eb_factor[p_ind]
			#t_meanbeam[ind] = meanbeam[p_ind]
			##t_nexp[ind] = nexp[p_ind]
			##t_locdens[ind] = locdens[p_ind]
	#except ValueError:
		#pass
	
#main_table['Reliability'] = t_rclass
#main_table['N_det'] = t_ndet
#main_table['EB_factor'] = t_eb
#main_table['Mean_Beam'] = t_meanbeam.round(4)
##main_table['N_exp'] = t_nexp
##main_table['Local_Density'] = t_locdens


#names = list(main_table['KGS_ID'])

###Create arrays to hold data


patti_tab = Table.read('KATALOGSS_RELEASE_extrastats_pc.vot')
pc_names = patti_tab["KGS_ID"]
reli = patti_tab['Reliability']
ndet = patti_tab['N_det']
eb_fac = patti_tab['EB_factor']
mean_beam = patti_tab['Mean_Beam']

pc_dict = {}

for i in xrange(len(pc_names)):
	pc_dict[pc_names[i]] = [reli[i],ndet[i],eb_fac[i],mean_beam[i]]

names = list(main_table['KGS_ID'])	
t_rclass = ones(len(names))*-1.0
t_eb = ones(len(names))*-1.0
t_ndet = ones(len(names))*-1.0
t_meanbeam = ones(len(names))*-1.0

failed_names = []
	
for i in xrange(len(names)):
	
	try:
		if '.0' in names[i]:
			name = names[i].split('.')[0]+'_v10'
		else:
			name = names[i]
		
		rel, nd, eb, mb = pc_dict[name]
		t_rclass[i] = rel
		t_eb[i] = eb
		t_ndet[i] = nd
		t_meanbeam[i] = mb
	except:
		failed_names.append(name)
		
		
#Old way of doing it
###Read in Pattis variables
import pickle

##titles ['sig_flux', 'sig_dec', 'R_class', 'snr_mean', 'ra_corr', 'N_det', 'EB_factor', 'flux', 'sig_ra', 'snr', 'ra', 'dec', 'id', 'dec_corr']

v11_details = pickle.load( open( "source_cand_v11.p", "rb" ) )

p_names = list(v11_details['id'])
rclass = v11_details['R_class']
ndet = v11_details['N_det']
eb_factor =  v11_details['EB_factor']
##nexp = v11_details['N_exp']
##meanbeam = v11_details['Mean_Beam']
##locdens = v11_details['Local_Density']
	
	
p_names = [str(float(name)) for name in p_names]

for name in failed_names:
	try:
		ind = names.index(name)
		if 'v10' in name:
			name = str(float(name.split('_')[0]))
		p_ind = p_names.index(name)
		t_rclass[ind] = rclass[p_ind]
		t_ndet[ind] = ndet[p_ind]
		t_eb[ind] = eb_factor[p_ind]
	except ValueError:
		pass
	
	
	
main_table['Reliability'] = t_rclass
main_table['N_det'] = t_ndet
main_table['EB_factor'] = t_eb
main_table['Mean_Beam'] = t_meanbeam.round(4)



##Get rid of positional offset with Patti's offset fixer
#import numpy as np
from astropy.modeling.models import Polynomial2D
from astropy.io import fits
from astropy.coordinates import SkyCoord
import warnings
warnings.filterwarnings('ignore')

def debias(ra,dec):
    ra = array(ra)
    dec = array(dec)
    ra[ra<180.]+=360.
    dra =  Polynomial2D(2, c0_0=-3.943863933015328, c1_0=0.021853307123963286, 
                        c2_0=-3.0533748559944215e-05, c0_1=-0.004094847465089692, 
                        c0_2=-5.476491790816323e-05, c1_1=3.7775645493556417e-06)(ra,dec)
    ddec = Polynomial2D(2, c0_0=-0.21274304736295002, c1_0=0.0012743212743134297, 
                        c2_0=-1.8989855539056251e-06, c0_1=-0.0005832444732804748, 
                        c0_2=7.253445630101684e-07, c1_1=1.3672088836566055e-06)(ra,dec)
    ra+=dra
    dec+=ddec
    ra[ra>=360.]-=360.
    return (ra,dec)


ra_corr, dec_corr = debias(main_table['KGS_RAJ2000'],main_table['KGS_DECJ2000'])

main_table['KGS_RAJ2000'] = ra_corr
main_table['KGS_DECJ2000'] = dec_corr

##Fix the names to represent the fixed positions
from sys import path
path.append('/home/jline/Documents/cataloguing/PUMA/scripts')
import make_table_lib as mkl

def make_name(ra,dec,ext):
	ra_str = mkl.deg_to_hour(ra)[:6]
	dec_str = mkl.deg_to_degmins(dec)[:7]
	if ext=='':
		return 'KGS J'+ra_str+dec_str
	else:
		return 'KGS J'+ra_str+dec_str+ext
	
names = []

for i in xrange(len(main_table['KGS_ID'])):
	name = make_name(ra_corr[i],dec_corr[i],'')
	names.append(name)

print names[0]

main_table['Name'] = array(names,dtype='a18')

print main_table['Name'][0]

##Whhhhhhhhhhy do these still exist??
#del main_table['N_exp']
#del main_table['Local_Density']
#del main_table['Number_cats']

main_table.write("%s_extrastats.vot" %options.output_name,format='votable',overwrite=True)

##Remove the KGS_ID tag,Outcomes and Chi_red column
del main_table['KGS_ID']
del main_table['Outcomes']
del main_table['Chi_sq_red']
del main_table['Intercept']

main_table.write("%s.vot" %options.output_name,format='votable',overwrite=True)
main_table.write("%s.fits" %options.output_name,format='fits',overwrite=True)

subprocess.call('rm KATALOGSS_PUMA_final_1.vot',shell=True)