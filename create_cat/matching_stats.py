import sys
from astropy.table import Table
import numpy as np

data = Table.read('KGS_compv5-v-m-s-n_comb_fix.vot')
names = data["KGS_ID"]
match_types = data['Match_Type']
outcomes = data['Outcomes']

v5_names = [1054.0,1134.0,2261.0,2290.0,235.0,2472.0,2508.0,2522.0,2710.0,2751.0,3189.0,3251.0,3462.0,3620.0,4032.0,4429.0,5080.0,59.0,6283.0]
v5_types = []
v5_outcomes = []
print len(v5_names), len(v5_types)

for i in xrange(len(names)):
	if float(names[i]) in v5_names:
		v5_types.append(match_types[i])
		v5_outcomes.append(outcomes[i])

##Add this is for the left lobe source we added in by hand
v5_types.append('multiple')
v5_outcomes.append('eyeball')
		
##These are sources replaced with v5 sources
data = Table.read('KGS_compv11-v-m-s-n_comb_fix.vot')
names = data["KGS_ID"]
match_types = data['Match_Type']
outcomes = list(data['Outcomes'])

v11_deletes = [24.0,1020.0,1378.0,1716.0,502.0,1142.0,2674.0,2694.0,334.0,112.0,2174.0]

##Need to explicitly look for these guys as they are in the -eyeball.txt file
v11_types = ['isolated','mulitple','mulitple','mulitple','mulitple','isolated','isolated','isolated','isolated','isolated','isolated']
v11_outcomes = ['accept','eyeball','eyeball','eyeball','eyeball','reject','reject','reject','reject','accept','accept']


#for i in xrange(len(names)):
	#if float(names[i]) in v11_deletes:
		#v11_types.append(match_types[i])
		#v11_outcomes.append(outcomes[i])
		
#print v11_types
#print v11_outcomes
		
def total_change(accept,match_type):
	all_accept = sum([1 for i in xrange(len(outcomes)) if outcomes[i]==accept and match_types[i]==match_type])
	v11_dels = sum([-1 for i in xrange(len(v11_outcomes)) if v11_outcomes[i]==accept and v11_types[i]==match_type])
	v5_dels = sum([1 for i in xrange(len(v5_outcomes)) if v5_outcomes[i]==accept and v5_types[i]==match_type])
	#print all_accept,v11_dels,v5_dels
	return all_accept + v11_dels + v5_dels

print 'ORIGINAL CATALOGUE===================================================='
print "Total number of sources: ",len(names) + len(v5_types) - len(v11_outcomes)
print "All sources accepted: " , sum([1 for out in outcomes if outcomes[i]=='accept']) + sum([1 for i in xrange(len(v5_outcomes)) if v5_outcomes[i]=='accept']) + sum([-1 for i in xrange(len(v11_outcomes)) if v11_outcomes[i]=='accept'])
print "\taccepted by isolated: ", total_change('accept','isolated')
print "\taccepted by dominant: ", total_change('accept','dominant')
print "\taccepted by multiple: ", total_change('accept','multiple')
print "All sources rejected: " , sum([1 for out in outcomes if outcomes[i]=='reject']) + sum([1 for i in xrange(len(v5_outcomes)) if v5_outcomes[i]=='reject']) + sum([-1 for i in xrange(len(v11_outcomes)) if v11_outcomes[i]=='reject'])
print "\trejected by isolated: ", total_change('reject','isolated')
print "\trejected by dominant: ", total_change('reject','dominant')
print "All sources retained to eyeball: " , total_change('eyeball','multiple')

final_table = Table.read('KATALOGSS_RELEASE_extrastats.vot')

f_typees = list(final_table['Match_Type'])
f_outcomes = list(final_table['Outcomes'])
f_num_cats = final_table['Number_cats']
f_inspecs = list(final_table['Inspected'])

def type_num(accept,type_name):
	accept_ls = sum([1 for i in xrange(len(f_inspecs)) if f_outcomes[i]==accept and f_typees[i]==type_name])
	mods = sum([1 for i in xrange(len(f_inspecs)) if f_outcomes[i]==accept and f_typees[i]==type_name and f_inspecs[i]==3])
	#test = [1 for i in xrange(len(inspecs)) if outcomes[i]==accept and typees[i]==type_name]
	#print len(test)
	return accept_ls,mods




print "FINAL TABLE======================================================"
print 'Accepted: %d' %sum([1 for out in f_outcomes if out=='accept'])
print '\tby isolated: %d  modified: %d' %(type_num('accept','isolated'))
print '\tby dominant: %d  modified: %d' %(type_num('accept','dominant'))
print '\tby multiple: %d  modified: %d' %(type_num('accept','multiple'))
print 'Rejected: %d' %sum([1 for out in f_outcomes if out=='reject'])
print '\tby isolated: %d  modified; %d' %(type_num('reject','isolated'))
print '\tby dominant: %d  modified: %d' %(type_num('reject','dominant'))