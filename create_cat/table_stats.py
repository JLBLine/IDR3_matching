
from astropy.io.votable import parse as vot_parse
table = vot_parse('KATALOGSS_RELEASE_extrastats.vot',pedantic=False).get_first_table()
data = table.array

outcomes = data['Outcomes']
inspecs = data['Inspected']
type_matchs = data['Match_Type']

zip_info = zip(outcomes,inspecs)

eyeballed = [outcome for outcome,inspec in zip_info if inspec > 0]
modded = [outcome for outcome,inspec in zip_info if inspec == 3]
modded_changed = [outcome for outcome,inspec in zip_info if inspec == 3 and outcome=='accept']
reject = [outcome for outcome,inspec in zip_info if outcome=='reject']
eyeball = [outcome for outcome,inspec in zip_info if outcome=='eyeball']

iso_mod = [i for i in xrange(len(outcomes)) if outcomes[i]=='accept' and inspecs[i]==3 and type_matchs[i]=='isolated']
dom_mod = [i for i in xrange(len(outcomes)) if outcomes[i]=='accept' and inspecs[i]==3 and type_matchs[i]=='dominant']
mul_mod = [i for i in xrange(len(outcomes)) if outcomes[i]=='accept' and inspecs[i]==3 and type_matchs[i]=='multiple']

isos = [i for i in xrange(len(outcomes)) if outcomes[i]=='accept' and type_matchs[i]=='isolated']
doms = [i for i in xrange(len(outcomes)) if outcomes[i]=='accept' and type_matchs[i]=='dominant']
muls = [i for i in xrange(len(outcomes)) if outcomes[i]=='accept' and type_matchs[i]=='multiple']

print 'Num eyeballed: %d' %len(eyeballed)
print 'Num modded: %d' %len(modded)
print 'Num changed: %d' %len(modded_changed)
print 'Num reject: %d' %len(reject)
print 'Num eyeball: %d' %len(eyeball)

print 'Num iso, changed, percentage overall: %d %d %.2f' %(len(isos), len(iso_mod),((len(iso_mod) + 7)/5618.0)*100.0)  ##7 sources were discarded
print 'Num dom, changed, percentage overall: %d %d %.2f' %(len(doms), len(dom_mod),(len(dom_mod)/float(len(doms)))*100.0)
print 'Num mul, changed, percentage overall: %d %d %.2f' %(len(muls), len(mul_mod),(len(mul_mod)/float(len(muls)))*100.0)


