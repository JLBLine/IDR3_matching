#!/usr/bin/python
from astropy.io.votable import parse as vot_parse
try:
	import pyfits as fits
except ImportError:
	from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import optparse
from statsmodels.robust.scale import mad
import statsmodels.api as sm


def open_table(name):
	if '.vot' in name:
		table = vot_parse(name,pedantic=False).get_first_table()
		return table.array
	elif '.fits' or '.FITS' in name:
		table = fits.open(name,pedantic=False)
		return table[1].data
	else:
		sys.exit('Entered table must either be VOTable or FITS')
		
eyeb_tab = open_table('eyeball_choices.vot')
rej_tab = open_table('reject_choices.vot')
puma_cat = open_table('/home/jline/Documents/cataloguing/puma_GLEAM/IDR2/internal_release/puma_gleammulti-v-m-s-n-a.fits')
	
SI_reject = rej_tab['SI']
SI_eyeball = eyeb_tab['SI']
SI_puma = puma_cat['SI']

fig = plt.figure(figsize=(5,6))
ax = fig.add_subplot(111)

colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']


from numpy import *
width = 0.15

bins = arange(-2.3 - width, 0.3 + width, width)

from statsmodels.robust.scale import mad

#ax.hist(SI_puma,color='k',,histtype='stepfilled',linewidth=3,normed=True,alpha=0.2,bins=100)
ax.hist(SI_puma,color='k',histtype='step',linewidth=1,normed=True,bins=bins,hatch='\\',label='PUMA accept\n (%.2f $\pm$ %.2f)' %(median(SI_puma),mad(SI_puma)))
ax.hist(SI_reject,color=colours[0],histtype='stepfilled',linewidth=1,normed=True,alpha=0.2,bins=bins)
ax.hist(SI_reject,color=colours[0],label='PUMA reject\n (%.2f $\pm$ %.2f)' %(median(SI_reject),mad(SI_reject)),histtype='step',linewidth=4,normed=True,bins=bins)
#ax.hist(SI_eyeball,color=colours[1],histtype='stepfilled',linewidth=1,normed=True,alpha=0.2,bins=bins)
ax.hist(SI_eyeball,color=colours[1],label='PUMA eyeball\n (%.2f $\pm$ %.2f)' %(median(SI_eyeball),mad(SI_eyeball)),histtype='step',linewidth=4,normed=True,bins=bins)


##ax.hist(SI_puma,color='k',,histtype='stepfilled',linewidth=3,normed=True,alpha=0.2,bins=100)
#ax.hist(SI_puma,color='k',histtype='step',linewidth=1,normed=True,bins=140,hatch='\\',label='PUMA accept')
#ax.hist(SI_reject,color=colours[0],histtype='stepfilled',linewidth=1,normed=True,alpha=0.2,bins=15)
#ax.hist(SI_reject,color=colours[0],label='PUMA reject',histtype='step',linewidth=4,normed=True,bins=15)
##ax.hist(SI_eyeball,color=colours[2],label='PUMA eyeball',histtype='stepfilled',linewidth=1,normed=True,linestyle='dashed',alpha=0.2)
#ax.hist(SI_eyeball,color=colours[1],label='PUMA eyeball',histtype='step',linewidth=4,normed=True,bins=16)



ax.set_xlim(-2.4,0.2)

ax.legend(loc='upper left')

ax.set_xlabel('Spectral Index')
ax.set_ylabel('Number of source (normalised)')

plt.tight_layout()
fig.savefig('SI_comparison_Brian.png',bbox_inches='tight')
plt.close()