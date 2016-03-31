from subprocess import call
from glob import glob

v2_names = open('IDR3multi_reject-eyeball_sources_v2.txt').read().split('\n')
v2_names = [name for name in v2_names if name != '']

current_plots = glob('../plots/*extended*png')

print len(v2_names), len(current_plots)

current_names = [string.split('/')[-1].split('_')[0] for string in current_plots]
delete_names = [name for name in current_names if name not in v2_names]

print len(delete_names)

#for name in delete_names:
	#cmd = "rm ../plots/%s*png" %name
	#call(cmd, shell=True)
