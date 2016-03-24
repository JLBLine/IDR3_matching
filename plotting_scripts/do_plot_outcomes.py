from subprocess import call
from sys import argv
from glob import glob

def make_plots(bayes,version):
	cmd = "plot_outcomes.py --input_bayes=%s" %bayes
	if version == 'multi':
		cmd += " --cat_freqs=76~84~92~99~107~115~122~130~143~151~158~166~174~181~189~197~204~212~220~227,74,408,843,1400,5000~8000~20000,150"
		cmd += " --matched_cats=gleam_multi,vlssr,mrc,sumss,nvss,atg20,TGSS --pref_cats=nvss,sumss,TGSS,gleam_multi,atg20,mrc,vlssr"
	else:
		cmd += " --cat_freqs=200,74,408,843,1400,5000~8000~20000,150"
		cmd += " --matched_cats=gleam_deep,vlssr,mrc,sumss,nvss,atg20,TGSS --pref_cats=nvss,sumss,TGSS,gleam_deep,atg20,mrc,vlssr,TGSS"
	cmd += " --prob_thresh=0.8,0.95 --epsilon_thresh=0.1 --chi_thresh=10 --resolution=00:02:20 --write=all"

	names = open(argv[1],'r').read().split('\n')

	name_str = ''

	for name in names:
		if name!='':
			name_str += '%s,' %name
		
	name_str = name_str[:-1]
		
	cmd += ' --query=%s' %name_str

	call(cmd,shell=True)

	output_dir = argv[2]

	if output_dir[-1] == '/':
		pass
	else:
		output_dir += '/'

	for name in names:
		if name == '':
			pass
		else:
			cmd = 'mv %s-pumaplot.png %s%s-%s_pumaplot.png' %(name,output_dir,name,version)
			call(cmd,shell=True)

make_plots('../bayes_IDR3multi-v-m-s-n-a-t.txt','multi')
#make_plots('../../deep/bayes_gleamdeep-v-m-s-n-a.txt','deep')
