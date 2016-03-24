#!/bin/sh
##This does the edited changed options
# ./plot_edits.py --matched_cats=gleam_multi,vlssr,mrc,sumss,nvss,atg20 \
# 	--pref_cats=nvss,sumss,gleam_multi,atg20,mrc,vlssr \
# 	--input_bayes=puma_gleammulti-v-m-s-n-a_outcomes.txt \
# 	--cat_freqs=76~84~92~99~107~115~122~130~143~151~158~166~174~181~189~197~204~212~220~227,74,408,843,1400,5000~8000~20000 \
# 	--prob_thresh=0.8,0.95 \
# 	--epsilon_thresh=0.1 --chi_thresh=10 \
# 	--resolution=00:02:20 \
# 	--choices=reject_reject_choices.txt --tag=reject
##This does the immediate combine all "c" used by use_choices.py
./plot_edits.py --matched_cats=gleam_multi,vlssr,mrc,sumss,nvss,atg20 \
	--pref_cats=nvss,sumss,gleam_multi,atg20,mrc,vlssr \
	--input_bayes=puma_gleammulti-v-m-s-n-a_outcomes.txt \
	--cat_freqs=76~84~92~99~107~115~122~130~143~151~158~166~174~181~189~197~204~212~220~227,74,408,843,1400,5000~8000~20000 \
	--prob_thresh=0.8,0.95 \
	--epsilon_thresh=0.1 --chi_thresh=10 \
	--resolution=00:02:20 \
	--choices=IDR2_reject_sources.txt --tag=reject --combine_all_only
