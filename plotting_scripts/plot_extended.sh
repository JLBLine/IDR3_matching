#!/bin/sh
python plot_extended.py	--matched_cats=gleam_multi,vlssr,mrc,sumss,nvss,atg20,TGSS \
	--input_bayes=../puma_IDR3multi-v-m-s-n-a-t-eyeball.txt \
	--cat_freqs=76~84~92~99~107~115~122~130~143~151~158~166~174~181~189~197~204~212~220~227,74,408,843,1400,5000~8000~20000,150 \
	--prob_thresh=0.8,0.95 \
	--epsilon_thresh=0.1 --chi_thresh=10 \
	--resolution=00:02:20 \
	--extended_names=../IDR3multi_reject-eyeball_sources.txt \
	--output_dir=../plots
