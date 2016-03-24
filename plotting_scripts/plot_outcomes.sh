#!/bin/sh
 plot_outcomes.py --matched_cats=gleam_multi,vlssr,mrc,sumss,nvss,atg20,TGSS \
 	--pref_cats=nvss,sumss,TGSS,gleam_multi,atg20,mrc,vlssr \
 	--input_bayes=../bayes_IDR3multi-v-m-s-n-a-t.txt \
 	--cat_freqs=76~84~92~99~107~115~122~130~143~151~158~166~174~181~189~197~204~212~220~227,74,408,843,1400,5000~8000~20000,150 \
 	--prob_thresh=0.8,0.95 \
 	--epsilon_thresh=0.1 --chi_thresh=10 \
 	--resolution=00:02:20 \
 	--query=GLM082231+0556 \
  	--write=all
