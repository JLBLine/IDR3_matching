./make_initial.py --matched_cats=comp_v11,vlssr,mrc,sumss,nvss \
	--pref_cats=nvss,sumss,comp_v11,mrc,vlssr \
	--input_bayes=/home/jline/Documents/cataloguing/FHD/carroll/paper/complete_v11+v10/bayes_v11+v10extras-v-m-s-n.txt \
	--cat_freqs=182.435,74,408,843,1400 \
	--prob_thresh=0.8,0.95 \
	--epsilon_thresh=0.1 --chi_thresh=10 \
	--resolution=00:02:18 \
	--output_name=KGS_compv10+v11-v-m-s-n --verbose
	
##Open up make_table_function.py and edit pref_cats etc to have comp_v5 in
# ./make_initial.py --matched_cats=comp_v5,vlssr,mrc,sumss,nvss \
# 	--pref_cats=nvss,sumss,comp_v5,mrc,vlssr \
# 	--input_bayes=bayes_v5-v-m-s-n.txt \
# 	--cat_freqs=182.435,74,408,843,1400 \
# 	--prob_thresh=0.8,0.95 \
# 	--epsilon_thresh=0.1 --chi_thresh=10 \
# 	--resolution=00:02:18 \
# 	--output_name=KGS_compv5-v-m-s-n --verbose #--split=00:02:00
