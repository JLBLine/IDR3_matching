
python create_master.py --base_table=KGS_compv10+v11-v-m-s-n.vot \
	--tags=extreme,multiple --delete_names=delete_choices.txt \
	--extras=eyeball_choices.vot,eyeball_edited.vot,\
reject_choices.vot,reject_edited.vot,no_matches_PUMA.vot,\
./v5_plots/v5_extreme_choices.vot,./v5_plots/v5_extreme_edited.vot,\
./v5_plots/v5_eyeball_choices.vot,./v5_plots/v5_eyeball_edited.vot,\
./v5_plots/v5_multiple_choices.vot,./v5_plots/v5_multiple_edited.vot,./v5_plots/v5_dominant_choices.vot,\
./v5_plots/v5_reject_choices.vot,./v5_plots/v5_reject_edited.vot,extras_edited.vot \
	--output_name=KATALOGSS_RELEASE

# python create_master.py --base_table=KGS_compv10+v11-v-m-s-n.vot \
# 	--tags=extreme --delete_names=delete_list.txt \
# 	--extras=eyeball_choices.vot,eyeball_edited.vot \
# 	--output_name=KATALOGSS_RELEASE
	
