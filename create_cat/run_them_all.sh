#!/bin/sh

##Convert the initial PUMA output into KGS table format
source make_initial.sh

##Create a table to contain the non-matches
python add_nomatches.py

##Create tables for inspected sources that required no
##changes or simple changes
source use_choices.sh

##Create tables for inspected sources that required
##more modification
source make_modifications.sh

##Do the same for the v5 replacements
cd ./v5_plots
source use_choices.sh
source make_modifications.sh
cd ..

# ##Create a list of sources that need deleting
# ##delete_list.txt is hand edited - bad sources or one off deletes
# ##Rest of sources are gathered from the v5 replacement lists in ./v5_mod
python create_delete.py

##Gather all the tables together and stitch into the
##uber final - also gather the KGS parameters.
source create_master.sh
