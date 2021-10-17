#!/bin/bash

if [ $# -lt 1 ]
then
	echo "./species_check.sh <speccies_name> "	
	echo "[note]: output Eukaryotic or Prokaryotic "
	exit 1
fi

# ----- input argument ----- #
species_name=$1


# ----- process ------- #
kingdom_name=`grep "$species_name" data/species_out_no_TaxTree | head -n1 | awk -F"\t" '{print $NF}' | awk '{print $NF}' `
check_name=`grep "$kingdom_name" data/Eukaryotic_or_Prokaryotic | awk -F"\t" '{print $NF}'`
echo "$check_name"

# ----- success exit ---- #
exit 0

