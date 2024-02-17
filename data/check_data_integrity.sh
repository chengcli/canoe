#! /bin/bash

#$1/fetch_hitran2012.sh $1
$1/fetch_hitran2020.sh $1
$1/fetch_jup_midir_H2broaden.sh $1
$1/fetch_jup_atm_moses_modelc.sh $1
$1/fetch_H2-He-cia.sh $1
$1/fetch_exogcm_opacity.sh $1
