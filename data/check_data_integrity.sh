#! /bin/bash

# arguments
# $1: path to the data directory
# $2: whether to fetch hitran data
# $3: whether to fetch jupiter data
# $4: whether to fetch exogcm data

if [ -z "$1" ]; then
    echo "Please provide the path to the data directory"
    exit 1
fi

if [ "$2" == "ON" ]; then
  #$1/fetch_hitran2012.sh $1
  $1/fetch_hitran2020.sh $1
fi

if [ "$3" == "ON" ]; then
  $1/fetch_hitran2020.sh $1
  $1/fetch_jup_midir_H2broaden.sh $1
  $1/fetch_jup_atm_moses_modelc.sh $1
  $1/fetch_H2-He-cia.sh $1
fi

if [ "$4" == "ON" ]; then
  $1/fetch_exogcm_opacity.sh $1
fi
