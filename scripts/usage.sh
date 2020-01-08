#!/usr/bin/env bash
# Get some stats on number of times prog used for tracking purposes.

user=$(whoami)
countfile="$(dirname $0)/../resource/count.txt"
now=$(date +'%F %H:%M:%S')

if [[ -e $countfile ]]; then
    last_count=($(tail -n1 $countfile | tr ',' ' '))
    new_tot=$(( ${last_count[-1]} + 1 ))
    echo "$now,$user,1,$new_tot" >> $countfile
else
    echo "Generating a new count file"
    echo "$now,$user,1,1" > $countfile
    exit
fi
