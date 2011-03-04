#!/bin/bash

dtau=0.0628318530718

i=1
while [ $i -le 100 ]
do
  tau=$(echo "$i*$dtau" |bc)
  echo "Calculating tau=$tau, $i"
  sed "s/tau = 0.01/tau = $tau/" src/fel.f > temp11_$i.f
  sed "s/tau_0.01.dat/tau_${i}.dat/" temp11_$i.f > temp12_$i.f
  g77 -o fel1 temp12_$i.f $cernlib
  ./fel1
  rm temp11_$i.f temp12_$i.f fel1
  i=$((i+1))
done

