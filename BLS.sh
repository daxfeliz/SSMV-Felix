#!/bin/bash


path=${PWD}
data=$1


qmin=$2
qmax=$3
Pmin=$4
Pmax=$5
Nfreq=$6
Nbins=$7
Npeaks=$8
TIC=$9

vartools -i $data -ascii -oneline -BLS q $qmin $qmax $Pmin $Pmax $Nfreq $Nbins 0 $Npeaks 1 $path/outdir 1 $path/outdir 0 nobinnedrms > TESSdata/${TIC}_BLS_output.txt

echo " "
echo " " #just extra spaces for aesthetic reasons 


