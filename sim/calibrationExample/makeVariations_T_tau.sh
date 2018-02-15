#!/bin/sh

for T in 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0; do
  for tau in 10 12 14 16 18 20 22 24 26 28 30; do
    projName=tmp_T${T}_tau${tau}
    cpinp kalibRef $projName
    perl -i -p -e "s/^.*savety time.*$/${T}/g" $projName.inp
    perl -i -p -e "s/^.*tau0.*$/${tau}/g" $projName.inp
    master $projName
    for xdet in 10000 12000 14000 16000; do
      mv $projName.x$xdet kalib_T${T}_tau${tau}.x$xdet
    done
    rm $projName.*
  done
done
