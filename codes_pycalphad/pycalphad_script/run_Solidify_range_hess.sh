#!/bin/bash

echo "Start"

# Generates .txt files with equilibrium composition and double derivatives at equilibrium compositions for each phase over the entire solidification range foran alloy composition 
# INPUT :  1 : tdb file , 2 : Solid Phase , 3 : Liquid Phase, 4 : Component for calculation, 5 :  Alloy composition

python Liquidus_Temp.py alzn_mey.tdb FCC_A1 LIQUID ZN 0.4 
python Thermo_write_Equi_Compositions.py alzn_mey.tdb FCC_A1 LIQUID ZN 0.4 

echo "Done "



