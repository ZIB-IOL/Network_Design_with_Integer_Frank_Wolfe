#!/bin/bash

#ds_names="Berlin-Friedrichshain Berlin-Tiergarten Anaheim"

penalties="1000.0"
powers="1.5"
methods="IFW IFW-I IFW-B BNDLMO-PR" 
frac_rem="0.01"
tot_scen="1 2"
solvers="SCIP HiGHS"


for ds_name in $ds_names; do
    for penalty in $penalties; do
        for pv in $powers; do
            for alg in $methods; do
                for fr in $frac_rem; do
                    for ts in $tot_scen; do
                        if [ "$alg" == "IFW" ]; then
                            for slvr in $solvers; do
                                julia run_stochastic_blmo.jl $ds_name $penalty $pv $alg $fr $ts $slvr
                            done
                        else
                            slvr="SCIP"
                            echo "Dataset $ds_name penalty $penalty power $pv algorithm $alg fraction removed $fr scenarios $ts solver $slvr"
                            julia run_stochastic_blmo.jl $ds_name $penalty $pv $alg $fr $ts $slvr
                        fi
                    done
                done
            done
        done
    done
done