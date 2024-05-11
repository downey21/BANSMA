#!/bin/sh

for i in {1..12}; do
    sed "s/\b\(pp = 1\)\b/pp = ${i}/" ./template/fit_outcome_Fulv_template.R > fit_outcome_Fulv_pp${i}.R
    sed "s/\b\(pp = 1\)\b/pp = ${i}/" ./template/fit_outcome_Palbo_template.R > fit_outcome_Palbo_pp${i}.R
    sed "s/\b\(pp = 1\)\b/pp = ${i}/" ./template/fit_outcome_Tamo_template.R > fit_outcome_Tamo_pp${i}.R
done

for i in {1..12}; do
    Rscript fit_outcome_Fulv_pp${i}.R > ./log/output_fit_outcome_Fulv_pp${i}.log 2>&1 &
    Rscript fit_outcome_Palbo_pp${i}.R > ./log/output_fit_outcome_Palbo_pp${i}.log 2>&1 &
    Rscript fit_outcome_Tamo_pp${i}.R > ./log/output_fit_outcome_Tamo_pp${i}.log 2>&1 &
done

wait

for i in {1..12}; do
    rm fit_outcome_Fulv_pp${i}.R
    rm fit_outcome_Palbo_pp${i}.R
    rm fit_outcome_Tamo_pp${i}.R
done