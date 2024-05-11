#!/bin/sh

for i in {1..12}; do
    sed "s/\b\(pp = 1\)\b/pp = ${i}/" ./template/fit_x_template.R > fit_x_pp${i}.R
done

for i in {1..12}; do
    Rscript fit_x_pp${i}.R > ./log/output_fit_x_pp${i}.log 2>&1
done

for i in {1..12}; do
    rm fit_x_pp${i}.R
done