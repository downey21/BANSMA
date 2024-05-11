#!/bin/sh

for i in {1..25}; do
    sed "s/\b\(r = 1\)\b/r = ${i}/" ./template/ER_template.R > ER_r${i}.R
    Rscript ER_r${i}.R > ./log/output_ER_r${i}.log 2>&1
    rm ER_r${i}.R
done

for i in {1..25}; do
    sed "s/\b\(r = 1\)\b/r = ${i}/" ./template/ER_probit_template.R > ER_probit_r${i}.R
    Rscript ER_probit_r${i}.R > ./log/output_ER_probit_r${i}.log 2>&1
    rm ER_probit_r${i}.R
done