#!/bin/sh

for i in {1..25}; do
    sed "s/\b\(r = 1\)\b/r = ${i}/" ./template/BA_template.R > BA_r${i}.R
    Rscript BA_r${i}.R > ./log/output_BA_r${i}.log 2>&1
    rm BA_r${i}.R
done

for i in {1..25}; do
    sed "s/\b\(r = 1\)\b/r = ${i}/" ./template/BA_probit_template.R > BA_probit_r${i}.R
    Rscript BA_probit_r${i}.R > ./log/output_BA_probit_r${i}.log 2>&1
    rm BA_probit_r${i}.R
done