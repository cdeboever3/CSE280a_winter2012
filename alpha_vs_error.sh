#!/bin/bash

# Compares alpha vs. error


# alpha: 0.01 error: 0.001


python generate_sample_tsv_complex.tsv a_0.01_e_0.001_truth.txt -o a_0.01_e_0.001.tsv -sl 1000000 -er 0.001 -sL 0.01 0.99

python mixed_variant_calling.py a_0.01_e_0.001.tsv -er 0.001 > a_0.01_e_0.001_genotypes.txt
wait



# alpha: 0.01 error: 0.01


python generate_sample_tsv_complex.tsv a_0.01_e_0.01_truth.txt -o a_0.01_e_0.01.tsv -sl 1000000 -er 0.01 -sL 0.01 0.99

python mixed_variant_calling.py a_0.01_e_0.01.tsv -er 0.01 > a_0.01_e_0.01_genotypes.txt
wait



# alpha: 0.01 error: 0.02


python generate_sample_tsv_complex.tsv a_0.01_e_0.02_truth.txt -o a_0.01_e_0.02.tsv -sl 1000000 -er 0.02 -sL 0.01 0.99

python mixed_variant_calling.py a_0.01_e_0.02.tsv -er 0.02 > a_0.01_e_0.02_genotypes.txt
wait



# alpha: 0.05 error: 0.001


python generate_sample_tsv_complex.tsv a_0.05_e_0.001_truth.txt -o a_0.05_e_0.001.tsv -sl 1000000 -er 0.001 -sL 0.05 0.95

python mixed_variant_calling.py a_0.05_e_0.001.tsv -er 0.001 > a_0.05_e_0.001_genotypes.txt
wait



# alpha: 0.05 error: 0.01


python generate_sample_tsv_complex.tsv a_0.05_e_0.01_truth.txt -o a_0.05_e_0.01.tsv -sl 1000000 -er 0.01 -sL 0.05 0.95

python mixed_variant_calling.py a_0.05_e_0.01.tsv -er 0.01 > a_0.05_e_0.01_genotypes.txt
wait



# alpha: 0.05 error: 0.02


python generate_sample_tsv_complex.tsv a_0.05_e_0.02_truth.txt -o a_0.05_e_0.02.tsv -sl 1000000 -er 0.02 -sL 0.05 0.95

python mixed_variant_calling.py a_0.05_e_0.02.tsv -er 0.02 > a_0.05_e_0.02_genotypes.txt
wait



# alpha: 0.1 error: 0.001


python generate_sample_tsv_complex.tsv a_0.1_e_0.001_truth.txt -o a_0.1_e_0.001.tsv -sl 1000000 -er 0.001 -sL 0.1 0.9

python mixed_variant_calling.py a_0.1_e_0.001.tsv -er 0.001 > a_0.1_e_0.001_genotypes.txt
wait



# alpha: 0.1 error: 0.01


python generate_sample_tsv_complex.tsv a_0.1_e_0.01_truth.txt -o a_0.1_e_0.01.tsv -sl 1000000 -er 0.01 -sL 0.1 0.9

python mixed_variant_calling.py a_0.1_e_0.01.tsv -er 0.01 > a_0.1_e_0.01_genotypes.txt
wait



# alpha: 0.1 error: 0.02


python generate_sample_tsv_complex.tsv a_0.1_e_0.02_truth.txt -o a_0.1_e_0.02.tsv -sl 1000000 -er 0.02 -sL 0.1 0.9

python mixed_variant_calling.py a_0.1_e_0.02.tsv -er 0.02 > a_0.1_e_0.02_genotypes.txt
wait



# alpha: 0.25 error: 0.001


python generate_sample_tsv_complex.tsv a_0.25_e_0.001_truth.txt -o a_0.25_e_0.001.tsv -sl 1000000 -er 0.001 -sL 0.25 0.75

python mixed_variant_calling.py a_0.25_e_0.001.tsv -er 0.001 > a_0.25_e_0.001_genotypes.txt
wait



# alpha: 0.25 error: 0.01


python generate_sample_tsv_complex.tsv a_0.25_e_0.01_truth.txt -o a_0.25_e_0.01.tsv -sl 1000000 -er 0.01 -sL 0.25 0.75

python mixed_variant_calling.py a_0.25_e_0.01.tsv -er 0.01 > a_0.25_e_0.01_genotypes.txt
wait



# alpha: 0.25 error: 0.02


python generate_sample_tsv_complex.tsv a_0.25_e_0.02_truth.txt -o a_0.25_e_0.02.tsv -sl 1000000 -er 0.02 -sL 0.25 0.75

python mixed_variant_calling.py a_0.25_e_0.02.tsv -er 0.02 > a_0.25_e_0.02_genotypes.txt
wait

