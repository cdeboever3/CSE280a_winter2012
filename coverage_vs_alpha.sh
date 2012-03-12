#!/bin/bash

# Compares coverage vs. alpha


# coverage: 25 alpha: 0.01


python generate_sample_tsv_complex.tsv c_25_a_0.01_truth.txt -o c_25_a_0.01.tsv -ac 25 -sl 1000000 -sL 0.01 0.99

python mixed_variant_calling.py c_25_a_0.01.tsv -er 0.001 > c_25_a_0.01_genotypes.txt
wait



# coverage: 25 alpha: 0.05


python generate_sample_tsv_complex.tsv c_25_a_0.05_truth.txt -o c_25_a_0.05.tsv -ac 25 -sl 1000000 -sL 0.05 0.95

python mixed_variant_calling.py c_25_a_0.05.tsv -er 0.001 > c_25_a_0.05_genotypes.txt
wait



# coverage: 25 alpha: 0.1


python generate_sample_tsv_complex.tsv c_25_a_0.1_truth.txt -o c_25_a_0.1.tsv -ac 25 -sl 1000000 -sL 0.1 0.9

python mixed_variant_calling.py c_25_a_0.1.tsv -er 0.001 > c_25_a_0.1_genotypes.txt
wait



# coverage: 25 alpha: 0.25


python generate_sample_tsv_complex.tsv c_25_a_0.25_truth.txt -o c_25_a_0.25.tsv -ac 25 -sl 1000000 -sL 0.25 0.75

python mixed_variant_calling.py c_25_a_0.25.tsv -er 0.001 > c_25_a_0.25_genotypes.txt
wait



# coverage: 50 alpha: 0.01


python generate_sample_tsv_complex.tsv c_50_a_0.01_truth.txt -o c_50_a_0.01.tsv -ac 50 -sl 1000000 -sL 0.01 0.99

python mixed_variant_calling.py c_50_a_0.01.tsv -er 0.001 > c_50_a_0.01_genotypes.txt
wait



# coverage: 50 alpha: 0.05


python generate_sample_tsv_complex.tsv c_50_a_0.05_truth.txt -o c_50_a_0.05.tsv -ac 50 -sl 1000000 -sL 0.05 0.95

python mixed_variant_calling.py c_50_a_0.05.tsv -er 0.001 > c_50_a_0.05_genotypes.txt
wait



# coverage: 50 alpha: 0.1


python generate_sample_tsv_complex.tsv c_50_a_0.1_truth.txt -o c_50_a_0.1.tsv -ac 50 -sl 1000000 -sL 0.1 0.9

python mixed_variant_calling.py c_50_a_0.1.tsv -er 0.001 > c_50_a_0.1_genotypes.txt
wait



# coverage: 50 alpha: 0.25


python generate_sample_tsv_complex.tsv c_50_a_0.25_truth.txt -o c_50_a_0.25.tsv -ac 50 -sl 1000000 -sL 0.25 0.75

python mixed_variant_calling.py c_50_a_0.25.tsv -er 0.001 > c_50_a_0.25_genotypes.txt
wait



# coverage: 75 alpha: 0.01


python generate_sample_tsv_complex.tsv c_75_a_0.01_truth.txt -o c_75_a_0.01.tsv -ac 75 -sl 1000000 -sL 0.01 0.99

python mixed_variant_calling.py c_75_a_0.01.tsv -er 0.001 > c_75_a_0.01_genotypes.txt
wait



# coverage: 75 alpha: 0.05


python generate_sample_tsv_complex.tsv c_75_a_0.05_truth.txt -o c_75_a_0.05.tsv -ac 75 -sl 1000000 -sL 0.05 0.95

python mixed_variant_calling.py c_75_a_0.05.tsv -er 0.001 > c_75_a_0.05_genotypes.txt
wait



# coverage: 75 alpha: 0.1


python generate_sample_tsv_complex.tsv c_75_a_0.1_truth.txt -o c_75_a_0.1.tsv -ac 75 -sl 1000000 -sL 0.1 0.9

python mixed_variant_calling.py c_75_a_0.1.tsv -er 0.001 > c_75_a_0.1_genotypes.txt
wait



# coverage: 75 alpha: 0.25


python generate_sample_tsv_complex.tsv c_75_a_0.25_truth.txt -o c_75_a_0.25.tsv -ac 75 -sl 1000000 -sL 0.25 0.75

python mixed_variant_calling.py c_75_a_0.25.tsv -er 0.001 > c_75_a_0.25_genotypes.txt
wait

