#! /bin/bash

### decompose ###
bcftools norm -m-both -o [vcf_decompose_out] [vcf_input]

### left-norm ###
bcftools norm -f [reffasta] -o [vcf_decompose_norm_out] [vcf_decompose_out]
