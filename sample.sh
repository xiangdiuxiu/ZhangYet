#!/bin/bash
#./mplink -bfile test/TEST -phenofile phenotype_txt -phenotypes auto_r_se_bend master_al_od_bend -o outputtest
./mplink  -bfile test/TEST -phenofile test/phenotype_txt  -phenotypes auto_r_se_bend  master_al_od_bend -pheno_type 1 1 -covar_name exam_age_bend sex  -o outputtest 
#./mplink  -bfile test/TEST -phenofile phenotype_txt  -phenotypes auto_r_se_bend   -pheno_type 1 1 1 -covar_name exam_age_bend sex master_al_od_bend  -o outputtest 