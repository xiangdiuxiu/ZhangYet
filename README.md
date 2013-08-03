mplink
======

#INSTALL
Just run the install.sh
    ./install.sh
If it is successfully installed, you will get a binary file mplink.
You can run the test.sh to check that.
    ./test.sh

#USAGE
    ./mplink -bfile [binary file] -phenofile [phenofile] -phenotypes [phenotypes] -o [output file]
    -bfile : include .bed, .fam, .bim file;
    -phenofile : 
    -phenotypes : if you have more than one phenotypes, seperate them by space;
    -o : all the results will record in the output file;

#EXAMPLE
    ./mplink -bfile TEST1 -phenofile phenotype_txt -phenotypes auto_r_se_bend master_al_od_bend -o ouputtest