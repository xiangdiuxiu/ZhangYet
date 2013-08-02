#ifndef __READBINDATA_H__
#define __READBINDATA_H__

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <bitset>
#include "Individual.h"
#include "Locus.h"
using namespace std;

class Parameter{
 public:
  string famfile;
  string bitfile;
  string snpfile; //respose to bitfilename_map;
  
  string phenofile;

  bool qt;
  bool bt;


  string missing_phenotype;

  bool SNP_major;
  
  vector<string> pheno_names;
  int MAX_LINE_LENGTH;
  bool all_pheno;
  bool coding01;

  string outputfile;

  Parameter()
    {
      famfile = "";
      bitfile = "";
      snpfile = "";
      phenofile = "";
      
      qt = false;
      bt = false;
      
      coding01 = false;
      
      missing_phenotype = "-9";
      
      SNP_major = false;
      
      pheno_names.resize(0);
      
      MAX_LINE_LENGTH = 1000;
      all_pheno = false;
      coding01 = false;

    }
  

  void initPar()
  {
  }
};

void extern error(string);
#endif
