#ifndef __MODEL_H__
#define __MODEL_H__

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <bitset>
#include <cmath>
#include "Individual.h"
#include "Locus.h"
#include "matrix.h"
#include "readBinData.h"

using namespace std;

class ModelPar 
{
 public:
  vector<bool> chr_sex;
  vector<bool> chr_haploid;
  
  ModelPar()
    {
      chr_sex.resize(0);
      chr_haploid.resize(0);
    }
  void defineHumanChromsomes()
  {
    chr_sex.resize(22+2+2+1);
    chr_haploid.resize(22+2+2+1);
    for(int i=0; i<=22; i++)
      {
	chr_haploid[i] = chr_sex[i] = false;
      }
    chr_sex[23] = true;
    chr_haploid[23] = false;
    chr_sex[24] = false;
    chr_haploid[24] = true;
    chr_sex[25] = false;
    chr_haploid[25] = false;

    chr_sex[26] = false;
    chr_haploid[26] = true;
  }

};


class Model
{
 public:
  enum terms {ADDITIVE};
  vector<int> additive;;
  int mAA;
  int mAB;
  int mBB;

  double mA, mB;
  
  vector< vector<double> > X;
  vector< vector<double> > Y; //This needs to be modified;
  vector<int> order;
  vector<int> type;
  
  vector<bool> xchr;
  vector<bool> haploid;
  vector<bool> miss;

  bool skip;

  Model();
  void addAdditiveSNP(int, vector<Locus*>&, vector<bool>&, vector<bool>&);
  double buildAdditive(Individual*, int snp);
  void buildDesignMatrix(vector<Individual*>&, vector<Locus*>&, vector<CSNP*>&);
  void setMissing(vector<Individual*>&);
  void setY(vector<Individual*>&);
  void displayResults(ofstream&, Locus*);
};



#endif
