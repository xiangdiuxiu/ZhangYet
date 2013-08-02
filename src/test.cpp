#include <iostream>
#include <fstream>
#include <algorithm>
#include <bitset>
#include <sstream>
#include <string>

#include "Individual.h"
#include "Locus.h"
#include "readBinData.h"
#include "model.h"
#include "matrix.h"

#include <gsl/gsl_cdf.h>

using namespace std;
void extern readBinData(vector<Individual*>&, vector<Locus*>&, vector<CSNP*>&);
void extern SNP2Ind(vector<Individual*>&, vector<CSNP*>&);
void extern readPhenoFile(vector<Individual*>&, vector<Locus*>&, vector<CSNP*>&);
void extern sizeMatrix(vector<vector<double> >&, int, int);
void extern multMatrix(vector<vector<double> >&, vector<vector<double> >&, vector<vector<double> >&); 
vector<vector<double> > svd_inverse(vector<vector<double> >&, bool &);
vector<double> var(vector<vector<double> >&);
vector<double> extern glmAssoc(vector<Individual*>&, vector<Locus*>&, vector<CSNP*>&);

int main()
{
  double res = gsl_cdf_chisq_P(4,2);
  cout<<"P 4 2 "<<res<<endl;
  vector<vector<double> > matA;
  matA.resize(4);
  int test[8] = {1,2,3,4,5,6,7,8};
  for(int i=0; i<4; i++)
    matA[i].resize(2,0);
  for(int i=0; i<4; i++)
    {
      for(int j=0; j<2; j++)
	matA[i][j] = test[i*2+j];
    }
  vector< vector<double> > tmatA;
  tmatA = trans(matA);
  multMatrix(matA,tmatA,tmatA);
  vector< vector<double> >::iterator it = tmatA.begin();
  while(it != tmatA.end())
    {
      vector<double>::iterator it2 = (*it).begin();
      while(it2!=(*it).end())
	{
	  cout<<(*it2)<<'\t';
	  it2++;
	}
      cout<<endl;
      it++;
    }

  /*
  for(int i=0; i<3; i++)
    {
      matA[i].resize(3);
      for(int j=0; j<3; j++)
	matA[i][j] = test[3*i+j];
    }
  cout<<"A:"<<endl;
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	cout<<matA[i][j]<<'\t';
      cout<<endl;
    }
  vector<vector<double> > matAA;
  
  bool flag = true;
  vector<vector<double> > result;
  result = svd_inverse(matA, flag);
  
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	cout<<result[i][j]<<'\t';
      cout<<endl;
    }

  multMatrix(result, matA, matAA);
  cout<<"matAA:"<<endl;
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	cout<<matAA[i][j]<<'\t';
      cout<<endl;
    }
  */
  /*
  vector<Individual*> sample;
  vector<Locus*> locus;
  vector<CSNP*> snp;
  readBinData(sample, locus, snp);
  readPhenoFile(sample, locus, snp);
  */
  /*
  vector<CSNP*>::iterator it = snp.begin();
  int i = 1;
  
  while(it!=snp.end())
    {
      vector<bool>::iterator it1 = (*it)->one.begin();
      vector<bool>::iterator it2 = (*it)->two.begin();
      cout<<i<<"th snp"<<endl;
      cout<<"one: "<<endl;
      while(it1!=(*it)->one.end())
	{
	  cout<<(*it1)<<'\t';
	  it1++;
	}
      cout<<"two: "<<endl;
      while(it2!=(*it)->two.end())
	{
	  cout<<(*it2)<<'\t';
	  it2++;
	}
      i++;
      it++;
    }
  */
  
  //SNP2Ind(sample, snp);
  /*
  for(int i=0; i<sample.size(); i++)
    {
      cout<<i<<"th sample: "<<endl;
      vector<bool>::iterator it1 = sample[i]->one.begin();
      while(it1!=sample[i]->one.end())
	{
	  cout<<(*it1)<<'\t';
	  it1++;
	}
    }
  */
  
  //ModelPar mpar;
  //mpar.defineHumanChromsomes();
  //int test2 = 0;
  //Model *lm = new Model;
  //lm->addAdditiveSNP(test2, locus, mpar.chr_sex, mpar.chr_haploid);
  //lm->buildDesignMatrix(sample, locus, snp);
  //lm->setY(sample);
  /*
  vector<double> test_res = glmAssoc(sample, locus, snp);
  cout<<"result: "<<endl;
  for(int i=0; i<test_res.size(); i++)
    cout<<test_res[i]<<'\t';
  cout<<endl;
  */
  /*
  vector<vector<double> > res;
  vector<vector<double> > tY ;
  tY = trans(lm->Y);
  cout<<"Y: r: "<<lm->Y.size()<<"c:"<<lm->Y[0].size()<<endl;
  cout<<"X: r: "<<lm->X.size()<<"c:"<<lm->X[0].size()<<endl;
  multMatrix(tY, lm->X, res);
  cout<<res[0][0]<<endl;
  cout<<res[1][0]<<endl;

  vector<double> res1 = var(lm->Y);
  vector<double> res2 = var(tY);
  
  vector<double>::iterator it;
  cout<<"lm-Y: "<<endl;
  for(it=res1.begin(); it!=res1.end(); it++)
    cout<<(*it)<<'\t';
  cout<<endl;
  cout<<"tY: "<<endl;
  for(it=res2.begin(); it!=res2.end(); it++)
    cout<<(*it)<<'\t';
  cout<<endl;
  */
  /*
  //Final test!

  vector<vector<double> > YG;
  multMatrix(tY, lm->X, YG);
  vector<vector<double> > tYG;
  tYG = trans(YG);
  vector<vector<double> > YY;
  multMatrix(tY, lm->Y, YY);
  vector<vector<double> > YY_inverse;
  bool flag = true;
  YY_inverse = svd_inverse(YY, flag);
  vector<vector<double> > t_res;
  multMatrix(tYG, YY_inverse, t_res);
  vector<vector<double> > f_res;
  multMatrix(t_res, YG, f_res);
  cout<<"Final test: "<<endl;
  cout<<f_res[0][0]<<endl;
  */
  /*
  vector<vector<double> >::iterator itY;
  vector<double>::iterator itY1;
  for(itY=lm->Y.begin(); itY!=lm->Y.end(); itY++)
    {
      for(itY1=(*itY).begin(); itY1!=(*itY).end(); itY1++)
	{
	  cout<<(*itY1)<<'\t';
	}
      cout<<endl;
    }
  */
  /*
  cout<<"lm->X:"<<endl;
  for(int i=0; i<lm->X.size(); i++)
    {
      for(int j=0; j<lm->X[i].size(); j++)
	cout<<lm->X[i][j]<<'\t';
      cout<<endl;
    }
  cout<<"Y:"<<endl;
  cout<<lm->Y[0][0]<<endl;
  
  */

  /*
    vector<Locus> test;
    string filename = "TEST.bim";
    readSNP(filename, test);
    string FamFile = "TEST.fam";
    
    vector<Individual*> sample;
    
    readFamFile(FamFile, sample);
    /*
    vector<Locus>::iterator it;
    for(it=test.begin(); it!=test.end(); it++)
        cout<<it->chr<<'\t'<<it->name<<'\t'<<it->pos<<it->freq<<endl;
    */
    /*
    vector<Individual*>::iterator it_ind;
    for(it_ind=sample.begin(); it_ind!=sample.end(); it_ind++)
        cout<<(*it_ind)->fid<<'\t'<<(*it_ind)->phenotypes[0]<<endl;
    */
    /*
    string PhenoFile = "phenotype_txt";
    readPhenoFile(PhenoFile, sample);
    vector<Individual*>::iterator it_ind;
    for(it_ind=sample.begin(); it_ind!=sample.end(); it_ind++)
      if((*it_ind)->missing == false)
	cout<<(*it_ind)->fid<<'\t'<<(*it_ind)->phenotypes[0]<<'\t'<<(*it_ind)->phenotypes[1]<<endl;
    */
  
  
    return 0;
}
