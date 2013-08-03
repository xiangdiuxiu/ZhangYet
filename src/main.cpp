#include "main.h"

int main(int argc, char *argv[])
{
  //parameters
  vector<string> args(argv+1, argv+argc);
  string output = getOutput(args);
  
  //Data
  vector<Individual*> sample;
  vector<Locus*> locus;
  vector<CSNP*> SNP;
  
  //read data from file
  readBinData(sample, locus, SNP, args);
  cout<<"sample:"<<sample.size()<<endl;
  cout<<"locus:"<<locus.size()<<endl;
  cout<<"SNP:"<<SNP.size()<<endl;
  readPhenoFile(sample, locus, SNP);
  cout<<"Finishing reading phenofile!"<<endl;
  //compute the statisc;
  vector<double> res;
  //output test
  //string output = argv[1];
  res = glmAssoc(sample, locus, SNP, output);
  for(vector<double>::iterator it=res.begin(); it!=res.end(); it++)
    cout<<(*it)<<'\t';
  cout<<endl;
  
  return 0;
}
