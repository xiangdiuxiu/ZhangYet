#include "model.h"

void extern error(string);

void extern SNP2Ind(vector<Individual*>& sample, vector<CSNP*>& snp)
{
  //cout<<"Converting data ot Individual major format\n";
  
  vector<Individual*>::iterator person = sample.begin();
  
  while(person!=sample.end())
    {
      (*person)->one.clear();
      (*person)->two.clear();
      person++;
    }
  
  vector<CSNP*>::iterator s = snp.begin();
  while(s!=snp.end())
    {
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      vector<Individual*>::iterator person = sample.begin();
      
      while(person!=sample.end())
	{
	  (*person)->one.push_back(*i1);
	  (*person)->two.push_back(*i2);
	  i1++;
	  i2++;
	  person++;
	}
      delete(*s);
      s++;
    }
  snp.clear();
  //par.SNP_major = false;  //PQ mark for parameter question
  cout<<"Finishing converting SNP to individual"<<endl;
}

void Model::addAdditiveSNP(int a, vector<Locus*>& locus, vector<bool>& chr_sex, vector<bool>& chr_haploid)
{
  additive.push_back(a);
  if(chr_sex[locus[a]->chr])
    xchr.push_back(true);
  else
    xchr.push_back(false);
  
  if(chr_haploid[locus[a]->chr])
    haploid.push_back(true);
  else
    haploid.push_back(false);

  type.push_back(ADDITIVE);
  order.push_back(additive.size()-1);
}


void Model::setMissing(vector<Individual*>& sample)
{
  miss.clear();
  miss.resize(sample.size(), false);
  for(int i=0; i<sample.size(); i++)
    if(sample[i]->missing)
      miss[i] = true;
}

double Model::buildAdditive(Individual * person , int snp )
{

  // Additive effects (assuming individual-major mode)
      
  int s = additive[snp];  

  bool i1 = person->one[s];
  bool i2 = person->two[s];
      
  if ( xchr[snp] )
    {

      /////////////////////////
      // X chromosome coding
      
      if ( person->sex ) // male
	{
	  if ( i1 ) 
	    {
	      if ( ! i2 ) 
		{
		  skip = true;
		  return 0;
		}
	      else
		return mA; 	      
	    }
	  else
	    {
	      if ( i2 )
		{
		  // This should not happen...
		  skip = true;
		  return 0;
		}		  
	      else
		return mB; 
	    }
	    }
      else // female x-chromosome
	{
	  if ( i1 ) 
	    {
	      if ( ! i2 ) 
		{
		  skip = true;
		  return 0;
		}
	      else
		return mAA; 	      
	    }
	  else
	    {
	      if ( i2 )
		return mAB; // het 
	      else
		return mBB; // hom
	    }		  
	}
      
    }
  else if ( haploid[snp] )
    {
      
      ///////////////////
      // Haploid coding
      
      if ( i1 ) 
	{
	  if ( ! i2 ) 
	    {
	      skip = true;
	      return 0;
	    }
	  else
	    return 0; 	      
	}
      else
	{
	  if ( i2 )
	    {
	      // haploid het
	      skip = true;
	      return 0;
	    }		  
	  else
	    return 1;
	    }
      
    }
  else 
    {
      
      ///////////////////////
      // Autosomal coding
      
      if ( i1 ) 
	{
	  if ( ! i2 ) 
	    {
	      skip = true;
	      return 0;
	    }
	  else
	    return mAA;
	}
      else
	{
	  if ( i2 )
	    return mAB; // het 
	  else
	    return mBB; // hom
	}
      
    }
  
}


void Model::setY(vector<Individual*>& sample)
{
  Y.clear();
  Y.resize(0);
  int nc;
  
  for(int i=0; i<sample.size(); i++)
    {
      if(!miss[i])
	nc = sample[i]->phenotypes.size();
    }
  vector<double> trow(nc,0);
  for(int i=0; i<sample.size(); i++)
    {
      if(!miss[i])
	{
	  for(int j=0; j<nc; j++)
	    trow[j] = sample[i]->phenotypes[j];
	  Y.push_back(trow);
	}
    }
  cout<<"Finishing setting Y"<<endl;
}

void Model::buildDesignMatrix(vector<Individual*>& sample, vector<Locus*>& locus, vector<CSNP*>& snp)
{
  int np = additive.size();
  
  
  for(int i=0; i<sample.size(); i++)
    {
      Individual * person = sample[i];
      if(miss[i])
	continue;

      skip = false;
      vector<double> trow(np);

      for(int p=0; p<np; p++)
	{
	  int pTtpye = type[p];
	  switch(pTtpye)
	    {
	    case ADDITIVE:
	      trow[p] = buildAdditive(person, order[p]);
	    }
	} // end of or p=0 p<np
      
      if(skip)
	{
	  miss[i] = true;
	  skip = false;
	  continue;
	}
      
      X.push_back(trow);
    } // end of for i=0 i<sample.size
  cout<<"Finishing building design matrix!"<<endl;
 
}

vector<double> extern var(vector<vector<double> >& X)
{
  int nr = X.size();
  int nc = X[0].size();
  vector<double> res, ave;
  res.resize(nc, 0);
  ave.resize(nc, 0);
  for(int j=0; j<nc; j++)
    {
      for(int i=0; i<nr; i++)
	ave[j] += X[i][j];
      ave[j] /= nr;
    }
  for(int j=0; j<nc; j++)
    {
      for(int i=0; i<nr; i++)
	res[j] += (X[i][j] - ave[j])*(X[i][j] - ave[j]);
      res[j] /= nr;
    }
  return res;
}

//temp test function

vector<double> YG(vector<vector<double> >& Y, vector<vector<double> >& G)
{
  int nc = Y[0].size();
  int nr = Y.size();
  vector<double> res;
  res.resize(nc, 0);
  for(int i=0; i<nr; i++)
    for(int j=0; j<nc; j++)
      res[j] += Y[i][j]*G[i][0];
  return res;
}

double inner(vector<double> & X)
{
  double res = 0;
  for(int i=0; i<X.size(); i++)
    res += X[i]*X[i];
  return res;
}

double YY(vector<vector<double> >& Y)
{
  double res = 0;
  for(int i=0; i<Y.size(); i++)
    res += inner(Y[i]);
  return res;
  
}


vector<vector<double> > scale(vector<vector<double> >& y)
{
  int nc = y[0].size();
  int nr = y.size();

  vector<double> max;
  vector<double> min;
  vector<double> ave;
  max.resize(0);
  min.resize(0);
  ave.resize(nc, 0);

  vector<vector<double> > res;
  res.resize(nr);
  for(int i=0; i<nr; i++)
    res[i].resize(nc,0);

  for(int j=0; j<nc; j++)
    {
      double temp_max=y[0][j];
      double temp_min=y[0][j];
      ave[j] += y[0][j];
      int i = 1;
      while(i<nr)
	{
	  if(y[i][j]>temp_max)
	    temp_max = y[i][j];
	  if(y[i][j]<temp_min)
	    temp_min = y[i][j];
	  ave[j] += y[i][j];
	  i++;
	}
      max.push_back(temp_max);
      min.push_back(temp_min);
      ave[j] /= nr;
    }

  
  for(int j=0; j<nc; j++)
    {
      double d = max[j] - min[j];
      for(int i=0; i<nr; i++)
	res[i][j] = (y[i][j]-ave[j])/ave[j];
      
    }
  return res;
}

vector<double> extern glmAssoc(vector<Individual*>& sample, vector<Locus*>& locus, vector<CSNP*>& snp, string outputfile)
{
  int ntests = snp.size();
  SNP2Ind(sample, snp);
  vector<double> results(ntests);
  
  ModelPar mpar;
  mpar.defineHumanChromsomes();

  //output
  ofstream ASC;
  ASC.open(outputfile.c_str(), ios::out);
  ASC << setw(4) <<"CHR" << " "
      << setw(10) <<"SNP" << " "
      << setw(4) << "BP" << " "
      << setw(10) << "STAT"<< " "<<"\n";

  //compute each snp
  for(int l=0; l<ntests; l++)
    {
      bool X = false;
      bool automaticSex = false;
      
      Model *lm = new Model();
      
      bool genotypic = false;
      string mainEffect = "ADD";

      lm->setMissing(sample);
      lm->addAdditiveSNP(l, locus, mpar.chr_sex, mpar.chr_haploid);
      lm->buildDesignMatrix(sample, locus, snp);
      lm->setY(sample);
      

      /*
      vector<vector<double> > Y;
      Y = scale(lm->Y);
      vector<double> yg;
      yg = YG(Y, lm->X);
      double yy = YY(Y);
      vector<double> varX = var(lm->X);
      results[l] = inner(yg)/yy/varX[0];
      */
      
      vector<vector<double> > Y, YG, tYG, YY, invYY, sY;
      vector<vector<double> > t_res, f_res;
      vector<double> varX;
      vector<vector<double> > temp;

      // compute the staticstic 
      sY = scale(lm->Y);
      Y = trans(sY);
      multMatrix(Y, lm->X, YG);
      tYG = trans(YG);
      multMatrix(Y, sY, YY);
      bool flag = true;
      temp = svd_inverse(YY, flag);
      invYY = trans(temp);
      // temp test

      multMatrix(tYG, temp, t_res);
      multMatrix(t_res, YG, f_res);
      varX = var(lm->X);
      results[l] = f_res[0][0]/varX[0];
      locus[l]->stat = results[l];
      lm->displayResults(ASC, locus[l]);
      delete lm;
      /*test
      cout<<"snp "<<l<<endl;
      for(int i=0; i<lm->X.size(); i++)
	cout<<"sample "<<i<<'\t'<<lm->X[i][0]<<endl;
      cout<<endl;
      */
    } // end of for int l=0

  return results;
}

Model::Model()
{
  order.clear();
  type.push_back(ADDITIVE);
  order.push_back(0);
  
  xchr.resize(0);
  haploid.resize(0);
  
  X.resize(0);
  Y.resize(0);

  mAA = 0;
  mAB = 1;
  mBB = 2;
  
  mA = 0;
  mB = 1;
}

void Model::displayResults(ofstream & OUT, Locus * loc)
{
  OUT << setw(4) << loc->chr << " "
      << setw(10)<< loc->name << " "
      << setw(4) << loc->bp  << " "
      << setw(10)<< loc->stat<< " ";
  OUT << "\n";
}
