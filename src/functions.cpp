#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "Individual.h"
#include "Locus.h"
using namespace std;

void extern readSNP(string filename, vector<Locus> & ordered)
{
    ifstream MAP(filename.c_str(), ios::in);
    MAP.clear();
    int c=0;
    while(!MAP.eof())
    {
        Locus *loc = new Locus;
        MAP >> loc->chr
            >> loc->name
            >> loc->pos
            >> loc->bp
            >> loc->allele1
            >> loc->allele2;
        if ( MAP.eof() )
        {
            delete loc;
            continue;
        }
        else if ( MAP.fail() )
        {
            delete loc;
            cout<<"Problem reading BIM file"<<endl;
        }
        loc->freq = c++;
        
        if (loc->name!="")
        {
            ordered.push_back(*loc);
        }
        else
            delete loc;
    }
    
    MAP.clear();
    MAP.close();
    
    stable_sort(ordered.begin(), ordered.end());
    c = 0;
    for (int i=0; i<ordered.size(); i++)
    {
        ordered[i].bp = (int) ordered[i].freq;
        //ordered[i].chr = l;
        ordered[i].freq = c++;
    }
    stable_sort(ordered.begin(), ordered.end());
}

/*
extern string pheno_name; //maybe some other day we may create a class;
extern int MAX_LINE_LENGTH = 1000;
extern string miss_phenotype = "-9";
*/
template <class T>
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss>>f>>t).fail();
}

void readFamFile(string filename, vector<Individual*>& sample)
{
  bool qt = false;
  bool bt = true;
  ifstream PED;
  PED.open(filename.c_str());
  PED.clear();

  vector<Individual*> ambiguous;
  int nmale = 0;
  int nfemale = 0;
  int nambig = 0;

  int c=0;
  while(!PED.eof())
    {
      Individual * person = new Individual;
      string phenotype;
      double temp; // use to init sample->phenotypes
      PED >> person->fid
	  >> person->iid
	  >> person->pat
	  >> person->mat
	  >> person->sexcode
	  >> phenotype;

      // neglect the code01 situation
      
      if (person->fid=="")
	{
	  delete person;
	  break;
	}
      
      if (person->fid=="FID")
	cout<<"FID is a reserved ID... please select a different family ID"<<endl;
      
      if (person->sexcode=="1")
	{
	  person->sex = true; // male
	  nmale++;
	}
      else if (person->sexcode=="2")
	{
	  person->sex = false;
	  nfemale ++;
	}
      else
	{
	  ambiguous.push_back(person);
	  nambig++;
	  person->missing = true; //in fact, we have an parameter to judge if we need to set missing to true.
	  
	}
      person->founder = (person->pat=="0" && person->mat=="0") ? true : false;
      if (phenotype == "-9")
	{
	  person->missing = true; //par::missing_phenotype="-9";
	  temp = person->phenotype;
	  person->phenotypes.push_back(temp);
	}
      else
	{
	  if ( !from_string<double>(person->phenotype, phenotype, std::dec))
	    person->missing = true;
	  else
	    {
	      temp = person->phenotype;
	      person->phenotypes.push_back(temp);

	      if ( phenotype != "0" &&
		   phenotype != "1" &&
		   phenotype != "2" )
		{
		  qt = true;
		  bt = false;
		}
	  }
	}
      c++;
      sample.push_back(person);

    } // end of while !PED.eof
  PED.clear();
  PED.clear();

  if (bt)
    for (int i=0; i<sample.size(); i++)
      if ( sample[i]->phenotype == 0)
	sample[i]->missing = true;
  
  int nm=0;
  for (int i=0; i<sample.size();i++)
    if(!sample[i]->missing) 
      nm++;
  
  if(bt)
    {
      //we negelect bt again, and we negelect those printLog part.
    } 

}

string pheno_name = "master_al_od_bend";//I will find another method
int MAX_LINE_LENGTH = 1000;
vector<string> pheno_names;

void get_pheno_names()
{
  pheno_names.push_back(pheno_name);
  pheno_names.push_back("auto_l_se_bend");
}

std::string int2str(int n)
{
  std::ostringstream s2(std::stringstream::out);
  s2 << n;
  return s2.str();
}


bool readPhenoFile(string filename, vector<Individual*> & sample)
{
  bool qt = false;
  bool bt = true;
  ifstream PHE(filename.c_str(), ios::in);

  get_pheno_names();
  
  map<string, Individual*> uid;
  map<string, Individual*>::iterator ii;
  for (int i=0; i<sample.size(); i++)
    {
      uid.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid, sample[i]));
      sample[i]->phenotype = -9;
      sample[i]->phenotypes.resize(0);
      sample[i]->missing = true;
    }
  int pheno_idx;
  vector<int> pheno_idxs;
  
  if ( pheno_name != "" && pheno_names.size() > 0)
    {
      string pfid, piid, ph;
      char cline[MAX_LINE_LENGTH];
      PHE.getline(cline, MAX_LINE_LENGTH, '\n');
      string sline = cline;
      if (sline !="" )
	{
	  string buf;
	  stringstream ss(sline);
	  vector<string> tokens;
	  while ( ss >> buf )
	    tokens.push_back(buf);
	  
	  if ( tokens[0] != "FID")
	    cout<<"First header field must be FID"<<endl;
	  if ( tokens[1] != "IID")
	    cout<<"Second header field must be IID"<<endl;
	  
	  pheno_idx = -1;
	  pheno_idxs.resize(0);
	  for ( int i=2; i<tokens.size(); i++)
	    if ( tokens[i] == pheno_name )
	      pheno_idx = i - 1;
	  
	  for ( int i=0; i<pheno_names.size(); i++)
	    {
	      
	      for ( int j=2; j<tokens.size(); j++)
		if ( tokens[j]==pheno_names[i] )
		  {
		    pheno_idxs.push_back(j-1);
		  }
	    }
	  
	  if (pheno_idx == -1)
	    cout<<"Did not find phenotype" << endl;
	  
	} // end of if sline != ""
    } // end of if pheno_name != ""
  
  int ccount = -1;
  bool all_pheno = false;
  while (!PHE.eof())
    {
      string pfid, piid, ph;
      vector<string> phs;
      char cline[MAX_LINE_LENGTH];
      PHE.getline(cline, MAX_LINE_LENGTH, '\n');
      
      string buf;
      string sline = cline;
      if (sline=="") continue;
      stringstream ss(sline);
      vector<string> tokens;
      while (ss >> buf)
	tokens.push_back(buf);
      
      if ( ccount < 0 )
	ccount = tokens.size();
      else if ( ccount != tokens.size() )
	cout<<("Wrong number of columns in file");
      
      if (tokens.size() < 2+pheno_idx )
	{
	  if (all_pheno)
	    {
	      cout<<"Processed all phenotypes\n";
	      return false;
	    }
	  else
	    cout<<"Proble with "<<filename<<endl;
	}
      pfid = tokens[0];
      piid = tokens[1];
      
      if ( pfid == "FID" || piid == "IID" )
	{
	  //phenotype_name = tokens[1+pheno_idx];
	  continue;
	}
      
      ph = tokens[1+pheno_idx];
      phs.resize(0);
      for (int i=0; i<pheno_idxs.size(); i++)
	{
	  phs.push_back(tokens[1+pheno_idxs[i]]);
	  //cout<<"tokens "<<i<<": "<<tokens[1+pheno_idxs[i]]<<endl;
	}
      //negelect 0/1 coding

      ii = uid.find(pfid+"_"+piid);
      if (ii != uid.end())
	{
	  (ii->second)->missing = true;
	  
	  if (ph != "-9")
	    {
	      (ii->second)->missing = false;
	    }
	  
	  if ( !from_string<double>((ii->second)->phenotype, ph, std::dec) )
	    (ii->second)->missing = false;
	  
	  double temp;
	  for(int i=0; i<pheno_idxs.size(); i++)
	    if( from_string<double>(temp, phs[i], std::dec) )
	      {
		(ii->second)->phenotypes.push_back(temp);
		(ii->second)->missing = false;
	      }
	    else
	      {
		(ii->second)->phenotypes.push_back(-9.0);
		(ii->second)->missing = true;
	      }
	  
	  if (ph != "-9" && ph != "0" &&
	      ph != "1" && ph != "2")
	    {
	      qt = true;
	      bt = false;
	    }
	} // end of if ii != uid.end
    } // end of while !PHE.eof()
  PHE.close();
  
  //negelect binary trait

  int new_nmissing = 0;
  for (int i=0; i<sample.size(); i++)
    if (!sample[i]->missing) new_nmissing++;
  cout <<int2str(new_nmissing)<<" individual with non_missing alternate phenotype\n";
  
  //negelect binary trait

  return true;
}
