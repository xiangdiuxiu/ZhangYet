#include "readBinData.h"

//init parameters
Parameter par; 

string extern getOutput(vector<string>& args)
{
  string res = "";
  for(vector<string>::iterator it=args.begin(); it!=args.end(); it++)
    if((*it)=="-o")
      {
	it++;
	res = (*it);
      }
  if(res=="")
    {
      error("No output file");
    }
  return res;
}

void extern initParameters(vector<string>& args)
{
  vector<string>::iterator it = args.begin();
	//par.pheno_names.push_back("auto_r_se_bend");
	//par.pheno_names.push_back("master_al_od_bend");
  bool lockon = false;

  while(it!=args.end())
    {
      string temp = (*it);
      if(temp=="-bfile")
	{
	  it++; lockon = false;
	  par.famfile = (*it)+".fam";
	  par.snpfile = (*it)+".bim";
	  par.bitfile = (*it)+".bed";
	  it++;
	}
      else if(temp=="-phenofile")
	{
	  it++; lockon = false;
	  par.phenofile = (*it);
	  it++;
	}
      else if(temp=="-phenotypes")
	{
	  lockon = true;
	  it++;
	}
      else if(temp=="-o")
	{
	  lockon = false;
	  it++;
	}
      else
	{
	  if(lockon)
	    {
	      par.pheno_names.push_back(*it);
	      it++;
	    }
	  else
	    it++;
	}
      
    }
  cout<<"phenotypes: "<<par.pheno_names.size()<<endl;
}

void extern error(string msg)
{
  cerr <<"\nERROR: " <<msg <<"\n";
  exit(1);
}

void checkFileExists(string f)
{
  ifstream inp;
  inp.open(f.c_str(), ifstream::in);
  if(inp.fail())
    {
      inp.clear(ios::failbit);
      inp.close();
      string msg = "No file [ "+f+" ] exists ";
      error(msg);
    }
  inp.close();
  return;
  
}

template <class T>
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss>>f>>t).fail();
}

string int2str(int n)
{
  ostringstream s2(stringstream::out);
  s2<<n;
  return s2.str();
}

void readFamFile(string filename, vector<Individual*> & sample)
{
  par.qt = false;
  par.bt = true;

  checkFileExists(filename);

  ifstream PED;
  PED.open(filename.c_str());
  PED.clear();

  int c=0;
  while(!PED.eof())
    {
      Individual * person = new Individual;
      string phenotype;

      PED >> person->fid
	  >> person->iid
	  >> person->pat
	  >> person->mat
	  >> person->sexcode
	  >> phenotype;

      if(par.coding01)
	{
	  if(phenotype=="1")
	    phenotype = "2";
	  else if (phenotype=="0")
	    phenotype = "1";
	  else
	    phenotype = "0";
	}

      if(person->fid=="")
	{
	  delete person;
	  break;
	}
      if(person->fid=="FID")
	error("FID is a reserved ID... please select a different family ID.");
      if(person->sexcode=="1")
	{
	  person->sex = true; //male
	}
      else if (person->sexcode=="2")
	{
	  person->sex = false; //make
	}
      else
	{
	  person->missing = true;
	}

      person->founder = (person->pat=="0" && person->mat=="0")?true:false;
      
      if (phenotype==par.missing_phenotype)
	{
	  person->missing = true;
	  person->phenotypes.push_back(-9);
	}
      else
	{
	  if(!from_string<double>(person->phenotype, phenotype, std::dec))
	    person->missing = true;
	  else
	    {
	      double temp = person->phenotype;
	      person->phenotypes.push_back(temp);
	    
	      if(phenotype!="0"&&phenotype!="1"&&phenotype!="2"&&phenotype!=par.missing_phenotype)
		{
		  par.qt = true;
		  par.bt = false;
		}
	    } // end of else from_string
	} // end of else phenotype== ...
      c++;
      sample.push_back(person);
    } // end of while !PED.eof
  PED.clear();
  PED.close();

  if(par.bt)
    for(int i=0; i<sample.size(); i++)
      if(sample[i]->phenotype==0)
	sample[i]->missing = true;
}

bool openBinaryFile(string s, ifstream& BIT)
{
  BIT.open(s.c_str(), ios::in|ios::binary);

  char ch[1];
  BIT.read(ch,1);
  bitset<8> b;
  
  b = ch[0];
  bool bfile_SNP_major = false;
  bool v1_bfile = true;
  
  if((b[2]&&b[3]&&b[5]&&b[6]) &&
     !(b[0]||b[1]||b[4]||b[7]))
    {
      BIT.read(ch,1);
      b = ch[0];
      if((b[0]&&b[1]&&b[3]&&b[4])&&
	 !(b[2]||b[5]||b[6]||b[7]))
	{
	  BIT.read(ch,1);
	  b = ch[0];
	  if(b[0]) bfile_SNP_major = true;
	  else bfile_SNP_major = false;
	  
	} else v1_bfile = false;
    } else v1_bfile = false; //end of if b....
  
  if(!v1_bfile)
    {
      BIT.close();
      BIT.clear();
      BIT.open(s.c_str(), ios::in|ios::binary);
      BIT.read(ch,1);
      b = ch[0];
    }

  if((!v1_bfile)&&(b[1]||b[2]||b[3]||b[4]||b[5]||b[6]||b[7]))
    {
      bfile_SNP_major = false;
      BIT.close();
      BIT.clear();
      BIT.open(s.c_str(), ios::in|ios::binary);
    }
  else if(!v1_bfile)
    {
      if(b[0]) bfile_SNP_major = true;
      else bfile_SNP_major = false;
    }

  return bfile_SNP_major;
}

void extern readBinData(vector<Individual*> & sample, vector<Locus*> & locus, vector<CSNP*>& SNP, vector<string>& args)
{
  initParameters(args);
  checkFileExists(par.bitfile);
  checkFileExists(par.famfile);
  checkFileExists(par.snpfile);

  vector<Locus> ordered;
  //read SNP;
  ifstream MAP(par.snpfile.c_str(), ios::in);
  MAP.clear();

  int c=0;
  while(!MAP.eof())
    {
      Locus* loc = new Locus;
      MAP >> loc->chr
	  >> loc->name
	  >> loc->pos
	  >> loc->bp
	  >> loc->allele1
	  >> loc->allele2;
      
      if(MAP.eof())
	{
	  delete loc;
	  continue;
	} // end of if MAP.eof
      else if (MAP.fail())
	{
	  delete loc;
	  error("Problem reading BIM file, line " + int2str(c+1) + "\n");
	}
      
      loc->freq = c++;
      
      if(loc->name!="")
	{
	  locus.push_back(loc);
	  ordered.push_back(*loc);
	}
      else
	delete loc;
    } // end of while !MAP.eof
  MAP.clear();
  MAP.close();

  if(locus.size()==0)
    error("No SNPs");
  
  //stable_sort(locus.begin(), locus.end(), less<Locus*>());
  stable_sort(ordered.begin(), ordered.end());
  
  

  c = 0;
  for(int i=0; i<locus.size(); i++)
    {
      ordered[i].bp = (int)ordered[i].freq;
      ordered[i].chr = 1;
      ordered[i].freq = c++;
    }
  //resort to get lookup table
  stable_sort(ordered.begin(), ordered.end());

  //Do we want to look at all the data;
  vector<int> include(0);
  int nl_actual = locus.size();
  for(int j=0; j<ordered.size(); j++)
    include.push_back((int)ordered[j].freq); //mind 
  //include mean the snp we need to check
  
  readFamFile(par.famfile, sample);
  
  ifstream BIT;
  bool bfile_SNP_major = openBinaryFile(par.bitfile, BIT);
  
  if(bfile_SNP_major)
    {
      for(int i=0; i<nl_actual; i++)
	{
	  CSNP * newlocus = new CSNP;
	  newlocus->one.resize(sample.size());
	  newlocus->two.resize(sample.size());
	  SNP.push_back(newlocus);
	} // end of for i<nl_actual
    }
  else
    {
      vector<Individual*>::iterator person = sample.begin();
      while(person != sample.end())
	{
	  (*person)->one.resize(nl_actual);
	  (*person)->one.resize(nl_actual);
	  person++;
	}
    } // end of if bfile_SNP_major
  
  if(bfile_SNP_major)
    {
      CSNP * snp;
      int s = 0;
      while(s<locus.size())
	{
	  if(include[s]>-1)
	    snp = SNP[include[s]];
	  else
	    snp = NULL;

	  //inner loop for individuals
	  int indx = 0;
	  int ss = sample.size();
	  while(indx<ss)
	    {
	      bitset<8> b;
	      char ch[1];
	      BIT.read(ch,1);
	      if(!BIT)
		error("Problem with BED file");
	      
	      b = ch[0];
	      int c = 0;
	      while(c<7 && indx<ss)
		{
		  if(snp)
		    {
		      snp->one[indx] = b[c++];
		      snp->two[indx] = b[c++];
		    }
		  else
		    {
		      c += 2;
		    }// end of if else snp
		  ++indx;
		}//end of while c<7 && index<ss
	    }//end of while indx<ss
	  s++;
	}// end of while s < locus.size
      par.SNP_major = true;
    }// end of if bfile_SNP_major
  else
    {
      vector<Individual*>::iterator person = sample.begin();
      while(person != sample.end())
	{
	  int s = 0;
	  while(s<locus.size())
	    {
	      char ch[1];
	      BIT.read(ch,1);
	      if(!BIT)
		error("Problem with the BED file");
	      
	      bitset<8> b;
	      b = ch[0];
      
	      int c = 0;
	      while(c<7 && s<locus.size())
		{
		  if(include[2]>-1)
		    {
		      (*person)->one[include[s]] = b[c++];
		      (*person)->two[include[s]] = b[c++];

		    }
		  else
		    {
		      c += 2;
		    }
		  s++;
		}// end of while c<7 &&
	    }// end of while s<locus.size
	  person++;
	}// end of while person != sample.end
      par.SNP_major = false;
    }// end of else bfile_SNP_major
  char ch[1];
  BIT.read(ch, 1);
  if(BIT)
    error("Problem with the BED file");

}

void extern readPhenoFile(vector<Individual*>& sample, vector<Locus*>& locus, vector<CSNP*>& snp)
{
  par.qt = false;
  par.bt = true;
  
  checkFileExists(par.phenofile);


  ifstream PHE(par.phenofile.c_str(), ios::in);
  
  map<string, Individual*> uid;
  map<string, Individual*>::iterator ii;
  
  for(int i=0; i<sample.size(); i++)
    {
      uid.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid, sample[i]));
      sample[i]->phenotype = -9;
      sample[i]->phenotypes.resize(0);
      sample[i]->missing = true;
    }
  
  vector<int> pheno_idx;
  pheno_idx.resize(0);
  if(par.pheno_names.size()>0)
    {
      string pfid, piid, ph;
      char cline[par.MAX_LINE_LENGTH];
      PHE.getline(cline, par.MAX_LINE_LENGTH, '\n');
      string sline = cline;
      if(sline!="")
	{
	  string buf;
	  stringstream ss(sline);
	  vector<string> tokens;
	  while(ss>>buf)
	    tokens.push_back(buf);

	  if(tokens[0] != "FID")
	    error("First header field must be FID");
	  if(tokens[1] != "IID")
	    error("Second header field must be IID");
	  
	  for(int j=0; j<par.pheno_names.size(); j++)
	    for(int i=2; i<tokens.size(); i++)
	      if(tokens[i] == par.pheno_names[j])
		{
		  pheno_idx.push_back(i-1);
		
		}

	  if(pheno_idx.size()<par.pheno_names.size())
	    error("Did not find phenotypes");

	  
	} // end of if sline != ""
    } // if end of par.pheno_names.size()>0

  int ccount = -1;
  int SAMPLE_SIZE = 0;
  while(!PHE.eof())
    {
      string pfid, piid;
      vector<string> ph;
      ph.resize(0);
      char cline[par.MAX_LINE_LENGTH];
      PHE.getline(cline, par.MAX_LINE_LENGTH, '\n');
      
      string sline = cline;
      if(sline=="") continue;
      
      string buf;
      stringstream ss(sline);
      vector<string> tokens;
      while(ss >> buf)
	tokens.push_back(buf);
      
      if(ccount<0)
	ccount = tokens.size();
      else if(ccount!=tokens.size())
	error("Wrong number of columns in file [ "+ par.phenofile+" ] line :\n"+sline);
      
      //skip tokens.size() < 2+par::mult_pheno

      pfid = tokens[0];
      piid = tokens[1];

      for(int i=0; i<pheno_idx.size(); i++)
	ph.push_back(tokens[1+pheno_idx[i]]);

      if(par.coding01)
	{
	  for(int i=0; i<ph.size(); i++)
	    if(ph[i]=="1")
	      ph[i] = "2";
	    else if(ph[i]=="2")
	      ph[i] = "1";
	    else 
	      ph[i] = "0";
	}

      ii = uid.find(pfid+"_"+piid);
      if(ii!=uid.end())
	{
	  (ii->second)->missing = true;
	  bool flag = true;
	  for(int i=0; i<ph.size(); i++)
	    if(ph[i] != par.missing_phenotype)
	      {
		flag = false;
	      }
	    else { flag = true;} // only all phenotypes != missing_phenotype, flag == false 
	  (ii->second)->missing = flag;
	  
	  //Convert to double, checking for illegal value
	  double temp;
	  (ii->second)->phenotypes.resize(ph.size(), 0);
	  vector<bool> flag2;
	  flag2.resize(0);
	  for(int i=0; i<ph.size(); i++)
	    if(from_string<double>(temp, ph[i], std::dec))
	      {
		(ii->second)->phenotypes[i] = temp;
		flag2.push_back(false);
	      }
	    else
	      {
		flag2.push_back(true);
	      }
	  flag = false;
	  for(int i=0; i<flag2.size(); i++)
	    {
	      if(flag2[i])
		{
		  flag = true;
		  break;
		}
	    }

	  (ii->second)->missing = flag;
	  
	  if(!flag)
	    SAMPLE_SIZE++;

	  for(int i=0; i<ph.size(); i++)
	    if (ph[i]!=par.missing_phenotype &&
		ph[i]!="0" &&
		ph[i]!="1" &&
		ph[i]!="2")
	      {
		par.qt = true;
		par.bt = false;
	      }
	} // end of if ii!=uid.end
      
      
    } // end of while !PHE.eof
  PHE.close();
  //in fact we treat every bt as qt

}


