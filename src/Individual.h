#ifndef __INDIVIDUAL_H__
#define __INDIVIDUAL_H__

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;
typedef vector<double> vector_t;
class WMLocus;
class GVariant;
class Family;
class Individual;

class Individual{
    public:
        string fid;
        string iid;
        string pat;
        string mat;
        
        Individual *pp;
        Individual *pm;
        
        int ip;
        int im;

        Individual *pperson;
        vector<Individual*> kids;
        vector<int> ikids;
        
        bool sex;
        string sexcode;
        double phenotype;
	vector<double> phenotypes;
        bool aff;
        double covar;
        bool bcovar;
        
        vector_t clist; //multiple covariates
        vector<bool> clistMissing;
        vector_t plist; //multiple phenotypes
        vector<bool> plistMissing;

        bool missing;
        bool missing2;
        bool flag;

        int sol;
        bool founder;
        Family *family;

        vector<bool> one;
        vector<bool> two;
        
        vector<GVariant*> gvar;

        //WMLocus wmlocus;
        double T;
        double B;
        double W;
	Individual()
	  {
	    fid=iid=pat=mat="";
	    ip=im=1;
	    sex=false;
	    phenotype = -9;
	    phenotypes.resize(0);
	    sexcode="";
	    aff=false;
	    covar=-9;
	    bcovar=false;
	    clist.resize(0);
	    clistMissing.resize(0);
	    plist.resize(0);
	    plistMissing.resize(0);
	    missing=false;
	    missing2=false;
	    flag=true;
	    one.resize(0);
	    two.resize(0);
	    sol=0;
	    founder=true;
	    pp=pm=NULL;
	    family=NULL;
	    kids.resize(0);
	    pperson=this;
	    T=W=B=0;
	    gvar.resize(0);

	  }
};

class WMLocus{
    public:
        int chr;
        string name;
        
        vector<string> allele;
        vector<double> weight;

        WMLocus() { chr=0; name="";  }
        void reset()
        {
            allele.clear();
            weight.clear();
        }
};


class Family
{
    public:
        bool include;
        bool parents;
        bool sibship;
        bool discordant_parents;
        bool singleton;
        bool TDT;

        Individual *pat;
        Individual *mat;
        vector<Individual*> kid;
        
        double B;
        Family()
        {
            include = false;
            parents = false;
            discordant_parents = false;
            singleton = false;
            sibship = false;
            TDT = false;
            pat = mat = NULL;
            kid.clear();
        }    
        void copy(const Family &rhs);
        Family &operator = (const Family &rhs)
	  {
	    copy(rhs);
	    return *this;
	  }
        
};


class GVariant {
    public:
        bool missing;
        int allele1;
        int allele2;
        float dosage1;
        float dosage2;
        GVariant()
        {
            missing = true;
            allele1 = allele2 = -1;
            dosage1 = dosage2 = 0;
        }
};

#endif
