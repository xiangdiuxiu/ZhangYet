#ifndef __LOCUS_H__
#define __LOCUS_H__

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <functional>
#include <new>

using namespace std;

class Locus {
    public:
        int chr;
        string name;
        string allele1;
        string allele2;

        double freq;    //of allele1
        double pos; //cM map positions
        int bp; //base-pair position
        int nm; // number of non-missing alleles

	double stat; //statisics

    public:
        Locus() { chr=0; name=""; allele1=""; allele2=""; freq=0; pos=0; bp=0; nm=0; }
        Locus(const Locus& h1) { copy(h1); }
        Locus & operator= (const Locus& h1) { copy(h1); return *this;}
        void copy(const Locus &h1);
        bool operator< (const Locus & p2) const
        {
            return (chr < p2.chr || (chr == p2.chr && bp < p2.bp));
        }
        bool operator== (const Locus & p2) const
        {
            return (name == p2.name);
        }   
};

class CSNP
{
 public:
  vector<bool> one;
  vector<bool> two;
  CSNP()
    {
      one.resize(0);
      two.resize(0);
    }
};

#endif
