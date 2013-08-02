#include "Locus.h"

void Locus::copy(const Locus &h1)
{
    chr = h1.chr;
    name = h1.name;
    allele1 = h1.allele1;
    allele2 = h1.allele2;
    freq = h1.freq;
    pos = h1.pos;
    bp = h1.bp;
    nm = h1.nm;
}


