#ifndef __MAIN_H__
#define __MAIN_H__

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <bitset>

#include "Individual.h"
#include "Locus.h"
#include "model.h"
#include "matrix.h"
#include "readBinData.h"

using namespace std;

//function table
string extern getOutput(vector<string>&);
void extern readBinData(vector<Individual*>&, vector<Locus*>&, vector<CSNP*>&, vector<string>&);
void extern readPhenoFile(vector<Individual*>&, vector<Locus*>&, vector<CSNP*>&);
vector<double> extern glmAssoc(vector<Individual*>&, vector<Locus*>&, vector<CSNP*>&, string);



#endif



