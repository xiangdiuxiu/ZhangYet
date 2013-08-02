#ifndef __STATS_H__
#define __STATS_H__

#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include "readBinData.h"
using namespace std;

void lubksb(vector<vector<double> > &a, vector<int> &indx, vector<double> &b);
void ludcmp(vector<vector<double> > &a, vector<int> &indx, double &d);
vector< vector<double> > inverse(vector< vector<double> > & m );

void extern error(string);

void sizeMatrix(vector< vector<double> >&, int, int);
void multMatrix(vector<vector<double> > &,
		vector<vector<double> > &,
		vector<vector<double> > & );

double SQR(double);
double pythag(const double, const double);
bool svdcmp(vector<vector<double> > &, 
	    vector<double> &, 
	    vector<vector<double> > &);
vector< vector<double> > svd_inverse(vector< vector<double> > & , bool & );
vector< vector<double> > trans(vector< vector<double> >&);
#endif
