#include "matrix.h"


vector< vector<double> > inverse(vector< vector<double> > & m )
{
  double d;
  int i, j;
  
  if (m.size() == 0) error("Internal error: matrix with no rows (inverse function)");
  if (m.size() != m[0].size() ) error("Internal error: cannot invert non-square matrix");
  int n = m.size();

  // indx is an integer array
  vector<int> indx(n);

  vector<double> col(n);
  vector<vector<double> > y(n);
  for (int i=0; i<n; i++) y[i].resize(n); 
  vector<vector<double> > tm;
  tm = m;
  
  ludcmp(tm,indx,d);
  
  for (j=0; j<n; j++)
    {
      for (i=0; i<n; i++) col[i]=0;
      col[j]=1;
      lubksb(tm,indx,col);
      for (i=0; i<n; i++) y[i][j]=col[i];
    }
  
  return y;
  
}

void ludcmp(vector<vector<double> > &a, vector<int> &indx, double &d)
{
  int i, imax = 0, j, k;
  double big, dum, sum, temp;
  int n = a.size();
  vector<double> vv(n);
  d=1;

  for (i=0; i<n; i++)
    {
      big=0;
      for (j=0; j<n; j++)
	if ((temp=fabs(a[i][j])) > big) big=temp;
      if (big==0) error("singular matrix in ludcmp");
      vv[i]=1/big;
    }
  
  for (j=0; j<n; j++)
    {
      for (i=0; i<j; i++)
	{
	  sum = a[i][j];
	  for (k=0; k<i; k++) sum -= a[i][k] * a[k][j];
	  a[i][j]=sum;
	}
      big=0;
      for (i=j; i<n; i++)
	{
	  sum=a[i][j];
	  for (k=0; k<j; k++)
	    sum -= a[i][k] * a[k][j];
	  a[i][j]=sum;
	  if ((dum=vv[i]*fabs(sum)) >= big)
	    {
	      big = dum;
	      imax = i;
	    }
	}
      if (j != imax)
	{
	  for (k=0; k<n; k++)
	    {
	      dum=a[imax][k];
	      a[imax][k]=a[j][k];
	      a[j][k]=dum;
	    }
	  d = -d;
	  vv[imax]=vv[j];
	}
      indx[j]=imax;
      if (a[j][j] == 0) a[j][j] = 1.0e-20;

      if (j != n-1)
	{
	  dum = 1/(a[j][j]);
	  for (i=j+1; i<n; i++) a[i][j] *= dum;
	}
    }
}

void lubksb(vector<vector<double> > &a, vector<int> &indx, vector<double> &b)
{

  int i, ii=0, ip, j;
  double sum;

  int n = a.size();

  for (i=0; i<n; i++)
    {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii != 0)
	for (j=ii-1; j<i; j++) sum -= a[i][j]*b[j];
      else if (sum != 0.0) ii=i+1;
      b[i]=sum;
    }
  for (i=n-1; i>=0; i--)
    {
      sum=b[i];
      for (j=i+1; j<n; j++) sum -= a[i][j]*b[j];
      b[i]=sum/a[i][i];
    }
}

void sizeMatrix(vector< vector<double> >& m, int r, int c)
{
  m.clear();
  m.resize(r);
  for(int i=0; i<r; i++)
    m[i].resize(c,0);
}

void multMatrix(vector<vector<double> > & a,
		vector<vector<double> > & b,
		vector<vector<double> > & c)
{

  int ar = a.size();
  int br = b.size();
  if (ar == 0 || br == 0)
    error("Internal error: multiplying 0-sized matrices");
  
  int ac = a[0].size();
  int bc = b[0].size();
  if ( ac != br )
    error("Internal error: non-conformable matrices in multMatrix()"); 
  
  int cr = ar;
  int cc = bc;

  c.clear();
  sizeMatrix(c,cr,cc);
  
  for (int i=0; i<ar; i++)
    for (int j=0; j<bc; j++)
      for (int k=0; k<ac; k++)
	c[i][j] += a[i][k] * b[k][j];

}

vector<vector<double> > trans(vector<vector<double> >& mat)
{
  vector<vector<double> > tmat;
  int nc = mat.size();
  int nr = mat[0].size();
  sizeMatrix(tmat, nr, nc);
  for(int i=0; i<nr; i++)
    for(int j=0; j<nc; j++)
      tmat[i][j] = mat[j][i];
  return tmat;
}

double pythag(const double a, const double b)
{
  double absa,absb;
 
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

double SQR(double a)
{
  return a*a;
}

template<class T>
inline const T SIGN(const T &a, const T &b)
{return b>=0?(a>=0?a:-a):(a>=0?-a:a);}

template<class T>
inline const T MAX(const T &a, const T &b)
{return b>a?(b):(a);}

template<class T>
inline const T MIN(const T &a, const T &b)
{return b<a?(b):(a);}

vector< vector<double> > svd_inverse(vector< vector<double> > & u , bool & flag )
{
  
  const double eps = 1e-24; 
  
  if (u.size() == 0) 
    error("Internal problem: matrix with no rows (inverse function)");
  if (u.size() != u[0].size() ) 
    error("Internal problem: Cannot invert non-square matrix");
  int n = u.size();
  
  vector<double> w(n,0);
  
  vector<vector<double> > v(n);
  for (int i=0; i<n; i++) 
    v[i].resize(n,0);

  flag = svdcmp(u,w,v); 
  
  // Look for singular values
  double wmax = 0;
  for (int i=0; i<n; i++)
    wmax = w[i] > wmax ? w[i] : wmax;
  double wmin = wmax * eps;
  for (int i=0; i<n; i++)
    {
      w[i] = w[i] < wmin ? 0 : 1/w[i];
    }
  

  // u w t(v)

  // row U * 1/w
  
  // results matrix
  vector<vector<double> > r(n);
  for (int i=0; i<n; i++)
    {
      r[i].resize(n,0);
      for (int j=0; j<n; j++)
	u[i][j] = u[i][j] * w[j];
    }

  // [nxn].[t(v)] 
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      for (int k=0; k<n; k++)
	r[i][j] += u[i][k] * v[j][k];
    
  return r;
}

bool svdcmp(vector<vector<double> > & a, 
	    vector<double> & w, 
	    vector<vector<double> > &v)
{
  bool flag;
  int i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z;
  double volatile temp;

  int m=a.size();
  if (m==0) error("Internal problem in SVD function (no observations left?)");
  int n=a[0].size();

  vector<double> rv1(n);
  g=scale=anorm=0.0;
  for (i=0;i<n;i++) {
    l=i+2;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a[k][i]);
      if (scale != 0.0) {
	for (k=i;k<m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l-1;j<n;j++) {
	  for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i+1 <= m && i+1 != n) {
      for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
      if (scale != 0.0) {
	for (k=l-1;k<n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l-1];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l-1]=f-g;
	for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
	for (j=l-1;j<m;j++) {
	  for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
	  for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l-1;k<n;k++) a[i][k] *= scale;
      }
    }
    anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g != 0.0) {
	for (j=l;j<n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=MIN(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) a[i][j]=0.0;
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j<n;j++) {
	for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<m;j++) a[j][i] *= g;
    } else for (j=i;j<m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its<30;its++) {
      flag=true;
      for (l=k;l>=0;l--) {
	nm=l-1;
	temp=fabs(rv1[l])+anorm;
	if (temp == anorm) {
	  flag=false;
	  break;
	}
	temp=fabs(w[nm])+anorm;
	if (temp == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<k+1;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  temp = fabs(f)+anorm;
	  if (temp == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == 29) 
	return false; // cannot converge: multi-collinearity?
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g=g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=0;jj<m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  return true;
}




