#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double func(DoubleVector vec1, DoubleVector vec2, DoubleVector vec3) {
  int n = vec1.size();
  long double sum=0;
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
    {
      for(int k=0;k<n;k++)
      {
        if(i!=j && j!=k && k!=i) sum += vec1[i]*vec2[j]*vec3[k];
      }
    }
  }
  sum /= n*(n-1)*(n-2)*1.0;
  return sum;
}


// [[Rcpp::export]]
NumericVector double_dist_sum(NumericMatrix  X)
{
  NumericMatrix Xm= clone(X);
  int n,p;
  n=Xm.nrow();
  p=Xm.ncol();
  NumericVector Z(n);
	
	for (int i = 0; i < n; i++)
	{
	  for (int j1 = 0; j1 < p; j1++)
	  {
	    for (int j2 = 0; j2 < p; j2++)
	    {
	      if (j1 < j2)
	      {
	        Z[i] += Xm(i,j1)*Xm(i,j2);
	        
	      }
	    }
	    
	  }
	  
	  
	}
	

	return Z;
}




// [[Rcpp::export]]
double func22(NumericMatrix X, DoubleVector Y) {
  int n = X.nrow();
  int p = X.ncol();
   long double sum=0;
  for(int j1=0; j1<p; j1++)
  {
    for(int i1=0; i1<n; i1++)
    {
	 for (int i2=0; i2<n; i2++ )
		{
		 for (int j2=0; j2<p; j2++)
		    {
			 for (int j3=0; j3<p; j3++)
			 {
					 if(i1<i2 && j2<j3) sum += X(i1,j1)*X(i2,j1)*X(i2,j2)*X(i2,j3)*Y[i1]*Y[i2];
			 }
			}	
		}
	  }
  }
  sum *= 1.0;
  return sum;
}




