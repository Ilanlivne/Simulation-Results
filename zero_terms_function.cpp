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