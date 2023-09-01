#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double pdDepthDists(NumericVector new_dists, NumericMatrix dists, int N1){
  double sum = 0.0;
  for(int i = 0; i < N1; i++){
    for(int j = 0; j < N1; j++){
      sum += ((new_dists[i] < dists(i,j)) || (new_dists[j] < dists(i,j)));
    }
  }
  return sum;
}

// [[Rcpp::export]]
double lcdVariantDepthDists(NumericVector new_dists, NumericMatrix dists, int N1){
  Rcpp::NumericMatrix weights(N1, N1);
  for(int i = 0; i < N1; i++){
    for(int j = 0; j < N1; j++){
      for(int k = 0; k < N1; k++){
        weights(i,j) += ((dists(i,k) <= dists(i,j)) || (dists(j,k) <= dists(i,j)));
      }
    }
  }
  
  double sum = 0.0;
  for(int i = 0; i < N1; i++){
    for(int j = 0; j < N1; j++){
      sum += ((new_dists[i] < dists(i,j)) || (new_dists[j] < dists(i,j))) * (1.0 - 1.0*weights(i,j)/N1);
    }
  }
  return sum;
}

// [[Rcpp::export]]
double lcdDepthDists(NumericVector new_dists, NumericMatrix dists, int N1){
  Rcpp::NumericVector weights(N1, 1.0);
  for(int i = 0; i < N1; i++){
    for(int j = 0; j < N1; j++){
      weights[i] += ((new_dists[j] <= new_dists[i]) || (dists(i,j) <= new_dists[i]));
    }
  }
  
  double sum = 0.0;
  for(int i = 0; i < N1; i++){
    for(int j = 0; j < N1; j++){
      sum += ((new_dists[j] < dists(i,j)) && (new_dists[j] < new_dists[i])) * 1.0/weights[i];
    }
  }
  return sum;
}

