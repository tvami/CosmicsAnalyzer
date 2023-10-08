#ifndef CosmicsAnalyzer_EarthAsDMAnalyzer_CommonFunction_h
#define CosmicsAnalyzer_EarthAsDMAnalyzer_CommonFunction_h

#include "TMath.h"
// #include <math.h>

// compute deltaR between two point (eta,phi) (eta,phi)
float deltaR(float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  while (dphi > M_PI)
    dphi -= 2 * M_PI;
  while (dphi <= -M_PI)
    dphi += 2 * M_PI;
  return sqrt(deta * deta + dphi * dphi);
}

#endif