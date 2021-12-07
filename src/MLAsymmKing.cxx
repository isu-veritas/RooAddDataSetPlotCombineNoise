#include "RooFit.h"

#include "Riostream.h"
#include <math.h>
#include <iostream>
#include <fstream>

#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include <TMath.h>
#include "RooAbsReal.h"
#include <TF1.h>

#include "MLAsymmKing.h"

MLAsymmKing::MLAsymmKing(): num(0) { }

//_____________________________________________________________________________
MLAsymmKing::MLAsymmKing(const char *name, const char *title,
                         RooAbsReal& _x, RooAbsReal& _y,
                         RooAbsReal& _xmean, RooAbsReal& _ymean,
                         RooAbsReal& _lambda, RooAbsReal& _sigmax, RooAbsReal& _sigmay) :
  RooAbsPdf(name,title),
  x("x","X-Observable",this,_x),
  y("y","Y-Observable",this,_y),
  xmean("xmean","X-Mean",this,_xmean),
  ymean("ymean","Y-Mean",this,_ymean),
  lambda("lambda","Slope",this,_lambda),
  sigmax("sigmax","Widthx",this,_sigmax),
  sigmay("sigmay","Widthy",this,_sigmay),
  is_symmetric(false)
{
}


//_____________________________________________________________________________                                                             
MLAsymmKing::MLAsymmKing(const MLAsymmKing& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  y("y",this,other.y),
  xmean("xmean",this,other.xmean),
  ymean("ymean",this,other.ymean),
  lambda("lambdax",this,other.lambda),
  sigmax("sigmax",this,other.sigmax),
  sigmay("sigmay",this,other.sigmay),
  is_symmetric(other.is_symmetric)
{
}


void MLAsymmKing::SetTheta(double thetaval) {
  theta=thetaval;
}

//_____________________________________________________________________________
Double_t MLAsymmKing::evaluate() const
{
  // Compute the rotation of the PSF
  //double theta = -TMath::ATan2(ymean, xmean) ;
  //double theta = -TMath::ATan2(0.0, 0.0) ;  
  //std::cout << "Compute the rotation of the PSF xmean " << xmean << " ymean " << ymean << " theta " << theta << std::endl;
  // Compute the transformed coordinates
  double xarg = x - xmean;
  double yarg = y - ymean;
  //double xarg = x - 0.0;
  //double yarg = y - 0.0;
  double x_prime = xarg*TMath::Cos(theta) - yarg*TMath::Sin(theta) ;
  double y_prime = xarg*TMath::Sin(theta) + yarg*TMath::Cos(theta) ;
  //std::cout << "Compute the transformed coordinates" << std::endl; 
  //std::cout << "xarg " << xarg << " yarg " << yarg << std::endl;
  //std::cout << " x_prime " << x_prime << " y_prime " << y_prime << std::endl;
  // If sigma and lambda are calculated via an MLPointSpreadFunctionVar, then
  // by caching the values here it removes the need for several calls to those
  // functions and thus many computations, expecially if they're being both
  // interpolated and not cached.
  double sigX = sigmax ;
  double sigY = sigmay ;
  double lamarg = lambda ;
    
  double norm1 = 0.5/(TMath::Pi()*sigX*sigY);
  double norm2 = (1.0 - (1.0/lamarg));
    
  double ret = norm1*norm2*std::pow(1.0 + ( (0.5/lamarg) * ((x_prime*x_prime)/(sigX*sigX)+(y_prime*y_prime)/(sigY*sigY)) ), -lamarg);
    
  if (ret < 0.0) ret = 0.0 ;
    
  return ret ;
}

//_____________________________________________________________________________
TF2* MLAsymmKing::getFormula()
{
  return nullptr ;
}

//_____________________________________________________________________________
void MLAsymmKing::snapshot(double new_lambda, double new_sigmax, double new_sigmay)
{
  cached_lambda = lambda ;
  cached_sigmax = sigmax ;
  cached_sigmay = sigmay ;
  cached_xmean  = xmean ;
  cached_ymean  = ymean ;
  std::cout << "SNAPSHOT" << std::endl;
  if (new_lambda > 0.0) cached_lambda = new_lambda ;
  if (new_sigmax > 0.0) cached_sigmax = new_sigmax ;
  if (new_sigmay > 0.0) cached_sigmay = new_sigmay ;
    
  double theta = -TMath::ATan2(cached_ymean, cached_xmean) ;
  double cosTheta = TMath::Cos(theta) ;
  cached_cosTheta_xnorm = cosTheta/cached_sigmax ;
  cached_cosTheta_ynorm = cosTheta/cached_sigmay ;
  double sinTheta = TMath::Sin(theta) ;
  cached_sinTheta_xnorm = sinTheta/cached_sigmax ;
  cached_sinTheta_ynorm = sinTheta/cached_sigmay ;
    
  return ;
}

//_____________________________________________________________________________
double MLAsymmKing::getFormula(double x, double y,
                               double newlambda,
                               double newsigmax, double newsigmay)
{
  // Return a TF2 representing a snapshot of this object
  std::cout << "Return a TF2 representing a snapshot of this object " << std::endl;
  // Precompute some things
  if (newlambda > 0.0) cached_lambda = newlambda ;
  if (newsigmax > 0.0) cached_sigmax = newsigmax ;
  if (newsigmay > 0.0) cached_sigmay = newsigmay ;
    
  // Compute the transformed coordinates
  double xarg = x - cached_xmean;
  double yarg = y - cached_ymean;
  double x_prime = xarg*cached_cosTheta_xnorm - yarg*cached_sinTheta_xnorm ;
  double y_prime = xarg*cached_sinTheta_ynorm + yarg*cached_cosTheta_ynorm ;
    
  // If sigma and lambda are calculated via an MLPointSpreadFunctionVar, then
  // by caching the values here it removes the need for several calls to those
  // functions and thus many computations, expecially if they're being both
  // interpolated and not cached.
  double norm1 = 0.5/(TMath::Pi()*cached_sigmax*cached_sigmay);
  double norm2 = (1.0 - (1.0/cached_lambda));

  // Evaluate and return
  return norm1*norm2*std::pow(1.0 + ( (0.5/cached_lambda) * ((x_prime*x_prime)+(y_prime*y_prime)) ), -cached_lambda);
}
/*
//_____________________________________________________________________________                                                           
double MLAsymmKing::getFormulaOneD(double x, double newlambda, double newsigma)
{
  // Return a TF2 representing a snapshot of this object                                                                                  
  // Precompute some things                                                                                                               
  if (newlambda > 0.0) cached_lambda = newlambda ;
  if (newsigma > 0.0) cached_sigmax = newsigma ;

  // Compute the transformed coordinates                                                                                                  
  double xarg = x - cached_xmean;
  //double x_prime = xarg*cached_cosTheta_xnorm - yarg*cached_sinTheta_xnorm ;
  // If sigma and lambda are calculated via an MLPointSpreadFunctionVar, then                                                             
  // by caching the values here it removes the need for several calls to those                                                            
  // functions and thus many computations, expecially if they're being both                                                               
  // interpolated and not cached.                                                                                                         
  double norm1 = 0.5/(TMath::Pi()*cached_sigmax*cached_sigmax);
  double norm2 = (1.0 - (1.0/cached_lambda));

  // Evaluate and return                                                                                                                  
  return norm1*norm2*std::pow(1.0 + ( (0.5/cached_lambda) * (xarg*xarg) ), -cached_lambda);
}
*/

//_____________________________________________________________________________
double MLAsymmKing::getFormulaUnnormalized(double x, double y,
                                           double newlambda, double newsigmax,
                                           double newsigmay)
{
  // Return a TF2 representing a snapshot of this object
  // Note that you need to have called MLAsymmKing::Snapshot() before this method.
  // Precompute some things
  if (newlambda > 1.0) cached_lambda = newlambda ;
  if (newsigmax > 0.0) cached_sigmax = newsigmax ;
  if (newsigmay > 0.0) cached_sigmay = newsigmay ;
    
  std::cout << "getFormulaUnnormalized" << std::endl;

  // Compute the transformed coordinates
  double xarg = x - cached_xmean;
  double yarg = y - cached_ymean;
  double x_prime = xarg*cached_cosTheta_xnorm - yarg*cached_sinTheta_xnorm ;
  double y_prime = xarg*cached_sinTheta_ynorm + yarg*cached_cosTheta_ynorm ;

    
  // Evaluate and return
  return TMath::Power(1.0 + ( (0.5/cached_lambda) * ((x_prime*x_prime)+(y_prime*y_prime)) ), -cached_lambda) ;
}
