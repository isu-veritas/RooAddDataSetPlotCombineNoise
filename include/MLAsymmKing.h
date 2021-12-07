//
#include <cmath>

// ROOT HEADERS
#include <TF2.h>
#include <TF1.h>

// ROOFIT HEADERS
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"


class MLAsymmKing : public RooAbsPdf {

 public:
  int num;///
  MLAsymmKing();
  // Symmetric King function                                                                                                               
  MLAsymmKing(const char *name, const char *title,
	      RooAbsReal& _x, RooAbsReal& _y,
	      RooAbsReal& _xmean, RooAbsReal& _ymean,
	      RooAbsReal& _lambda, RooAbsReal& _sigmax, RooAbsReal& _sigmay);

  // Copy constructor                                                                                                                      
  MLAsymmKing(const MLAsymmKing& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new MLAsymmKing(*this,newname); }
  inline virtual ~MLAsymmKing() {} ;

  // Methods for defining the x,y position of the PSF
  void setXMean(double newXMean) { xmean = newXMean ; }
  void setYMean(double newYMean) { ymean = newYMean ; }
  double getXMean() { return xmean ; }
  double getYMean() { return ymean ; }
  Double_t evaluate() const;
  // Method for getting a TF2 representing a snapshot of this object with the current                                                       
  // values of sigma and lambda                                                                                                             
  TF2* getFormula() ;

  void SetTheta(double thetaval) ;
  // Return the properly normalized version of this function
  double getFormula(double x, double y,
		    double newlambda=-1.0, double newsigmax=-1.0,
		    double newsigmay=-1.0) ;
  //double getFormulaOneD(double x, double newlambda=-1.0, double newsigma=-1.0);
  // Return the value stripped of all of its normalization constants to simplify
  // the number of computations needed in the extended source computation
  double getFormulaUnnormalized(double x, double y,
				double newlambda=-1.0, double newsigmax=-1.0,
				double newsigmay=-1.0) ;


  void snapshot(double new_lambda=-1.0, double new_sigmax=-1.0, double new_sigmay=-1.0) ;

  double getCurLambda() {return lambda ;}
  double getCurSigmaX() {return sigmax ;}
  double getCurSigmaY() {return sigmay ;}
  double getCurMeanX()  {return xmean ;}
  double getCurMeanY()  {return ymean ;}    
  bool IsSymmetric() {return is_symmetric ;}


 protected:
    
  // Proxy variables for referencing the values of the passed RooRealVars
  RooRealProxy x ;
  RooRealProxy y ;
  RooRealProxy xmean ;
  RooRealProxy ymean ;
  RooRealProxy lambda ;
  RooRealProxy sigmax, sigmay ;
    
  double cached_sigmax, cached_sigmay ;
  double cached_lambda ;
  double cached_xmean ;
  double cached_ymean ;
  double cached_cosTheta_xnorm=1.0 ;
  double cached_cosTheta_ynorm=1.0 ;
  double cached_sinTheta_xnorm=1.0 ;
  double cached_sinTheta_ynorm=1.0 ;
    
  bool is_symmetric = false ;

  double theta;
  //Double_t evaluate() const;

};
