//A class that gets a RooDataSet and TTree from a root file                                                                             
// C++ HEADERS                                                                                                                          
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <getopt.h>
#include <regex>

// ROOT INCLUDES                                                                                                                        
#include "TFile.h"
#include "TObject.h"
#include "TString.h"
#include "TTree.h"
#include "TMath.h"
#include "TNamed.h"
#include <THn.h>
#include "TF1.h"
#include "TF2.h"

// ROOFIT INCLUDES                                                                                                                      
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooRandom.h"
#include "RooDataHist.h"
#include "RooThresholdCategory.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include <RooFitResult.h>
#include "RooPlot.h"
#include <RooDataHist.h>

#include "MLAsymmKing.h"
//#include "MLAsymmKingOneD.h" 
//#include "MLBinning.h"

class AddRooDataSet {

 public:
  // Default constructor                                                                                                             
  AddRooDataSet();
  // Constructor from a single datafile
  AddRooDataSet(std::string ss1, std::string ss2, std::ofstream& myfile, std::ofstream& intfile, int numTels);
  virtual ~AddRooDataSet();


  // Some constants to help us out                                                                                                                        
  const double r2d = TMath::RadToDeg();
  const double d2r = TMath::DegToRad();

  //MLBinning* AxisBinning;
  // Some constants to help us out                                                                                                   
  //==================================
  // Public Data Members
  //==================================

  // TODO : change these variables to RooRealProxy
  
  RooRealVar* xCoord ;        // x-position
  RooRealVar* yCoord ;        // y-position
  RooRealVar* msw ;           // mean scaled width
  RooRealVar* cosZ ;          // zenith angle (degrees)
  RooRealVar* azimuth ;       // azimuth angle (degrees)
  RooRealVar* Erec ;          // reconstructed energy (TeV)
  RooRealVar* noise ;         // noise level
  RooRealVar* Tels ;          // telescope multiplicity
  RooRealVar* Etrue ;         // true energy (TeV, only valid for simulations)
  
  //bool fillFromFile(std::string s5, std::string rooname);
  bool addDataSetFromFile(std::string ss1, std::string ss2, std::ofstream& myfile, std::ofstream& intfile, int numTels);
  bool cleanupDataSet();
  Double_t psfeq(Double_t *x, Double_t *par);
  void PrintInfo(std::string rootname1,std::string rootname2);
  bool setEtrueBinning(std::vector<double> binbounds);

  // Boolean for whether this file is a simulation or not                                                                               
  bool isSimData() {return _isSim ;}

  // Return the epoch, season, or run number of this dataset                                                                            
  int getEpoch() {return _epoch ;}
  int getSeason() {return _ATM ;}
  //virtual long long getRunNumber() {return _runNum ;}

  // Return the tracking coordinates as a 2 element vector in radians                                                                   
  std::vector<double> getTrackingCoordinates() {return _trackPos ;}
  // Return the tracking RA of this dataset in degrees                                                                                  
  double getTrackingRA_Deg() {return _trackPos[0];}
  //double getTrackingRA_Deg() {return _trackPos[0]*TMath::RadToDeg() ;}
  //double getTrackingRA_Rad() {return _trackPos[0] ;}
  // Return the tracking Dec of this dataset                                                                                            
  double getTrackingDec_Deg() {return _trackPos[1];} 
  //double getTrackingDec_Deg() {return _trackPos[1]*TMath::RadToDeg() ;}
  //double getTrackingDec_Rad() {return _trackPos[1] ;}

  // Return the source coordinates as a 2 element vector in radians                                                                     
  std::vector<double> getSourceCoordinates() {return _sourcePos ;}
  // Return the source RA from this file                                                                                                
  double getSourceRA_Deg() {return _sourcePos[0];} 
  //double getSourceRA_Deg() {return _sourcePos[0]*TMath::RadToDeg() ;}
  //double getSourceRA_Rad() {return _sourcePos[0] ;}
  // Return the source Dec from this file                                                                                               
  double getSourceDec_Deg() {return _sourcePos[1];} 
  //double getSourceDec_Deg() {return _sourcePos[1]*TMath::RadToDeg() ;}
  //double getSourceDec_Rad() {return _sourcePos[1] ;}

  double getSourceOffset() {return _sourceOffset ;}
  std::string getSourceName() {return _SourceName ;}

  std::vector<double> ebinsbounds;
  int tel_multiplicity ;          // Holds a bitset representing which telescope multiplicities                                                           
  double _sourceOffset; // Angular distance from tracking position to source (degrees)                                                                    
  int _epoch;
  int _ATM;
  double _livetime;   //livetime in seconds                                                                                                               
  double _avgZenith;  //not cosZenith                                                                                                                    
  double _avgAzimuth;
  double _avgNoise;
  double _offset;
  double _sigmax;
  double _sigmay;
  double _lambda;



 protected:

  RooDataSet* _dataset[2]; //pointer to stored dataset
  RooDataSet* _totaldataset;
  TTree* _datatree[2]; //pointer to stored ttree

  //vectors
  std::vector<double> *input_ObsInfoVect[2];
  //std::vector<double> *input_ObsInfoVect;
  std::vector<std::string> *input_NameVect[2];
  //std::vector<std::string> *input_NameVect; 
  std::vector<bool> *input_SimDataVect[2];
  //std::vector<bool> *input_SimDataVect;


  //VECTOR VARIABLES//
  std::string _SourceName ;       // Source name for this run 
  bool _isSim; // true if simulation header exists 
  std::vector<double> _sourcePos ;// Holds the RA/Dec of the source position                                                             
  std::vector<double> _trackPos ; // Holds the RA/Dec of the tracking position 
  /*
  int tel_multiplicity ;          // Holds a bitset representing which telescope multiplicities
  double _sourceOffset; // Angular distance from tracking position to source (degrees) 
  int _epoch;
  int _ATM;
  double _livetime;   //livetime in seconds                                                                                              
  double _avgZenith;  //not cosZenith                                                                                                    
  double _avgAzimuth;
  double _avgNoise;
  double _offset;
  double sigma;
  double ll;
  */
  //SetName("name") ;
  //_isSim = false ;

  //_trackPos[0] = 0.0 ;
  //_trackPos[1] = 0.0 ;
  //_sourcePos[0] = 0.0 ;
  //_sourcePos[1] = 0.0 ;


};
