//A class that gets a RooDataSet and TTree from a root file
// C++ HEADERS
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <memory>
#include <getopt.h>
#include <regex>

//EIGEN
#include <Eigen/Dense>
#include <Eigen/Geometry>

// ROOT INCLUDES
#include "TFile.h"
#include "TObject.h"
#include "TString.h"
#include "TTree.h"
#include "TMath.h"
#include "TNamed.h"
#include "TStyle.h"
#include "TPad.h"
#include "TFrame.h"
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPave.h>
#include <TPavesText.h>
#include <TLine.h>
#include <RooHist.h>
#include "TPrincipal.h"
#include "TNtuple.h"
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TMatrixDSymEigen.h>

// ROOFIT INCLUDES
#include "RooAbsArg.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRandom.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include <RooFitResult.h>
#include "RooPlot.h"
#include <RooDataHist.h>

#include "AddRooDataSet.h"
//#include "MLAsymmKing.h"


// Default constructor, this doesnt do much of anything
AddRooDataSet::AddRooDataSet(void){
  std::cout << "Object for reading RooDataSet and TTree created " << std::endl;
  // initialize the source and tracking positions to some generic location                                                              
  _trackPos  = std::vector<double>(2,0.0) ;
  _sourcePos = std::vector<double>(2,0.0) ;

}

// Constructor with a real data file 
AddRooDataSet::AddRooDataSet(std::string ss1, std::string ss2, std::ofstream& myfile, std::ofstream& intfile, int numTels){
 //bool err;
  _trackPos  = std::vector<double>(2,0.0) ;
  _sourcePos = std::vector<double>(2,0.0) ;
  addDataSetFromFile(ss1, ss2, myfile, intfile, numTels);
  return ;
}
// Destructor
AddRooDataSet::~AddRooDataSet(void)
{
  // What do I need to delete so that memory management
  // doesnt become a problem?
  std::cout << "Object is being deleted" << std::endl;
}

bool AddRooDataSet::addDataSetFromFile(std::string ss1, std::string ss2, std::ofstream& myfile, std::ofstream& intfile, int numTels) {

  myfile.open("kingparams.txt", std::ios_base::app);
  intfile.open("integralfrac.txt", std::ios_base::app);

  char out_consts[150],out_integral[150];
  const int NoAzbins = 4;
  //const int NoAzbins = 2;
  std::string Abins[NoAzbins] = {"0-90","90-180","180-270","270-360"};
  //std::string Abins[NoAzbins] = {"0-180","180-360"};
  std::string MSWbin[2] = {"0p8-1p1","1p1-1p3"};
  float MSWLo[2] = {0.8,1.1};
  float MSWHi[2] = {1.1,1.3};
  int AzRng1LowerBin[NoAzbins] = {0,90,180,270};
  int AzRng1UpperBin[NoAzbins] = {90,180,270,359};
  int AzRngBin[NoAzbins] = {45,135,225,315};
  //int AzRng1LowerBin[NoAzbins] = {0,180};                                                                                                
  //int AzRng1UpperBin[NoAzbins] = {180,359};                                                                                              
  //int AzRngBin[NoAzbins] = {90,270};
  std::string Offsetstring[9] = {"0p0","0p25","0p5","0p75","1p0","1p25","1p5","1p75","2p0"};
  double OffsetBin[9] = {0.00,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00};
  bool FitStatus; //bool FitStatusBinned[NofEbins][NoAzbins];


  //open root file
  TFile* rfinRoo[2];

  rfinRoo[0] = new TFile(ss1.c_str(), "read");                                                                                        
  _dataset[0] = (RooDataSet*)rfinRoo[0]->Get("dataout");
  _datatree[0] = (TTree*)rfinRoo[0]->Get("datatree");

  rfinRoo[1] = new TFile(ss2.c_str(), "read");
  _dataset[1] = (RooDataSet*)rfinRoo[1]->Get("dataout");
  _datatree[1] = (TTree*)rfinRoo[1]->Get("datatree");

  rfinRoo[0]->GetObject("NameVect", input_NameVect[0]);                                                                                      
  rfinRoo[0]->GetObject("ObsInfoVect", input_ObsInfoVect[0]);
  rfinRoo[0]->GetObject("Sim_Data", input_SimDataVect[0]);
  rfinRoo[1]->GetObject("NameVect", input_NameVect[1]);
  rfinRoo[1]->GetObject("ObsInfoVect", input_ObsInfoVect[1]);
  rfinRoo[1]->GetObject("Sim_Data", input_SimDataVect[1]);

  //DIAGNOSTICS
  if ((_datatree[0]->GetEntries() == 0)||(_datatree[1]->GetEntries() == 0)) {
    std::cout << "Number of events recorded in either TTree: " 
	      << _datatree[0]->GetEntries() << " or " << _datatree[1]->GetEntries() << std::endl;
    delete _datatree[0]; delete _datatree[1];
    delete _dataset[0]; delete _dataset[1];
    rfinRoo[0]->Close(); rfinRoo[1]->Close();
    delete rfinRoo[0]; delete rfinRoo[1];
    return 0;
  } else {

    std::cout << "Number of events recorded in first TTree: " << _datatree[0]->GetEntries() << std::endl;
    std::cout << "Number of events recorded in second TTree: " << _datatree[1]->GetEntries() << std::endl;

    if( (abs(input_ObsInfoVect[0]->at(7)-input_ObsInfoVect[1]->at(7))>1.5) ){
      std::cout << "Both root files have large differences in zenith!" << std::endl;
      return 0;
    }

    else if( (abs(input_ObsInfoVect[0]->at(10)-input_ObsInfoVect[1]->at(10))>0.1) ){
      std::cout << "Both root files have large differences in offset!" << std::endl;
      return 0;
    }

    else if( (abs(input_ObsInfoVect[0]->at(9)-input_ObsInfoVect[1]->at(9))>1.6) ){
      std::cout << "Both root files have large differences in noise!" << std::endl;
      return 0;
    }

    else{
      std::cout << "Source Name: " << input_NameVect[0]->at(0) << std::endl;
      std::cout << "Epoch: " << input_ObsInfoVect[0]->at(0) << std::endl;
    }

    _epoch = input_ObsInfoVect[0]->at(0);
    _livetime = input_ObsInfoVect[0]->at(1);
    _ATM = input_ObsInfoVect[0]->at(2);
    _trackPos = {(input_ObsInfoVect[0]->at(3)), (input_ObsInfoVect[0]->at(4))};
    _sourcePos = {(input_ObsInfoVect[0]->at(5)), (input_ObsInfoVect[0]->at(6))};
    _avgZenith = input_ObsInfoVect[0]->at(7);
    if(_avgZenith<2.0) {_avgZenith=0.0;}
    _avgAzimuth = input_ObsInfoVect[0]->at(8);
    _avgNoise = (input_ObsInfoVect[0]->at(9)+input_ObsInfoVect[1]->at(9))/2;
    std::cout << "Average Noise File1: " << input_ObsInfoVect[0]->at(9) << "Average Noise File2: " << input_ObsInfoVect[1]->at(9) << std::endl;
    std::cout << "Average Noise: " << _avgNoise << std::endl;
    //_avgNoise = input_ObsInfoVect->at(9);
    _offset = input_ObsInfoVect[0]->at(10);
    _isSim = input_SimDataVect[0]->at(0);
    _SourceName = input_NameVect[0]->at(0);

    std::cout << "Number of elements in first dataset " << _dataset[0]->numEntries() << std::endl;
    std::cout << "Number of elements in second dataset " << _dataset[1]->numEntries() << std::endl;
    _dataset[0]->append(*_dataset[1]);
    std::cout << "Number of elements in first dataset(after merging) " << _dataset[0]->numEntries() << std::endl;
    //_totaldataset = _dataset[0];
    //delete _dataset[0]; delete _dataset[1];
    //_totaldataset->numEntries()
    //std::cout << "Number of elements in total dataset " << _totaldataset->numEntries() << std::endl;

    const RooArgSet* row = _dataset[0]->get();
    xCoord = (RooRealVar*)row->find("xpos");
    yCoord = (RooRealVar*)row->find("ypos");
    msw = (RooRealVar*)row->find("msw");
    cosZ = (RooRealVar*)row->find("cosZ");       // zenith angle (degrees)                                                                 
    azimuth = (RooRealVar*)row->find("azimuth");       // azimuth angle (degrees)                                                          
    //Erec =;          // reconstructed energy (TeV)                                                                                       
    noise = (RooRealVar*)row->find("noise");         // noise level                                                                        
    Tels = (RooRealVar*)row->find("Tels");          // telescope multiplicity                                                              
    Etrue = (RooRealVar*)row->find("Etrue");         // true energy (TeV, only valid for simulations)  

    std::cout<< "Mean Noise: " << _dataset[0]->mean(*noise) << " " << std::endl;
    std::cout<< "Mean Zenith: " << std::acos(_dataset[0]->mean(*cosZ))*TMath::RadToDeg() << " " << std::endl;
    std::cout<< "Mean Xpos: " << _dataset[0]->mean(*xCoord) << " " << std::endl;
    std::cout<< "Mean Ypos: " << _dataset[0]->mean(*yCoord) << " " << std::endl;
    std::cout<< "Offset: " << _offset << " " << std::endl;
    std::cout<< "Mean Etrue: " << _dataset[0]->mean(*Etrue) << " " << std::endl;
    std::cout<< "Get Entries " << _dataset[0]->numEntries() << " " << std::endl;
    
    // Get the x,y, MSW, azimuth, and noise variables for the dataset                                                                      
    RooRealVar* tmpX = xCoord ;
    RooRealVar* tmpY = yCoord ;
    // The following parameters will need to be cut on                                                                                     
    RooRealVar* tmpTels  = Tels ;
    RooRealVar* tmpAz    = azimuth ;
    RooRealVar* tmpEtrue = Etrue ;
    RooRealVar* tmpMSW = msw;
    // Set the mean of the source position to the mean x & y of the dataset                                                                
    RooRealVar xmean("xmean","xmean", _dataset[0]->mean(*xCoord),
		     tmpX->getMin(), tmpX->getMax()) ;
    RooRealVar ymean("ymean","ymean", _dataset[0]->mean(*yCoord),
		     tmpY->getMin(), tmpY->getMax()) ;

    // Create the model that will be used to fit to the data                                                                               
    // The following parameters will be passed to the RooKing2D function                                                                                       
    bool doAsymmPSF=1;
    const double sigmaStart = 0.04 ;
    const double lambdaStart = 2.5 ;
    RooRealVar* sigmax = new RooRealVar("sigmax","sigma X (degrees)", sigmaStart, 0.001, 0.3) ;
    sigmax->setAttribute("StoreError") ;
    //RooRealVar* sigmay = sigmax ;   // Default for the symetric case                                                                     
    RooRealVar* sigmay;
    if (doAsymmPSF) {
      // If an asymmetric PSF is requested, initialize sigmay to it's own value                                                         
      std::cout << "Running asymmetric case" << std::endl; 
      sigmay = new RooRealVar("sigmay","sigma Y (degrees)", sigmaStart, 0.001, 0.3) ;
      sigmay->setAttribute("StoreError") ;
    }
    else { std::cout << "Running symmetric case" <<std::endl; sigmay = sigmax ;}

    RooRealVar lambda("lambda","lambda", lambdaStart, 1.0, 8.0) ;
    std::cout << "Lambda name: " << lambda.GetName() << std::endl;

    // Check if we're fixing the lambda value                                                                                       
    const double fixed_lambda=2.5;
    if (fixed_lambda >= 1.0) {
      // Set the range so that lambda is actually a valid value (just in case)                                                             
      // otherwise if 'fixed_lambda' is outside the range it will be set to                                                                
      // the closest boundary value.                                                                                                       
      lambda.setRange(1.0, 2.0*fixed_lambda) ;
      // Set lambda and fix its value so that it wont change                                                                                                      
      lambda.setVal(fixed_lambda) ;
      lambda.setConstant(kTRUE) ;
    }
    
    //asymmetric
    MLAsymmKing king("king","Asymmetric King function", *tmpX, *tmpY, xmean, ymean, lambda, *sigmax, *sigmay) ;
    //symmetric
    //MLAsymmKing king("king","Asymmetric King function", *tmpX, *tmpY, xmean, ymean, lambda, *sigmax);
    // Some variables that we'll be needing
    double AzVal, AzRng1Lower, AzRng1Upper, AzRng2Lower, AzRng2Upper ;
    double EtrueLo, EtrueHi ;
    float mswLo=0.8; float mswGamma=1.1; float mswHi=1.3;

    // Create a cut string        
    RooDataSet* mswCutDataSet;
    RooDataSet* azCutDataSet;
    //    RooDataSet* energyCutDataSet;
    std::ostringstream mswgammacutstring,mswbkgcutstring,energycutstring ;
    //std::ostringstream mswbkgcutstring ;
    mswgammacutstring << "(" << tmpMSW->GetName() << ">=" << mswLo << ")&&"                                                                                     
          << "(" << tmpMSW->GetName() << "<=" << mswGamma << ")" ;                                                                                             
    mswbkgcutstring << "(" << tmpMSW->GetName() << ">" << mswGamma << ")&&"
		    << "(" << tmpMSW->GetName() << "<=" << mswHi << ")" ;


    //azcutsring << "(" << tmpAz->GetName() << ">=" << AzRng1Lower << ")&&"
    //	       << "(" << tmpAz->GetName() << "<" << AzRng1Upper << ")" ;
    energycutstring << "(" << tmpEtrue->GetName() << ">=" << EtrueLo << ")&&"
		    << "(" << tmpEtrue->GetName() << "<" << EtrueHi << ")" ;
    // Reduce the data based on the telescope cut variable. This will hopefully save time later on                                                                  
    // so that we wont have to cut on this same parameter over and over and over and over again.                                                                    
    std::ostringstream telCutString ;
    RooDataSet* telCutData;
    if(numTels<34) {
      std::cout << "Num of Tels: " << numTels << std::endl;
      telCutString << "(" << tmpTels->GetName() << ">" << numTels-0.5 << ")&&"
		   << "(" << tmpTels->GetName() << "<" << numTels+0.5 << ")" ;
      telCutData = dynamic_cast<RooDataSet*>(_dataset[0]->reduce(telCutString.str().c_str()) ) ;
    }
    else {
      std::cout << "Num of Tels: " << numTels << std::endl;
    }
    setEtrueBinning(ebinsbounds);
    std::cout << "Size of Energy Vector " << ebinsbounds.size()  << std::endl;
    bool FitStatusBinned[ebinsbounds.size()-1][NoAzbins];
    std::cout << "First Bin of Energy Vector " << ebinsbounds.at(0)  << std::endl;
    std::cout << "Second Bin of Energy Vector " << ebinsbounds.at(1)  << std::endl;

    char canvasnamex[100],canvasnamey[100],savepsfx[100],savepsfy[100];
    char histname[100],title[100],outnamex[100],outnamey[100],testname[100]; 
    char subtitle[50], cdfname[100], savecdfname[100]; TText T;
    char nameking[100], namedata[100];
    //TCanvas * plotskyMSWbins[azbin][ebin][mswbin];
    TCanvas * psfcanvasx[4][60][2]; TCanvas * psfcanvasy[4][60][2]; TCanvas * cdfcanvasx[4][60][2];
    TCanvas * test[4][60][2]; TH1 * hist[4][60][2]; TH2F * hist2d[4][60][2];
    RooPlot * mkplotsx[4][60][2]; RooPlot * mkplotsy[4][60][2];
    RooPlot * mkplotsxpull[4][60][2]; RooPlot * mkplotsypull[4][60][2];
    RooPlot * cdfframe[4][60][2];
    RooHist *hpull;
    double xxyy, tmeanX, tmeanY;

    // Loop over each point in azimuth                                                                                                    
    for (int aIndx=0; aIndx<4; aIndx++) { 
      // The azimuth value we are interested in is the lower bound on this bin                                                             
      AzVal = AzRng1LowerBin[aIndx] ;
      AzRng1Lower = AzVal ;
      AzRng1Upper = AzRng1UpperBin[aIndx] ;
      //  int AzRng1LowerBin[NoAzbins] = {0,90,180,270};
      //int AzRng1UpperBin[NoAzbins] = {90,180,270,359};
      std::ostringstream azcutstring;
      // If we're observing near zenith, include all values of azimuth to prevent limiting statistics                                     
      if(_avgZenith<45.0) {
	AzRng1Lower = 0.0 ;
	AzRng1Upper = 360.0 ;
	azcutstring << "(" << tmpAz->GetName() << ">=" << AzRng1Lower << ")&&"
		    << "(" << tmpAz->GetName() << "<" << AzRng1Upper << ")" ;
      }
      else {
	//select north (270-90)
	if((AzRng1Upper==90.0)||(AzRng1Upper==359.0)) {
	  AzRng1Lower = 0.0 ;
	  AzRng1Upper = 90.0 ;
	  AzRng2Lower = 270.0 ;
	  AzRng2Upper = 359.0 ;
	  azcutstring << "(" << "(" << tmpAz->GetName() << ">=" << AzRng1Lower << ")&&"
		      << "(" << tmpAz->GetName() << "<" << AzRng1Upper << ")" << ")||("
		      << "(" << tmpAz->GetName() << ">=" << AzRng2Lower << ")&&"
                      << "(" << tmpAz->GetName() << "<" << AzRng2Upper << ")" << ")";
	}
	//select south (90-270)
	if((AzRng1Upper==180.0)||(AzRng1Upper==270.0)) {
	  AzRng1Lower = 90.0 ;
	  AzRng1Upper = 270.0 ;
	  azcutstring << "(" << tmpAz->GetName() << ">=" << AzRng1Lower << ")&&"
		      << "(" << tmpAz->GetName() << "<" << AzRng1Upper << ")" ;
	}
      }
      if(numTels<34){
	azCutDataSet = dynamic_cast<RooDataSet*>(telCutData->reduce(azcutstring.str().c_str()) ) ;
      }
      else {
	azCutDataSet = dynamic_cast<RooDataSet*>(_dataset[0]->reduce(azcutstring.str().c_str()) ) ;
      }

      for (int mswIndx=0; mswIndx<2; mswIndx++) {
	if(mswIndx<1.0) {mswCutDataSet = dynamic_cast<RooDataSet*>(azCutDataSet->reduce(mswgammacutstring.str().c_str()) ) ;}
	if(mswIndx>0.0) {mswCutDataSet = dynamic_cast<RooDataSet*>(azCutDataSet->reduce(mswbkgcutstring.str().c_str()) ) ;}
	
	//for (int eIndx=0; eIndx<NofEbins; eIndx++) {                                                                              
	for (int eIndx=0; eIndx<(ebinsbounds.size()-1); eIndx++) {
	  // Get the lower and upper values of this bin in actual true energy                                                                                       
	  EtrueLo = ebinsbounds.at(eIndx);
	  EtrueHi = ebinsbounds.at(eIndx+1);
	  std::cout << "Bin " << eIndx << " ETrueLo " << EtrueLo << " ETrueHi " << EtrueHi << std::endl;
	  
	  // Create a cut string                                                                                                                     
	  std::ostringstream energycutstring ;
	  energycutstring << "(" << tmpEtrue->GetName() << ">=" << EtrueLo << ")&&"
			  << "(" << tmpEtrue->GetName() << "<" << EtrueHi << ")" ;

	  // Get a reduced version of the dataset by reducing on the cuts string                             
	  RooDataSet* energyCutDataSet = dynamic_cast<RooDataSet*>(mswCutDataSet->reduce(energycutstring.str().c_str()) ) ;

	  if ((EtrueLo>-0.55)&&(EtrueLo<1.00)&&(energyCutDataSet->numEntries() >=10)) {
	    //Create example rotation matrices
	    /*
	    //float theta = 3.14159/4; //45 degrees
	    //Eigen::Rotation2Df t(theta);
	    //std::cout << "THETA TEST: " << t.toRotationMatrix() << std::endl;
	    //Eigen::Vector2f v; v << 1.0f,0.0f; //(1,0)
	    //Eigen::Vector2f rotatedVect = t.toRotationMatrix()*v;
	    //std::cout << "rotatedVect: " << rotatedVect << std::endl;

	    //positive X direction
	    Eigen::VectorXcd xd(2); xd[0] = 1.0; xd[1] = 0.0;
	    //positive Y direction
	    Eigen::VectorXcd yd(2); yd[0] = 0.0; yd[1] = 1.0;
	    Eigen::MatrixXd A(5,2);
	    Eigen::MatrixXd B(5,2);
	    Eigen::MatrixXd C(5,2);
	    Eigen::MatrixXd D(5,2);
	    for (int row = 0; row < 5; ++row)
	      {
		A(row,0) = -1*(row+1); //x
		A(row,1) = -1*(row+1);
		B(row,0) = 1*(row+1); //x
                B(row,1) = 1*(row+1);
                C(row,0) = -1*(row+1); //x
                C(row,1) = 1*(row+1);
                D(row,0) = 1*(row+1); //x
                D(row,1) = -1*(row+1);

	      }
	    //std::cout << "matrix= " << A << std::endl;
	    Eigen::MatrixXd centered_A = A.rowwise()-A.colwise().mean();
	    Eigen::Matrix2d cov_A = (centered_A.adjoint() * centered_A);
	    Eigen::MatrixXd centered_B = B.rowwise()-B.colwise().mean();
	    Eigen::Matrix2d cov_B = (centered_B.adjoint() * centered_B);
	    Eigen::MatrixXd centered_C = C.rowwise()-C.colwise().mean();
	    Eigen::Matrix2d cov_C = (centered_C.adjoint() * centered_C);
	    Eigen::MatrixXd centered_D = D.rowwise()-D.colwise().mean();
	    Eigen::Matrix2d cov_D = (centered_D.adjoint() * centered_D);

	    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver_A(cov_A);
	    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver_B(cov_B);
	    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver_C(cov_C);
	    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver_D(cov_D);

	    Eigen::Rotation2D<double> rot2dA(cov_A);
	    Eigen::Rotation2D<double> rot2dB(cov_B);
	    Eigen::Rotation2D<double> rot2dC(cov_C);
	    Eigen::Rotation2D<double> rot2dD(cov_D);
	    //std::cout << "angle: " << rot2d.fromRotationMatrix() <<  std::endl;
	    std::cout << "matrix = " << std::endl;
	    std::cout << A << std::endl;
	    std::cout << "The rotation angle: " << rot2dA.angle() << " " << rot2dA.angle()*r2d <<  std::endl;
	    std::cout << "The smallest positive angle in [0,2pi]: " << rot2dA.smallestPositiveAngle()  << " " << rot2dA.smallestPositiveAngle()*r2d <<  std::endl;
	    
	    std::cout << "matrix = " << std::endl;
	    std::cout << B << std::endl;
	    std::cout << "The rotation angle: " << rot2dB.angle() << " " << rot2dB.angle()*r2d  <<  std::endl;
	    std::cout << "The smallest positive angle in [0,2pi]: " << rot2dB.smallestPositiveAngle()  << " " << rot2dB.smallestPositiveAngle()*r2d <<  std::endl;

	    std::cout << "matrix = " << std::endl;
	    std::cout << C << std::endl;
	    std::cout << "The rotation angle: " << rot2dC.angle() << " " << rot2dC.angle()*r2d  <<  std::endl;
	    std::cout << "The smallest positive angle in [0,2pi]: " << rot2dC.smallestPositiveAngle()  << " " << rot2dC.smallestPositiveAngle()*r2d <<  std::endl;

	    std::cout << "matrix = " << std::endl;
	    std::cout << D << std::endl;
	    std::cout << "The rotation angle: " << rot2dD.angle() << " " << rot2dD.angle()*r2d  <<  std::endl;
	    std::cout << "The smallest positive angle in [0,2pi]: " << rot2dD.smallestPositiveAngle()  << " " << rot2dD.smallestPositiveAngle()*r2d <<  std::endl;
	    */
	    int esize = energyCutDataSet->numEntries();
	    std::cout << "Number of events in the sub-roodataset: " << esize << std::endl;
	    Eigen::MatrixXd AA(energyCutDataSet->numEntries(),2);
	    std::cout << "Created Matrix for x and y positions " <<std::endl;
	    for (int row = 0; row < energyCutDataSet->numEntries(); ++row)
              {
		const RooArgSet * looprow = energyCutDataSet->get(row);
		RooRealVar * xloop = (RooRealVar*)looprow->find("xpos");
		RooRealVar * yloop = (RooRealVar*)looprow->find("ypos"); 
		//std::cout << "Event " << row << " X: " << xloop->getVal() << " Y: " << yloop->getVal() << std::endl;
		AA(row,0) = (xloop->getVal()); //x
		AA(row,1) = yloop->getVal(); //y
              }
	    //std::cout << "matrix= " << std::endl;
	    //std::cout << AA << std::endl;
	    Eigen::MatrixXd centeredAA = AA.rowwise()-AA.colwise().mean();
	    Eigen::Matrix2d covAA = (centeredAA.adjoint() * centeredAA);
            //Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(A.rows() - 1);                                                 
	    std::cout << "cov matrix= " << std::endl;
	    std::cout << covAA << std::endl;

	    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolverAA(covAA);
            //if (eigensolver.info() != Success) abort();                                                                                    
	    std::cout << "The eigenvalues of A are:\n" << eigensolverAA.eigenvalues() << std::endl;
	    std::cout << "Here's a matrix whose columns are eigenvectors of A \n"
		      << "corresponding to these eigenvalues:\n"
                      << eigensolverAA.eigenvectors() << std::endl;
	    
	    std::cout << "The first eigenvector of the covar matrix of ones is:"
                 << std::endl << eigensolverAA.eigenvectors().col(0) << std::endl;
	    std::cout << "The second eigenvector of the covar matrix of ones is:"
		      << std::endl << eigensolverAA.eigenvectors().col(1) << std::endl;


            //if (eigensolver.info() != Success) abort();                                                                                   
            //std::cout << "The eigenvalues of A are:\n" << eigensolver_A.eigenvalues() << std::endl;                                       
            //std::cout << "Here's a matrix whose columns are eigenvectors of A \n"                                                         
            //   << "corresponding to these eigenvalues:\n"                                                                                 
            //        << eigensolver_A.eigenvectors() <
	    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolvertest(covAA);
	    Eigen::Matrix2d testmatrix = eigensolvertest.eigenvectors();

	    Eigen::Rotation2D<double> rot2d(testmatrix);
	    std::cout << "The rotation angle: " << rot2d.angle() << " " << rot2d.angle()*r2d << std::endl;
	    std::cout << "The smallest positive angle in [0,2pi]: " << rot2d.smallestPositiveAngle()  << " " << rot2d.smallestPositiveAngle()*r2d <<  std::endl;
	    double passtheta = rot2d.smallestPositiveAngle()-(90*d2r);
	    king.SetTheta(passtheta);

	    // Fit the king function to the reduced dataset
	    int trials(0) ;
	    RooFitResult* fitresult ;
	    // Try the fit 10 times before giving it up as a lost cause
	    while ( trials++ < 10) {
	      fitresult = king.fitTo(*energyCutDataSet, RooFit::PrintLevel(-1), RooFit::Offset(true), RooFit::Save(true)) ;
	    }         
	    // Assess whether the fit has actually successfully converged
	    if ((fitresult->covQual() == 3) && (fitresult->status() == 0)) {
	      //FitStatus = true ;
	      FitStatusBinned[eIndx][aIndx] = true;
	    }                
	    std::cout << "CovQual " << fitresult->covQual() << " FitStatus " << fitresult->status() << " MSWbin " << MSWbin[mswIndx] << " Abin " << Abins[aIndx] << " Ebin " << eIndx << " : " << EtrueLo << " " << EtrueHi << std::endl;
	    std::printf("Indices (%d) ", int(FitStatusBinned[eIndx][aIndx])) ;

	    // Save the fitted value of sigma and lambda from the king function 
	    //READOUT
	    std::printf(" MSWLo  |  MSWHi  |  ELo  |  Ehi  |  SigmaX  |  SigmaY  | Lambda | Counts\n" ) ;
	    std::printf(" %1.1f  |  %1.1f  |  %1.3f  |  %1.3f  |   %6.5f  |  %6.5f  |  %6.5f   |   %d\n\n", MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi, sigmax->getVal(), sigmay->getVal(), lambda.getVal(), energyCutDataSet->numEntries()) ;

	    std::printf(" MSWLo  |  MSWHi  |  ETrueLo | ETrueHi | Zenith | Avg. Azimuth | Avg. Noise | Offset | SigmaX | SigmaY | Lambda |  Counts\n" ) ;
	    std::printf(" %1.1f  |  %1.1f  |  %1.3f  |  %1.3f  |   %2.0f   |  %f   |    %2.2f    |  %1.2f  |   %6.5f   |  %6.5f  |  %6.5f   |   %d\n", MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi,_avgZenith, _avgAzimuth, _avgNoise, _offset, sigmax->getVal(), sigmay->getVal(), lambda.getVal(), energyCutDataSet->numEntries()) ;

	    //FILE OUTPUT
	    sprintf( out_consts, "%2.0f\t%1.2f\t%2.2f\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.2f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%1.3f\t%d\n", _avgZenith, _offset, _avgNoise, AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi, sigmax->getVal(), sigmax->getError(), sigmay->getVal(), sigmay->getError(), lambda.getVal(), passtheta, energyCutDataSet->numEntries());
	    if((EtrueLo>-0.650)&&(EtrueLo<1.05)) {myfile << out_consts;}

	    sprintf(histname, "PSF_az_%d_msw_%1.1f-%1.1f_energy_%1.3f-%1.3f", AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);      
	    sprintf(title, "PSF for Azimuth %d MSW %1.1f-%1.1f Energy %1.3f-%1.3f", AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);

	    TLine *ll; TLine *l1; TLine *l2; TLine *l1int; TLine *l2int;
            l1 = new TLine(0.1,0.0,0.1,0.35); l1->SetLineColor(kRed); l1->SetLineStyle(10);
            l2 = new TLine(-0.1,0.0,-0.1,0.35); l2->SetLineColor(kRed); l2->SetLineStyle(10);
            l1int = new TLine(0.15,0.0,0.15,0.35); l1int->SetLineColor(kRed); l1int->SetLineStyle(1);
            l2int = new TLine(-0.15,0.0,-0.15,0.35); l2int->SetLineColor(kRed); l2int->SetLineStyle(1);

	    float scale, intdata, intking, intdiff, intdatax, intdatay, intpdfy, intpdfx, intxdiff, intydiff;
            TPaveText* paveText = new TPaveText(0.64,0.70,0.82,0.88,"NDC");
            paveText->SetBorderSize(0.0);
            paveText->SetFillColor(kWhite);
            paveText->SetFillStyle(0);
            paveText->SetTextSize(0.038);
            paveText->SetTextColor(kRed); float chisqr;

	    const RooArgSet* rowplot = energyCutDataSet->get();
	    RooRealVar * xplot = (RooRealVar*)rowplot->find("xpos");
	    RooRealVar * yplot = (RooRealVar*)rowplot->find("ypos");//
	    //RooArgSet nset(*(energyCutDataSet->get()));
	    
	    //// Create and fill ROOT 2D histogram (20x20 bins) with contents of dataset
	    TH2F *hh_data = dynamic_cast<TH2F*>(energyCutDataSet->createHistogram("xpos,ypos",100,100));
	    scale = 1/(hh_data->Integral());
	    hh_data->Scale(scale);
	    //This histogram is already scaled
	    TH2F *hh_pdf = dynamic_cast<TH2F*>(king.createHistogram("xpos,ypos",100,100));

	    intking=hh_pdf->Integral(48,52,48,52);
	    intdata=hh_data->Integral(48,52,48,52);
	    intdiff=(intdata/intking);

	    TH1D *hh_data1dx = hh_data->ProjectionX(); TH1D *hh_data1dy = hh_data->ProjectionY();
	    TH1D *hh_pdf1dx = hh_pdf->ProjectionX(); TH1D *hh_pdf1dy = hh_pdf->ProjectionY();

	    hh_data1dx->SetLineColor(kBlack); hh_data1dy->SetLineColor(kBlack);
	    hh_data1dx->SetMarkerStyle(20); hh_data1dy->SetMarkerStyle(20);
	    hh_pdf1dx->SetMarkerStyle(kCircle); hh_pdf1dy->SetMarkerStyle(kCircle);
	    hh_pdf1dx->SetLineColor(kBlue); hh_pdf1dy->SetLineColor(kBlue);
	    hh_pdf1dx->SetMarkerColor(kBlue); hh_pdf1dx->SetMarkerColor(kBlue);

	    intdatax=hh_data1dx->Integral(48,52); intdatay=hh_data1dy->Integral(48,52);
	    intpdfx=hh_pdf1dx->Integral(48,52); intpdfy=hh_pdf1dy->Integral(48,52);
	    intxdiff=(intdatax/intpdfx); intydiff=(intdatay/intpdfy);

	    std::cout << "Integral of king function: " << intking << " Integral of 2d dist: " << intdata << " intdata/intking " << intdiff << std::endl;
	    sprintf( out_integral, " %2.0f\t%1.2f\t%2.2f\t%d\t%1.1f\t%1.1f\t%2.2f\t%1.2f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%d\n", _avgZenith, _offset, _avgNoise, AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi, intdatax, intdatay, intpdfx, intpdfy, intdata, intking, passtheta, energyCutDataSet->numEntries());
	    if((EtrueLo>-0.650)&&(EtrueLo<1.05)) {intfile << out_integral;}

	    //Create cumulative distributions
	    //Return value of king function normalized over x in range [-2,2]
	    //std::cout << "king function unnormalized integral " << king.getVal() << std::endl;
	    RooAbsReal* kingcdf = king.createCdf(*(energyCutDataSet->get()));
	    //std::cout << "king function normalized integral " << king.getVal(&nset) << std::endl;
	    //error: no matching function for call to ???MLAsymmKing::createIntegral(const RooArgSet*)???
	    //RooAbsReal* createIntegral(const RooArgSet& iset, const char* rangeName)
	    RooAbsReal* kingint = king.createIntegral(*rowplot,*rowplot,"xpos");
	    std::cout << "king int = " << kingint->getVal() << std::endl;
	    //RooAbsReal* kingcdf = king.createCdf(rowplot);
	    //Plot cdf of king versus x
	    cdfframe[aIndx][eIndx][mswIndx] = xplot->frame(-0.52,0.52,25);
	    kingcdf->plotOn(cdfframe[aIndx][eIndx][mswIndx]);

	    gStyle->SetOptStat(0);//suppress events/stddev/entries tables 
	    //DRAW PLOT ON CANVAS
	    sprintf(cdfname, "CDF for Azimuth %d MSW %1.1f-%1.1f Energy %1.3f-%1.3f", AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);
	    sprintf(namedata, "XY Data for Azimuth %d MSW %1.1f-%1.1f Energy %1.3f-%1.3f", AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);
	    sprintf(nameking, "King Function for Azimuth %d MSW %1.1f-%1.1f Energy %1.3f-%1.3f", AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);
	    hh_data->SetNameTitle(namedata,namedata);
	    hh_pdf->SetNameTitle(nameking,nameking);
	    sprintf(testname, "plots/test_zenith_%2.0f_noise_%1.2f_offset_%1.2f_az_%d_msw_%1.1f-%1.1f_energy_%1.3f-%1.3f.pdf", _avgZenith, _avgNoise, _offset, AzRngBin[aIndx], MSWLo[mswIndx],MSWHi[mswIndx], EtrueLo, EtrueHi);
	    sprintf(savecdfname, "plots/cdfx_zenith_%2.0f_noise_%1.2f_offset_%1.2f_az_%d_msw_%1.1f-%1.1f_energy_%1.3f-%1.3f.pdf", _avgZenith, _avgNoise, _offset, AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);
	    test[aIndx][eIndx][mswIndx] = new TCanvas("test","test", 1000, 1000);
	    test[aIndx][eIndx][mswIndx]->Divide(2,2);

	    test[aIndx][eIndx][mswIndx]->cd(1);
	    hh_data->GetXaxis()->SetRangeUser(-0.5, 0.5);
	    hh_data->GetYaxis()->SetRangeUser(-0.5, 0.5);
	    hh_data->Draw("colz");

	    test[aIndx][eIndx][mswIndx]->cd(2);
	    hh_pdf->GetXaxis()->SetRangeUser(-0.5, 0.5);
            hh_pdf->GetYaxis()->SetRangeUser(-0.5, 0.5);
	    hh_pdf->Draw("colz");

	    test[aIndx][eIndx][mswIndx]->cd(3);
	    hh_data1dx->SetNameTitle("Xproj","Xproj");
	    hh_data1dx->GetXaxis()->SetRangeUser(-0.5, 0.5);
	    hh_data1dx->GetYaxis()->SetRangeUser(0.0, 0.35);
	    hh_data1dx->Draw();
	    hh_pdf1dx->Draw("same");
	    l1->Draw(); l2->Draw(); l1int->Draw(); l2int->Draw();

	    test[aIndx][eIndx][mswIndx]->cd(4);
	    hh_data1dy->SetNameTitle("Yproj","Yproj");
	    hh_data1dy->GetXaxis()->SetRangeUser(-0.5, 0.5);
	    hh_data1dy->GetYaxis()->SetRangeUser(0.0, 0.35);
	    hh_data1dy->Draw();
	    hh_pdf1dy->Draw("same");
	    l1->Draw(); l2->Draw(); l1int->Draw(); l2int->Draw();

	    test[aIndx][eIndx][mswIndx]->SaveAs(testname);

	    cdfcanvasx[aIndx][eIndx][mswIndx] = new TCanvas(cdfname,cdfname, 200, 10, 700, 900);
	    cdfcanvasx[aIndx][eIndx][mswIndx]->cd();
	    cdfframe[aIndx][eIndx][mswIndx]->Draw();
	    //cdfcanvasx[aIndx][eIndx][mswIndx]->SaveAs(savecdfname);
	    
	    //X
	    ll = new TLine(-0.5,0.0,0.5,0.0); ll->SetLineColor(kRed); ll->SetLineStyle(10);
	    mkplotsx[aIndx][eIndx][mswIndx] = xplot->frame(-0.52,0.52,25);
	    energyCutDataSet->plotOn(mkplotsx[aIndx][eIndx][mswIndx]);                                                                  
            king.plotOn(mkplotsx[aIndx][eIndx][mswIndx]);                                                                               
            mkplotsx[aIndx][eIndx][mswIndx]->SetXTitle("Projected X Axis");                                                             
            mkplotsx[aIndx][eIndx][mswIndx]->SetYTitle("events");   

	    hpull = mkplotsx[aIndx][eIndx][mswIndx]->pullHist();
	    mkplotsxpull[aIndx][eIndx][mswIndx] = xplot->frame(-0.52,0.52,25);
	    mkplotsxpull[aIndx][eIndx][mswIndx]->addPlotable(hpull,"P");
	    paveText->AddText(Form("CHI-SQ = %1.3f ", (mkplotsx[aIndx][eIndx][mswIndx]->chiSquare())));
            chisqr=(mkplotsx[aIndx][eIndx][mswIndx]->chiSquare());
	    sprintf(canvasnamex, "X Proj PSF for Azimuth %d MSW %1.1f-%1.1f Energy %1.3f-%1.3f", AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);
	    sprintf( savepsfx, "plots/psfx_zenith_%2.0f_noise_%1.2f_offset_%1.2f_az_%d_msw_%1.1f-%1.1f_energy_%1.3f-%1.3f.pdf", _avgZenith, _avgNoise, _offset, AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);
	    psfcanvasx[aIndx][eIndx][mswIndx] = new TCanvas(canvasnamex, canvasnamex, 200, 10, 700, 900);
            mkplotsx[aIndx][eIndx][mswIndx]->SetTitle(canvasnamex);
            psfcanvasx[aIndx][eIndx][mswIndx]->SetFillColor(10);
	    psfcanvasx[aIndx][eIndx][mswIndx]->Divide(1,2);
	    psfcanvasx[aIndx][eIndx][mswIndx]->cd(1);
            mkplotsx[aIndx][eIndx][mswIndx]->Draw();
            paveText->Draw();
            psfcanvasx[aIndx][eIndx][mswIndx]->cd(2);
            mkplotsxpull[aIndx][eIndx][mswIndx]->SetName(canvasnamex);
            mkplotsxpull[aIndx][eIndx][mswIndx]->Draw();
	    ll->Draw();
	    //psfcanvasx[aIndx][eIndx][mswIndx]->SaveAs(savepsfx);
	    sprintf( outnamex, "plots/psfx_zenith_%2.0f_noise_%1.2f_offset_%1.2f_az_%d_msw_%1.1f-%1.1f_energy_%1.3f-%1.3f.root", _avgZenith, _avgNoise, _offset, AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);
	    TFile *rfoutx = new TFile( outnamex, "recreate");
            psfcanvasx[aIndx][eIndx][mswIndx]->Write();
            rfoutx->Write();
            rfoutx->Close();
	    
	    /*
	    //Y
	    ll = new TLine(energyCutDataSet->mean(*tmpY)-0.5,0.0,energyCutDataSet->mean(*tmpY)+0.5,0.0);
	    ll->SetLineColor(kRed); ll->SetLineStyle(10);
	    mkplotsy[aIndx][eIndx][mswIndx] = yplot->frame(energyCutDataSet->mean(*tmpY)-0.52,energyCutDataSet->mean(*tmpY)+0.52,25);// 

	    energyCutDataSet->plotOn(mkplotsy[aIndx][eIndx][mswIndx]);//                                                                  
            king.plotOn(mkplotsy[aIndx][eIndx][mswIndx]);//                                                                               
            mkplotsy[aIndx][eIndx][mswIndx]->SetXTitle("Projected Y Axis");//                                                             
            mkplotsy[aIndx][eIndx][mswIndx]->SetYTitle("events");//  
	    
	    // Construct a histogram with the pulls of the data w.r.t the curve
	    hpull = mkplotsy[aIndx][eIndx][mswIndx]->pullHist();
	    // Create a new frame to draw the pull distribution and add the distribution to the frame
	    mkplotsypull[aIndx][eIndx][mswIndx] = yplot->frame(energyCutDataSet->mean(*tmpY)-0.52,energyCutDataSet->mean(*tmpY)+0.52,25);
	    mkplotsypull[aIndx][eIndx][mswIndx]->addPlotable(hpull,"P");

	    paveText->AddText(Form("CHI-SQ = %1.3f ", (mkplotsy[aIndx][eIndx][mswIndx]->chiSquare()))); 
	    chisqr=(mkplotsy[aIndx][eIndx][mswIndx]->chiSquare());

	    sprintf(canvasnamey, "Y Proj PSF for Azimuth %d MSW %1.1f-%1.1f Energy %1.3f-%1.3f", AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);//
	    psfcanvasy[aIndx][eIndx][mswIndx] = new TCanvas(canvasnamey, canvasnamey, 0,0,600,600);//
	    mkplotsy[aIndx][eIndx][mswIndx]->SetTitle(canvasnamey);//
	    psfcanvasy[aIndx][eIndx][mswIndx]->Divide(1,2);
	    psfcanvasy[aIndx][eIndx][mswIndx]->cd(1);
	    mkplotsy[aIndx][eIndx][mswIndx]->Draw();
	    paveText->Draw();
	    psfcanvasy[aIndx][eIndx][mswIndx]->cd(2);                                                                                   
	    mkplotsypull[aIndx][eIndx][mswIndx]->SetName(canvasnamex);
            mkplotsypull[aIndx][eIndx][mswIndx]->Draw();
	    ll->Draw();
	    sprintf( savepsfy, "plots/psfy_az_%d_msw_%1.1f-%1.1f_energy_%1.3f-%1.3f.pdf", AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);//
	    psfcanvasy[aIndx][eIndx][mswIndx]->SaveAs(savepsfy);//

	    sprintf( outnamey, "plots/psfy_az_%d_msw_%1.1f-%1.1f_energy_%1.3f-%1.3f.root", AzRngBin[aIndx], MSWLo[mswIndx], MSWHi[mswIndx], EtrueLo, EtrueHi);//
	    TFile *rfouty = new TFile( outnamey, "recreate");//                                                                           
            psfcanvasy[aIndx][eIndx][mswIndx]->Write();// 
	    rfouty->Write();//
	    rfouty->Close();//
	    */
	  }
	}  
      }
    }
    std::cout << "Tracking RA: " << _trackPos[0] << std::endl;                                                               
    std::cout << "Tracking Dec: " << _trackPos[1] << std::endl;                                                              
    std::cout << "Source RA: " << _sourcePos[0] << std::endl;                                                                 
    std::cout << "Source Dec: " << _sourcePos[1] << std::endl;
    //std::cout << "Name: " << xCoord->GetName() << " " << xCoord->getVal() << std::endl;
    
    PrintInfo(ss1,ss2);
    //delete _datatree;
    //delete _dataset;
    rfinRoo[0]->Close();
    rfinRoo[1]->Close();
    myfile.close();
    intfile.close();
    //delete rfinRoo;
  
    return 1;
  }
  return 0;
}


bool AddRooDataSet::cleanupDataSet()
{
  // The goal of this method is to basically reinitialize all                                                                           
  // the variables associated with this dataset as well as                                                                              
  // reset the dataset                                                                                                                  

  // set the name to something generic and reset to not a sim                                                                           
  //SetName("name") ;
  _isSim = false ;

  // start with reinitializing all the data parameters                                                                                  
  _ATM      = 0 ;
  _epoch    = 0 ;
  _livetime = 0 ;

  // Move the source and tracking VACoordinatePair positions to                                                                         
  // some generic location                                                                                                              
  _trackPos[0] = 0.0 ;
  _trackPos[1] = 0.0 ;
  _sourcePos[0] = 0.0 ;
  _sourcePos[1] = 0.0 ;

  // Now actually cleanup the dataset, but only if we've                                                                                
  // actually set it to be something                                                                                                    
  //if ((_dataset[0] != nullptr)||(_dataset[1] != nullptr)||(_dataset[0] != nullptr)) {
    //        _dataset->reset() ;                                                                                                     
  //}

  return 0 ;
}


Double_t AddRooDataSet::psfeq(Double_t *x, Double_t *par)
{
  // Return a TF2 representing a snapshot of this object                                                                                  
  // Precompute some things                                                                                                               
  Double_t ll = par[0] ;
  Double_t ss = par[1] ;

  // Compute the transformed coordinates                                                                                                  
  Double_t xarg = x[0] - par[2];
  Double_t yarg = x[1] - par[3];
  //double x_prime = xarg*cached_cosTheta_xnorm - yarg*cached_sinTheta_xnorm ;
  // If sigma and lambda are calculated via an MLPointSpreadFunctionVar, then                                                             
  // by caching the values here it removes the need for several calls to those                                                            
  // functions and thus many computations, expecially if they're being both                                                               
  // interpolated and not cached.                                                                                                         
  double norm1 = 0.5/(TMath::Pi()*ss*ss);
  double norm2 = (1.0 - (1.0/ll));

  // Evaluate and return                                                                                                                  
  return norm1*norm2*std::pow(1.0 + ( (0.5/ll) * (((xarg*xarg)/(ss*ss))+((yarg*yarg)/(ss*ss))) ), -ll);
}

void AddRooDataSet::PrintInfo(std::string ss1, std::string ss2)
{
  // Print out some useful information about this run                                                                                   
  std::cout << "----------------------------------------------------" << std::endl ;
  std::cout << "INFORMATION STORED IN RUNS    : " << ss1 << " " << ss2 << std::endl ;
  std::cout << "     Saved Name        : " << _SourceName << std::endl ;
  std::cout << "     Is Simulation     : " << _isSim  << std::endl ;
  std::cout << "     Livetime (Sec)    : " << _livetime << std::endl ;
  std::cout << "     Tracking RA       : " << getTrackingRA_Deg() << std::endl ;
  std::cout << "     Tracking Dec      : " << getTrackingDec_Deg() << std::endl ;
  std::cout << "     Source RA         : " << getSourceRA_Deg() << std::endl;
  std::cout << "     Source Dec        : " << getSourceDec_Deg() << std::endl;
  //std::cout << "     Energy Range [TeV]: " << EnergyLowerBound << "-" << EnergyUpperBound << std::endl;
  std::cout << "     Events Stored     : " << _dataset[0]->numEntries() << std::endl ;
  std::cout << "     Epoch             : " << _epoch << std::endl ;
  std::cout << "     Season (ATM)      : " << _ATM << std::endl ;
  //double avg_cosz = getAvgCosZ() ;
  std::cout << "     Avg. Zenith       : " << _avgZenith << std::endl ;
  std::cout << "     Avg. Azimuth      : " << _avgAzimuth << std::endl ;
  std::cout << "     Avg. Noise        : " << _avgNoise << std::endl ;
  std::cout << "     Offset            : " << _offset << std::endl;

  return ;
}

//_____________________________________________________________________________                                                       
// E-true [TeV]                                                                                                                       
bool AddRooDataSet::setEtrueBinning(std::vector<double> binbounds)
{
  // binbounds is a vector holding the bin boundaries which define the ranges                                                         
  // for bins in Etrue. This will also include the bins at which we do not                                                           
  // have the simulations to fill them.                                                                                               
  // NOTE : For this function, binbounds should be in units of [TeV]                                                                  
  // Handle the case for binbounds being an empty vector (this is the default case)                                                   
  if (binbounds.empty()) {
    // True Energy is broken into bins associated with incremental values                                                             
    // in E_true. The bounds are separated by 0.05 in log(E_true). True Energy for this                                               
    // class is in units of TeV                                                                                                       
    double binwidth = 0.05 ;   // bin width                                                                                           
    double log10EtrueLo = -1.5 ;     // 30 GeV or 0.01 TeV                                                                            
    double log10EtrueHi = 2.0  ;   // 100000 GeV or 100 TeV^S                                                                         
    int numBins = (log10EtrueHi-log10EtrueLo) / binwidth ;
    std::cout << "Number of Energy bins " << numBins << std::endl;

    // Loop over each value in the vector and set it appropriately                                                                    
    for (int boundIndx=0; boundIndx<=numBins; boundIndx++) {
      // set the appropriate bin boundary value as Etrue/[TeV]                                                                        
      //binbounds.emplace_back(std::pow(10.0,log10EtrueLo + (binwidth*boundIndx))) ;
      binbounds.emplace_back(log10EtrueLo + (binwidth*boundIndx));
    }
  }
  // Check that there are at least 2 bounds in "binbounds" (this is the minimum                                                       
  // necessary to establish a single bin)                                                                                             
  if (binbounds.size() < 2) {
    // Well, that's a problem, so throw an error.                                                                                     
    std::cout << "[ERROR] MLBinning::setEtrueBinning() : Number of bin boundaries should be >= 2."  << std::endl;
    //std::string errMsg("[ERROR] MLBinning::setEtrueBinning() : Number of bin boundaries should be >= 2.") ;
    //std::cerr << errMsg << std::endl;
    //throw MLException(errMsg) ;
    return false ;
  }

  // Store the vector in the map of vectors                                                                                           
  ebinsbounds = binbounds ;
  return true ;

}
