//A class that gets a RooDataSet and TTree from a root file                                                                             
// C++ HEADERS                                                                                                                                    
#include <ctime>
#include <time.h>
#include <iostream>
#include <fstream>
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
#include "TList.h"
#include "TMath.h"
#include "TNamed.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// ROOFIT INCLUDES       
#include <TF1.h>                                                                                                                         
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include <RooFitResult.h>

#include "AddRooDataSet.h"
//#include "MLBinning.h"
//#include "MLAsymmKing.h"
//#include "MLPointSpreadFunction.h"

int main(int argc, char * argv[])
{
  
  //Set binning
  //MLBinning* _AxisBinning;

  std::cout << "Name of file list: " << std::string(argv[1]) << std::endl;
  //Determine number of files in list                                                                                                
  std::string fileslist=std::string(argv[1]);
  int Tels = atoi(argv[2]);
  int linenum=0; std::string line,line1,line2;
  std::ifstream readfiles(fileslist);
  std::ifstream loopfiles(fileslist);

  std::ofstream mf,intfile;
  mf.open("kingparams.txt");
  intfile.open("integralfrac.txt");
  //intfile.open("integralfrac.txt");
  //0.8   1.1   20   0.50   4.28   45   -1.100   -1.050   0.06258   2.50000
  //mf << "MSWLow | MSWHi | Zenith | Avg.Noise | Azimuth | Offset |  EnergyLow  |  EnergyHi  |  Sigma  |  Lambda\n";
  mf << "Zenith | Offset | Avg.Noise | Azimuth |  MSWLow  |  MSWHi  |  EnergyLow  |  EnergyHi  |  SigmaX  |  SigmaXError  |  SigmaY  |  SigmaYError  |  Lambda  |  Tilt Angle  |  Num.Counts\n";
  intfile << "Zenith | Offset | Avg.Noise | Azimuth |  MSWLow  |  MSWHi  |  EnergyLow  |  EnergyHi  |  intdatax  | intdatay  |  intkingx  |  intkingy  |  intdata  |  intking  |  Tilt Angle  |  NumCounts\n";
  mf.close();
  intfile.close();
  //char out_consts[150];

  while(std::getline(readfiles, line)){
    linenum++;
  }
  std::cout << "Number of root files in " << fileslist.c_str() << ": " << linenum << std::endl;

  //AddRooDataSet testfile[linenum];

  while(std::getline(loopfiles, line1)){
    std::getline(loopfiles, line2);
    AddRooDataSet* ds = new AddRooDataSet(line1, line2, mf, intfile, Tels);

    //AddRooDataSet* ds3 = new AddRooDataSet(line2, 3);
    //AddRooDataSet* ds4 = new AddRooDataSet(line2, 4);
    //std::cout << "min x " << dataset->xCoord->getMin() << std::endl;
    //std::cout << "min xpos " << ds->xCoord->getMin() << std::endl; 
    //ii++;
    // We're done with the dataset, so delete it.                                                                                 

    delete ds;
  }

  //mf.close();

  return 0;
}
