#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>
#include <random>
//root
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TProfile2D.h"
//anita
#include "AnitaConventions.h"
#include "AnitaGeomTool.h"
#include "RawAnitaHeader.h"
#include "CalibratedAnitaEvent.h"
#include "UsefulAnitaEvent.h"

#include "sineUtils.h"

#include "/Users/brotter/Downloads/bestFit/BestFit.h"

using namespace std;

class ellipseResult {
  
public:
  double centerX;
  double centerY;
  double majorAxis;
  double minorAxis;
  double rotation;
};


void ellipseFitting(TGraph* inGraph,const int numPoints,ellipseResult *result) {

  double points[numPoints*2];

  double x,y;
  for (int pt=0; pt<numPoints; pt++) {
    inGraph->GetPoint(pt,x,y);
    points[pt*2] = x;
    points[pt*2+1] = y;
  }

  cout << "filled input" << endl;

  BestFitIO input;
  input.numPoints = numPoints;
  input.points = points;

  cout << "BestFitIO" << endl;

  BestFitIO output;

  int type = BestFitFactory::Ellipse;
  BestFit *fit = BestFitFactory::Create(type, std::cout);
  fit->Compute(input,output);

  cout << "computed" << endl;

  result->centerX = output.outputFields[BestFitIO::EllipseCentreX];
  result->centerY = output.outputFields[BestFitIO::EllipseCentreY];
  result->majorAxis = output.outputFields[BestFitIO::EllipseMajor];
  result->minorAxis = output.outputFields[BestFitIO::EllipseMinor];
  result->rotation = output.outputFields[BestFitIO::EllipseRotation];
  
  delete fit;

  return;
}

int getBinIndex(int surf, int lab, int rco, int bin) {

  return (surf*8 + lab*2 + rco)*260 + bin;

}

int main(int argc, char** argv) {


  if (argc != 2) {
    cout << "Usage: " << argv[0] << " [out base fileame] " << endl;
    return -1;
  }
  else {
    cout << "Using outfile base" << argv[1] << endl;
  }

  int run = 10105; //432.1MHz
  //  int run = 10024; //567.8MHz (less events)
  char* dataDir;
  if (run < 1000) dataDir = getenv("ANITA3_DATA");
  else            dataDir = getenv("ANITA3_CALDATA");

  
  stringstream name,title;  
  //Events Waveforms
  TChain *rawEventTree = new TChain("eventTree","");
  name.str("");
  name << dataDir << "run" << run << "/eventFile" << run << ".root";
  rawEventTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  cout << "rawEventTree Entries: " << rawEventTree->GetEntries() << endl;

  //Event Headers
  TChain *headTree = new TChain("headTree","");
  name.str("");
  name << dataDir << "run" << run << "/headFile" << run << ".root";
  headTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  cout << "headTree Entries: " << headTree->GetEntries() << endl;


  //Set Branch Addresses
  RawAnitaEvent *rawEvent = NULL;
  rawEventTree->SetBranchAddress("event",&rawEvent);

  RawAnitaHeader *header = NULL;
  headTree->SetBranchAddress("header",&header);

    

  //storage hists
  TH2D *hEllipse[96*260]; //12 surfs, 4 labs/surf, 2rco/lab = 96 * 260 bins
  TGraph *allPoints[96*260];


  for (int surf=0; surf<12; surf++) {
    for (int lab=0; lab<4; lab++) {
      for (int rco=0; rco<2; rco++) {
	for (int bin=0; bin<260; bin++) {
	  int index = getBinIndex(surf,lab,rco,bin);
	  name.str("");
	  name << "Surf" << surf << " Lab" << lab << " Rco" << rco;
	  name << " Bin" << bin << " vs Bin" << bin+1 <<  "; Sum (mV); Difference (mV)";
	  title.str("");
	  title << "hEllipse_s" << surf << "l" << lab << "r" << rco << "b" << bin;
	  hEllipse[index] = new TH2D(title.str().c_str(),name.str().c_str(),
				     201,-200.5,200.5, 201,-200.5,200.5);
	  allPoints[index] = new TGraph();
	  allPoints[index]->SetName(name.str().c_str());
	  allPoints[index]->SetTitle(name.str().c_str());
	  
	  

	}
      }
    }
  }
  cout << "made storage histograms and TGraphs" << endl;


  //where the sine waves are (starting at 1 so subtract 1 later)
  int chans[12] = {4,2,4,2,4,2,4,4,2, 2, 2, 2};

  //loop through entries
  int lenEntries = 10000; // can't do too many, since I'm filling a billion TGraphs
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%10 ==0) {
      cout << entry << "/" << lenEntries << "\r";
      fflush(stdout);
    }
    headTree->GetEntry(entry);
    rawEventTree->GetEntry(entry);

    UsefulAnitaEvent *useful = new UsefulAnitaEvent(rawEvent,WaveCalType::kJustUnwrap,header);
    //I don't want to resample the alfa so I have to turn it's filtering off
    useful->setAlfaFilterFlag(false);

    for (int surf=0; surf<12; surf++) {
      int chan = chans[surf] - 1;

      int usefulIndex = surf*9 + chan;

      //get lab number and rco
      int lab = rawEvent->getLabChip(usefulIndex);
      int rco = rawEvent->getRCO(usefulIndex);

      TGraph *gr = useful->getGraphFromSurfAndChan(surf,chan);
      for (int pt=0; pt<gr->GetN()-1; pt++) {
	int cap = useful->fCapacitorNum[usefulIndex][pt];
	double binLow = gr->GetY()[pt];
	double binHigh = gr->GetY()[pt+1];
	double X = binLow+binHigh;
	double Y = binLow-binHigh;
	int index = getBinIndex(surf,lab,rco,cap);
	hEllipse[index]->Fill(X,Y);
	allPoints[index]->SetPoint(allPoints[index]->GetN(),X,Y);
      }
      
      delete gr;
    }
    delete useful;
  }

  cout << "Finished filling all the graphs and histograms" << endl;

  name.str("");
  name << argv[1] << ".txt";
  ofstream outFile(name.str().c_str());

  name.str("");
  name << argv[1] << ".root";
  TFile *outRootFile = TFile::Open(name.str().c_str(),"recreate");

  for (int surf=0; surf<12; surf++) {
    for (int lab=0; lab<4; lab++) {
      for (int rco=0; rco<2; rco++) {
	for (int bin=0; bin<260; bin++) {
	  int index = getBinIndex(surf,lab,rco,bin);
	  outRootFile->cd();
	  hEllipse[index]->Write();
	  allPoints[index]->Write();
	  
	  outFile << index << " " << surf << " " << lab << " " << rco << " " << bin << " ";

	  ellipseResult *result = new ellipseResult;
	  ellipseFitting(allPoints[index],allPoints[index]->GetN(),result);
	  
	  outFile << result->centerX   << " ";
	  outFile << result->centerY   << " ";
	  outFile << result->majorAxis << " ";
	  outFile << result->minorAxis << " ";
	  outFile << result->rotation  << endl;

	  delete result;
	}
      }
    }
  }

  outRootFile->Close();
  outFile.close();

  return 1;
  
}

  
