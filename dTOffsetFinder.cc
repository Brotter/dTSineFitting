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

using namespace std;

/*

  Ben Rotter November 2016 - University of Hawaii at Manoa

  The dTs probably aren't right.

  So now a good thing to try would be to determine the __X__ offset, instead of the Y offset between
  the fit and the measured point.  If I make a histogram of that vs pt I might be able to find a systematic
  offset which can be applied as a correction.

  I wrote a lot of the equations for doing this in the book, and there are comments throughout.


 */

//Run10105 is a longer version with 540k in it!  Maybe this one works.
#define SINE_RUN 10105						




//lets try to make this a bit easier to read

//GLOBALS:
//I need the ANITA data to be visible to any subroutine, so those will be global
TChain *rawEventTree;
RawAnitaEvent *rawEvent = NULL;
TChain *headTree;
RawAnitaHeader *header  = NULL;

//also I guess I'll make the whole sine fit tree stuff too
TTree *fitTree = new TTree("fitTree","fitTree");
double amp,freq,phase,offset,resid;
int chanIndex,rco,lab,surf,chan,eventNumber;

//and the ped corrections
double pedCorrections[12*8*4*260];


void makeFitTree() {
  fitTree = new TTree("fitTree","fitTree");
  fitTree->Branch("eventNumber",&eventNumber);
  fitTree->Branch("amp",&amp);
  fitTree->Branch("phase",&phase);
  fitTree->Branch("residual",&resid);
  fitTree->Branch("chanIndex",&chanIndex);
  fitTree->Branch("rco",&rco);
  fitTree->Branch("lab",&lab);
  fitTree->Branch("surf",&surf);
  fitTree->Branch("chan",&chan);

  return;
}
  


void openAnitaData() {

  int run = SINE_RUN;

  stringstream name;
  //Events Waveforms
  rawEventTree = new TChain("eventTree","");
  name.str("");
  name << "/Volumes/ANITA3Data/antarctica14/root/run" << run << "/eventFile" << run << ".root";
  rawEventTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  int lenRawEventTree = rawEventTree->GetEntries();
  cout << "rawEventTree Entries: " << lenRawEventTree << endl;

  //Event Headers
  headTree = new TChain("headTree","");
  name.str("");
  name << "/Volumes/ANITA3Data/antarctica14/root/run" << run << "/headFile" << run << ".root";
  headTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  int lenHeadTree = headTree->GetEntries();
  cout << "headTree Entries: " << lenHeadTree << endl;

  if (lenRawEventTree != lenHeadTree) {
    cout << "the header and event trees have different lenghts so BYE THIS ISN'T CORRECT" << endl;
    return;
  }


  //Set Branch Addresses
  rawEventTree->SetBranchAddress("event",&rawEvent);
  headTree->SetBranchAddress("header",&header);

  return;

}


TF1* sineWaveFitter(TGraph *graphToFit) {
  
  //Lets define the sine fit function (I don't know what I should do for the range...)
  //Since the phase is the biggest generator of error, and the freq is very constant, lets try
  /// hard coding those?
  TF1 *sinFit = new TF1("sinFit","[0]*sin(x*0.4321*2.*3.14159+[1])",20.,80.);
  sinFit->SetParName(0,"Amplitude");
  sinFit->SetParameter(0,200.);
  
  sinFit->SetParName(1,"Phase");
  sinFit->SetParameter(1,0);

  int fitStatus = graphToFit->Fit(sinFit,"NQ"); //N=no draw, Q=quiet, M=more (better) fit?
  if (fitStatus != 0) {   //if the fit is no good tell me about it
    cout << "fitStatus != 0" << endl;
    return NULL;
  }
  else {   //otherwise if the fit is good, fill the histograms
    amp = sinFit->GetParameter(0);
    phase = sinFit->GetParameter(1);
	
    double absResidSum = 0;
    for (int pt = 0; pt<graphToFit->GetN(); pt++) {
      absResidSum += pow(graphToFit->GetY()[pt] - sinFit->Eval(graphToFit->GetX()[pt]),2);
    }
    absResidSum /= graphToFit->GetN();
    resid = TMath::Sqrt(absResidSum);

    fitTree->Fill();
  }
  return sinFit;
  
}



int makeOffsetHists(TFile* outFile,int startEntry, int stopEntry) {

  //Channels with the sine wave in them (All Surfs, 1->12):
  //4,2,4,2,4,2,4,4,2,2,2,2
  int surfs[12] = {1,2,3,4,5,6,7,8,9,10,11,12};
  int chans[12] = {4,2,4,2,4,2,4,4,2, 2, 2, 2};

  
  //output storage histograms (need 96 2D histograms I guess... (bin# vs offset))
  stringstream title,name;
  TH2D *storageHists[96];
  TH2D *diffHists[96];
  for (int surfi=0; surfi<12; surfi++) {
    for (int labi=0; labi<4; labi++) {
      for (int rcoi=0; rcoi<2; rcoi++) {
	name.str("");
	name << "surf" << surfi << "_lab" << labi << "_rco" << rcoi;
	title.str("");
	title << name.str() << ";Storage bin;Corrected dT (ns)";
	int dTArrayIndex = surfi*8 + labi*2 + rcoi;
	storageHists[dTArrayIndex] = new TH2D(name.str().c_str(),title.str().c_str(),
					      260,-0.5,259.5, 1001,0.0,1.0);

	name.str("");
	name << "Correction_surf" << surfi << "_lab" << labi << "_rco" << rcoi;
	title.str("");
	title << name.str() << ";Storage bin;Corrected dT (ns)";
	diffHists[dTArrayIndex] = new TH2D(name.str().c_str(),title.str().c_str(),
					      260,-0.5,259.5, 1001,-0.5,0.5);
      }
    }
  }

  TH2D *devVsVolt = new TH2D("devVsVolt","Deviation vs Voltage;Voltage (mV);Deviation (nS)",
			     250,0,3.,  251,-1.,1.);

  TGraph *gWeirdBin = new TGraph();
  gWeirdBin->SetName("gWeirdBin");

  cout << "made storage stuff" << endl;


  //fitTree will have 8 repeated eventnumbers in a row most likely, so I don't have to make all
  // those UsefulAnitaEvents.
  int prevEventNumber = -1;

  //LOOP THROUGH EVENTS
  int lenEntries = headTree->GetEntries();
  if (stopEntry == -1) stopEntry = lenEntries;
  for (int entry=startEntry; entry<stopEntry; entry++) {
    if (entry%10 == 0) {
      cout << entry << "/" << stopEntry-startEntry << "\r";
      fflush(stdout);
      }

    rawEventTree->GetEntry(entry);
    headTree->GetEntry(entry);

    UsefulAnitaEvent *useful = new UsefulAnitaEvent(rawEvent,WaveCalType::kFull,header);

    for (surf=0; surf<12; surf++) {
      chan = chans[surf]-1;

      TGraph *waveform = getCorrectedPedsGraph(useful,surf,chan,pedCorrections);

      TF1* sineFit = sineWaveFitter(waveform);
      if (sineFit == NULL) {
	cout << "Entry" << entry << " ... Bad Fit ... Skipping ..." << endl;
	break;
      }

      //lets make a cut on a sort of "fit goodness" as define by me
      //    if (!(freq<0.44 && freq>0.43 && phase<4 && phase > -4 && amp > 80 && amp < 400
      //	  && residual<100)) {
      //	continue;
      //    }
    
      
      //I need to make a new array of x values and to save time I'll encode the cap num and rco phase in it
      //x point is the "corrected" time point where it matches with the sine wave
      //y point is the cap number.  negative is RCO1, positive is RCO0 ( pow(-1,rco) )
      TGraph *newXArray = new TGraph();
      
      //storage values for x,y and cap number data
      double x,y;
      int prevCapNum = -1;
      
      //see what anitaEventCalibrator thinks the "first" rco phase is
      int rcoFlipper = pow(-1,useful->fRcoArray[surf]);
      
      //loop to create the new, corrected, xArray
      for (int pt=0; pt<waveform->GetN(); pt++) {
	waveform->GetPoint(pt,x,y);
	//if we are within 10% of the maximum, set the correction to -999 and I wont use those
	//(because the fit deviations will be dominated by vNoise)
	//also there is some bug where I can't do things with x<0 because I am dumb so lets just fuck that for a sec
	double xCorrected; //need to define it outside the if/else statements

	int usefulIndex = surf*9+chan;
	int capNum = useful->fCapacitorNum[usefulIndex][pt];
	lab = useful->getLabChip(usefulIndex);
	rco = useful->getRCO(chanIndex);


	if (abs(y) > amp*0.9 || x<1) {
	  xCorrected = -999;
	}
	else {
	  xCorrected = x - findXOffset(amp,0.4321,phase,0,x,y); 
	}

	if (pt!=0 && (capNum-prevCapNum < 0) ) { //this should find the rollaround point
	  rcoFlipper *= -1; }
	prevCapNum = capNum;
	newXArray->SetPoint(newXArray->GetN(),xCorrected,capNum*rcoFlipper);
	
      }

      //Loop to store dT values from the xArray:
      //now that we have the new X array, we have to store the dT values (x0->x1 = dT0)
      for (int pt=0; pt<newXArray->GetN() - 1; pt++) {
	double x1 = newXArray->GetX()[pt+1];
	double x0 = newXArray->GetX()[pt];
	if (x1==-999 || x0==-999) { //if either of those points were in the 10% cut region don't do anything
	  continue;
	}
	double dT = x1-x0;
	int capNum = newXArray->GetY()[pt];
	int rcoToStore = 0;            //if the capNum is positive, it is rco0
	if (capNum<0) rcoToStore = 1;  //else it is rco1 (from previous loop)
	int storageIndex = surf*8 + lab*2 + rcoToStore;
	storageHists[storageIndex]->Fill(abs(capNum),dT);
	double origdT = waveform->GetX()[pt+1] - waveform->GetX()[pt];
	diffHists[storageIndex]->Fill(abs(capNum),origdT-dT);
	devVsVolt->Fill(abs(y)/amp,origdT-dT);
	if (capNum == -67 && surf==11 && lab==3) { //this cap has a two peaked distribution
	  gWeirdBin->SetPoint(gWeirdBin->GetN(),eventNumber,dT);
	}
      }
      delete waveform;
      delete newXArray;
      
    }//end loop through surfs


  } //END LOOP THROUGH EVENTSx
    

  //write out the storage hists and diff hists
  outFile->cd();
  for (int i=0; i<96; i++) {
    storageHists[i]->Write();
    delete storageHists[i];
    diffHists[i]->Write();
    delete diffHists[i];
  }
  //write out the misc stuff
  devVsVolt->Write();
  gWeirdBin->Write();

  //write out the sine fit tree
  fitTree->Write();




  return 1;
}
  


int main(int argc, char** argv) {

  
  string outFileName;
  int startEntry,stopEntry;
  if (argc == 2) {
    outFileName = argv[1];
    cout << "Using " << outFileName << " as output file" << endl;
    cout << "Using all events" << endl;
    startEntry = 0;
    stopEntry = -1;
  }
  if (argc == 4) {
    outFileName = argv[1];
    cout << "Using " << outFileName << " as output file" << endl;
    startEntry = atoi(argv[2]);
    stopEntry  = atoi(argv[3]);
    cout << "Only looking at entries " << startEntry << " through " << stopEntry << endl;
  }
  else {
    cout << "Usage: " << argv[0] << " [Output file name] optional: [entry start #] [entry stop #]" << endl;
    return -1;
  }

  cout << "Starting Physics!" << endl;
  
  TFile *outFile = TFile::Open(outFileName.c_str(),"recreate");
  
  openAnitaData();
  loadPedCorrections(pedCorrections);
  makeFitTree();
  makeOffsetHists(outFile,startEntry,stopEntry);

  outFile->Close();

  return 1;

}
