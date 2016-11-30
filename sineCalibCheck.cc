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

using namespace std;

/*

  Ben Rotter November 2016 - University of Hawaii at Manoa

  So the varience in dTs is three times larger in anita2 and anita3 than it was in anita1

  This code it to do a few things:
  A) Just do a fit on the sine waves measured in the calibration runs (10023 is a good one)
        to see what the residuals are.  It should be really good.
  B) Actually re-calculate the dTs by measuring the fraction of zero crossings in every time bin

  BenS did this, but if I do it too maybe we can get these back to a normal <10% variance.  Otherwise
   we are going to have to make some sort of hit-bus dependant transfer function I guess.

 */

//Run10023 seems okay so I'll just define it for now (432.1MHz for 650k entries, overnight)
// It actually has a bunch of missing header files and the constant "Clock problem found" errors :(
//#define SINE_RUN 10023
//Run10010 is shorter (6656 events), but the same configuration
//#define SINE_RUN 10010
//Run10105 is a longer version with 540k in it!  Maybe this one works.
#define SINE_RUN 10105						

void sineFitting(TFile *outFile){
  
  int run = SINE_RUN;
  
  stringstream name;  
  
  //Events Waveforms
  TChain *rawEventTree = new TChain("eventTree","");
  name.str("");
  name << "/Volumes/ANITA3Data/antarctica14/root/run" << run << "/eventFile" << run << ".root";
  rawEventTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  cout << "rawEventTree Entries: " << rawEventTree->GetEntries() << endl;

  //Event Headers
  TChain *headTree = new TChain("headTree","");
  name.str("");
  name << "/Volumes/ANITA3Data/antarctica14/root/run" << run << "/headFile" << run << ".root";
  headTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  cout << "headTree Entries: " << headTree->GetEntries() << endl;

  //something is wrong with some head files, and they are missing events SUPER ANNOYINGLY
  headTree->BuildIndex("eventNumber");

  //Set Branch Addresses
  RawAnitaEvent *rawEvent = NULL;
  rawEventTree->SetBranchAddress("event",&rawEvent);

  RawAnitaHeader *header = NULL;
  headTree->SetBranchAddress("header",&header);
  
  //Channels with the sine wave in them (All Surfs, 1->12):
  //4,2,4,2,4,2,4,4,2,2,2,2
  int surfs[12] = {1,2,3,4,5,6,7,8,9,10,11,12};
  int chans[12] = {4,2,4,2,4,2,4,4,2, 2, 2, 2};




  //then lets make a histogram for each of those parameters, plus the residuals (which I'll need to calc)
  TH1D *hAmp = new TH1D("hAmp","Fit Amplitudes",100,0,500);
  TH1D *hFreq = new TH1D("hFreq","Fit Frequencies",1000,0.400,0.500);
  TH1D *hPhase = new TH1D("hPhase","Fit Phases",200,-1*M_PI-1,M_PI+1);
  TH1D *hOffset = new TH1D("hOffset","Fit Offsets",100,-10,10);
  TH1D *hResid = new TH1D("hResid","Fit Residual RMS",1000,0,50);

  //Actually lets make a whole damn tree!
  TTree *fitTree = new TTree("fitTree","fitTree");
  //and add some of those things to it
  double amp,freq,phase,offset,resid;
  int chanIndex,rco,lab,surf,chan,eventNumber;
  fitTree->Branch("eventNumber",&eventNumber);
  fitTree->Branch("amp",&amp);
  fitTree->Branch("freq",&freq);
  fitTree->Branch("phase",&phase);
  fitTree->Branch("offset",&offset);
  fitTree->Branch("residual",&resid);
  fitTree->Branch("chanIndex",&chanIndex);
  fitTree->Branch("rco",&rco);
  fitTree->Branch("lab",&lab);
  fitTree->Branch("surf",&surf);


  //get the noise RMS to normalize (function to generate in this file)
  TFile *noiseFile = TFile::Open("noise.root");
  TH2D *allNoiseRMS = (TH2D*)noiseFile->Get("rmsHist");

  double chanNoise[12][4][2];
  for (int surf=0; surf<12; surf++) {
    for (int lab=0; lab<4; lab++) {
      for (int rco=0; rco<2; rco++) {
	int chanIndex = surf*8 + lab*2 + rco + 1;
	chanNoise[surf][lab][rco] = allNoiseRMS->ProjectionY("temp",chanIndex,chanIndex+1)->GetMean();
	cout << surf << " " << lab << " " << rco << " : " << chanNoise[surf][lab][rco] << endl;
      }
    }
  }



  int lenEntries = headTree->GetEntries();
  //  int lenEntries = 5000;
  for (int entry=0; entry<lenEntries; entry++) {
    cout << "Entry:" << entry << "/" << lenEntries << "\r";
    fflush(stdout);
    rawEventTree->GetEntry(entry);
    headTree->GetEntryWithIndex(rawEvent->eventNumber);
    if (rawEvent->eventNumber != header->eventNumber) {
      cout << "Event Number Mismatch! " << rawEvent->eventNumber << "!=" <<  header->eventNumber << endl;
      continue;
    }
    eventNumber = rawEvent->eventNumber;

    UsefulAnitaEvent *usefulRawEvent = new UsefulAnitaEvent(rawEvent,WaveCalType::kFull,header);
    //I don't want to resample the alfa so I have to turn it's filtering off
    usefulRawEvent->setAlfaFilterFlag(false);
    
    for (int surfi=0; surfi<12; surfi++) {
      surf = surfs[surfi]-1;
      chan = chans[surfi]-1;
      TGraph *gr = usefulRawEvent->getGraphFromSurfAndChan(surf,chan);

      //from simpleStructs.h -> chanIndex = chan+9*surf
      chanIndex = chan+9*surf;
      rco = rawEvent->getRCO(chanIndex);
      lab = rawEvent->getLabChip(chanIndex);


      //Lets define the sine fit function (I don't know what I should do for the range...)
      TF1 *sinFit = new TF1("sinFit","[0]*sin(x*[1]+[2])+[3]",20.,80.);
      sinFit->SetParName(0,"Amplitude");
      sinFit->SetParameter(0,200.);

      //      sinFit->SetParName(1,"FrequencyX2Pi");
      //      sinFit->SetParameter(1,2.*M_PI*0.4321); //in GHz since x is in ns

      sinFit->SetParName(2,"Phase");
      sinFit->SetParameter(2,0);

      sinFit->SetParName(3,"Offset");
      sinFit->SetParameter(3,0);

      int fitStatus = gr->Fit(sinFit,"NQ"); //N=no draw, Q=quiet, M=more (better) fit?
      //if the fit is no good tell me
      if (fitStatus != 0) {
	cout << "fitStatus!=0 for entry=" << entry << " surfi=" << surfi << endl;
	//	for (int param=0; param<4; param++) {
	//	  cout << sinFit->GetParameter(param) << " ";
	//	}
	//	cout << endl;
      }
      //if the fit is good, fill the histograms
      else {
	amp = sinFit->GetParameter(0);
	freq = sinFit->GetParameter(1)/(2.*M_PI);
	phase = sinFit->GetParameter(2);
	offset = sinFit->GetParameter(3);
	

	hAmp->Fill(amp);
	hFreq->Fill(freq);
	hPhase->Fill(phase);
	hOffset->Fill(offset);
	
	double absResidSum = 0;
	for (int pt = 0; pt<gr->GetN(); pt++) {
	  absResidSum += pow(gr->GetY()[pt] - sinFit->Eval(gr->GetX()[pt]),2);
	}
	absResidSum /= gr->GetN();
	resid = TMath::Sqrt(absResidSum);
	resid /= chanNoise[surf][lab][rco];
	hResid->Fill(resid);

	fitTree->Fill();
      }
      delete(gr);   
      delete(sinFit);
    }

    delete(usefulRawEvent);

  }//entry

  
  outFile->cd();
  hAmp->Write();
  hFreq->Write();
  hPhase->Write();
  hOffset->Write();
  hResid->Write();
  fitTree->Write();
  return;

}


void noiseRMS(TFile* outFile) {
  /*
    i need to see what the normal noise RMS is without fitting

    the difference between this and the residual RMS is going to be the actual "wrongness" of the fit

   */

  

  int run = 10438; //this is a pretty long run (17k events) where all the amps are off
  
  stringstream name;  
  
  //Events Waveforms
  TChain *rawEventTree = new TChain("eventTree","");
  name.str("");
  name << "/Volumes/ANITA3Data/antarctica14/root/run" << run << "/eventFile" << run << ".root";
  rawEventTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  cout << "rawEventTree Entries: " << rawEventTree->GetEntries() << endl;

  //Event Headers
  TChain *headTree = new TChain("headTree","");
  name.str("");
  name << "/Volumes/ANITA3Data/antarctica14/root/run" << run << "/headFile" << run << ".root";
  headTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  cout << "headTree Entries: " << headTree->GetEntries() << endl;

  //Set Branch Addresses
  RawAnitaEvent *rawEvent = NULL;
  rawEventTree->SetBranchAddress("event",&rawEvent);

  RawAnitaHeader *header = NULL;
  headTree->SetBranchAddress("header",&header);
  
  //Channels with the sine wave in them (All Surfs, 1->12):
  //4,2,4,2,4,2,4,4,2,2,2,2
  int surfs[12] = {1,2,3,4,5,6,7,8,9,10,11,12};
  int chans[12] = {4,2,4,2,4,2,4,4,2, 2, 2, 2};

  TH2D* rmsHist = new TH2D("rmsHist","Channel Noise RMS;Noise RMS (mV);(Surf*8)+(Chip*2)+RCO",
			   96,-0.5,95.5,
			   100,0,10);


  int lenEntries = headTree->GetEntries();
  //  int lenEntries = 1000;
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%100 == 0) {
      cout << "Entry:" << entry << "/" << lenEntries << "\r";
      fflush(stdout); }
    rawEventTree->GetEntry(entry);
    headTree->GetEntry(entry);
    if (rawEvent->eventNumber != header->eventNumber) {
      cout << "Event Number Mismatch! " << rawEvent->eventNumber << "!=" <<  header->eventNumber << endl;
      continue;
    }

    UsefulAnitaEvent *usefulRawEvent = new UsefulAnitaEvent(rawEvent,WaveCalType::kFull,header);
    //I don't want to resample the alfa so I have to turn it's filtering off
    usefulRawEvent->setAlfaFilterFlag(false);
    
    
    for (int surfi=0; surfi<12; surfi++) {
      int surf = surfs[surfi]-1;
      int chan = chans[surfi]-1;
      TGraph *gr = usefulRawEvent->getGraphFromSurfAndChan(surf,chan);

      //from simpleStructs.h -> chanIndex = chan+9*surf
      int chanIndex = chan+9*surf;
      int rco = rawEvent->getRCO(chanIndex);
      int lab = rawEvent->getLabChip(chanIndex);
      
      rmsHist->Fill(surf*8+lab*2+rco,TMath::RMS(gr->GetN(),gr->GetY()));
      delete gr;
    }

    delete usefulRawEvent;

}

  outFile->cd();
  rmsHist->Write();

  return;

}
      

    
int main(int argc, char** argv) {

  string outFileName;
  if (argc == 2) {
    outFileName = argv[1];
    cout << "Using " << outFileName << " as output file" << endl;
  }
  else {
    cout << "Usage: " << argv[0] << " [Output file name]" << endl;
    return -1;
  }

  cout << "Starting Physics!" << endl;
  
  TFile *outFile = TFile::Open(outFileName.c_str(),"recreate");
  
  sineFitting(outFile);
  //  noiseRMS(outFile);

  outFile->Close();
  
  return 1;
  
}
