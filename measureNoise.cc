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
  
  Ben Rotter - December 2016 - University of Hawaii at Manoa



  I am decoupling the code used to determine the digitizer noise of each bin into this code


 */







void fillStorageHists(WaveCalType::WaveCalType_t waveCalType, TH2D** storageHists, int run=-1) {
  /*
    i need to see what the normal noise RMS is without fitting

    the difference between this and the residual RMS is going to be the actual "wrongness" of the fit

   */

  
  //default run
  if (run==-1) {
      run = 10438; //this is a pretty long run (17k events) where all the amps are off
    }  

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

  int lenEntries = headTree->GetEntries();
  //  int lenEntries = 1000;
  for (int entry=0; entry<lenEntries; entry++) {
    if (entry%100 == 0) {
      cout << "Entry:" << entry << "/" << lenEntries << "\r";
      fflush(stdout); }
    rawEventTree->GetEntry(entry);
    headTree->GetEntry(entry);

    //some runs are weird, so just check to see if there is a mismatch
    if (rawEvent->eventNumber != header->eventNumber) {
      cout << "Event Number Mismatch! " << rawEvent->eventNumber << "!=" <<  header->eventNumber << endl;
      continue;
    }

    //make the useful event
    UsefulAnitaEvent *usefulRawEvent = new UsefulAnitaEvent(rawEvent,waveCalType,header);
    //I don't want to resample the alfa so I have to turn it's filtering off
    usefulRawEvent->setAlfaFilterFlag(false);
    
    
    for (int surf=0; surf<NUM_SURF; surf++) {
      for (int chan=0; chan<CHANNELS_PER_SURF; chan++) {

	TGraph *gr = usefulRawEvent->getGraphFromSurfAndChan(surf,chan);

	//from simpleStructs.h -> chanIndex = chan+9*surf
	int chanIndex = surf*9 + chan;
	int rco = rawEvent->getRCO(chanIndex);
	int lab = rawEvent->getLabChip(chanIndex);
      
	//I am doing a different indexing system because I'm special (and I don't care about the clock)
	int index = surf*8 + chan;
	storageHists[index]->Fill(surf*8+lab*2+rco,TMath::RMS(gr->GetN(),gr->GetY()));
      delete gr;
      }
    }
    
    delete usefulRawEvent;
      
  }
  cout << endl;

  delete headTree;
  delete rawEventTree;

  cout << "Filled storage hists" << endl;
  return;

}

void makeStorageHists(WaveCalType::WaveCalType_t waveCalType,TH2D** storageHists) {
  //Lets separate this into its own function!
  // I THINK I can call it like this...

  
  //lets make some storage histograms?  I have to divide them by wavecaltype later I guess....
  //I'll put that in their names

  //this is going to be different based on which wavecaltype, I'm also going to do the full range because why not
  //I'll always do 4096 bins tho
  double voltsMin,voltsMax;
  if (waveCalType == WaveCalType::kOnlyTiming) {
    //adc counts, 12 bits 2**12=4096
    voltsMin = 0;
    voltsMax = 4096; }
  if (waveCalType == WaveCalType::kFull) {
    //mV, I think we have like a -1V to 1V range?
    voltsMin = -1.0;
    voltsMax =  1.0;
  }
   
  //make the array to store them

  stringstream title,name;
  for (int surf=0; surf<NUM_SURF; surf++) {
    for (int chan=0; chan<CHANNELS_PER_SURF ;chan++) {
      title.str("");
      title << "Channel Noise RMS --";
      title << "waveCalType::" << WaveCalType::calTypeAsString(waveCalType);
      title << " Surf" << surf << "Chan" << chan;
      title << "; SCA Bin Number; Noise Value (mV/adc)";
      name.str("");
      name << "caltype" << WaveCalType::calTypeAsString(waveCalType);
      name << "_surf" << surf << "_chan" << chan;

      int index = surf*8 + chan;
      storageHists[index] = new TH2D(name.str().c_str(),title.str().c_str(),
				     260,-0.5,259.5,4096,voltsMin,voltsMax);
    }
  }
			       
  cout << "Created storage hists" << endl;

  return;
  
}


void deleteStorageHists(TH2D** storageHists) {
  for (int surf=0; surf<NUM_SURF; surf++) {
    for (int chan=0; chan<CHANNELS_PER_SURF ;chan++) {
      int index=surf*8+chan;
      delete storageHists[index];
    }
  }
  cout << "Deleted storage hists" << endl;
  return;
}
      



void saveStorageHists(TFile *outFile, TH2D** storageHists){

  outFile->cd();
  for (int surf=0; surf<NUM_SURF; surf++) {
    for (int chan=0; chan<CHANNELS_PER_SURF ;chan++) {
      int index=surf*8+chan;
      storageHists[index]->Write();
    }
  }
  cout << "Saved storage hists" << endl;
  return;
}
  





int main(int argc, char** argv) {

  cout << "Lets do some physics!" << endl;

  
  string outFileName;
  if (argc == 2) {
    outFileName = argv[1];
    cout << "Using " << outFileName << " as output file" << endl;
  }
  else {
    cout << "Usage: " << argv[0] << " [Output file name]" << endl;
    return -1;
  }

  TFile *outFile = TFile::Open(outFileName.c_str(),"recreate");


  WaveCalType::WaveCalType_t waveCalType;

  TH2D* storageHists[96];

  cout << "Doing kFull" << endl;
  waveCalType = WaveCalType::kFull;
  makeStorageHists(waveCalType,storageHists);
  fillStorageHists(waveCalType,storageHists);
  saveStorageHists(outFile,storageHists);
  deleteStorageHists(storageHists);
  
  cout << "Doing kOnlyTiming" << endl;
  waveCalType = WaveCalType::kOnlyTiming;
  makeStorageHists(waveCalType,storageHists);
  fillStorageHists(waveCalType,storageHists);
  saveStorageHists(outFile,storageHists);
  deleteStorageHists(storageHists);


  outFile->Close();
  
 

  return 1;
}
