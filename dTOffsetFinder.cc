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

  The dTs probably aren't right.

  So now a good thing to try would be to determine the __X__ offset, instead of the Y offset between
  the fit and the measured point.  If I make a histogram of that vs pt I might be able to find a systematic
  offset which can be applied as a correction.

  I wrote a lot of the equations for doing this in the book, and there are comments throughout.


 */

//Run10105 is a longer version with 540k in it!  Maybe this one works.
#define SINE_RUN 10105						



double findOffset(double amp, double freq, double phase, double offset, double xValue, double yValue) {

  //  cout << "params " <<  amp << " " << freq << " " << phase << " " << offset << endl;


  //first lets determine which half-period of the sine wave this is (from x I guess)
  double period = 1./freq;
  double xAngle = fmod((xValue*freq*2*M_PI + phase),(2.*M_PI)); //this should be 0->2pi


  //I have to figure out if it is in quadrant 2 or 3
  int quadrant = 0;
  if (xAngle > M_PI/2. && xAngle <= M_PI) {
    quadrant = 2;
  }
  else if (xAngle > M_PI && xAngle <= (3./2.)*M_PI) {
    quadrant = 3;
  }

  //  cout << "x=" << xValue << "|" << fmod(xValue,period) << " xAngle " << xAngle << " " << quadrant << endl;

  //determine what the angle should be of something with that y
  double temp0 = (yValue - offset)/amp;
  //thats the normalized value, so we can figure out what angle that is;
  double temp1 = asin( temp0 );
  //if the x value is in the "odd" region not returned by asin, need to shift it into that region
  double temp2;
  if (quadrant == 2) {
    temp2 = M_PI - temp1;
  }
  else if (quadrant == 3) {
    temp2 = 3.*M_PI - temp1;
  }
  else {
    temp2 = temp1;
  }
  double temp3 = temp2 - phase;


  //  cout << "y=" << yValue << " temps: " << temp0 << " " << temp1  << " " << temp2 << " " << temp3 << "<-----------------------------" << endl;

  //  double x = asin(( (yValue*pow(-1,perN) )-offset)/amp)-phase+(M_PI*perN)/freq;
  //  cout << x << endl;

  //so then the "correction" should be the difference betwen that and the initial value!
  double xCorrectionAngle = temp2-xAngle;
  //this is circular again, so we need to correct for things that are 2pi separated
  if (xCorrectionAngle > M_PI) {
    xCorrectionAngle -= 2.*M_PI;
  }
  else if (xCorrectionAngle < -M_PI) {
    xCorrectionAngle += 2.*M_PI;
  }
  //and the x is that angle scaled by frequency
  double xCorrection = xCorrectionAngle/(freq*2.*M_PI);
      
  //  cout << "Final: " << xCorrectionAngle << " " << xCorrection << endl;
  
  return xCorrection;

}



int makeOffsetHists(TFile* outFile) {

  int run = SINE_RUN;

  stringstream name;
  //Events Waveforms
  TChain *rawEventTree = new TChain("eventTree","");
  name.str("");
  name << "/Volumes/ANITA3Data/antarctica14/root/run" << run << "/eventFile" << run << ".root";
  rawEventTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  int lenRawEventTree = rawEventTree->GetEntries();
  cout << "rawEventTree Entries: " << lenRawEventTree << endl;

  //Event Headers
  TChain *headTree = new TChain("headTree","");
  name.str("");
  name << "/Volumes/ANITA3Data/antarctica14/root/run" << run << "/headFile" << run << ".root";
  headTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  int lenHeadTree = headTree->GetEntries();
  cout << "headTree Entries: " << lenHeadTree << endl;

  if (lenRawEventTree != lenHeadTree) {
    cout << "the header and event trees have different lenghts so BYE THIS ISN'T CORRECT" << endl;
    return -1;
  }

  //make the indexes
  rawEventTree->BuildIndex("eventNumber");
  headTree->BuildIndex("eventNumber");

  //Set Branch Addresses
  RawAnitaEvent *rawEvent = NULL;
  rawEventTree->SetBranchAddress("event",&rawEvent);

  RawAnitaHeader *header = NULL;
  headTree->SetBranchAddress("header",&header);

  

  TFile *fitFile = TFile::Open("/Volumes/ANITA3Data/bigAnalysisFiles/sineCalibCheck_all10105_AmpPhase.root");
  TTree *fitTree = (TTree*)fitFile->Get("fitTree");
  double amp,freq,phase,residual,offset;
  int eventNumber,chanIndex,rco,lab,surf;
  fitTree->SetBranchAddress("eventNumber",&eventNumber);
  fitTree->SetBranchAddress("amp",&amp);
  //  fitTree->SetBranchAddress("freq",&freq);
  freq = 0.4321;
  fitTree->SetBranchAddress("phase",&phase);
  //  fitTree->SetBranchAddress("offset",&offset);
  offset = 0.0;
  fitTree->SetBranchAddress("residual",&residual);
  fitTree->SetBranchAddress("chanIndex",&chanIndex);
  fitTree->SetBranchAddress("rco",&rco);
  fitTree->SetBranchAddress("lab",&lab);
  fitTree->SetBranchAddress("surf",&surf);
  int lenFitTree = fitTree->GetEntries();
  cout << "fitTree->GetEntres()=" << lenFitTree << endl;

  //Channels with the sine wave in them (All Surfs, 1->12):
  //4,2,4,2,4,2,4,4,2,2,2,2
  int surfs[12] = {1,2,3,4,5,6,7,8,9,10,11,12};
  int chans[12] = {4,2,4,2,4,2,4,4,2, 2, 2, 2};

  
  //output storage histograms (need 96 2D histograms I guess... (bin# vs offset))
  stringstream title;
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

  cout << "made storage stuff" << endl;


  //fitTree will have 8 repeated eventnumbers in a row most likely, so I don't have to make all
  // those UsefulAnitaEvents.
  int prevEventNumber = -1;
  UsefulAnitaEvent *useful = NULL;


  //LOOP THROUGH ALL EVENTS
  lenFitTree = 250000; //or just some of the events
  for (int entry=0; entry<lenFitTree; entry++) {
    if (entry%1000 == 0) {
	cout << entry << "/" << lenFitTree << "\r";
	fflush(stdout);
      }


    fitTree->GetEntry(entry);
    
    //lets make a cut on a sort of "fit goodness" as define by me
    //    if (!(freq<0.44 && freq>0.43 && phase<4 && phase > -4 && amp > 80 && amp < 400
    //	  && residual<100)) {
    //	continue;
    //    }
    

    //When we move from one event number to another
    if (eventNumber != prevEventNumber) {
      rawEventTree->GetEntryWithIndex(eventNumber);
      headTree->GetEntryWithIndex(eventNumber);
      if (useful != 0x0)  delete useful; 
      useful = new UsefulAnitaEvent(rawEvent,WaveCalType::kFull,header);
    }

    prevEventNumber = eventNumber;

    //figure out what channel that should be (index starts at 0 for getGraph)
    int channel = chans[surf]-1;
    
    //index for usefulAnitaEvent->fCapacitorNum
    int usefulIndex = (surf)*9+channel;

    //    cout << "ev#:" << eventNumber << " " << surf << "," << channel << endl;
    TGraph *waveGraph = useful->getGraphFromSurfAndChan(surf,channel);


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
    for (int pt=0; pt<waveGraph->GetN(); pt++) {
      waveGraph->GetPoint(pt,x,y);
      //if we are within 10% of the maximum, set the correction to -999 and I wont use those
      //(because the fit deviations will be dominated by vNoise)
      //also there is some bug where I can't do things with x<0 because I am dumb so lets just fuck that for a sec
      double xCorrected; //need to define it outside the if/else statements
      if (abs(y) > amp*0.9 || x<1) {
	xCorrected = -999;
      }
      else {
	xCorrected = x - findOffset(amp,freq,phase,offset,x,y); 
      }
      int capNum = useful->fCapacitorNum[usefulIndex][pt];
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
      double origdT = waveGraph->GetX()[pt+1] - waveGraph->GetX()[pt];
      diffHists[storageIndex]->Fill(abs(capNum),origdT-dT);
      devVsVolt->Fill(abs(y)/amp,origdT-dT);
    }
    delete waveGraph;
    delete newXArray;

  } //END LOOP THROUGH EVENTS

  outFile->cd();
  for (int i=0; i<96; i++) {
    storageHists[i]->Write();
    delete storageHists[i];
    diffHists[i]->Write();
    delete diffHists[i];
  }
  devVsVolt->Write();

  return 1;
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

  makeOffsetHists(outFile);

  outFile->Close();

  return 1;

}
