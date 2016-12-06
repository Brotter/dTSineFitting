#ifndef SINEUTILS
#define SINEUTILS

#include <math.h>
#include <iostream>

//#include "UsefulAnitaEvent.h"

using namespace std;



int storageIndex(int surf, int chan, int lab) {
  return surf*32 + chan*4 + lab;
}

int pedIndex(int surf, int chan, int lab, int sample) {
  if (sample==0) cout << "pedIndex(): WARNING!  You're indexing me wrong, start at 1" << endl;
  return (surf*8*4*259) + (chan*4*259) + (lab*259) + sample;
}


//and the ped corrections
double pedCorrections[12*8*4*260];
void loadPedCorrections() {

  ifstream inFile("pedCorrections.txt");

  int surf,chan,lab,sample;
  double pedCorr;
  while (inFile >> surf >> chan >> lab >> sample >> pedCorr) {
    pedCorrections[pedIndex(surf,chan,lab,sample)] = pedCorr;
  }

  return;

}


TGraph *getCorrectedPedsGraph(UsefulAnitaEvent *useful,int surf,int chan,double* pedCorrections) {
  TGraph *gr = useful->getGraphFromSurfAndChan(surf,chan);
  TGraph *grCorr = new TGraph();
  int usefulIndex = surf*9 + chan;
  int lab = useful->getLabChip(usefulIndex);
  double x,y;
  for (int pt=0; pt<gr->GetN(); pt++) {
    gr->GetPoint(pt,x,y);
    double pedCorr = pedCorrections[pedIndex(surf,chan,lab,useful->fCapacitorNum[usefulIndex][pt])];
    grCorr->SetPoint(pt,x,y-pedCorr);
  }

  delete gr;
  return grCorr;
}


double findXOffset(double amp, double freq, double phase, double offset, double xValue, double yValue) {

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




TF1* sineWaveFitter(TGraph *graphToFit) {
  
  //Lets define the sine fit function (I don't know what I should do for the range...)
  //Since the phase is the biggest generator of error, and the freq is very constant, lets try
  /// hard coding those?
  double xMin = graphToFit->GetX()[0];
  double xMax = graphToFit->GetX()[graphToFit->GetN()-1];
  TF1 *sinFit = new TF1("sinFit","[0]*sin(x*0.4321*2.*3.14159+[1])",xMin,xMax);
  sinFit->SetParName(0,"Amplitude");
  sinFit->SetParameter(0,200.);
  
  sinFit->SetParName(1,"Phase");
  sinFit->SetParameter(1,0);

  int fitStatus = graphToFit->Fit(sinFit,"NQ"); //N=no draw, Q=quiet, M=more (better) fit?
  if (fitStatus != 0) {   //if the fit is no good tell me about it
    cout << "fitStatus != 0" << endl;
    return NULL;
  }

  return sinFit;
  
}


//also I guess I'll make the whole sine fit tree stuff too
TTree *fitTree = NULL;
double amp,freq,phase,offset,resid;
int chanIndex,rco,lab,surf,chan,eventNumber;

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

void setFitTreeBranches(TTree *inTree){
  inTree->SetBranchAddress("eventNumber",&eventNumber);
  inTree->SetBranchAddress("amp",&amp);
  inTree->SetBranchAddress("phase",&phase);
  inTree->SetBranchAddress("residual",&resid);
  inTree->SetBranchAddress("chanIndex",&chanIndex);
  inTree->SetBranchAddress("rco",&rco);
  inTree->SetBranchAddress("lab",&lab);
  inTree->SetBranchAddress("surf",&surf);
  inTree->SetBranchAddress("chan",&chan);
}


#endif

