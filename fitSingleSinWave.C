#include "AnitaConventions.h"



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
      
  cout << "Final: " << xCorrectionAngle << " " << xCorrection << endl;
  
  return xCorrection;

}







void fitSingleSinWave(int eventNumber=-1){
  /*
    I need to be able to plot single sine waves to see why the residuals are so weird
   */
  int run = 10105;
  
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

  //I am calling event numbers so I need to build some indexes
  rawEventTree->BuildIndex("eventNumber");
  headTree->BuildIndex("eventNumber");

  //Set Branch Addresses
  RawAnitaEvent *rawEvent = NULL;
  rawEventTree->SetBranchAddress("event",&rawEvent);

  RawAnitaHeader *header = NULL;
  headTree->SetBranchAddress("header",&header);


  if (eventNumber==-1) {
    headTree->GetEntry(0);
    eventNumber = header->eventNumber;
  }
  rawEventTree->GetEntryWithIndex(eventNumber);
  headTree->GetEntryWithIndex(eventNumber);

  UsefulAnitaEvent *usefulRawEvent = new UsefulAnitaEvent(rawEvent,WaveCalType::kFull,header);
  //I don't want to resample the alfa so I have to turn it's filtering off
  usefulRawEvent->setAlfaFilterFlag(false);
    
  int surf = 0;
  int chan = 3;
  TGraph *gr = usefulRawEvent->getGraphFromSurfAndChan(surf,chan);

  //Lets define the sine fit function (I don't know what I should do for the range...)
  TF1 *sinFit = new TF1("sinFit","[0]*sin(x*[1]+[2])+[3]",20.,80.);
  sinFit->SetParName(0,"Amplitude");
  sinFit->SetParameter(0,200.);

  sinFit->SetParName(1,"FrequencyX2Pi");
  sinFit->SetParameter(1,2.*M_PI*0.4321); //in GHz since x is in ns
  
  sinFit->SetParName(2,"Phase");
  sinFit->SetParameter(2,0);
  
  sinFit->SetParName(3,"Offset");
  sinFit->SetParameter(3,0);
  
  int fitStatus = gr->Fit(sinFit,"NQ"); //N=no draw, Q=quiet, M=more (better) fit?
  //if the fit is no good tell me
  if (fitStatus != 0) {
    cout << "THIS ONE DOESN'T WORK" << endl;
	//	for (int param=0; param<4; param++) {
	//	  cout << sinFit->GetParameter(param) << " ";
	//	}
	//	cout << endl;
      }
  //if the fit is good, fill the histograms
  else {
    double amp = sinFit->GetParameter(0);
    double freq = sinFit->GetParameter(1)/(2.*M_PI);
    double phase = sinFit->GetParameter(2);
    double offset = sinFit->GetParameter(3);
    
    double absResidSum = 0;
    for (int pt = 0; pt<gr->GetN(); pt++) {
      absResidSum += pow(gr->GetY()[pt] - sinFit->Eval(gr->GetX()[pt]),2);      
    }
    absResidSum /= gr->GetN();
    double resid = TMath::Sqrt(absResidSum);

    cout << "eventNumber " << eventNumber << endl;
    cout << "amp | freq | phase | offset | resid" << endl;
    cout << amp << " | " << freq << " | " << phase << " | " << offset << " | " << resid << endl;


    TGraph *grFit = new TGraph();
    grFit->SetMarkerColor(kRed);
    grFit->SetLineColor(kRed);
    double grX = 0;
    double grY = 0;
    for (int pt=0; pt<1000; pt++) {
      grX = pt*0.1;
      grFit->SetPoint(grFit->GetN(),grX,sinFit->Eval(grX));
    }
    
    TGraph *grResid = new TGraph();
    for (int pt=0; pt<gr->GetN(); pt++) {
      grX = gr->GetX()[pt];
      grResid->SetPoint(grResid->GetN(),grX,gr->GetY()[pt] - sinFit->Eval(grX));
    }




    TCanvas *c1 = new TCanvas("canvas","canvas",1024,800);
    c1->Divide(1,2);

    c1->cd(1);
    gr->SetTitle("Sine Wave Fitting");
    gr->Draw("alp");
    grFit->Draw("lpSame");
    TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry(gr,"Calibrated Anita Event","lp");
    leg->AddEntry(grFit,"Sine Wave Fit","lp");
    leg->Draw("same");

    //also draw the "x" residuals
    for (int pt=0; pt<gr->GetN(); pt++) {
      gr->GetPoint(pt,grX,grY);
      if (abs(grY) > 0.8*amp) {
	continue;
      }

      double xCorrection = findOffset(amp,freq,phase,offset,grX,grY);
      TArrow *temp = new TArrow(grX,grY,grX+xCorrection,grY,0.005,"|>");
      temp->SetLineColor(kBlue);
      temp->SetFillColor(kBlue);
      cout << grX << " " << grX-xCorrection << endl;
      temp->Draw();
    }




    c1->cd(2);
    grResid->Draw("alp");
  }

}
