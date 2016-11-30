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
      
  //  cout << "Final: " << xCorrectionAngle << " " << xCorrection << endl;
  
  return xCorrection;

}







void plotMultipleFits(){
  /*
    I need to be able to plot single sine waves to see why the residuals are so weird
   */
  int run = 10105;
  
  stringstream name;  


  //get the fit file
  //TFile *fitFile = TFile::Open("/Volumes/ANITA3Data/bigAnalysisFiles/sineCalibCheck_all10105_normalized.root");
  TFile *fitFile = TFile::Open("sineCalibCheck_adc.root");
  TTree *fitTree = (TTree*)fitFile->Get("fitTree");
  double amp,freq,phase,residual,offset;
  int eventNumber,chanIndex,rco,lab,surf;
  fitTree->SetBranchAddress("eventNumber",&eventNumber);
  fitTree->SetBranchAddress("amp",&amp);
  fitTree->SetBranchAddress("freq",&freq);
  fitTree->SetBranchAddress("phase",&phase);
  fitTree->SetBranchAddress("offset",&offset);
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

  
  //make the canvas
  TCanvas *c1 = new TCanvas("canvas","canvas",1024,800);


  //Lets define the sine fit function (I don't know what I should do for the range...)
  TF1 *sinFit = new TF1("sinFit","[0]*sin(x*[1]+[2])+[3]",20.,80.);
    

  for (int entry=0; entry<fitTree->GetEntries(); entry++) {
    fitTree->GetEntry(entry);
    rawEventTree->GetEntryWithIndex(eventNumber);
    headTree->GetEntryWithIndex(eventNumber);

    UsefulAnitaEvent *usefulRawEvent = new UsefulAnitaEvent(rawEvent,WaveCalType::kFull,header);
    //I don't want to resample the alfa so I have to turn it's filtering off
    usefulRawEvent->setAlfaFilterFlag(false);
    
    //get the graph (and the chan)
    int chan = chans[surf]-1;
    TGraph *gr = usefulRawEvent->getGraphFromSurfAndChan(surf,chan);
      
    //create the fit graph
    sinFit->SetParameter(0,amp);
    sinFit->SetParameter(1,freq*2.*M_PI);
    sinFit->SetParameter(2,phase);
    sinFit->SetParameter(3,offset);
    TGraph *grFit = new TGraph();
    grFit->SetMarkerColor(kRed);
    grFit->SetLineColor(kRed);
    double grX = 0;
    double grY = 0;
    for (int pt=0; pt<1000; pt++) {
      grX = (pt-50)*0.1;
      grFit->SetPoint(grFit->GetN(),grX,sinFit->Eval(grX));
    }
     
    //create the graph of X residuals
    TGraph *grResid = new TGraph();
    for (int pt=0; pt<gr->GetN(); pt++) {
      grX = gr->GetX()[pt];
      grResid->SetPoint(grResid->GetN(),grX,gr->GetY()[pt] - sinFit->Eval(grX));
    }
      
    //clear the canvas and set it up
    c1->Clear();
    c1->Divide(1,2);

    //plot the event with a title
    c1->cd(1);
    name.str("");
    name << "surf" << surf+1 << " lab" << lab << " rco" << rco << ";Time(ns);Voltage(mV)";
    gr->SetTitle(name.str().c_str());
    gr->Draw("alp");
    grFit->Draw("lpSame");
    TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry(gr,"Calibrated Anita Event","lp");
    leg->AddEntry(grFit,"Sine Wave Fit","lp");
    leg->Draw("same");
      
    //also draw the "x" residuals and the arrow from one to another
    TGraph *xResiduals = new TGraph();
    xResiduals->SetTitle("X Residuals;Time(ns);Time Residual(ns)");
    for (int pt=0; pt<gr->GetN(); pt++) {
      gr->GetPoint(pt,grX,grY);
      if (abs(grY) > 0.8*amp || grX<1) { //I can't seem to do stuff below 0 correctly
	continue;
      }
      double xCorrection = findOffset(amp,freq,phase,offset,grX,grY);
      TArrow *temp = new TArrow(grX,grY,grX+xCorrection,grY,0.005,"|>");
      xResiduals->SetPoint(xResiduals->GetN(),grX,xCorrection);
      temp->SetLineColor(kBlue);
      temp->SetFillColor(kBlue);
      //	cout << grX << " " << grX-xCorrection << endl;
      temp->Draw();
    }
    
    c1->cd(2);
    //    grResid->SetTitle("Residuals;Time (ns);Voltage(mV)");
    //    grResid->Draw("alp");
    xResiduals->Draw("alp");
    
    c1->Update();
    usleep(500000);

    //delete all the crap I made
    delete xResiduals;
    delete leg;
    delete grResid;
    delete grFit;
    delete gr;
    delete usefulRawEvent;

  }

  return;
  
}
