#include "AnitaConventions.h"
#include "sineUtils.h"


void plotMultipleFits(){
  /*
    I need to be able to plot single sine waves to see why the residuals are so weird
   */
  int run = 10105;
  
  stringstream name;  


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

  //Set Branch Addresses
  RawAnitaEvent *rawEvent = NULL;
  rawEventTree->SetBranchAddress("event",&rawEvent);

  RawAnitaHeader *header = NULL;
  headTree->SetBranchAddress("header",&header);


  //get the pedestal correction
  double pedCorrections[12*8*4*260];
  loadPedCorrections(pedCorrections);


  //and I'm only looking at one channel and one lab
  int surf=0;
  int chan=3;
  int usefulIndex = surf*9+chan;
  int lab=0;
  int rco=0;


  //I want to keep the residuals by cap updating
  TH2D *hXResidCap = new TH2D("hXResidCap","X Residuals By Cap;Cap Number; Time Residual(ns)",
			      260,-0.5,259.5,  101,-0.1,0.1);
  hXResidCap->SetStats(0);

  TFile *outFile = TFile::Open("plotMultipleFits.root","recreate");

  int lenEntries = headTree->GetEntries();
  lenEntries = 1000;
  for (int entry=0; entry<10; entry++) {
    rawEventTree->GetEntry(entry);
    headTree->GetEntry(entry);

    UsefulAnitaEvent *usefulRawEvent = new UsefulAnitaEvent(rawEvent,WaveCalType::kFull,header);
    //I don't want to resample the alfa so I have to turn it's filtering off
    usefulRawEvent->setAlfaFilterFlag(false);

    //only interested in the first lab!
    if (usefulRawEvent->getLabChip(usefulIndex) != lab) {
      delete usefulRawEvent;
      continue;
    }
    //only interested when we are looking at RCO 0
    if (usefulRawEvent->getRCO(usefulIndex) != rco) {
      delete usefulRawEvent;
      continue;
    }

  
    //make the canvas
    name.str("");
    name << "cEntry" << entry;
    TCanvas *c1 = new TCanvas(name.str().c_str(),name.str().c_str(),1024,800);

    //get the waveform
    TGraph *rawWaveform = getCorrectedPedsGraph(usefulRawEvent,surf,chan,pedCorrections);

    //only do one RCO phase at a time I guess (so make a new waveform with only that many points
    TGraph *waveform = new TGraph();
    int capNumber = usefulRawEvent->fCapacitorNum[usefulIndex][0];
    for (int pt=0; pt<rawWaveform->GetN(); pt++) {
      int newCapNumber = usefulRawEvent->fCapacitorNum[usefulIndex][pt];
      if (abs(capNumber - newCapNumber) > 1) {
	break; //stop if you find the place it wraps around;
      }
      else {
	waveform->SetPoint(waveform->GetN(),rawWaveform->GetX()[pt],rawWaveform->GetY()[pt]);
	capNumber = newCapNumber;
      }
    }

    cout << waveform->GetN() << endl;

    //fit a sine wave to it
    TF1* sineFit = sineWaveFitter(waveform);
    if (sineFit == NULL) {
      cout << "Entry" << entry << " ... Bad Fit ... Skipping ..." << endl;
      delete usefulRawEvent;
      delete waveform;
      continue;
    }

    double amp = sineFit->GetParameter(0);
    double phase = sineFit->GetParameter(1);
    double freq = 0.4321;
    double offset = 0.0;

    //make a TGraph from that fit
    TGraph *grFit = new TGraph();
    grFit->SetMarkerColor(kRed);
    grFit->SetLineColor(kRed);
    double grX = 0;
    double grY = 0;
    for (int pt=0; pt<1000; pt++) {
      grX = (pt-50)*0.1;
      if (grX-10 > waveform->GetX()[waveform->GetN()-1]){
	break;
      }
      grFit->SetPoint(grFit->GetN(),grX,sineFit->Eval(grX));
    }
     
    //create the graph of Y residuals
    TGraph *grYResid = new TGraph();
    for (int pt=0; pt<waveform->GetN(); pt++) {
      grX = waveform->GetX()[pt];
      grYResid->SetPoint(grYResid->GetN(),grX,waveform->GetY()[pt] - sineFit->Eval(grX));
    }
      

    //clear the canvas and set it up
    c1->Clear();
    c1->Divide(1,3);

    //plot the event with a title
    c1->cd(1);
    name.str("");
    name << "eventNumber: " << header->eventNumber << ";Time(ns);Voltage(mV)";
    waveform->SetTitle(name.str().c_str());
    waveform->Draw("ap");
    grFit->Draw("lSame");
    TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry(waveform,"Calibrated Anita Event","p");
    leg->AddEntry(grFit,"Sine Wave Fit","l");
    leg->Draw("same");
      
    //also draw the "x" residuals and the arrow from one to another
    // as well as the residuals AS A FUNCTION OF CAPACITOR NUMBER
    TGraph *grXResid = new TGraph();
    grXResid->SetTitle("X Residuals;Time(ns);Time Residual(ns)");
    for (int pt=0; pt<waveform->GetN(); pt++) {
      waveform->GetPoint(pt,grX,grY);
      if (abs(grY) > 0.8*abs(amp) || grX<1) { //I can't seem to do stuff below 0 correctly
	continue;
      }
      double xCorrection = findXOffset(amp,freq,phase,offset,grX,grY);
      TArrow *temp = new TArrow(grX,grY,grX+xCorrection,grY,0.005,"|>");
      grXResid->SetPoint(grXResid->GetN(),grX,xCorrection);
      if (abs(xCorrection) < 0.2 ) { //only store them if they aren't crazy to keep the scale good
	hXResidCap->Fill(usefulRawEvent->fCapacitorNum[usefulIndex][pt],xCorrection);
      }
      temp->SetLineColor(kBlue);
      temp->SetFillColor(kBlue);
      //	cout << grX << " " << grX-xCorrection << endl;
      temp->Draw();
    }



    c1->cd(3);
    hXResidCap->Draw("colz");

    
    c1->cd(2);
    //    grResid->SetTitle("Residuals;Time (ns);Voltage(mV)");
    //    grResid->Draw("alp");
    if (grXResid->GetN()==0) {
      cout << "no x resids?" << endl;
      c1->Update();
      sleep(2);
    }
    else {
      grXResid->Draw("alp");
      c1->Update();
    }

    outFile->cd();
    c1->Write();
    delete c1;

    usleep(500000);

    //delete all the crap I made
    delete leg;
    delete grXResid;
    delete grYResid;
    delete grFit;
    delete waveform;
    delete rawWaveform;
    delete usefulRawEvent;
    delete sineFit;
    //    delete grXResidCap;


  }

  return;
  
}
