void linkSinToEvents() {

  /*
    There is a multi-peak distribution in the fit residuals that doesn't correspond to any of the fit
    parameters, so maybe it corresponds to like hit-bus location or something.

    To link the results of the fits to the actual events I can't just use TTree->Draw() because
    the entries aren't identical, so I have to make an index and stuff

    Also added in noise normalization

   */

  TFile *outFile = TFile::Open("linkSinToEvents.root","recreate");


  TFile *fitFile = TFile::Open("/Volumes/ANITA3Data/bigAnalysisFiles/sineCalibCheck_all10105_fixed.root");
  if (fitFile == 0) {
    cout << "can't find the fit file :(" << endl;
    return -1;}
  TTree *fitTree = (TTree*)fitFile->Get("fitTree");
  double amp,freq,phase,offset,resid;
  int chanIndex,rco,lab,surf,chan,eventNumber;
  fitTree->SetBranchAddress("eventNumber",&eventNumber);
  fitTree->SetBranchAddress("amp",&amp);
  fitTree->SetBranchAddress("freq",&freq);
  fitTree->SetBranchAddress("phase",&phase);
  fitTree->SetBranchAddress("offset",&offset);
  fitTree->SetBranchAddress("residual",&resid);
  fitTree->SetBranchAddress("chanIndex",&chanIndex);
  fitTree->SetBranchAddress("rco",&rco);
  fitTree->SetBranchAddress("lab",&lab);
  fitTree->SetBranchAddress("surf",&surf);


  TFile *noiseFile = TFile::Open("noise.root");
  if (noiseFile == 0) {
    cout << "can't find the noise file :(" << endl;
    return -1;
  }
  
  TH2D *allNoiseRMS = (TH2D*)noiseFile->("rmsHist");

  double chanNoise[12][4][2];
  for (int surf=0; surf<12; surf++) {
    for (int lab=0; lab<4; lab++) {
      for (int rco=0; rco<2; rco++) {
	int chanIndex = surf*8 + lab*2 + rco + 1;
	chanNoise[surf][lab][rco] = allNoiseRMS->ProjectionY("temp",chanIndex,chanIndex+1)->GetMean();
      }
    }



  TFile *eventFile = TFile::Open("/Volumes/ANITA3Data/antarctica14/root/run10105/eventFile10105.root");
  if (eventFile == 0) {
    cout << "can't find the event file :(" << endl;
    return -1;}
  
  TTree *eventTree = (TTree*)eventFile->Get("eventTree");
  eventTree->BuildIndex("eventNumber");
  RawAnitaEvent *event = NULL;
  eventTree->SetBranchAddress("event",&event);


  TH2D *hLastHitBus = new TH2D("hLastHitBus","hLastHitBus;residual RMS(mV);last hit bus",
			       100,0,20,260,0.5,260.5);
  TH2D *hWrappedHitBus = new TH2D("hWrappedHitBus","hWrappedHitBus;residual RMS(mV);wrapped hit bus",
			       100,0,20,260,0.5,260.5);
  


  stringstream name;
  stringstream chanNam;
  TH1D *outResid[12][4][2];
  TH2D *hFirstHitBus[12][4][2];
  TH2D *hLastHitBis[12][4][2];
  TH2D *hWrappedHitBus[12][4][2];
  for (int surf=0; surf<12; surf++) {
    for (int lab=0; lab<4; lab++) {
      for (int rco=0; rco<2; rco++) {
	chanName.str("");
	chanName << "Surf" << surf+1 << "_Lab" << lab << "_rco" << rco;
	name.str("");
	name << "Residuals_" << chanName.str();
	outResid[surf][lab][rco] = new TH1D(name.str().c_str(),name.str().c_str(),500,0,10);
	name.str("");
	name << "hFirstHitBus" << chanName.str();
	hFirstHitBus[surf][lab][rco] = new TH2D(name,"hFirstHitBus;residual RMS(mV);first hit bus",
						200,0,20,260,0.5,260.5);
      }
    }
  }


  int lenEntries = fitTree->GetEntries();
  for (int entry=0; entry<lenEntries; entry++) {
 

    if (entry%10000 == 0) {
      cout << 100.*entry / lenEntries << "%\r";
      fflush(stdout);}

    fitTree->GetEntry(entry);

    if (freq<0.44 && freq>0.43 && phase<4 && phase > -4 && resid < 25 ) {
      eventTree->GetEntryWithIndex(eventNumber);
      hFirstHitBus->Fill(resid,event->getFirstHitBus(3));
      hLastHitBus->Fill(resid,event->getLastHitBus(3));
      hWrappedHitBus->Fill(resid,event->getWrappedHitBus(3));

    }
  }

  TCanvas *c1 = new TCanvas("hi","a canvas",1024,800);
  c1->Divide(3,1);
  c1->cd(1);
  hFirstHitBus->Draw("colz");
  c1->cd(2);
  hLastHitBus->Draw("colz");
  c1->cd(3);
  hWrappedHitBus->Draw("colz");

}
  

  
