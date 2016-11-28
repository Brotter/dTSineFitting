void comparePhases() {

  /*

    The calibration attempts to align the waveforms, so the phases should be related
    is that true?

   */


  //get the fit file
  TFile *fitFile = TFile::Open("/Volumes/ANITA3Data/bigAnalysisFiles/sineCalibCheck_all10105_normalized.root");
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



  //make a histogram for the comparison

  TH2D *hPhaseComp = new TH2D("hPhaseComp","Phase difference between N and Surf0;Surf;Phase Offset (ns)",
			      11,0.5,11.5, 401, -10.0,10.0);

  double surfZeroPhase = -999;
  for (int entry=0; entry<fitTree->GetEntries(); entry++) {
    fitTree->GetEntry(entry);
    if (surf == 0) {
      surfZeroPhase = phase;
    }
    else {
      double diff = surfZeroPhase-phase;
      if (diff < - M_PI/2.) diff += M_PI;
      if (diff >   M_PI/2.) diff -= M_PI;
      //period = 1/2.6GHz in 2pi radians
      hPhaseComp->Fill(surf,diff*2.6);
    }
  }


  hPhaseComp->Draw("colz");

}
