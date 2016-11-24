{

  /*

    The residual RMS should be normalized by the digitizer noise (which is stored in a different file)

   */


  TFile *outFile = TFile::Open("normalizedFitRMS.root","recreate");

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
  
  TH2D *allNoiseRMS = (TH2D*)noiseFile->Get("rmsHist");

  TH2D *normalHist = new TH2D("normResid","Normalized Fit Residuals;surf*8 + lab*2 + rco; residual RMS",
			      96,-0.5,95.5, 250,0,10);

  stringstream name;
  stringstream cut;
  for (int surf=0; surf<12; surf++) {
    for (int lab=0; lab<4; lab++) {
      for (int rco=0; rco<2; rco++) {
	int chanIndex = surf*8 + lab*2 + rco + 1;
	cout << chanIndex << endl;
	fflush(stdout);
	double noiseRMS = allNoiseRMS->ProjectionY("temp",chanIndex,chanIndex+1)->GetMean();
	name.str("");
	name << "Surf" << surf+1 << "_Lab" << lab << "_rco" << rco;
	currHist = new TH1D(name.str().c_str(),name.str().c_str(),500,0,10);
	name.str("");
	name << "residual/" << noiseRMS << " >> " << name.str();
	cut.str("");
	cut << "freq<0.44 && freq>0.43 && phase<4 && phase > -4 && amp > 80 && amp < 140 && residual<100 ";
	cut << "&& surf == " << surf << " && lab == " << lab << " && rco == " << rco;
	fitTree->Draw(name.str().c_str(),cut.str().c_str());
	outFile->cd();
	currHist->Write();
	delete currHist;
      }
    }
  }

  outFile->Close();

}
  
  
