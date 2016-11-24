void normalizeResidToNoise(){

  /*

    The residual RMS should be normalized by the digitizer noise (which is stored in a different file)

    Basically this just plots it now since I integrated it into the actual sine wave fitting

   */

  //open some default out file
  TFile *outFile = TFile::Open("normalizedFitRMS.root","recreate");


  //open the fit file for all the sine wave fits
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


  //open the noise file with the noise rmses in it
  TFile *noiseFile = TFile::Open("noise.root");
  if (noiseFile == 0) {
    cout << "can't find the noise file :(" << endl;
    return -1;
  }
  
  TH2D *allNoiseRMS = (TH2D*)noiseFile->Get("rmsHist");
  allNoiseRMS->SetTitle("Digitizer Noise RMS; ;RMS (mV)");

  TH2D *fitResid = new TH2D("fitResid","Absolute Fit Residuals;(Surf*8) +(Lab*2) + Rco; Residual RMS (mV)",
			      96,-0.5,95.5, 250,0,25);


  TCanvas *cRMS2D = new TCanvas("cRMS2D","cRMS2D",1000,800);
  cRMS2D->Divide(1,2);
  cRMS2D->cd(1);
  allNoiseRMS->Draw("colz");
  cRMS2D->cd(2);
  fitTree->Draw("residual:surf*8 + lab*2 + rco >> fitResid","residual<20","colz");
  

  TH1D *rmsRatio = new TH1D("rmsRatio","Ratio Between Noise and Fit Residual RMS;ratio;count",250,0,10);
  for (int surfi=0; surfi<12; surfi++) {
    for (int labi=0; labi<4; labi++) {
      for (int rcoi=0; rcoi<2; rcoi++) {
	int index = surfi*8 + labi*2 + rcoi;
	double meanNoise = allNoiseRMS->ProjectionY("temp",index,index+1)->GetMean();
	double meanFit   = fitResid->ProjectionY("temp",index,index+1)->GetMean();
	rmsRatio->Fill(meanFit/meanNoise);
      }
    }
  }

  TCanvas *cComp = new TCanvas("cComp","cComp",1000,800);
  rmsRatio->Draw();



  outFile->Close();

  return;

}
  
  
