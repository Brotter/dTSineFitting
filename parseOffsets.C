int parseOffsets(bool savePlots=false){


  /*

    Once the c++ code generates the histograms of offsets, it would be nice to be able to easily manipulate them!

    I can probably even use this to generate a "correction table"

   */

  //get the offset file
  TFile *offsetFile = TFile::Open("quick.root");

  //histogram to look at distribution of dT "incorrectness"
  TH1D *hMeanOffsets = new TH1D("meanOffsets","meanOffsets;count;offset(ns)",250,-0.1,0.1);
  
  //saving the values of each mean
  double meanOffsets[96][260];
  
  stringstream name,saveName;
  TCanvas *c1 = new TCanvas("offsetCanvas","offsetCanvas",1200,800);
  for (int surfi=0; surfi<12; surfi++) {
    for (int labi=0; labi<4; labi++) {
      for (int rcoi=0; rcoi<2; rcoi++) {
	name.str("");
	name << "surf" << surfi << "_lab" << labi << "_rco" << rcoi;
	int dTArrayIndex = surfi*8 + labi*2 + rcoi;
	TH2D *storageHist = (TH2D*)offsetFile->Get(name.str().c_str());
	if (savePlots == true) {
	  saveName.str("");
	  saveName << "plots/" << name.str() << ".png";
	  c1->cd();
	  storageHist->Draw("colz");
	  c1->SaveAs(saveName.str().c_str());
	}
	for (int bin=0; bin<260; bin++) {
	  double mean = storageHist->ProjectionY("temp",bin,bin+1)->GetMean();
	  hMeanOffsets->Fill(mean);
	  meanOffsets[dTArrayIndex][bin] = mean;
	}
      }
    }
  }

  //  offsetFile->Close();

  //hMeanOffsets->Draw();
  TFile *outFile = TFile::Open("temp.root","recreate");
  TTree *dTcalTree = new TTree("dTCalTree","dTCalTree");
  dTcalTree->ReadFile("/Users/brotter/anita16/local/share/anitaCalib/justBinByBin.dat",
		      "surf/I:lab/I:rco/I:bin/I:dT/F");
  float dT;
  int surf,lab,rco,bin;
  dTcalTree->SetBranchAddress("surf",&surf);
  dTcalTree->SetBranchAddress("lab",&lab);
  dTcalTree->SetBranchAddress("rco",&rco);
  dTcalTree->SetBranchAddress("bin",&bin);
  dTcalTree->SetBranchAddress("dT",&dT);
    

  TH1D *origVar = new TH1D("origVar","Original Variance;dT(ns);count",200,0,1);
  TH1D *corrVar = new TH1D("corrVar","Corrected Variance;dT(ns);count",200,0,1);
  TH1D *var = new TH1D("var","dT correction Offset;dT(ns);count",200,-0.02,0.02);

  TH2D *corr2D = new TH2D("corr2D","corr2D",200,0,1,200,-0.1,-0.1);

  for (int entry=0; entry<dTcalTree->GetEntries(); entry++) {
    dTcalTree->GetEntry(entry);
    int dTArrayIndex = surf*8 + lab*2 + rco;
    //    cout << meanOffsets[dTArrayIndex][bin] << " " << dT << endl;
    corrVar->Fill(meanOffsets[dTArrayIndex][bin]+dT);
    origVar->Fill(dT);
    var->Fill(meanOffsets[dTArrayIndex][bin]);
    corr2D->Fill(meanOffsets[dTArrayIndex][bin],dT);
  }
  
  c1->Divide(1,2);
  c1->cd(1);
  corrVar->Draw();
  origVar->SetLineColor(kRed);
  origVar->Draw("same");
  c1->cd(2);
  var->Draw();
  
    //corr2D->Draw("colz");

  return 1;
}
