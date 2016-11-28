int parseOffsets(bool savePlots=false){


  /*

    Once the c++ code generates the histograms of offsets, it would be nice to be able to easily manipulate them!

    I can probably even use this to generate a "correction table"

    Also a bunch of plots, plots are great

   */

  //get the offset file
  TFile *offsetFile = TFile::Open("dTOffsetFinder_3.root");

  //histogram to look at distribution of dT "incorrectness"
  TH1D *hMeanOffsets = new TH1D("meanOffsets","meanOffsets;count;offset(ns)",200,0,1);
  
  //saving the values of each mean
  double corrdTArray[96][260];
  double stdDevArray[96][260];

  stringstream name,saveName;
  TCanvas *c1 = new TCanvas("offsetCanvas","offsetCanvas",1200,800);
  c1->SetLogz();
  for (int surfi=0; surfi<12; surfi++) {
    for (int labi=0; labi<4; labi++) {
      for (int rcoi=0; rcoi<2; rcoi++) {
	name.str("");
	name << "surf" << surfi << "_lab" << labi << "_rco" << rcoi;
	int dTArrayIndex = surfi*8 + labi*2 + rcoi;
	TH2D *storageHist = (TH2D*)offsetFile->Get(name.str().c_str());
	if (0) {
	  saveName.str("");
	  saveName << "plots/" << name.str() << ".png";
	  c1->cd();
	  storageHist->Draw("colz");
	  c1->SaveAs(saveName.str().c_str());
	}
	for (int bin=1; bin<260; bin++) { //bins start at 1
	  double mean = storageHist->ProjectionY("temp",bin,bin)->GetMean();
	  stdDevArray[dTArrayIndex][bin] = storageHist->ProjectionY("temp",bin,bin)->GetStdDev();
	  hMeanOffsets->Fill(mean);
	  corrdTArray[dTArrayIndex][bin] = mean;
	}
      }
    }
  }

  if (savePlots == true) {
    TH2D *storageHist = (TH2D*)offsetFile->Get("surf11_lab3_rco0");
    for (int bin=1;bin<260;bin++) {
      saveName.str("");
      saveName << "plots/bin" << bin << ".png";
      name.str("");
      name << "Bin " << bin << ";Fit offset(ns);count";
      TH1D* tempProj = storageHist->ProjectionY("temp",bin,bin+1);
      c1->Clear();
      c1->cd();
      tempProj->SetTitle(name.str().c_str());
      tempProj->Draw();
      c1->SaveAs(saveName.str().c_str());
      delete tempProj;
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
  TH1D *correction = new TH1D("correction","correction;dTorig-dTcorr(ns);count",200,-0.1,0.1);
  TH1D *hSigma = new TH1D("hSigma","Standard Deviation; Standard Deviation(ns);count",200,0,0.5);

  TH2D *corr2D = new TH2D("corr2D","corr2D",200,0,1,200,-0.1,-0.1);

  for (int entry=0; entry<dTcalTree->GetEntries(); entry++) {
    dTcalTree->GetEntry(entry);
    int dTArrayIndex = surf*8 + lab*2 + rco;
    if (dT==0) continue; //these are the "last dT" and shouldn't be included
    //    cout << meanOffsets[dTArrayIndex][bin] << " " << dT << endl;
    corrVar->Fill(corrdTArray[dTArrayIndex][bin]);
    origVar->Fill(dT);
    corr2D->Fill(corrdTArray[dTArrayIndex][bin],dT);
    hSigma->Fill(stdDevArray[dTArrayIndex][bin]);
  }
  

  c1->Clear();
  c1->Divide(2,1);
  c1->cd(1);
  corrVar->Draw();
  origVar->SetLineColor(kRed);
  origVar->Draw("sames");
  TLegend *leg = new TLegend(0.1,0.8,0.3,0.9);
  leg->AddEntry(origVar,"Original dT Distribution","l");
  leg->AddEntry(corrVar,"Corrected dT Distribution","l");
  leg->Draw("same");
  c1->cd(2);
  hSigma->Draw();


    //corr2D->Draw("colz");

  return 1;
}
