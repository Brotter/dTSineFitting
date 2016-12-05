{


  /*
    The cluster nicely separates the run into 256 different time ordered histograms, so
    lets look at a single bin over time
  */


  
  //first lets get the same values the cluster starter gets
  int numEntries=540800;
  int entriesPerCluster=numEntries/4;
  int entriesPerCore=numEntries/256;
  int clusterNum=0;


  string clusterDir = "/home/brotter/nfsShared/results/dTSineFitting/";

  //this bin has a weird distribution
  int surf=0;
  int lab=0;
  int rco=0;
  int binX=167;

  //make the time hist
  TH2D *timeHist = new TH2D("timeHist","Time Dependance -- SURF0 LAB0 RCO0 BIN167; time bin; dT(ns)",
			    256,-0.5,255.5,   1001,0,1.0);

  stringstream name;
  for (int cluster=0; cluster<4; cluster++) {
    for (int core=0; core<64; core++) {
      cout << "cluster" << cluster << " core" << core << endl;
      int timeBin = cluster*64 + core;
      int startEv = entriesPerCluster*cluster + core*entriesPerCore;
      name.str("");
      name << clusterDir << "dTOffsetFinder_stEv" << startEv << ".root";
      TFile *inFile = TFile::Open(name.str().c_str());
      name.str("");
      name << "surf" << surf << "_lab" << lab << "_rco" << rco;
      TH2D *inputHist = (TH2D*)inFile->Get(name.str().c_str());
      int dTArrayIndex = surf*8 + lab*2 + rco;
      for (int binY=0; binY<1001;binY++) {
	double binCenterY = inputHist->GetYaxis()->GetBinCenter(binY);
	timeHist->Fill(timeBin,binCenterY,inputHist->GetBinContent(binX,binY));
      }
      inFile->Close();
    }
  }

  timeHist->Draw("colz");
  
  return;


}

