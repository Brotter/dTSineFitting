{


  /*

    The cluster spits out a bunch of histograms that need to be smooshed together, so this smooshes

    Also there is a TTree that needs to be smooshed

   */


  //first lets get the same values the cluster starter gets
  int numEntries=540800;
  int entriesPerCluster=numEntries/4;
  int entriesPerCore=numEntries/256;
  int clusterNum=0;


  string clusterDir = "/home/brotter/nfsShared/results/dTSineFitting/";


  TFile *outFile = TFile::Open("dTGlobbedClusterResults.root","recreate");


  //output storage histograms (need 96 2D histograms I guess... (bin# vs offset))
  stringstream title,name;
  TH2D *storageHists[96];
  TH2D *diffHists[96];
  for (int surfi=0; surfi<12; surfi++) {
    for (int labi=0; labi<4; labi++) {
      for (int rcoi=0; rcoi<2; rcoi++) {
	name.str("");
	name << "surf" << surfi << "_lab" << labi << "_rco" << rcoi;
	title.str("");
	title << name.str() << ";Storage bin;Corrected dT (ns)";
	int dTArrayIndex = surfi*8 + labi*2 + rcoi;
	storageHists[dTArrayIndex] = new TH2D(name.str().c_str(),title.str().c_str(),
					      260,-0.5,259.5, 1001,0.0,1.0);

	name.str("");
	name << "Correction_surf" << surfi << "_lab" << labi << "_rco" << rcoi;
	title.str("");
	title << name.str() << ";Storage bin;Corrected dT (ns)";
	diffHists[dTArrayIndex] = new TH2D(name.str().c_str(),title.str().c_str(),
					      260,-0.5,259.5, 1001,-0.5,0.5);
      }
    }
  }


  //output Fit tree
  outFile->cd();
  makeFitTree();


  for (int cluster=0; cluster<4; cluster++) {
    for (int core=0; core<64; core++) {
      cout << "cluster" << cluster << " core" << core << endl;
      int startEv = entriesPerCluster*cluster + core*entriesPerCore;
      name.str("");
      name << clusterDir << "dTOffsetFinder_stEv" << startEv << ".root";
      TFile *inFile = TFile::Open(name.str().c_str());

      TTree *inTree = (TTree*)inFile->Get("fitTree");
      for (int entry=0; entry<inTree->GetN(); entry++) {
	inTree->GetEntry(entry);
	fitTree->Fill();
      }

    for (int surf=0; surf<12; surf++) {
      for (int lab=0; lab<4; lab++) {
	for (int rco=0; rco<2; rco++) {
	  name.str("");
	  name << "surf" << surf << "_lab" << lab << "_rco" << rco;
	  TH2D *inputHist = (TH2D*)inFile->Get(name.str().c_str());
	  int dTArrayIndex = surf*8 + lab*2 + rco;
	  for (int binX=0; binX<260; binX++) {
	    for (int binY=0; binY<1001;binY++) {
	      double binCenterX = inputHist->GetXaxis()->GetBinCenter(binX);
	      double binCenterY = inputHist->GetYaxis()->GetBinCenter(binY);
	      storageHists[dTArrayIndex]->Fill(binCenterX,binCenterY,inputHist->GetBinContent(binX,binY));
	    }
	  }
	}
      }
    }
    inFile->Close();
    }
  }


  
    for (int surf=0; surf<12; surf++) {
      for (int lab=0; lab<4; lab++) {
	for (int rco=0; rco<2; rco++) {
	  int dTArrayIndex = surf*8 + lab*2 + rco;
	  outFile->cd();
	  storageHists[dTArrayIndex]->Write();
	}
      }
    }
    outFile->cd();
    fitTree->Write();


    return;
}
  
