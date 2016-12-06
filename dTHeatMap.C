{


  /*
    For some crazy stupid reason, drawing a ttree to a TH2D messes up the bins.

    Who knows
  */


  TTree *dTTree = new TTree("dTTree","dTTree");
  dTTree->ReadFile("/Users/brotter/anita16/local/share/anitaCalib/justBinByBin.dat",
		   "surf/I:chip/I:rco/I:sample/I:binSize/D");


  TH2D* hbinSize = new TH2D("hbinSize","ANITA3 Sampling Bin Size;Sample Number;(Surf*8) + (Chip*2) + RCO",
			    259,-0.5,258.5,       96,-0.5,95.5);

  TH2D* hbinDiff = new TH2D("hbinDiff","Bin Deviation from Nominal (%);Sample Number;(Surf*8)+(Chip*2)+RCO",
			    259,-0.5,258.5,       96,-0.5,95.5);

  TProfile2D* hbinSizePerSurf = new TProfile2D("hbinSizePerSurf","ANITA3 Sampling Bin Size;Sample Number;Surf",
					       259,-0.5,258.5,       12,-0.5,11.5);

  int surf,chip,rco,sample;
  double binSize;

  dTTree->SetBranchAddress("surf",&surf);
  dTTree->SetBranchAddress("chip",&chip);
  dTTree->SetBranchAddress("rco",&rco);
  dTTree->SetBranchAddress("sample",&sample);
  dTTree->SetBranchAddress("binSize",&binSize);


  cout << dTTree->GetEntries() << endl;
  for (int entry=0; entry<dTTree->GetEntries(); entry++) {
    dTTree->GetEntry(entry);
    int index = surf*8 + chip*2 + rco;
    hbinSize->Fill(sample,index,binSize);
    hbinDiff->Fill(sample,index,100*(abs(binSize-(1./2.6))*2.6));
    hbinSizePerSurf->Fill(sample,surf,binSize);
  }

  hbinSizePerSurf->Draw("colz");

}
    

