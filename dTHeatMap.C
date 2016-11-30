{


  /*
    For some crazy stupid reason, drawing a ttree to a TH2D messes up the bins.

    Who knows
  */


  TTree *dTTree = new TTree("dTTree","dTTree");
  dTTree->ReadFile("/Users/brotter/anita16/local/share/anitaCalib/justBinByBin.dat",
	       "surf/I:chip/I:rco/I:sample/I:binSize/D");


  TH2D* hbinSize = new TH2D("hbinSize","Bin Size Larger than 500ps;Sample Number;(Surf*8) + (Chip*2) + RCO",
			    259,-0.5,258.5,       96,-0.5,95.5);

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
    if ( binSize < 0.5) {
      continue; 
    }
    int index = surf*8 + chip*2 + rco;
    hbinSize->Fill(sample,index);
  }

  hbinSize->Draw("col");

}
    

