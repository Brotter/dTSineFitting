void elipsePlottingPlot() {


  TTree *myTree = new TTree("myTree","myTree");
  myTree->ReadFile("elipsePlotting_new.txt",
		   "index/I:surf/I:lab/I:rco/I:bin/I:centerX/F:centerY/F:major/F:minor/F:rotation/F");
  float major,minor,rotation;
  int surf,lab,bin,rco;
  myTree->SetBranchAddress("major",&major);
  myTree->SetBranchAddress("minor",&minor);
  myTree->SetBranchAddress("rotation",&rotation);
  myTree->SetBranchAddress("surf",&surf);
  myTree->SetBranchAddress("lab",&lab);
  myTree->SetBranchAddress("bin",&bin);
  myTree->SetBranchAddress("rco",&rco);


  TTree *benTree = new TTree("benTree","benTree");
  benTree->ReadFile("/Users/brotter/anita16/local/share/anitaCalib/justBinByBin.dat",
		   "surf/I:lab/I:rco/I:bin/I:width/F");
  float width;
  benTree->SetBranchAddress("width",&width);


  TH2D *diff = new TH2D("diff","Ellipse Method vs Nominal; Nominal dT (ns); Ellipse dT (ns)",
			201,0,1,   201,0,1);

  
  TProfile2D *newVolts = new TProfile2D("newVolts", "Voltage Calibration (Ellipse);surf; lab",
					12,-0.5,11.5,  4,-0.5,3.5);

  for (int entry=1; entry<24960; entry++) {
    benTree->GetEntry(entry-1);
    myTree->GetEntry(entry);
    double eWidth;
    //weird rotations
    if ( ( (rotation > 1) && (rotation < 2) ) || ( (rotation > 4) && (rotation < 5) ) ) {
      eWidth = TMath::ATan(major/minor) / (3.14159 * 0.4321);}
    else {
      eWidth = TMath::ATan(minor/major) / (3.14159 * 0.4321); }
    //    cout << entry << " " << width << " " << eWidth << endl;

    
    //    if (abs(width - eWidth) > 0.1) {
    if (width < 0.06) {
      cout << surf << " " << lab << " " << rco << " " << bin;
      cout << " | " << width << " " << eWidth << " " << rotation << endl;
    }
    else     diff->Fill(width,eWidth);

    double eAmp = TMath::Sqrt(pow(major,2) + pow(minor,2))/2.;
    newVolts->Fill(surf,lab,eAmp/10.);
    
  }

 
  diff->Draw("colz");


  TTree *voltTree = new TTree("harmTree","harmTree");
  voltTree->ReadFile("/Users/brotter/anita16/local/share/anitaCalib/simpleVoltageCalibrationHarm.txt",
		   "surf/I:chan/I:lab/I:volts/D");
  int surf2,lab2,chan2;
  double volts;
  voltTree->SetBranchAddress("surf",&surf2);
  voltTree->SetBranchAddress("lab",&lab2);
  voltTree->SetBranchAddress("chan",&chan2);
  voltTree->SetBranchAddress("volts",&volts);



  int chans[12] = {4,2,4,2,4,2,4,4,2, 2, 2, 2};

  TH2D *oldVolts = new TH2D("oldVolts", "Voltage Calibration (Harm's);surf; lab",
			    12,-0.5,11.5,  4,-0.5,3.5);

  for (int entry=0; entry<voltTree->GetEntries(); entry++) {
    voltTree->GetEntry(entry);
    if ( chan2 == chans[surf2]-1 ) oldVolts->Fill(surf2,lab2,volts);
  }

  /*
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  newVolts->Draw("colz");
  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  oldVolts->Draw("colz");
  */

  return;

}




void elipsePlottingGifMaker(){


  stringstream name;
  for (int i=0; i<100; i++) {
    cout << i << endl;

    name.str("");
    name << "elipsePlotting_s5c1l3n" << i << ".root";
    TFile *openFile = TFile::Open(name.str().c_str());
    TH2D *hEllipse = (TH2D*)openFile->Get("hEllipse");

    TCanvas *c1 = new TCanvas("c1","c1",1000,800);
    hEllipse->Draw("colz");

    if (i==0) c1->SaveAs("elipsePlottingPlot_s5c1l3.gif");
    else c1->SaveAs("elipsePlottingPlot_s5c1l3.gif+");

    openFile->Close();
    delete c1;
  }


  return;

}
