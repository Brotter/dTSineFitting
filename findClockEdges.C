#include "AnitaConventions.h"

void findClockEdges(int startRun, int stopRun,char* outFilename){

  /*

    It is useful to know where the sync clock edges are, or if they are evenly distributed

   */



  //make storage histograms (one per surf I guess, so just a TH2D!)
  TH2D * hDownEdgesPt  = new TH2D("hDownEdgesPt","Sync Clock Location; Surf; Sync Clock Upgoing Edge Ordered Bin",
				 12,-0.5,11.5,   260,-0.5,259.5);
  TH2D * hDownEdgesCap = new TH2D("hDownEdgesCap","Sync Clock Location; Surf; Sync Clock Upgoing Edge Capacitor",
				 12,-0.5,11.5,   260,-0.5,259.5);


  
  //import the sine wave run
  //  const int run = 10105;
  string dataDir = getenv("ANITA3_DATA");

  stringstream name;
  TChain *rawEventTree = new TChain("eventTree","");
  TChain *headTree = new TChain("headTree","");
  for (int run=startRun; run<stopRun; run++) {
    if (run==11) continue;  //this run is broken;

    //Events Waveforms       
    name.str("");
    name << dataDir << "/run" << run << "/eventFile" << run << ".root";
    rawEventTree->Add(name.str().c_str());
    cout << "Adding: " << name.str() << endl;
    int lenRawEventTree = rawEventTree->GetEntries();
    cout << "rawEventTree Entries: " << lenRawEventTree << endl;
  
    //Event Headers
    name.str("");
    name << dataDir << "/run" << run << "/headFile" << run << ".root";
    headTree->Add(name.str().c_str());
    cout << "Adding: " << name.str() << endl;
    int lenHeadTree = headTree->GetEntries();
    cout << "headTree Entries: " << lenHeadTree << endl;
    
    if (lenRawEventTree != lenHeadTree) {
      cout << "the header and event trees have different lenghts so BYE THIS ISN'T CORRECT" << endl;
      return;
    }
  }

  //Set Branch Addresses
  RawAnitaEvent *rawEvent = NULL;
  rawEventTree->SetBranchAddress("event",&rawEvent);
  RawAnitaHeader *header = NULL;
  headTree->SetBranchAddress("header",&header);


  
  int lenEntries = headTree->GetEntries();;
  if (lenEntries == 0) {
    cout << "There aren't any events, I quit goodbye!" << endl;
    return;
  }

  //  lenEntries = 50000;
  //Loop Through Entries
  for (int entry=0; entry<lenEntries; entry++) {
    headTree->GetEntry(entry);
    rawEventTree->GetEntry(entry);

    if (entry%100 == 0) {
      cout << entry << "/" << lenEntries << "\r"; 
      fflush(stdout); }

    UsefulAnitaEvent *useful = new UsefulAnitaEvent(rawEvent,WaveCalType::kFull,header);

    for (int surf=0; surf<12; surf++) {
      TGraph *gSyncClock = useful->getGraphFromSurfAndChan(surf,8);
      double x,y;
      double prevY = gSyncClock->GetY()[0];
      for (int pt=0; pt<gSyncClock->GetN(); pt++) {
	gSyncClock->GetPoint(pt,x,y);
	if ( (prevY>0) && (y<0) ) { //downgoing edge
	  hDownEdgesPt->Fill(surf,pt);
	  hDownEdgesCap->Fill(surf,useful->fCapacitorNum[surf*9+8][pt]);
	}
	prevY = y;
      }
      delete gSyncClock;
    }
    delete useful;
  }


  hDownEdgesCap->Draw("colz");


  TFile *outFile = TFile::Open(outFilename,"recreate");
  hDownEdgesPt->Write();
  hDownEdgesCap->Write();
  outFile->Close();


  return;

}
