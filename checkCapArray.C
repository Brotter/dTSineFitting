#include "AnitaConventions.h"

void checkCapArray(){



  /* The capacitor array (which cap for each time ordered sampled) is supposed to be saved in UsefulAnitaEvent
     Does that actually work?  Lets find out!
  */


  int run = 10105;
  
  stringstream name;  
  
  //Events Waveforms
  TChain *rawEventTree = new TChain("eventTree","");
  name.str("");
  name << "/Volumes/ANITA3Data/antarctica14/root/run" << run << "/eventFile" << run << ".root";
  rawEventTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  cout << "rawEventTree Entries: " << rawEventTree->GetEntries() << endl;

  //Event Headers
  TChain *headTree = new TChain("headTree","");
  name.str("");
  name << "/Volumes/ANITA3Data/antarctica14/root/run" << run << "/headFile" << run << ".root";
  headTree->Add(name.str().c_str());
  cout << "Adding: " << name.str() << endl;
  cout << "headTree Entries: " << headTree->GetEntries() << endl;

  //I am calling event numbers so I need to build some indexes

  //Set Branch Addresses
  RawAnitaEvent *rawEvent = NULL;
  rawEventTree->SetBranchAddress("event",&rawEvent);

  RawAnitaHeader *header = NULL;
  headTree->SetBranchAddress("header",&header);

  rawEventTree->GetEntry(0);
  headTree->GetEntry(0);

  UsefulAnitaEvent *usefulRawEvent = new UsefulAnitaEvent(rawEvent,WaveCalType::kFull,header);
  //I don't want to resample the alfa so I have to turn it's filtering off
  usefulRawEvent->setAlfaFilterFlag(false);
    
  int surf = 0;
  int chan = 3;

  for (int chan=0; chan<NUM_DIGITZED_CHANNELS; chan++) {
    cout << endl << chan << " | ";
    for (int samp=0; samp<NUM_SAMP; samp++){
      cout << usefulRawEvent->fCapacitorNum[chan][samp] << " ";
    }
  }


  return 1;
}
