#include <math.h>
#include <iostream>

using namespace std;

int pedIndex(int surf, int chan, int lab, int sample) {
  return (surf*8*4*259) + (chan*4*259) + (lab*259) + sample;
}

void loadPedCorrections(double* pedCorrections) {

  ifstream inFile("pedCorrections.txt");

  int surf,chan,lab,sample;
  double pedCorr;
  while (inFile >> surf >> chan >> lab >> sample >> pedCorr) {
    pedCorrections[pedIndex(surf,chan,lab,sample)] = pedCorr;
  }

  return;


}


double findXOffset(double amp, double freq, double phase, double offset, double xValue, double yValue) {

  //  cout << "params " <<  amp << " " << freq << " " << phase << " " << offset << endl;


  //first lets determine which half-period of the sine wave this is (from x I guess)
  double period = 1./freq;
  double xAngle = fmod((xValue*freq*2*M_PI + phase),(2.*M_PI)); //this should be 0->2pi


  //I have to figure out if it is in quadrant 2 or 3
  int quadrant = 0;
  if (xAngle > M_PI/2. && xAngle <= M_PI) {
    quadrant = 2;
  }
  else if (xAngle > M_PI && xAngle <= (3./2.)*M_PI) {
    quadrant = 3;
  }

  //  cout << "x=" << xValue << "|" << fmod(xValue,period) << " xAngle " << xAngle << " " << quadrant << endl;

  //determine what the angle should be of something with that y
  double temp0 = (yValue - offset)/amp;
  //thats the normalized value, so we can figure out what angle that is;
  double temp1 = asin( temp0 );
  //if the x value is in the "odd" region not returned by asin, need to shift it into that region
  double temp2;
  if (quadrant == 2) {
    temp2 = M_PI - temp1;
  }
  else if (quadrant == 3) {
    temp2 = 3.*M_PI - temp1;
  }
  else {
    temp2 = temp1;
  }
  double temp3 = temp2 - phase;
