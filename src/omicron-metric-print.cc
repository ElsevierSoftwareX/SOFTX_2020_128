//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include <TF1.h>
#include <TriggerMetric.h>
#include "OmicronUtils.h"


using namespace std;

void PrintUsage(void){
  cerr<<endl;
  cerr<<"Usage:"<<endl;
  cerr<<endl;
  cerr<<"omicron-metric-print trigger-gps-end=[GPS end] \\"<<endl;
  cerr<<"                     trigger-gps-start=[GPS start] \\"<<endl;
  cerr<<"                     trigger-function=[TFormula] \\"<<endl;
  cerr<<"                     trigger-m1=[first mass] \\"<<endl;
  cerr<<"                     trigger-m2=[second mass] \\"<<endl;
  cerr<<"                     trigger-mchirp=[chirp mass] \\"<<endl;
  cerr<<"                     channel=[channel name] \\"<<endl;
  cerr<<"                     file=[trigger file pattern] \\"<<endl;
  cerr<<endl;
  cerr<<"[GPS end]                 Stopping GPS time (required)"<<endl;
  cerr<<"[GPS start]               Starting GPS time (required)"<<endl;
  cerr<<"[TFormula]                Trigger waveform function in the time-frequency plane"<<endl;
  cerr<<"                          The function returns the time [s] when given a frequency [Hz]:"<<endl;
  cerr<<"                          t=t(f)"<<endl;
  cerr<<"                          The function must be given as a ROOT TFormula:"<<endl;
  cerr<<"                          https://root.cern.ch/doc/master/classTFormula.html"<<endl;
  cerr<<"                          By default trigger-function=\"pow(x, -8.0/3.0)/[0]+[1]\", corresponding"<<endl;
  cerr<<"                          to a CBC Newtonian waveform"<<endl;
  cerr<<"[first mass]              First mass m1 in solar masses. By default =1.4"<<endl;
  cerr<<"[second mass]             Second mass m2 in solar masses. By default =1.4"<<endl;
  cerr<<"[chirp mass]              Chirp mass mc in solar masses. By default it is computed using m1 and m2"<<endl;
  cerr<<"[channel name]            Channel name used to retrieve centralized Omicron triggers."<<endl;
  cerr<<"                          By default: V1:Hrec_hoft_16384Hz"<<endl;
  cerr<<"[trigger file pattern]    File pattern to ROOT omicron trigger files"<<endl;
  cerr<<endl;
  return;
}

int main (int argc, char* argv[]){

  // number of arguments
  if(argc<3){
    PrintUsage();
    return 1;
  }

  // list of parameters + default
  string chname="";     // channel name
  string tfile_pat="";  // trigger file pattern
  double gps_start=-1;  // GPS start
  double gps_end=-1;    // GPS end
  string sfunc="pow(x, -8.0/3.0)/[0]+[1]";// user func
  double m1=1.4;        // m1
  double m2=1.4;        // m2
  double mc=-1.0;       // mc

  // loop over arguments
  vector <string> sarg;
  for(int a=1; a<argc; a++){
    sarg=SplitString((string)argv[a],'=');
    if(sarg.size()!=2) continue;
    if(!sarg[0].compare("channel"))           chname=(string)sarg[1];
    if(!sarg[0].compare("file"))              tfile_pat=(string)sarg[1];
    if(!sarg[0].compare("trigger-gps-start")) gps_start=stod(sarg[1]);
    if(!sarg[0].compare("trigger-gps-end"))   gps_end=stod(sarg[1]);
    if(!sarg[0].compare("trigger-function"))  sfunc=(string)sarg[1];
    if(!sarg[0].compare("trigger-m1"))        m1=stod(sarg[1]);
    if(!sarg[0].compare("trigger-m2"))        m2=stod(sarg[1]);
    if(!sarg[0].compare("trigger-mc"))        mc=stod(sarg[1]);
  }

  // required time range
  if(gps_start<0||gps_end<0){
    cerr<<"A trigger time range must be specified"<<endl;
      cerr<<"Type omicron-metric-print for help"<<endl;
      return 1;
  }

  
  return 0;
}

