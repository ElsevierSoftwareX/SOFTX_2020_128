//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

using namespace std;

int main (int argc, char* argv[]){
  gErrorIgnoreLevel = 3000;

  // check the argument
  if(argc<3){
    cerr<<endl;
    cerr<<" usage:"<<endl; 
    cerr<<endl;
    cerr<<" omiscan [central time] [option file #1] [option file #2] ... [option file #N]"<<endl; 
    cerr<<endl;
    cerr<<" [central time]    =  GPS where the maps will be centred"<<endl; 
    cerr<<" [option file #i]  =  Omicron option file"<<endl; 
    cerr<<"                      For available options, please visit:"<<endl; 
    cerr<<"                      http://virgo.in2p3.fr/GWOLLUM/Friends/omicron.html"<<endl; 
    cerr<<endl;
    return 1;
  }
  
  ///////////////////////////////////////////////////////////////////////////////
  /////////                        COMMAND LINE                         /////////
  ///////////////////////////////////////////////////////////////////////////////
  
  // central GPS
  double gps=atof(argv[1]);

  // test GPS time
  if(gps<700000000){
    cerr<<"Omiscan ERROR: GPS time '"<<gps<<"' is not reasonable"<<endl;
    return 2;
  }

  // number of option files
  int nopt = argc-2;

  // get option files
  string *optionfile = new string [nopt];
  for(int o=0; o<nopt; o++){
    optionfile[o] =(string)argv[2+o];

    // test option file
    if(!IsTextFile(optionfile[o])){
      cerr<<"Omiscan ERROR: option file '"<<optionfile[o]<<"' cannot be found/read"<<endl;
      delete [] optionfile;
      return 2;
    }
  }


  // Omicron object
  Omicron *O;

  // loop over option files
  for(int o=0; o<nopt; o++){

    // init omicron for this option file
    O = new Omicron(optionfile[o]);
    if(!O->GetStatus()){
      cerr<<"Omiscan ERROR: input options of "<<optionfile[o]<<" are invalid"<<endl;
      if(o==nopt-1) O->ScanReport(); // report if this is the last round
      delete O;
      continue;
    }

    // scan it
    if(!O->Scan(gps)){
      cerr<<"Omiscan ERROR: scan failed for "<<optionfile[o]<<endl;
      if(o==nopt-1) O->ScanReport(); // report if this is the last round
      delete O;
      continue;
    }

    // report if this is the last round
    if(o==nopt-1) O->ScanReport();
    delete O;
  }

  delete [] optionfile;
  return 0;
}
  
