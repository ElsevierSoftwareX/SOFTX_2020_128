//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

using namespace std;

int main (int argc, char* argv[]){
  gErrorIgnoreLevel = 3000;

  // check the argument
  if(argc!=3){
    cerr<<argv[0]<<" usage:"<<endl<<endl; 
    cerr<<argv[0]<<" [central time] [option file]"<<endl; 
    return 1;
  }
  
  ///////////////////////////////////////////////////////////////////////////////
  /////////                        COMMAND LINE                         /////////
  ///////////////////////////////////////////////////////////////////////////////
  // command line parameters
  double gps=atof(argv[1]);
  string optionfile=(string)argv[2];

  // check command line parameters
  if(!IsTextFile(optionfile)){
    cerr<<"Omiscan ERROR: option file '"<<optionfile<<"' cannot be found"<<endl;
    return 2;
  }
  if(gps<700000000){
    cerr<<"Omiscan ERROR: GPS time '"<<gps<<"' is not reasonable"<<endl;
    return 2;
  }

  // init Omicron
  Omicron *O = new Omicron(optionfile);
  if(!O->GetStatus()){
    cerr<<"Omiscan ERROR: input options are invalid"<<endl;
    delete O;
    return 3;
  }

  // scan
  if(!O->Scan(gps)){
    cerr<<"Omiscan ERROR: scan failed"<<endl;
    delete O;
    return 4;
  }


  delete O;
  return 0;
}
  
