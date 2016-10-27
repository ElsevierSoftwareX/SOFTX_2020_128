//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "OmicronUtils.h"

using namespace std;

void printhelp(void){
    cerr<<endl;
    cerr<<"This program prints the list of Omicron trigger files."<<endl;
    cerr<<endl;
    cerr<<"Usage:"<<endl; 
    cerr<<endl;
    cerr<<"CASE1: omicron-listfile [channel name] [GPS start time] [GPS stop time]"<<endl; 
    cerr<<"|____ prints the list of files associated to given channel between 2 GPS times."<<endl; 
    cerr<<endl;
    cerr<<"CASE2: omicron-listfile [channel name] [GPS time]"<<endl; 
    cerr<<"|____ prints the file associated to a given channel containing a GPS time."<<endl; 
    cerr<<endl;
    return;
}
int main (int argc, char* argv[]){

  // check the command line
  if(argc<3){
    printhelp();
    return -1;
  }

  // channel name
  string channelname = (string)argv[1];

  // timing
  int start = atoi(argv[2]);
  int stop=start;  
  if(argc>3) stop = atoi(argv[3]);

  cout<<GetOmicronFilePattern(channelname,start,stop)<<endl;
  
  return 0;
}

