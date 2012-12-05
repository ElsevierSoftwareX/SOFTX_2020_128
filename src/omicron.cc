//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string.h>
#include <stdio.h>
#include "Omicron.h"
#include "Segments.h"

using namespace std;

int main (int argc, char* argv[]){

  // check the argument
  if(argc!=4&&argc!=3){
    cerr<<argv[0]<<" usage:"<<endl; 
    cerr<<argv[0]<<" [start time] [stop time] [option file]"<<endl; 
    cerr<<"OR:"<<endl; 
    cerr<<argv[0]<<" [segment file] [option file]"<<endl; 
    return 1;
  }

  cout<<"\nRead parameters..."<<endl;
  int start;
  int stop;
  vector<double> segment_starts; vector<double> segment_ends;
  Segments *segments;
  string optionfile;
  string segmentfile;

  if(argc==4){// start and stop time
    start=atoi(argv[1]); segment_starts.push_back(start);
    stop=atoi(argv[2]);  segment_ends.push_back(stop);
    optionfile=(string)argv[3];
 
    // check the time is reasonable
    if(start<700000000||stop<=start){
      cerr<<"argv[0]: the start time is not reasonable!"<<endl; 
      return 1;
    }
    if(stop<4){
      cerr<<"argv[0]: the stop time is not reasonable!"<<endl; 
      return 1;
    }
    
    // segment object
    segments = new Segments(segment_starts,segment_ends);

  }
  else{// segment file
    segmentfile=(string)argv[1];
    optionfile=(string)argv[2];

    // check that the option file exists
    if(!IsTextFile(segmentfile)) return 1;

    // load segment file
    segments = new Segments(segmentfile);
  }

  // check that the option file exists
  if(!IsTextFile(optionfile)) return 1;

  // let's go!
  cout<<"\nCreate Omicron object..."<<endl;
  Omicron *O = new Omicron(segments,optionfile);

  // tiling
  cout<<"\nMake Q-F-T tiling..."<<endl;
  if(!O->MakeTiling()){
    cerr<<" failed"<<endl;
    delete O;
    delete segments;
    return 1;
  }
  
  // processing
  cout<<"\nStart Processing..."<<endl;
  if(!O->Process()){
    cerr<<" failed"<<endl;
    delete O;
    delete segments;
    return 1;
  }
  
  
  delete O;
  delete segments;
  return 0;
}

