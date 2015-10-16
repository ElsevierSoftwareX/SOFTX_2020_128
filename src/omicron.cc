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
    cerr<<endl;
    cerr<<argv[0]<<":"<<endl;
    cerr<<"This program runs an Omicron analysis."<<endl;
    cerr<<endl;
    cerr<<"Usage:"<<endl; 
    cerr<<endl;
    cerr<<"CASE1: omicron [GPS start time] [GPS stop time] [option file]"<<endl; 
    cerr<<"|__ runs the Omicron algorithm between 2 GPS times (integers)."<<endl; 
    cerr<<endl;
    cerr<<"CASE2: omicron [segment file] [option file]"<<endl; 
    cerr<<"|__ runs the Omicron algorithm over a segment list."<<endl; 
    cerr<<endl;
    cerr<<"CASE3: omicron [GPS time] [option file]"<<endl; 
    cerr<<"|__ runs the Omicron algorithm over one single chunk of data"<<endl;
    cerr<<"    centered on a GPS time."<<endl; 
    cerr<<endl;
    cerr<<">>> In all cases, the user must must provide an option file listing"<<endl;
    cerr<<"    the Omicron parameters."<<endl;
    cerr<<">>> In all cases, the GPS times must be integer values"<<endl;
    cerr<<endl;
   return 1;
  }

  int gps=0;
  int start;
  int stop;
  Segments *segments;
  string optionfile;
  string segmentfile;

  //**** CASE1: start and stop time
  if(argc==4){
    start=(int)round(atof(argv[1]));
    stop=(int)round(atof(argv[2]));
    optionfile=(string)argv[3];
 
    // segment object
    segments = new Segments(start,stop);
    if(!segments->GetLiveTime()){
      cerr<<"The input time segment to process does not make sense."<<endl;
      delete segments;
      return -1;
    }
  }

  //****** CASE2: segment file
  else if(atoi(argv[1])<700000000&&IsTextFile(argv[1])){
    segmentfile=(string)argv[1];
    optionfile=(string)argv[2];

    // load segment file
    segments = new Segments(segmentfile);
    if(!segments->GetLiveTime()){
      cerr<<"The input segments are corrupted."<<endl;
      delete segments;
      return -1;
    }
  }

  //****** CASE3: central time
  else if(atoi(argv[1])>700000000){
    gps=(int)round(atof(argv[1]));
    optionfile=(string)argv[2];
    // empty segment object
    segments = new Segments();
  }

  //******
  else{
    cerr<<"The option sequence cannot be interpreted."<<endl;
    cerr<<"Just type 'omicron' for help."<<endl;
    return -1;
  }

  // check that the option file exists
  if(!IsTextFile(optionfile)){
    cerr<<"omicron: The option file "<<optionfile<<" cannot be found."<<endl;
    delete segments;
    return -1;
  }

  // init omicron
  Omicron *O = new Omicron(optionfile);
  if(!O->GetStatus()){
    cerr<<"omicron: The Omicron object is corrupted."<<endl;
    delete segments;
    delete O;
    return -2;
  }

  O->PrintMessage("Omicron has been successfully initiated");

  // update for CASE3
  if(!segments->GetLiveTime())
    segments->AddSegment(gps-O->GetChunkDuration()/2,gps+O->GetChunkDuration()/2);

  // locals
  int dsize;
  double *dvector;
  int res;
  
  // channel list
  vector <string> Channels = O->GetChannels();

  // init segments
  if(!O->InitSegments(segments, gps)) return 1;

  // create trigger directories
  if(!O->MakeDirectories(gps)) return 2;

  O->PrintMessage("Start looping over chunks and channels");

  // loop over chunks
  while(O->NewChunk()){
      
    // new channels
    while(O->NewChannel()){
      
      // get data vector
      dvector=NULL; dsize=0;
      if(!O->LoadData(&dvector,&dsize)) continue;
      
      // condition data vector
      res=O->Condition(dsize, dvector);
      if(res<0){
	delete dvector;
	return 3;// fatal
      }
      if(res>0){
	delete dvector;
	continue;
      }
      delete dvector;// not needed anymore

      // project data
      if(!O->Project()) continue;

      // write chunk outputs
      if(!O->WriteOutput()) continue;
    }

  }
    
  O->PrintMessage("Omicron processing is over");

  // prints summary report
  O->PrintStatusInfo();
  
  // cleaning  
  delete O;
  delete segments;
  return 0;
}

