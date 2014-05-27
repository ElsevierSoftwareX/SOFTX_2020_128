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
    cerr<<endl;
    cerr<<"This command runs the omicron algorithm over a stretch of data."<<endl;
    cerr<<"The user must define the data to process with either a starting and stoping time"<<endl;
    cerr<<"or a segment file."<<endl;
    cerr<<"The user must provide an option file listing the Omicron parameters."<<endl;
   return 1;
  }

  int start;
  int stop;
  Segments *segments;
  string optionfile;
  string segmentfile;

  //****** start and stop time ******
  if(argc==4){
    start=atoi(argv[1]);
    stop=atoi(argv[2]);
    optionfile=(string)argv[3];
 
    // segment object
    segments = new Segments(start,stop);
    if(!segments->GetStatus()){
      cerr<<argv[0]<<": The input segment to process does not make sense."<<endl;
      delete segments;
      return 2;
    }
  }


  //****** segment file ******
  else{
    segmentfile=(string)argv[1];
    optionfile=(string)argv[2];

    // load segment file
    segments = new Segments(segmentfile);
    if(!segments->GetStatus()){
      cerr<<argv[0]<<": The input segments are corrupted."<<endl;
      delete segments;
      return 2;
    }
  }

  // check that the option file exists
  if(!IsTextFile(optionfile)){
    cerr<<argv[0]<<": The option file "<<optionfile<<" cannot be found."<<endl;
    delete segments;
    return 3;
  }

  // let's go!
  Omicron *O = new Omicron(optionfile);
  if(!O->GetStatus()){
    cerr<<argv[0]<<": The Omicron object is corrupted."<<endl;
    delete segments;
    delete O;
    return 2;
  }

  // processing
  if(!O->Process(segments)){
    cerr<<argv[0]<<": The Omicron processing failed."<<endl;
    O->PrintStatusInfo();
    delete O;
    delete segments;
    return 1;
  }

  // prints summary
  O->PrintStatusInfo();
  
  // cleaning  
  delete O;
  delete segments;
  return 0;
}

