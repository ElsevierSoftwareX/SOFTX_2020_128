//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string.h>
#include <stdio.h>
#include "Omicron.h"
#include "IO.h"

using namespace std;

int main (int argc, char* argv[]){

  // check the argument
  if(argc!=4){
    cerr<<argv[0]<<" usage:"<<endl<<endl; 
    cerr<<argv[0]<<" [channel file] [option file] [central time]"<<endl; 
    return 1;
  }

  ///////////////////////////////////////////////////////////////////////////////
  // read command line parameters
  string channelfile=(string)argv[1];
  string optionfile=(string)argv[2];
  double gps=atof(argv[3]);

  // check command line parameters
  if(!IsTextFile(channelfile)){
    cerr<<"Omiscan ERROR: channel file '"<<channelfile<<"' cannot be found"<<endl;
    return 2;
  }
  if(!IsTextFile(optionfile)){
    cerr<<"Omiscan ERROR: option file '"<<optionfile<<"' cannot be found"<<endl;
    return 2;
  }
  if(gps<700000000){
    cerr<<"Omiscan ERROR: GPS time '"<<gps<<"' is not reasonable"<<endl;
    return 2;
  }

  ///////////////////////////////////////////////////////////////////////////////
  // Parse option file
  IO *io = new IO(optionfile.c_str());

  //**** ffl file
  string fflfile;
  if(io->GetOpt("DATA","FFL", fflfile)){
    if(!IsTextFile(fflfile)){
      cerr<<"Omiscan ERROR: FFL file '"<<fflfile<<"' cannot be found"<<endl;
      return 3;
    }
  }
  else{
    cerr<<"Omiscan ERROR: FFL file is required"<<endl;
    return 3;
  }

  //**** windows
  //FIXME: power of 2?
  vector <int> windows;
  if(io->GetOpt("PARAMETER","WINDOW", windows)){
    if(!windows.size()){
      cerr<<"Omiscan ERROR: the window option is not correct"<<endl;
      return 3;
    }
  }
  else{// default windows
    windows.clear();
    windows.push_back(4); windows.push_back(16); windows.push_back(64);
  }

  delete io;

  ///////////////////////////////////////////////////////////////////////////////
  // Channels to scan
  vector<string> channels;
  ReadAscii *C = new ReadAscii(channelfile,"s");
  if(!C->GetCol(channels, 0)){
    cerr<<"Omiscan ERROR: channel file '"<<channelfile<<"' has not the right format"<<endl;
    delete C;
    return 3;
  }
  if(!channels.size()){
    cerr<<"Omiscan ERROR: no channel to scan"<<endl;
    delete C;
    return 3;
  }
  delete C;

  ///////////////////////////////////////////////////////////////////////////////
  FrFile *frfile = FrFileINew((char*)fflfile.c_str());
  
  FrVect *chanvect = NULL; // data vector

  // timing
  int timerange=windows[(int)windows.size()-1];
  double start=gps-(double)timerange/2.0;
  double stop=gps+(double)timerange/2.0;
  Segments *segment = new Segments();
  if(!segment->AddSegment(start,stop)){
    cerr<<"Omiscan ERROR: problem with segment"<<endl;
    delete segment;
    return 4;
  }
  
  int sampling, sampling_new;
  int powerof2;
  Odata *data;
  Otile *tiles;

  ///////////////////////////////////////////////////////////////////////////////
  // loop over channels
  for(int c=0; c<(int)channels.size(); c++){

    cout<<"Omiscan: prepare channel "<<channels[c]<<endl;

    // get data
    chanvect = NULL;
    chanvect = FrFileIGetVectDN(frfile,(char*)(channels[c].c_str()),start,stop-start);
 
    // check data vector
    if(chanvect==NULL){
      cout<<"Omiscan WARNING: missing channel --> skip"<<endl;
      continue;
    }
    if(!chanvect->nData){
      cout<<"Omiscan WARNING: no data --> skip"<<endl;
      FrVectFree(chanvect);
      continue;
    }
    if((int)(chanvect->nData)%(int)(timerange)){
      cout<<"Omiscan WARNING: weird sampling --> skip"<<endl;
      FrVectFree(chanvect);
      continue;
    }

    // sampling
    sampling=chanvect->nData/(timerange);
    powerof2=(int)floor(log(sampling)/log(2.0));
    sampling_new=TMath::Min((int)pow(2,powerof2),4096);// FIXME: hard-coded
    cout<<"         sampling = "<<sampling<<"Hz --> "<<sampling_new<<"Hz"<<endl;

    // init data
    data = new Odata("ONLINE", channels[c], sampling, sampling_new, segment, stop-start, stop-start, 8, 2);// FIXME: 2 and 8 is hard-coded

    // init tiles
    tiles = new Otile(timerange, 8, 3, 141, 2, sampling_new/2.0, sampling_new, 0.3, 8);

    

    delete tiles;
    delete data;

  }

  // cleaning
  delete segment;

  return 0;
}

