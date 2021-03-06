/**
 * @file 
 * @brief Program to run an Omicron analysis.
 * @details `omicron` is a command line program to perform a Q-transform analysis of data time-series.
 * The program must be given a timing and an option file.
 * @verbatim
omicron [GPS start time] [GPS stop time] [option file]
 @endverbatim
 * This command runs the Omicron algorithm between 2 GPS times.
 * The option file is a text file in which the parameters are listed using keywords.
 * For more details about the syntax, see Omicron::ReadOptions().
 *
 * @verbatim
omicron [segment file] [option file]
 @endverbatim
 * This command runs the Omicron algorithm over a time segment list.
 * The segment file is a text file with 2 columns listing the time segments (GPS start - GPS end).
 *
 * @verbatim
omicron [GPS time] [option file]
 @endverbatim
 * This command runs the Omicron algorithm over once single time segment centered on `[GPS time]`.
 *
 * @note all GPS time values must be integer values.
 *
 * Optional arguments can be provided.
 * - the output triggers can be time-selected using a list of segments:
 * @verbatim
omicron [some input timing] [option file] [GPS start time] [GPS stop time]
 @endverbatim
 * In that case, output triggers are only saved if they start between the 2 GPS times given in the command line.
 * @verbatim
omicron [some input timing] [option file] [segment file]
 @endverbatim
 * In that case, output triggers are only saved if they start in one of the time segment given in the command line.
 * - One last argument can be provided in the command line: `strict`.
 * In this mode, the program exits whenever an error is met.
 *
 * @author Florent Robinet - <a href="mailto:florent.robinet@ijclab.in2p3.fr">florent.robinet@ijclab.in2p3.fr</a>
 */
#include <iostream>
#include <string.h>
#include <stdio.h>
#include "Oomicron.h"
#include "Segments.h"

using namespace std;

/**
 * @brief Print the program usage message.
 */
void printhelp(void){
    cerr<<endl;
    cerr<<"This program runs an Omicron analysis."<<endl;
    cerr<<endl;
    cerr<<"Usage:"<<endl; 
    cerr<<endl;
    cerr<<"CASE1: omicron [GPS start time] [GPS stop time] [option file]"<<endl; 
    cerr<<"|____ runs the Omicron algorithm between 2 GPS times."<<endl; 
    cerr<<endl;
    cerr<<"CASE2: omicron [segment file] [option file]"<<endl; 
    cerr<<"|____ runs the Omicron algorithm over a segment list."<<endl; 
    cerr<<endl;
    cerr<<"CASE3: omicron [GPS time] [option file]"<<endl; 
    cerr<<"|____ runs the Omicron algorithm over one single chunk of data"<<endl;
    cerr<<"      centered on a GPS time."<<endl; 
    cerr<<endl;
    cerr<<">>> In all cases, the user must must provide an option file listing"<<endl;
    cerr<<"    the Omicron parameters."<<endl;
    cerr<<">>> In all cases, the GPS times must be integer values"<<endl;
    cerr<<endl;
    cerr<<endl;
    cerr<<"-----------------------------------------------------------------------------------"<<endl;
    cerr<<"Optionally, the output triggers can be time-selected using a list of segments"<<endl;
    cerr<<"This option is activated by adding extra arguments to the command line:"<<endl;
    cerr<<endl;
    cerr<<"OPTION1: omicron [some input timing] [option file] [GPS start time] [GPS stop time]"<<endl; 
    cerr<<"|____ output triggers are only saved if they start between 2 GPS times."<<endl;
    cerr<<endl;
    cerr<<"OPTION2: omicron [some input timing] [option file] [segment file]"<<endl; 
    cerr<<"|____ output triggers are only saved if they start in some list of time segments."<<endl;
    cerr<<"-----------------------------------------------------------------------------------"<<endl;
    cerr<<endl;
    cerr<<endl;
    cerr<<"-----------------------------------------------------------------------------------"<<endl;
    cerr<<"One last argument can be provided in the command line:"<<endl;	
    cerr<<"If the last argument is the character string \"strict\", the Omicron algorithm is run"<<endl;	
    cerr<<"in a failure mode. It means that the program exits whenever an error is met."<<endl;	
    cerr<<"-----------------------------------------------------------------------------------"<<endl;
    return;
}

/**
 * @brief Main program.
 */
int main (int argc, char* argv[]){

  if(argc>1&&!((string)argv[1]).compare("version")){
    PrintVersion();
    return 0;
  }

  // check the command line
  if(argc<3||argc>7){
    printhelp();
    return -1;
  }
 
  // first check the strict option
  bool strict=false;
  if(!strcmp(argv[argc-1],"strict")){
    strict=true;
    argc--;
  }
  
  int start=0;
  int stop=0;
  int sub;
  string segmentfile="none";
  string optionfile="none";
  int ostart=0;
  int ostop=0;
  string osegmentfile="none";

  // get first argument
  start=(int)floor(atof(argv[1]));
  sub=(int)((atof(argv[1])-start)*1000.);
  if(!start)// a segment file is provided
    segmentfile=(string)argv[1];

  // get second argument
  stop=(int)floor(atof(argv[2]));
  if(!stop)// a parameter file is provided
    optionfile=(string)argv[2];

  // get third argument
  if(argc>3){
    ostart=(int)floor(atof(argv[3]));
    if(!ostart){
      if(!optionfile.compare("none")) optionfile=(string)argv[3];
      else osegmentfile=(string)argv[3];
    }
  }

  // get fourth argument
  if(argc>4){
    ostop=(int)floor(atof(argv[4]));
    if(!ostop) osegmentfile=(string)argv[4];
    else if(!ostart){ ostart=ostop; ostop=0; }
    else;
  }

  // get fifth argument
  if(argc>5){
    ostop=(int)round(atof(argv[5]));
  }

  // check parameter file
  if(!IsTextFile(optionfile)){
    cerr<<"omicron: A valid parameter file must be provided."<<endl;
    return -2;
  }

  // init omicron
  Omicron *O;
  if(start) O = new Omicron(optionfile,start,strict);
  else if(stop) O = new Omicron(optionfile,start,strict);
  else O = new Omicron(optionfile,-1,strict);
  
  if(!O->GetStatus()){
    cerr<<"omicron: The Omicron object is corrupted."<<endl;
    delete O;
    return -3;
  }
  O->PrintMessage("Omicron has been successfully initiated");

  // input segments
  Segments *insegments;
  if(segmentfile.compare("none")) insegments = new Segments(segmentfile);
  else if(start&&stop) insegments = new Segments(start,stop);
  else if(start){
    insegments = new Segments(start-O->GetChunkDuration()/2, start+O->GetChunkDuration()/2);
    O->SetPlotTimeOffset((double)sub/1000.0);
  }
  else{
    cerr<<"omicron: A valid input timing must be provided."<<endl;
    delete O;
    return -4;
  }

  // output segments
  Segments *outsegments;
  if(osegmentfile.compare("none")) outsegments = new Segments(osegmentfile);
  else if(ostart&&ostop) outsegments = new Segments(ostart,ostop);
  else outsegments = new Segments(insegments->GetStarts(),insegments->GetEnds());

  // check segments
  if(!insegments->GetLiveTime()||!outsegments->GetLiveTime()){
    cerr<<"omicron: A valid timing must be provided."<<endl;
    delete O;
    delete insegments;
    delete outsegments;
    return -5;
  }

  
  // locals
  int dsize;
  double *dvector;
  int res;
  
  // init segments
  if(!O->InitSegments(insegments,outsegments)) return 1;
  delete insegments;
  delete outsegments;

  // create specific trigger directories
  if(!stop&&!segmentfile.compare("none")){
    if(!O->MakeDirectories((double)start+(double)sub/1000.0)) return 2;
  }
  else{
    if(!O->MakeDirectories()) return 2;
  }

  O->PrintMessage("Start looping over chunks and channels");

  // loop over chunks
  while(O->NewChunk()){
      
    // new channels
    while(O->NewChannel()){
      
      // get data vector
      dvector=NULL; dsize=0;
      if(!O->LoadData(&dvector,&dsize)){
	if(strict){ delete O; return 3; }
	else continue;
      }
	
      // condition data vector
      res=O->Condition(dsize, dvector);
      if(res<0){
	delete dvector;
	return 3;// fatal
      }
      if(res>0){
	delete dvector;
	if(strict){ delete O; return 4; }
	else continue;
      }
      delete dvector;// not needed anymore

      // project data
      if(O->Project()<0){
	if(strict){ delete O; return 5; }
	else continue;
      }
            
      // write chunk outputs
      if(!O->WriteOutput()){
	if(strict){ delete O; return 6; }
	else continue;
      }

    }
  }
  
  O->PrintMessage("Omicron processing is over");
  
  // prints summary report
  O->PrintStatusInfo();
  
  // cleaning  
  delete O;
  return 0;
}

