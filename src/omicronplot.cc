//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string.h>
#include <stdio.h>
#include "TriggerPlot.h"
#include "EventMap.h"
#include "Triggers.h"
#include "Segments.h"

using namespace std;

int main (int argc, char* argv[]){

  // check the argument
  if(argc!=6){
    cerr<<argv[0]<<" usage:"<<endl; 
    cerr<<argv[0]<<" [title] [outdir] [list of root files] [GPS start] [GPS stop]"<<endl; 
    cerr<<endl; 
    cerr<<" [title]: plot title prefix"<<endl; 
    cerr<<" [outdir]: output directory where to save the plots"<<endl; 
    cerr<<" [list of root files]: list of trigger ROOT files separated by spaces."<<endl; 
    cerr<<" This list can be composed of file path and/or file patterns"<<endl; 
    cerr<<" [GPS start]: GPS start"<<endl; 
    cerr<<" [GPS stop]: GPS stop"<<endl; 
    cerr<<endl; 
    return 1;
  }

  // outdir
  string title = (string)argv[1];
  
  // outdir
  string outdir = (string)argv[2];
  system(("mkdir -p "+outdir).c_str());

  // input files
  string infiles = (string)argv[3];
   
  // timing
  int start = atoi(argv[4]);
  int stop = atoi(argv[5]);
  string sstart = (string)argv[4];
  string sstop = (string)argv[5];
 
  // input segments
  ReadTriggerSegments *insegments = new ReadTriggerSegments(infiles,"",0);
  if(!insegments->GetNFiles()){
    cerr<<"No trigger files for this time segment"<<endl;
    delete insegments;
    return 1;
  }
  
  // select files of interest
  string trigfiles = insegments->GetTriggerFiles(start,stop);
  if(!trigfiles.compare("none")){
    cerr<<"No trigger files for this time segment"<<endl;
    delete insegments;
    return 1;
  }
  delete insegments;
  
  // input triggers
  TriggerPlot *triggers;
  if(trigfiles.size()<10000) triggers = new TriggerPlot(trigfiles,"",1);
  else triggers = new TriggerPlot(infiles,"",1);
  if(!triggers->GetNTrig()){
    cerr<<"No triggers to plot"<<endl;
    delete triggers;
    return 2;
  }
  trigfiles.clear();

  // livetime
  double live = triggers->GetSegments()->GetLiveTime(start,stop);

  // make plots
  triggers->SetTitlePrefix(title+": ");
  triggers->OmicronPlot("MAP",start,stop,outdir+"/map_"+sstart+"_"+sstop+".gif");
  triggers->OmicronPlot("FREQ",start,stop,outdir+"/freq_"+sstart+"_"+sstop+".gif");
  triggers->OmicronPlot("SNR",start,stop,outdir+"/snr_"+sstart+"_"+sstop+".gif");
  triggers->OmicronPlot("RATE",start,stop,outdir+"/rate_"+sstart+"_"+sstop+".gif");
  triggers->OmicronPlot("SNRFREQ",start,stop,outdir+"/snrfreq_"+sstart+"_"+sstop+".gif");
  triggers->OmicronPlot("SNRTIME",start,stop,outdir+"/snrtime_"+sstart+"_"+sstop+".gif");

  // segment plot
  triggers->ResizePlot(1000,100);
  triggers->OmicronPlot("SEG",start,stop,outdir+"/segments_"+sstart+"_"+sstop+".gif");

  // info file
  ofstream infofile((outdir+"/info_"+sstart+"_"+sstop+".txt").c_str());
  infofile.flags ( ios::fixed );
  infofile<<setprecision(5)<<endl;
  infofile<<"NTRIGGERS "<<triggers->GetNTrig()<<endl;
  infofile<<"NCLUSTERS "<<triggers->GetNCluster()<<endl;
  infofile<<"NCLUSTERS0 "<<triggers->GetNTrigPlot(1)<<endl;
  infofile<<"NCLUSTERS1 "<<triggers->GetNTrigPlot(2)<<endl;
  infofile<<"NCLUSTERS2 "<<triggers->GetNTrigPlot(4)<<endl;
  infofile<<"NCLUSTERS3 "<<triggers->GetNTrigPlot(8)<<endl;
  infofile<<"RATE0 "<<(double)triggers->GetNTrigPlot(1)/live<<endl;
  infofile<<"RATE1 "<<(double)triggers->GetNTrigPlot(2)/live<<endl;
  infofile<<"RATE2 "<<(double)triggers->GetNTrigPlot(4)/live<<endl;
  infofile<<"RATE3 "<<(double)triggers->GetNTrigPlot(8)/live<<endl;
  infofile<<"FREQMIN "<<triggers->GetFreqMinPlot(1)<<endl;
  infofile<<"FREQMAX "<<triggers->GetFreqMaxPlot(1)<<endl;
  infofile<<"SNRMIN "<<triggers->GetSNRMinPlot(1)<<endl;
  infofile<<"SNRMAX "<<triggers->GetSNRMaxPlot(1)<<endl;
  infofile<<"TIMESNRMAX "<<triggers->GetTimeSNRMaxPlot(1)<<endl;
  infofile<<"FREQSNRMAX "<<triggers->GetFreqSNRMaxPlot(1)<<endl;
  infofile<<"DELTAT "<<triggers->GetDeltaT()<<endl;
  infofile<<"LIVETIME "<<live<<endl;
  infofile.close();

  // event mapping of the loudest event
  EventMap *loudestevent= new EventMap(triggers->GetTriggerFile(triggers->GetTimeSNRMaxPlot(1)),"", 0);
  loudestevent->BuildMap(0,triggers->GetTimeSNRMaxPlot(1));
  loudestevent->PrintMap(0,outdir+"/loudest_"+sstart+"_"+sstop+".gif");
  delete loudestevent;
  
  // cleaning
  delete triggers;
  
  return 0;
}

