//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Oplot.h"

using namespace std;

int main (int argc, char* argv[]){

  // check the argument
  if(argc<5){
    cerr<<argv[0]<<" usage:"<<endl; 
    cerr<<argv[0]<<" [outdir] [trigger root file pattern] [GPS start] [GPS stop] [filename]"<<endl; 
    cerr<<endl; 
    cerr<<" [outdir]                   : output directory where to save the plots"<<endl; 
    cerr<<" [trigger root file pattern]: list of trigger ROOT files."<<endl; 
    cerr<<"                              This list can be composed of file paths"<<endl;
    cerr<<"                              and/or file patterns"<<endl; 
    cerr<<" [GPS start]                : GPS start"<<endl; 
    cerr<<" [GPS stop]                 : GPS stop"<<endl; 
    cerr<<" [filename]                 : this (optional) parameter force the output file naming to:"<<endl; 
    cerr<<"                              [filename]_type.gif"<<endl;
    cerr<<endl; 
    return 1;
  }

  // outdir
  string outdir = (string)argv[1];
  system(("mkdir -p "+outdir).c_str());

  // input files
  string infiles = (string)argv[2];
   
  // timing
  int start = atoi(argv[3]);
  int stop = atoi(argv[4]);
  string sstart = (string)argv[3];
  string sstop = (string)argv[4];

  // file name
  string filename = "";
  if(argc>5) filename=(string)argv[5];
 
  // input segments
  ReadTriggerMetaData *insegments = new ReadTriggerMetaData(infiles,"",0);
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
  Oplot *triggers;
  if(trigfiles.size()<10000) triggers = new Oplot(trigfiles,"",1);
  else triggers = new Oplot(infiles,"",1);
  trigfiles.clear();

  triggers->SetTimeRange(start,stop);
  triggers->MakeCollections();

  triggers->PrintPlot("snr");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_snr.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_snr_"+sstart+"_"+sstop+".gif");

  triggers->PrintCollectionPlot("rate");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_rate.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_rate_"+sstart+"_"+sstop+".gif");

  triggers->PrintCollectionPlot("frequency");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_frequency.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_frequency_"+sstart+"_"+sstop+".gif");

  triggers->PrintCollectionPlot("freqtime");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_freqtime.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_freqtime_"+sstart+"_"+sstop+".gif");
 
  triggers->PrintCollectionPlot("snrtime");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_snrtime.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_snrtime_"+sstart+"_"+sstop+".gif");

  triggers->PrintCollectionPlot("snrfreq");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_snrfreq.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_snrfreq_"+sstart+"_"+sstop+".gif");
   
  triggers->PrintCollectionPanel();
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_panel.gif",2);
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_panel_"+sstart+"_"+sstop+".gif",2);
 
  // cleaning
  delete triggers;
  
  return 0;
}

