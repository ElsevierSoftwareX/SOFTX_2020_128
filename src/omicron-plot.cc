//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "TriggerPlot.h"
#include "OmicronUtils.h"


using namespace std;

void PrintUsage(void){
  cerr<<endl;
  cerr<<"Usage:"<<endl;
  cerr<<endl;
  cerr<<"omicron-plot.exe channel=[channel name] \\"<<endl;
  cerr<<"                 file=[trigger file pattern] \\"<<endl;
  cerr<<"                 gps-start=[GPS start] \\"<<endl;
  cerr<<"                 gps-end=[GPS end] \\"<<endl;
  cerr<<"                 snr-thresholds=[list of SNR thresholds] \\"<<endl;
  cerr<<"                 use-cluster=[cluster flag] \\"<<endl;
  cerr<<"                 cluster-dt=[cluster time window] \\"<<endl;
  cerr<<"                 use-date=[date flag] \\"<<endl;
  cerr<<"                 outdir=[output directory] \\"<<endl;
  cerr<<"                 outformat=[output file format] \\"<<endl;
  cerr<<"                 style=[style]"<<endl;
  cerr<<endl;
  cerr<<"[channel name]            channel name used to retrieve centralized Omicron triggers"<<endl;
  cerr<<"[trigger file pattern]    file pattern to ROOT trigger files (GWOLLUM convention)"<<endl;
  cerr<<"[GPS start]               starting GPS time (integer only)"<<endl;
  cerr<<"[GPS end]                 stopping GPS time (integer only)"<<endl;
  cerr<<"[list of SNR thresholds]  list of SNR thresholds. By default, snr-thresholds=\"5;8;10;20\""<<endl;
  cerr<<"[cluster tag]             1 = use clusters, 0 = use triggers. By default, use-cluster=1"<<endl;
  cerr<<"[cluster time window]     cluster time window [s]. By default, cluster-dt=0.1"<<endl;
  cerr<<"[date tag]                1 = use date, 0 = use GPS time. By default, use-date=1"<<endl;
  cerr<<"[output directory]        output directory. By default, outdir=."<<endl;
  cerr<<"[output file format]      output file format. By default, outformat=png"<<endl;
  cerr<<"[style]                   GWOLLUM-supported style. By default, style=GWOLLUM"<<endl;
  cerr<<endl;
  return;
}

int main (int argc, char* argv[]){

  // number of arguments
  if(argc<2){
    PrintUsage();
    return 1;
  }

  // list of parameters + default
  string chname=""; // channel name
  string tfile_pat="";  // trigger file pattern
  string snrthrs="5;8;10;20";  // list of SNR thresholds
  int gps_start=-1;      // GPS start
  int gps_end=-1;  // GPS end
  string style="GWOLLUM"; // style
  bool usecluster=true; // use cluster
  double cluster_dt=0.1; // cluster time window
  bool usedate=true; // use date
  string outdir=".";
  string outformat="png";

  // loop over arguments
  vector <string> sarg;
  for(int a=1; a<argc; a++){
    sarg=SplitString((string)argv[a],'=');
    if(sarg.size()!=2) continue;
    if(!sarg[0].compare("channel"))        chname=(string)sarg[1];
    if(!sarg[0].compare("file"))           tfile_pat=sarg[1];
    if(!sarg[0].compare("gps-start"))      gps_start=atoi(sarg[1].c_str());
    if(!sarg[0].compare("gps-end"))        gps_end=atoi(sarg[1].c_str());
    if(!sarg[0].compare("snr-thresholds")) snrthrs=(string)sarg[1];
    if(!sarg[0].compare("style"))          style=(string)sarg[1];
    if(!sarg[0].compare("use-cluster"))    usecluster=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("cluster-dt"))     cluster_dt=atof(sarg[1].c_str());
    if(!sarg[0].compare("use-date"))       usedate=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("outdir"))         outdir=(string)sarg[1];
    if(!sarg[0].compare("outformat"))      outformat=(string)sarg[1];
  }

  // centralized trigger files
  if(!tfile_pat.compare("")){
    if(chname.compare("")&&gps_start>0&&gps_end>0) tfile_pat=GetOmicronFilePattern(chname,gps_start,gps_end);// get centralized trigger files
    else{
      cerr<<"Triggers must be specified, either with a channel name and a time range or a file pattern"<<endl;
      cerr<<"Type omicron-plot.exe for help"<<endl;
      return 1;
    }
  }

  // check file pattern
  if(!tfile_pat.compare("")){
    cerr<<"No trigger files"<<endl;
    return 2;
  }
  
  // snr thresholds
  sarg=SplitString(snrthrs,';');
  vector <double> snrthr;
  for(int s=0; s<(int)sarg.size(); s++) snrthr.push_back(atof(sarg[s].c_str()));

  stringstream tmpstream;

  // triggers
  TriggerPlot *TP = new TriggerPlot((int)snrthr.size(),tfile_pat,"",style);
  if(!TP->GetStatus()||!TP->GetSegments()->GetLiveTime()) return 2;

  // use trigger time range
  if(gps_start<0) gps_start=TP->GetTimeMin();
  if(gps_end<0)   gps_end=TP->GetTimeMax();

  // clusterize
  TP->SetClusterizeDt(cluster_dt);
  if(usecluster) TP->Clusterize(1);

  // SNR max
  double snrmax;
  if(TP->GetCollectionSelection(0)->GetSNRMax()<50.0)         snrmax=50.0;
  else if(TP->GetCollectionSelection(0)->GetSNRMax()<100.0)   snrmax=100.0;
  else if(TP->GetCollectionSelection(0)->GetSNRMax()<1000.0)  snrmax=1000.0;
  else if(TP->GetCollectionSelection(0)->GetSNRMax()<10000.0) snrmax=10000.0;
  else snrmax=100000.0;

  // time binning
  int ntbins;
  if(gps_end-gps_start<=3600)        ntbins = (gps_end-gps_start)/60+1;
  else if(gps_end-gps_start<=100000) ntbins = (gps_end-gps_start)/600+1;
  else                               ntbins = (gps_end-gps_start)/3600+1;


  // Apply plot selections
  for(int s=0; s<(int)snrthr.size(); s++){

    // snr thresholds
    TP->GetCollectionSelection(s)->SetSNRRange(snrthr[s],snrmax);
    tmpstream<<"SNR \\ge "<<fixed<<setprecision(1)<<snrthr[s];
    TP->SetCollectionLegend(s,tmpstream.str());
    tmpstream.str(""); tmpstream.clear();

    // time range
    TP->GetCollectionSelection(s)->SetTimeRange(gps_start, gps_end);

    // time resolution (FIXME)
    TP->GetCollectionSelection(s)->SetTimeResolution(ntbins);

    // use clusters
    if(usecluster) TP->SetCollectionUseClusters(s,1);

    // time format CHECKME: must be user-defined
    TP->SetDateFormat(usedate);

    // set collection marker
    TP->SetCollectionMarker(s,20, TMath::Max(0.3,(double)(s+1)/(double)(snrthr.size())));
		       
    // set collection color
    TP->SetCollectionColor(s,TP->GetColorPalette((int)((double)(s+1)/(double)(snrthr.size())*(double)(TP->GetNumberOfColors()))-2));
      
  }

  // make collections
  TP->MakeCollections();

  TP->PrintPlot("snr");
  TP->DrawLegend();
  TP->Print(outdir+"/snr."+outformat);

  TP->PrintCollectionPlot("rate");
  TP->DrawLegend();
  TP->Print(outdir+"/rate."+outformat);

  TP->PrintCollectionPlot("freqtime");
  TP->DrawLegend();
  TP->Print(outdir+"/freqtime."+outformat);

  TP->PrintCollectionPlot("snrtime");
  TP->DrawLegend();
  TP->Print(outdir+"/snrtime."+outformat);

  TP->PrintCollectionPlot("snrfreq");
  TP->DrawLegend();
  TP->Print(outdir+"/snrfreq."+outformat);

  delete TP;

  /*
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

  cout<<"omiplot: print SNR plot"<<endl;
  triggers->PrintPlot("snr");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_snr.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_snr_"+sstart+"_"+sstop+".gif");

  cout<<"omiplot: print rate plot"<<endl;
  triggers->PrintCollectionPlot("rate");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_rate.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_rate_"+sstart+"_"+sstop+".gif");

  cout<<"omiplot: print frequency plot"<<endl;
  triggers->PrintCollectionPlot("frequency");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_frequency.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_frequency_"+sstart+"_"+sstop+".gif");

  cout<<"omiplot: print frequency vs. time plot"<<endl;
  triggers->PrintCollectionPlot("freqtime");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_freqtime.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_freqtime_"+sstart+"_"+sstop+".gif");
 
  cout<<"omiplot: print SNR vs. time plot"<<endl;
  triggers->PrintCollectionPlot("snrtime");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_snrtime.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_snrtime_"+sstart+"_"+sstop+".gif");

  cout<<"omiplot: print SNR vs. frequency plot"<<endl;
  triggers->PrintCollectionPlot("snrfreq");
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_snrfreq.gif");
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_snrfreq_"+sstart+"_"+sstop+".gif");
   
  cout<<"omiplot: print plot panel"<<endl;
  triggers->PrintCollectionPanel();
  triggers->DrawLegend();
  if(filename.compare("")) triggers->Print(outdir+"/"+filename+"_panel.gif",2);
  else triggers->Print(outdir+"/"+triggers->GetStreamName()+"_panel_"+sstart+"_"+sstop+".gif",2);
 
  cout<<"omiplot: print loudest event map"<<endl;
  string fn;
  if(filename.compare("")) fn=outdir+"/"+filename+"_loudest.gif";
  else fn=outdir+"/"+triggers->GetStreamName()+"_loudest_"+sstart+"_"+sstop+".gif";
  triggers->PrintLoudestEventMap(fn);
 
  // cleaning
  delete triggers;
  */  
  return 0;
}

