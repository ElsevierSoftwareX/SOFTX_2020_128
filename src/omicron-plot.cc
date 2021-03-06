/**
 * @file 
 * @brief Program to plot omicron triggers.
 * @details `omicron-plot` is a command line program to plot omicron triggers from trigger files.
 * The program must be given a minimum set of options:
 * @verbatim
omicron-plot channel=[channel name] gps-start=[GPS start] gps-end=[GPS end] style=STANDARD
 @endverbatim
 * This command plots omicron triggers of a given channel between 2 GPS times.
 * This command assumes that omicron trigger root files are saved in a standard place pointed by the environement variable `$OMICRON_TRIGGERS`.
 *
 * One can also plot triggers from a list of root files with:
 * @verbatim
 omicron-plot file=[trigger file pattern] gps-start=[GPS start] gps-end=[GPS end] style=STANDARD
 @endverbatim
 * where `[trigger file pattern]` can contain wild cards. For example: `file="/path1/to/triggers/*.root /path2/to/triggers/*.root"`.
 *
 * The `omicron-plot` command comes with many additional options. Type `omicron-plot` to get the full list of options. In particular, triggers can be filtered in frequency, snr and so on.
 * @snippet this omicron-plot-usage 
 *
 * \anchor omicron_plot_freqtime
 * \image html plot_V1-LSC_DARM-1217376018-187200_freqtime.png "Triggers are plotted in the time-frequency plane." width=700
 * \anchor omicron_plot_snrtime
 * \image html plot_V1-LSC_DARM-1217376018-187200_snrtime.png "Triggers are plotted in the time-SNR plane." width=700
 * \anchor omicron_plot_snrfreq
 * \image html plot_V1-LSC_DARM-1217376018-187200_snrfreq.png "Triggers are plotted in the SNR-frequency plane." width=700
 * \anchor omicron_plot_rate
 * \image html plot_V1-LSC_DARM-1217376018-187200_rate.png "Trigger rates are plotted as a function of time." width=700
 * \anchor omicron_plot_snr
 * \image html plot_V1-LSC_DARM-1217376018-187200_snr.png "Trigger SNR distribution is plotted." width=700
 *
 *
 * @author Florent Robinet - <a href="mailto:florent.robinet@ijclab.in2p3.fr">florent.robinet@ijclab.in2p3.fr</a>
 */
#include "TriggerPlot.h"
#include "OmicronUtils.h"
#include "Oconfig.h"

using namespace std;

/**
 * @brief Print the program usage message.
 */
void PrintUsage(void){
  //! [omicron-plot-usage]
  cerr<<endl;
  cerr<<"Usage:"<<endl;
  cerr<<endl;
  cerr<<"omicron-plot.exe channel=[channel name] \\"<<endl;
  cerr<<"                 file=[trigger file pattern] \\"<<endl;
  cerr<<"                 gps-start=[GPS start] \\"<<endl;
  cerr<<"                 gps-end=[GPS end] \\"<<endl;
  cerr<<"                 segments=[segment file] \\"<<endl;
  cerr<<"                 snr-thresholds=[list of SNR thresholds] \\"<<endl;
  cerr<<"                 freq-min=[minimum frequency] \\"<<endl;
  cerr<<"                 freq-max=[maximum frequency] \\"<<endl;
  cerr<<"                 use-cluster=[cluster flag] \\"<<endl;
  cerr<<"                 cluster-dt=[cluster time window] \\"<<endl;
  cerr<<"                 use-date=[date flag] \\"<<endl;
  cerr<<"                 outdir=[output directory] \\"<<endl;
  cerr<<"                 outformat=[output file format] \\"<<endl;
  cerr<<"                 file-prefix=[file prefix] \\"<<endl;
  cerr<<"                 file-name=[file name] \\"<<endl;
  cerr<<"                 style=[style] \\"<<endl;
  cerr<<"                 drawtimeline=[GPS time] \\"<<endl;
  cerr<<"                 plot-width=[plot width] \\"<<endl;
  cerr<<"                 plot-height=[plot height]"<<endl;
  cerr<<"                 plot-star=[plot star flag]"<<endl;
  cerr<<endl;
  cerr<<"[channel name]            channel name used to retrieve centralized Omicron triggers"<<endl;
  cerr<<"[trigger file pattern]    file pattern to ROOT trigger files (GWOLLUM convention)"<<endl;
  cerr<<"[GPS start]               starting GPS time (integer only)"<<endl;
  cerr<<"[GPS end]                 stopping GPS time (integer only)"<<endl;
  cerr<<"[segment file]            path to a file with a list of segments"<<endl;
  cerr<<"[list of SNR thresholds]  list of SNR thresholds. By default, snr-thresholds=\"5;8;10;20\""<<endl;
  cerr<<"[minimum frequency]       minimum frequency value [Hz]"<<endl;
  cerr<<"[maximum frequency]       maximum frequency value [Hz]"<<endl;
  cerr<<"[cluster tag]             1 = use clusters, 0 = use triggers. By default, use-cluster=1"<<endl;
  cerr<<"[cluster time window]     cluster time window [s]. By default, cluster-dt=0.1"<<endl;
  cerr<<"[date tag]                1 = use date, 0 = use GPS time. By default, use-date=1"<<endl;
  cerr<<"[output directory]        output directory. By default, outdir=."<<endl;
  cerr<<"[output file format]      output file format. By default, outformat=png"<<endl;
  cerr<<"[file prefix]             file name prefix. By default, file-prefix=plot"<<endl;
  cerr<<"[file name]               file name. By default, file-name=default"<<endl;
  cerr<<"[style]                   GWOLLUM-supported style. By default, style=GWOLLUM"<<endl;
  cerr<<"[GPS time]                GPS time at which drawing a vertical line (time plots only)"<<endl;
  cerr<<"[plot width]              plot width in pixels. By default, plot-width=850"<<endl;
  cerr<<"[plot height]             plot height in pixels. By default, plot-height=500"<<endl;
  cerr<<"[plot star flag]          plot a star on the loudest event (=1, default). =0: do not plot the star"<<endl;
  cerr<<endl;
  //! [omicron-plot-usage]
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
  string segfile=""; // segment file
  double freqmin=-1; // freq min
  double freqmax=1e20; // freq max
  string style="GWOLLUM"; // style
  bool usecluster=true; // use cluster
  double cluster_dt=0.1; // cluster time window
  bool usedate=true; // use date
  string outdir="."; // output directory
  string outformat="png"; // file format
  string fileprefix="plot"; // file name prefix
  string filename="default"; // file name
  double vline=-1.0;// do not draw a line
  int plot_w=850;  // plot width
  int plot_h=500;  // plot height
  int plot_star=1; // plot star

  // loop over arguments
  vector <string> sarg;
  for(int a=1; a<argc; a++){
    sarg=SplitString((string)argv[a],'=');
    if(sarg.size()!=2) continue;
    if(!sarg[0].compare("channel"))        chname=(string)sarg[1];
    if(!sarg[0].compare("file"))           tfile_pat=sarg[1];
    if(!sarg[0].compare("gps-start"))      gps_start=atoi(sarg[1].c_str());
    if(!sarg[0].compare("gps-end"))        gps_end=atoi(sarg[1].c_str());
    if(!sarg[0].compare("segments"))       segfile=(string)sarg[1];
    if(!sarg[0].compare("snr-thresholds")) snrthrs=(string)sarg[1];
    if(!sarg[0].compare("freq-min"))       freqmin=atof(sarg[1].c_str());
    if(!sarg[0].compare("freq-max"))       freqmax=atof(sarg[1].c_str());
    if(!sarg[0].compare("style"))          style=(string)sarg[1];
    if(!sarg[0].compare("use-cluster"))    usecluster=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("cluster-dt"))     cluster_dt=atof(sarg[1].c_str());
    if(!sarg[0].compare("use-date"))       usedate=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("outdir"))         outdir=(string)sarg[1];
    if(!sarg[0].compare("outformat"))      outformat=(string)sarg[1];
    if(!sarg[0].compare("file-prefix"))    fileprefix=(string)sarg[1];
    if(!sarg[0].compare("file-name"))      filename=(string)sarg[1];
    if(!sarg[0].compare("drawtimeline"))   vline=atof(sarg[1].c_str());
    if(!sarg[0].compare("plot-width"))     plot_w=atoi(sarg[1].c_str());
    if(!sarg[0].compare("plot-height"))    plot_h=atoi(sarg[1].c_str());
    if(!sarg[0].compare("plot-star"))      plot_star=atoi(sarg[1].c_str());
  }

  // selection segment file
  Segments *SelSeg = NULL;
  if(segfile.compare("")){
    
    SelSeg = new Segments(segfile);
    if(!SelSeg->GetStatus() || !SelSeg->GetNsegments()){
      cerr<<"Invalid segment file: "<<segfile<<endl;
      cerr<<"Type omicron-plot.exe for help"<<endl;
      return 1;
    }

    // trim segments
    if(gps_start>0) SelSeg->Intersect(gps_start,20000000000);
    else gps_start=(int)SelSeg->GetStart(0);
    if(gps_end>0) SelSeg->Intersect(700000000, gps_end);
    else gps_end=(int)ceil(SelSeg->GetEnd(SelSeg->GetNsegments()-1));
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

  stringstream tmpstream;

  // check file pattern
  if(!tfile_pat.compare("")){
    cerr<<"No trigger files"<<endl;

    // make output name
    tmpstream<<"no_trigger-"<<gps_start<<"-"<<gps_end-gps_start;
    if(!filename.compare("default")) filename=tmpstream.str();

    // print no-trigger plots
    GwollumPlot *GP = new GwollumPlot("notrigger",style);
    GP->ResizePlot(plot_w,plot_h);
    GP->AddText("NO TRIGGER", 0.1, 0.1, 0.2);
    GP->Print(outdir+"/"+fileprefix+"_"+filename+"_snr."+outformat);
    GP->Print(outdir+"/"+fileprefix+"_"+filename+"_rate."+outformat);
    GP->Print(outdir+"/"+fileprefix+"_"+filename+"_freqtime."+outformat);
    GP->Print(outdir+"/"+fileprefix+"_"+filename+"_snrtime."+outformat);
    GP->Print(outdir+"/"+fileprefix+"_"+filename+"_snrfreq."+outformat);
    delete GP;
    return 2;
  }
  
  // snr thresholds
  sarg=SplitString(snrthrs,';');
  vector <double> snrthr;
  for(int s=0; s<(int)sarg.size(); s++) snrthr.push_back(atof(sarg[s].c_str()));

  // triggers
  TriggerPlot *TP = new TriggerPlot((int)snrthr.size(),tfile_pat,"",style);
  if(!TP->GetStatus()||!TP->GetSegments()->GetLiveTime()) return 2;
  TP->ResizePlot(plot_w,plot_h);

  // use trigger time range
  if(gps_start<0) gps_start=TP->GetTimeMin();
  if(gps_end<0)   gps_end=TP->GetTimeMax();

  // clusterize
  TP->SetClusterizeDt(cluster_dt);
  if(usecluster) TP->Clusterize(1);

  // star for the loudest event
  TP->PlotLoudestEvent((bool)plot_star);

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

  // selection segments
  TP->SetSelectionSegments(SelSeg);
  
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

    // time format
    TP->SetDateFormat(usedate);

    // set frequency range
    if(TP->GetFrequencyMin()<freqmin) TP->GetCollectionSelection(s)->SetFrequencyMin(freqmin);
    if(TP->GetFrequencyMax()>freqmax) TP->GetCollectionSelection(s)->SetFrequencyMax(freqmax);

    // set collection marker
    TP->SetCollectionMarker(s,20, TMath::Max(0.3,(double)(s+1)/(double)(snrthr.size())));
		       
    // set collection color
    TP->SetCollectionColor(s,TP->GetColorPalette((int)((double)(s+1)/(double)(snrthr.size())*(double)(TP->GetNumberOfColors()))-2));
  }

  // force minimum marker size for first collection
  if(gps_end-gps_start>40000) TP->SetCollectionMarker(0,1,1);

  // make collections
  TP->MakeCollections();

  // vertical line
  TLine *tvline = new TLine(vline,0,vline,1);
  tvline->SetLineColor(2);
  
  // make output name
  tmpstream<<TP->GetNamePrefix()<<"-"<<TP->GetNameSuffix()<<"-"<<gps_start<<"-"<<gps_end-gps_start;
  if(!filename.compare("default")) filename=tmpstream.str();

  // print plots
  TP->PrintPlot("snr");
  TP->DrawLegend();
  TP->Print(outdir+"/"+fileprefix+"_"+filename+"_snr."+outformat);

  TP->PrintCollectionPlot("rate");
  tvline->SetY1(TP->GetYmin("rate",0));
  tvline->SetY2(TP->GetYmax("rate",0));
  TP->Draw(tvline,"same");
  TP->DrawLegend();
  TP->Print(outdir+"/"+fileprefix+"_"+filename+"_rate."+outformat);

  TP->PrintCollectionPlot("freqtime");
  tvline->SetY1(TP->GetYmin("freqtime",0));
  tvline->SetY2(TP->GetYmax("freqtime",0));
  TP->Draw(tvline,"same");
  TP->DrawLegend();
  TP->Print(outdir+"/"+fileprefix+"_"+filename+"_freqtime."+outformat);

  TP->PrintCollectionPlot("snrtime");
  tvline->SetY1(TP->GetYmin("snrtime",0));
  tvline->SetY2(TP->GetYmax("snrtime",0));
  TP->Draw(tvline,"same");
  TP->DrawLegend();
  TP->Print(outdir+"/"+fileprefix+"_"+filename+"_snrtime."+outformat);

  TP->PrintCollectionPlot("snrfreq");
  TP->DrawLegend();
  TP->Print(outdir+"/"+fileprefix+"_"+filename+"_snrfreq."+outformat);
  
  delete TP;
  delete tvline;
  if(SelSeg!=NULL) delete SelSeg;
  
  return 0;
}

