//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string.h>
#include <stdio.h>
#include "GwollumPlot.h"
#include "Omicron.h"
#include "IO.h"

using namespace std;

int main (int argc, char* argv[]){
  gErrorIgnoreLevel = 3000;

  // check the argument
  if(argc!=4){
    cerr<<argv[0]<<" usage:"<<endl<<endl; 
    cerr<<argv[0]<<" [channel file] [option file] [central time]"<<endl; 
    return 1;
  }
  
  ///////////////////////////////////////////////////////////////////////////////
  /////////                        COMMAND LINE                         /////////
  ///////////////////////////////////////////////////////////////////////////////
  // command line parameters
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
  /////////                        OPTION FILE                          /////////
  ///////////////////////////////////////////////////////////////////////////////
  // Parse option file
  IO *io = new IO(optionfile.c_str());

  //**** output directory
  string outdir;
  if(io->GetOpt("OUTPUT","DIRECTORY", outdir)){
    if(!IsDirectory(outdir)){
      cerr<<"Omiscan ERROR: output directory '"<<outdir<<"' does not exist"<<endl;
      return 3;
    }
  }
  else outdir=".";// by default current directory
  system(("mkdir -p "+outdir+"/plots").c_str());
  
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
  vector <int> windows;
  if(io->GetOpt("PARAMETER","WINDOW", windows)){
    if(!windows.size()){
      cerr<<"Omiscan ERROR: the window option is not correct"<<endl;
      return 3;
    }
  }
  else{// default windows
    windows.clear();
    windows.push_back(2); windows.push_back(8); windows.push_back(32);
  }

  //**** frequency ranges
  vector <int> fsample_cat;
  vector <double> frange;
  if(io->GetOpt("PARAMETER","SAMPLING_CAT", fsample_cat)){
    if(!fsample_cat.size()){
      cerr<<"Omiscan ERROR: the sampling frequency category option is not correct"<<endl;
      return 3;
    }
    if(io->GetOpt("PARAMETER","FRANGE", frange)){
      if(frange.size()!=2*(fsample_cat.size()+1)){
	cerr<<"Omiscan ERROR: the frequency ranges option is incompatible with the sampling categories"<<endl;
	return 3;
      }
      for(int f=0; f<(int)frange.size(); f+=2){
	if(frange[f]>=frange[f+1]){
	  cerr<<"Omiscan ERROR: the frequency range "<<frange[f]<<" "<<frange[f+1]<<" is incorrect"<<endl;
	  return 3;
	}
      }
    }
    else{
      cerr<<"Omiscan ERROR: the frequency ranges option is missing"<<endl;
      return 3;
    }
  }
  else{// default
    if(io->GetOpt("PARAMETER","FRANGE", frange)){
      fsample_cat.clear();
      if(frange.size()!=2||frange[0]>=frange[1]){
	cerr<<"Omiscan ERROR: the frequency range "<<frange[0]<<" "<<frange[1]<<" is incorrect"<<endl;
	return 3;
      }
    }
    else{
      fsample_cat.clear(); fsample_cat.push_back(1025);
      frange.clear(); 
      frange.push_back(0.1); frange.push_back(64);
      frange.push_back(16); frange.push_back(1024);
    }
  }

  //**** SNR threshold
  double snrthr;
  if(io->GetOpt("PARAMETER","SNR_THRESHOLD", snrthr)){
    if(snrthr<0){
      cerr<<"Omiscan ERROR: the SNR threshold must be >=0"<<endl;
      return 3;
    }
  }
  else{// default SNR threshold
    snrthr=5;
  }
  
  //**** verbosity
  int verbose;
  if(io->GetOpt("VERBOSITY","LEVEL", verbose)){
    if(verbose<0){
      cerr<<"Omiscan ERROR: the verbosity option is not correct"<<endl;
      return 3;
    }
  }
  else verbose=0;// default
  
  delete io;
  
  if(verbose){
    cout<<"Omiscan: list of input options:"<<endl;
    cout<<"         output directory = "<<outdir<<endl;
    cout<<"         FFL file         = "<<fflfile<<endl;
    cout<<"         verbosity level  = "<<verbose<<endl;
    cout<<"         time ranges      = ";
    for(int w=0; w<(int)windows.size(); w++) cout<<windows[w]<<"s ";
    cout<<endl;
    cout<<"         SNR threshold    = "<<snrthr<<endl;
  }
  
  ///////////////////////////////////////////////////////////////////////////////
  /////////                        CHANNEL FILE                         /////////
  ///////////////////////////////////////////////////////////////////////////////
  cout<<"\nOmiscan: parsing channel file..."<<endl;
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
  if(verbose>1) cout<<"Omiscan: number of channels to process = "<<channels.size()<<endl;
  if(verbose>2){
    for(int c=0; c<(int)channels.size(); c++) cout<<channels[c]<<endl;
  }
  
  ///////////////////////////////////////////////////////////////////////////////
  /////////                       INITIALIZATION                        /////////
  ///////////////////////////////////////////////////////////////////////////////
  
  // FFL file  
  if(verbose) cout<<"Omiscan: load FFL file "<<fflfile<<"..."<<endl;
  FrFile *frfile = FrFileINew((char*)fflfile.c_str());
  
  // data vector
  FrVect *chanvect = NULL;

  // gwollum plots
  GwollumPlot *GPlot = new GwollumPlot ("omiscan");
  
  // declarations
  int timerange, pad;        // analysis time window and padding
  int start, stop;           // analysis starting/ending time
  Segments *segment = new Segments();// analysis segment
  int sampling, sampling_new;// sampling frequency before-after
  int psdsize;               // psd size
  int powerof2;              // closest power of 2 below
  double fmin, fmax;         // search frequency range
  Odata *data;               // data structure
  Otile *tiles;              // tiling structure
  double *c_data[2];         // condition data
  double *psd;               // psd vector - DO NOT DELETE!
  double *f_bin;             // frequency bins
  int nfbins;                // number of frequency bins
  TH2D *qmap;                // Q map - DO NOT DELETE!
  TGraph *graph;             // graph plots
  TH2D **map = new TH2D * [(int)windows.size()];// overall maps
  int bin_start, bin_stop, bin_start_t, bin_stop_t, bin_start_f, bin_stop_f, dummy;// bin indexes
  double content;            // tile content
  ostringstream tmpstream;   // stream
  int loudest_bin, loudest_bin_t, loudest_bin_f, loudest_bin_s;// loudest bin
  double loudest_t, loudest_f, loudest_s;// loudest tile
  
  ofstream summary((outdir+"/summary.txt").c_str());// outfile stream
  summary<<"# channel name"<<endl;
  summary<<"# native sampling [Hz]"<<endl;
  summary<<"# effective sampling [Hz]"<<endl;
  summary<<"# Loudest event GPS [s]"<<endl;
  summary<<"# Loudest event frequency [Hz]"<<endl;
  summary<<"# Loudest event SNR"<<endl;
  summary<<"# Number of Q planes"<<endl;
  summary<<"# List of Q values"<<endl;

  ///////////////////////////////////////////////////////////////////////////////
  /////////                      PROCESS CHANNELS                       /////////
  ///////////////////////////////////////////////////////////////////////////////

  // loop over channels
  for(int c=0; c<(int)channels.size(); c++){

    cout<<"\nOmiscan: process channel "<<channels[c]<<"..."<<endl;
    
    // get data
    if(verbose) cout<<"         Optimizing parameters..."<<endl;
    chanvect = NULL;
    chanvect = FrFileIGetVectDN(frfile,(char*)(channels[c].c_str()),(int)floor(gps),1);
 
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
    if(chanvect->dataD[0]==0.0&&chanvect->dataD[chanvect->nData-1]==0.0){
      cout<<"Omiscan WARNING: zero data --> skip"<<endl;
      FrVectFree(chanvect);
      continue;
    }

    // native sampling
    sampling=chanvect->nData;

    // get frequency range for this category
    fmin=frange[2*fsample_cat.size()];
    fmax=frange[2*fsample_cat.size()+1];
    for(int f=0; f<(int)fsample_cat.size(); f++){
      if(sampling<fsample_cat[f]){
	fmin=frange[2*f];
	fmax=frange[2*f+1];
	break;
      }
    }
    
    // optimize sampling for this range
    powerof2=(int)floor(log(TMath::Min(sampling,2*(int)fmax))/log(2.0));
    sampling_new=(double)pow(2,powerof2);// optimized sampling for this range
    fmax=TMath::Min(fmax,(double)sampling_new/2.0);

    // print sampling info
    if(verbose) cout<<"           native sampling    = "<<sampling<<"Hz"<<endl;
    if(verbose) cout<<"           re-sampling        = "<<sampling_new<<"Hz"<<endl;
    if(verbose) cout<<"           frequency range    = "<<setprecision(2)<<fmin<<"Hz --> "<<fmax<<"Hz"<<endl;
    if(fmax<=fmin){
      cout<<"Omiscan WARNING: frequency range is not adapted --> skip"<<endl;
      FrVectFree(chanvect);
      continue;
    }

    // optimized analysis time window
    pad=TMath::Max((int)ceil(8.0/fmin),4);// padding
    timerange=windows[(int)windows.size()-1]+2*pad;
    powerof2=(int)ceil(log(timerange)/log(2.0));
    timerange=(int)pow(2,powerof2);
    start=(int)floor(gps)-(double)timerange/2.0;
    stop=(int)floor(gps)+(double)timerange/2.0;
    if(verbose){
      cout<<"           analysis timerange = "<<timerange<<"s"<<endl;
      cout<<"           GPS start          = "<<start<<endl;
      cout<<"           GPS stop           = "<<stop<<endl;
      if(verbose>1) cout<<"           analysis pad       = "<<pad<<"s"<<endl;
    }

    // get data vector
    if(verbose) cout<<"         Loading data vector..."<<endl;
    FrVectFree(chanvect); chanvect = NULL;
    chanvect = FrFileIGetVectDN(frfile,(char*)(channels[c].c_str()),start,timerange);
 
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
    if(chanvect->dataD[0]==0.0&&chanvect->dataD[chanvect->nData-1]==0.0){
      cout<<"Omiscan WARNING: zero data --> skip"<<endl;
      FrVectFree(chanvect);
      continue;
    }

    // segment structure
    if(verbose) cout<<"         Create analysis window..."<<endl;
    if(!segment->Reset()||!segment->AddSegment(start,stop)){
      cerr<<"Omiscan WARNING: problem with generating segment window --> skip"<<endl;
      FrVectFree(chanvect);
      continue;
    }


    // init and load data
    if(verbose) cout<<"         Create Odata object..."<<endl;
    data = new Odata("ONLINE", channels[c], sampling, sampling_new, segment, timerange, timerange, 2*pad, fmin);
    if(!data->ReadVect(chanvect)){
      cout<<"Omiscan WARNING: Odata object is corrupted --> skip"<<endl;
      FrVectFree(chanvect);
      delete data;
      continue;
    }
    FrVectFree(chanvect);

    // init conditioned data containers
    c_data[0] = new double [timerange*sampling_new/2]; // real part
    c_data[1] = new double [timerange*sampling_new/2]; // imaginary part

    // condition data
    if(verbose) cout<<"         Get conditionned data..."<<endl;
    if(!data->GetConditionedData(0,c_data[0],c_data[1],&psd,psdsize)){
      cout<<"Omiscan WARNING: conditionned data are corrupted --> skip"<<endl;
      delete data;
      delete c_data[0];
      delete c_data[1];
      continue;
    }
    
    // make tiling
    if(verbose) cout<<"         Make tiling..."<<endl;
    tiles = new Otile(timerange, pad, 3, 141, fmin, fmax, sampling_new, 0.2);
    if(!tiles->SetPowerSpectrum(psd,psdsize)){
      cout<<"Omiscan WARNING: cannot normalize tiling --> skip"<<endl;
      delete tiles;
      delete data;
      delete c_data[0];
      delete c_data[1];
      continue;
    }


    //*****  TIME SERIES plots  ***************************************
    if(verbose) cout<<"         Make time-series..."<<endl;
    graph = data->GetTimeSeries();
    tmpstream<<channels[c]<<" raw time series";
    graph->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();

    // offset to be centered on 0
    for(int p=0; p<graph->GetN(); p++) graph->GetX()[p]-=gps;

    GPlot->SetLogy(0);
    GPlot->SetLogx(0);
    GPlot->SetGridx(1);
    GPlot->SetGridy(1);
    GPlot->Draw(graph,"APL");
    
    // cosmetics
    graph->GetXaxis()->SetTimeOffset(-gps);
    graph->GetXaxis()->SetTitle("time [s]");
    graph->GetXaxis()->SetTitleOffset(1.1);
    graph->GetXaxis()->SetLabelSize(0.045);
    graph->GetYaxis()->SetLabelSize(0.045);
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitle("Amplitude");
    graph->GetXaxis()->SetNoExponent();
    graph->SetLineWidth(1);
    graph->SetLineColor(7);
    graph->SetMarkerColor(7);
    
    // window resize
    for(int w=(int)windows.size()-1; w>=0; w--){
      
      // zoom
      graph->GetXaxis()->SetLimits(-(double)windows[w]/2.0,+(double)windows[w]/2.0);
 
      // save plot
      tmpstream<<outdir<<"/plots/"<<channels[c]<<"_raw_dt"<<windows[w]<<".gif";
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();
      //tmpstream<<outdir<<"/plots/th_"<<channels[c]<<"_raw_dt"<<windows[w]<<".gif";
      //GPlot->Print(tmpstream.str(),0.5);
      //tmpstream.str(""); tmpstream.clear();
    }
    delete graph;

    //*****  PSD plot  ***********************************************
    if(verbose) cout<<"         Make PSD..."<<endl;
    graph = new TGraph();
    tmpstream<<channels[c]<<" amplitude spectral density";
    graph->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
    for(int i=0; i<psdsize; i++)
      graph->SetPoint(i,i*(double)sampling_new/2.0/(double)psdsize,sqrt(psd[i]/(double)sampling_new/(double)psdsize/2.0));
    GPlot->SetLogy(1);
    GPlot->SetLogx(1);
    GPlot->SetGridx(1);
    GPlot->SetGridy(1);
    GPlot->Draw(graph,"APL");

    // cosmetics
    graph->GetXaxis()->SetTitle("Frequency [Hz]");
    graph->GetXaxis()->SetLimits(fmin,fmax);
    graph->SetLineWidth(2);
    graph->SetLineColor(7);
    graph->SetMarkerColor(7);
	
    // save PSD plot
    tmpstream<<outdir<<"/plots/"<<channels[c]<<"_psd.gif";
    GPlot->Print(tmpstream.str());
    tmpstream.str(""); tmpstream.clear();
    //tmpstream<<outdir<<"/plots/th_"<<channels[c]<<"_psd.gif";
    //GPlot->Print(tmpstream.str(),0.5);
    //tmpstream.str(""); tmpstream.clear();
    delete graph;

    // frequency binning
    if(fmin>=1) nfbins = fmax-fmin;
    else nfbins = 1000;
    f_bin = new double [nfbins+1];
    for(int f=0; f<nfbins+1; f++) f_bin[f]=fmin*pow(10,f*log10(fmax/fmin)/nfbins);

    // create overall map
    if(verbose) cout<<"         Create overall maps..."<<endl;
    for(int w=0; w<(int)windows.size(); w++){
      tmpstream<<"map_"<<windows[w];
      map[w] = new TH2D(tmpstream.str().c_str(),tmpstream.str().c_str(),1000,-(double)windows[w]/2.0,+(double)windows[w]/2.0,nfbins,f_bin);
      tmpstream.str(""); tmpstream.clear();
    }
    delete f_bin;

    // get map for each plane
    if(verbose) cout<<"         Populate Q maps..."<<endl;
    for(int q=0; q<tiles->GetNQPlanes(); q++){
      if(verbose>2) cout<<"           make map Q"<<q<<"..."<<endl;
      qmap=tiles->GetMap(q,c_data[0],c_data[1]);
      if(qmap==NULL){
	cout<<"Omiscan WARNING: missing map --> skip"<<endl;
	continue;
      }

      // fill overall maps
      if(verbose>2) cout<<"           update overall maps with map Q"<<q<<"..."<<endl;
      for(int bt=1; bt<=qmap->GetNbinsX(); bt++){
	for(int bf=1; bf<=qmap->GetNbinsY(); bf++){
	  content=qmap->GetBinContent(bt,bf);
	  for(int w=0; w<(int)windows.size(); w++){
	    if(qmap->GetXaxis()->GetBinLowEdge(bt)<-(double)windows[w]/2.0) continue;
	    if(qmap->GetXaxis()->GetBinUpEdge(bt)>(double)windows[w]/2.0) continue;
	    if(qmap->GetYaxis()->GetBinLowEdge(bf)<map[w]->GetYaxis()->GetBinLowEdge(1)) continue;
	    if(qmap->GetYaxis()->GetBinUpEdge(bf)>map[w]->GetYaxis()->GetBinUpEdge(nfbins)) continue;
	    bin_start = map[w]->FindBin(qmap->GetXaxis()->GetBinLowEdge(bt),qmap->GetYaxis()->GetBinLowEdge(bf));
	    bin_stop = map[w]->FindBin(qmap->GetXaxis()->GetBinUpEdge(bt),qmap->GetYaxis()->GetBinUpEdge(bf));
	    map[w]->GetBinXYZ(bin_start,bin_start_t,bin_start_f,dummy);
	    map[w]->GetBinXYZ(bin_stop,bin_stop_t,bin_stop_f,dummy);
	    for(int bbt=bin_start_t; bbt<=bin_stop_t; bbt++){// time-sweep the tile
	      for(int bbf=bin_start_f; bbf<=bin_stop_f; bbf++){// freq-sweep the tile
		if(content>map[w]->GetBinContent(bbt,bbf)) map[w]->SetBinContent(bbt,bbf,content);
	      }
	    }
	  }
	}
      }
      
      // draw q map
      if(verbose>2) cout<<"           draw map Q"<<q<<"..."<<endl;
      GPlot->SetLogx(0);
      GPlot->SetLogy(1);
      GPlot->SetLogz(1);
      GPlot->SetGridx(1);
      GPlot->SetGridy(1);
      GPlot->Draw(qmap,"COLZ");

      // cosmetics
      qmap->GetXaxis()->SetTitle("Time [s]");
      qmap->GetYaxis()->SetTitle("Frequency [Hz]");
      qmap->GetZaxis()->SetTitle("SNR");
      qmap->GetZaxis()->SetRangeUser(1,50);

      // window resize for Qmap
      for(int w=(int)windows.size()-1; w>=0; w--){
      
	// zoom
	qmap->GetXaxis()->SetRangeUser(-(double)windows[w]/2.0,(double)windows[w]/2.0);
 
	// get loudest tile
	loudest_bin=qmap->GetMaximumBin();
	qmap->GetBinXYZ(loudest_bin,loudest_bin_t,loudest_bin_f,loudest_bin_s);
	loudest_t=qmap->GetXaxis()->GetBinCenter(loudest_bin_t);
	loudest_f=qmap->GetYaxis()->GetBinCenter(loudest_bin_f);
	loudest_s=qmap->GetBinContent(loudest_bin);

	// plot title
	tmpstream<<"Q = "<<tiles->GetQ(q)<<" "<<channels[c]<<" loudest tile: GPS = "<<setprecision(12)<<loudest_t+gps<<" f = "<<setprecision(5)<<loudest_f<<"Hz SNR = "<<loudest_s;
	qmap->SetTitle(tmpstream.str().c_str());
	tmpstream.str(""); tmpstream.clear();

	// save qmap
	tmpstream<<outdir<<"/plots/"<<channels[c]<<"_map_Q"<<q<<"_dt"<<windows[w]<<".gif";
	GPlot->Print(tmpstream.str());
	tmpstream.str(""); tmpstream.clear();
	//tmpstream<<outdir<<"/plots/th_"<<channels[c]<<"_map_Q"<<q<<"_dt"<<windows[w]<<".gif";
	//GPlot->Print(tmpstream.str(),0.5);
	//tmpstream.str(""); tmpstream.clear();
      }


    }
    
    delete c_data[0];
    delete c_data[1];
    delete data;

    // draw overall map
    if(verbose) cout<<"         Draw overall maps..."<<endl;
    for(int w=(int)windows.size()-1; w>=0; w--){
      GPlot->SetLogx(0);
      GPlot->SetLogy(1);
      GPlot->SetLogz(1);
      GPlot->SetGridx(1);
      GPlot->SetGridy(1);
      GPlot->Draw(map[w],"COLZ");
       
      // cosmetics
      map[w]->GetXaxis()->SetTitle("Time [s]");
      map[w]->GetYaxis()->SetTitle("Frequency [Hz]");
      map[w]->GetZaxis()->SetTitle("SNR");
      map[w]->GetZaxis()->SetRangeUser(1,50);
      
      // get loudest tile
      loudest_bin=map[w]->GetMaximumBin();
      map[w]->GetBinXYZ(loudest_bin,loudest_bin_t,loudest_bin_f,loudest_bin_s);
      loudest_t=map[w]->GetXaxis()->GetBinCenter(loudest_bin_t);
      loudest_f=map[w]->GetYaxis()->GetBinCenter(loudest_bin_f);
      loudest_s=map[w]->GetBinContent(loudest_bin);
      
      // plot title
      tmpstream<<channels[c]<<" loudest tile: GPS = "<<setprecision(12)<<loudest_t+gps<<" f = "<<setprecision(5)<<loudest_f<<"Hz SNR = "<<loudest_s;
      map[w]->SetTitle(tmpstream.str().c_str());
      tmpstream.str(""); tmpstream.clear();

      // print summary
      // FIXME: the snr threshold is only applied to the summary!
      if(!w && loudest_s>snrthr){
	summary<<channels[c]<<" "<<sampling<<" "<<sampling_new<<" "<<setprecision(12)<<loudest_t+gps<<" "<<setprecision(5)<<loudest_f<<" "<<loudest_s<<" "<<tiles->GetNQPlanes()<<" ";
	for(int q=0; q<tiles->GetNQPlanes(); q++) summary<<setprecision(5)<<tiles->GetQ(q)<<" ";
	summary<<endl;
      }
      
      // save map
      tmpstream<<outdir<<"/plots/"<<channels[c]<<"_map_dt"<<windows[w]<<".gif";
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();
      //tmpstream<<outdir<<"/plots/th_"<<channels[c]<<"_map_dt"<<windows[w]<<".gif";
      //GPlot->Print(tmpstream.str(),0.5);
      //tmpstream.str(""); tmpstream.clear();
      delete map[w];
    }
        
    delete tiles;
  }

  summary<<"*** omiscan done ***"<<endl;
  summary.close();

  // cleaning
  delete segment;
  delete GPlot;
  delete map;

  return 0;
}

