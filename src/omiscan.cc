//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

using namespace std;

int main (int argc, char* argv[]){
  gErrorIgnoreLevel = 3000;

  // check the argument
  if(argc!=3){
    cerr<<argv[0]<<" usage:"<<endl<<endl; 
    cerr<<argv[0]<<" [central time] [option file]"<<endl; 
    return 1;
  }
  
  ///////////////////////////////////////////////////////////////////////////////
  /////////                        COMMAND LINE                         /////////
  ///////////////////////////////////////////////////////////////////////////////
  // command line parameters
  double gps=atof(argv[1]);
  string optionfile=(string)argv[2];

  // check command line parameters
  if(!IsTextFile(optionfile)){
    cerr<<"Omiscan ERROR: option file '"<<optionfile<<"' cannot be found"<<endl;
    return 2;
  }
  if(gps<700000000){
    cerr<<"Omiscan ERROR: GPS time '"<<gps<<"' is not reasonable"<<endl;
    return 2;
  }

  // init Omicron
  Omicron *O = new Omicron(optionfile);
  if(!O->GetStatus()){
    cerr<<"Omiscan ERROR: input options are invalid"<<endl;
    delete O;
    return 3;
  }

  // scan
  if(!O->Scan(gps)){
    cerr<<"Omiscan ERROR: scan failed"<<endl;
    delete O;
    return 4;
  }


  delete O;
  return 0;
}
  
  /*
  ///////////////////////////////////////////////////////////////////////////////
  /////////                        OPTION FILE                          /////////
  ///////////////////////////////////////////////////////////////////////////////


  // declarations
  int timerange, pad;        // analysis time window and padding
  int start, stop;           // analysis starting/ending time
  double toffset;            // time offset from working with int times
  Segments *segment = new Segments();// analysis segment
  Sample *sample = new Sample(4096,0);
  int sampling, sampling_new;// sampling frequency before-after
  int asdsize;               // asd size
  int powerof2;              // closest power of 2 below
  double fmin, fmax;         // search frequency range
  double *chanvect;          // data vector
  double *TukeyWindow;       // tukey window
  Otile *tiles;              // tiling structure
  double *c_data[2];         // condition data
  double *asd;               // asd vector - DO NOT DELETE!
  double *f_bin;             // frequency bins
  int nfbins;                // number of frequency bins
  TH2D **qmap;               // Q maps - DO NOT DELETE!
  TGraph *graph;             // graph plots
  TH2D **map = new TH2D * [(int)windows.size()];// overall maps
  TH1D *projt, *projf;       // Projections
  int bin_start, bin_stop, bin_start_t, bin_stop_t, bin_start_f, bin_stop_f, dummy;// bin indexes
  double content;            // tile content
  ostringstream tmpstream;   // stream
  int loudest_bin, loudest_bin_t, loudest_bin_f, loudest_bin_s;// loudest bin
  double loudest_t, loudest_f, loudest_s, loudest_q;// loudest tile
  TFile *rootfile;           // output rootfile
  
  ofstream summary((outdir+"/summary.txt").c_str());// outfile stream
  summary<<"# Channel name"<<endl;
  summary<<"# Native sampling [Hz]"<<endl;
  summary<<"# Effective sampling [Hz]"<<endl;
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
    cout<<"\nOmiscan: process channel "<<c+1<<"/"<<channels.size()<<": "<<channels[c]<<"..."<<endl;
    
    


      // print summary
      if(!w){
	summary<<channels[c]<<" "<<sampling<<" "<<sampling_new<<" "<<setprecision(12)<<loudest_t+gps<<" "<<setprecision(5)<<loudest_f<<" "<<loudest_s<<" "<<loudest_q<<" "<<tiles->GetNQPlanes()<<" ";
	for(int q=0; q<tiles->GetNQPlanes(); q++) summary<<setprecision(5)<<tiles->GetQ(q)<<" ";
	summary<<endl;
      }
      
      // save map
      tmpstream<<outdir<<"/plots/"<<channels[c]<<"_map_dt"<<windows[w]<<".gif";
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();

      // Make projections
      GPlot->SetLogx(0);
      GPlot->SetLogy(1);
      GPlot->SetGridx(1);
      GPlot->SetGridy(1);
      projt=map[w]->ProjectionX("_pfx",loudest_bin_f,loudest_bin_f);
      tmpstream<<"SNR vs time at f = "<<setprecision(5)<<loudest_f<<"Hz";
      projt->SetTitle(tmpstream.str().c_str());
      tmpstream.str(""); tmpstream.clear();
      projt->GetXaxis()->SetTitle("Time [s]");
      projt->GetYaxis()->SetTitle("SNR");
      projt->GetXaxis()->SetTitleOffset(1.1);
      projt->GetXaxis()->SetLabelSize(0.045);
      projt->GetYaxis()->SetLabelSize(0.045);
      projt->GetXaxis()->SetTitleSize(0.045);
      projt->GetYaxis()->SetTitleSize(0.045);
      GPlot->Draw(projt);
      tmpstream<<outdir<<"/plots/"<<channels[c]<<"_projt_dt"<<windows[w]<<".gif";
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();
      delete projt;
      GPlot->SetLogx(1);
      projf=map[w]->ProjectionY("_pfy",loudest_bin_t,loudest_bin_t);
      tmpstream<<"SNR vs frequency at GPS = "<<setprecision(12)<<loudest_t+gps;
      projf->SetTitle(tmpstream.str().c_str());
      tmpstream.str(""); tmpstream.clear();
      projf->GetXaxis()->SetTitle("Frequency [Hz]");
      projf->GetYaxis()->SetTitle("SNR");
      projf->GetXaxis()->SetTitleOffset(1.1);
      projf->GetXaxis()->SetLabelSize(0.045);
      projf->GetYaxis()->SetLabelSize(0.045);
      projf->GetXaxis()->SetTitleSize(0.045);
      projf->GetYaxis()->SetTitleSize(0.045);
      GPlot->Draw(projf);
      tmpstream<<outdir<<"/plots/"<<channels[c]<<"_projf_dt"<<windows[w]<<".gif";
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();
      delete projf;
      delete map[w];
    }
    delete tiles;
    delete qmap;
        
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
    if(writeroot){
      rootfile->cd();
      graph->Write();
    }

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
    }
    delete graph;

    //*****  ASD plot  ***********************************************
    if(verbose) cout<<"         Make ASD..."<<endl;
    graph = new TGraph();
    graph->SetName("ASD");
    tmpstream<<channels[c]<<" amplitude spectral density";
    graph->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
    for(int i=0; i<asdsize; i++)
      graph->SetPoint(i,i*(double)sampling_new/2.0/(double)asdsize,sqrt(asd[i]));
    delete data;

    GPlot->SetLogy(1);
    GPlot->SetLogx(1);
    GPlot->SetGridx(1);
    GPlot->SetGridy(1);
    GPlot->Draw(graph,"APL");
    if(writeroot){
      rootfile->cd();
      graph->Write();
    }

    // cosmetics
    graph->GetXaxis()->SetTitle("Frequency [Hz]");
    graph->GetYaxis()->SetTitle("Amp / #sqrt{Hz}");
    graph->GetXaxis()->SetLimits(fmin,fmax);
    graph->GetXaxis()->SetTitleOffset(1.1);
    graph->GetXaxis()->SetLabelSize(0.045);
    graph->GetYaxis()->SetLabelSize(0.045);
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetTitleSize(0.045);
    graph->SetLineWidth(2);
    graph->SetLineColor(7);
    graph->SetMarkerColor(7);
	
    // save ASD plot
    tmpstream<<outdir<<"/plots/"<<channels[c]<<"_asd.gif";
    GPlot->Print(tmpstream.str());
    tmpstream.str(""); tmpstream.clear();
    delete graph;

    if(writeroot) rootfile->Close();
  }

  summary<<"*** omiscan done ***"<<endl;
  summary.close();

  // cleaning
  delete segment;
  delete GPlot;
  delete map;
  */
