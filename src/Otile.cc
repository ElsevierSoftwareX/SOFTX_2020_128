//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Otile.h"

ClassImp(Otile)

////////////////////////////////////////////////////////////////////////////////////
Otile::Otile(const int aTimeRange,
	     const double aQMin, const double aQMax, 
	     const double aFrequencyMin, const double aFrequencyMax, 
	     const int aSampleFrequency, const double aMaximumMismatch, 
	     const int aVerbosity): GwollumPlot("otile","GWOLLUM"){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // Plots
  SetLogx(0); SetLogy(1); SetLogz(1);

  // save parameters
  int TimeRange=(int)fabs(aTimeRange);
  double QMin=fabs(aQMin);
  double QMax=fabs(aQMax);
  double FrequencyMin=fabs(aFrequencyMin);
  double FrequencyMax=fabs(aFrequencyMax);
  int SampleFrequency=(int)fabs(aSampleFrequency);
  double MaximumMismatch=TMath::Max(0.01,fabs(aMaximumMismatch));
  fVerbosity=aVerbosity;
  
  ////// Adjust parameters   ////////////////////////////////
  if(QMax<QMin){ QMin=fabs(aQMax); QMax=fabs(aQMin); }
  if(QMin<sqrt(11.0)) QMin=sqrt(11.0);
  if(QMax<sqrt(11.0)) QMax=sqrt(11.0);

  if(MaximumMismatch>0.5){
    MaximumMismatch=0.5;
    cerr<<"Otile::Otile: maximum mismatch is not reasonable --> set to 0.5"<<endl;
  }
  if(!IsPowerOfTwo(SampleFrequency)){
    SampleFrequency=NextPowerOfTwo(SampleFrequency);
    cerr<<"Otile::Otile: the sampling frequency must be a power of two --> set to "<<SampleFrequency<<endl;
  }
  ///////////////////////////////////////////////////////////

  double MismatchStep=2.0*sqrt(MaximumMismatch/3.0);

  // compute Q values
  vector <double> Qs = ComputeQs(QMin,QMax,MaximumMismatch);
  nq = (int)Qs.size();
 
  // create Q planes
  if(aVerbosity) cout<<"Otile::Otile: creating "<<nq<<" Q-planes"<<endl;
  qplanes = new Oqplane* [nq];
  for(int q=0; q<nq; q++){
    qplanes[q]=new Oqplane(Qs[q],SampleFrequency,TimeRange,FrequencyMin,FrequencyMax,MismatchStep);
    if(aVerbosity>1) qplanes[q]->PrintParameters();
  }
  Qs.clear();

  // update parameters  
  TimeRange=qplanes[0]->GetTimeRange();
  FrequencyMin=qplanes[0]->GetFrequencyMin();
  FrequencyMax=qplanes[nq-1]->GetFrequencyMax();

  // frequency binning
  int nfbins = 1000;// FIXME
  double *fbins = new double [nfbins+1];
  for(int b=0; b<=nfbins; b++) fbins[b] = FrequencyMin * exp((double)b/(double)nfbins*log(FrequencyMax/FrequencyMin));

  // create combined tiling
  tilemap = new TH2D("Otiling","Otiling",
		     qplanes[0]->GetBandNtiles(qplanes[0]->GetNBands()-1),-(double)TimeRange/2.0,(double)TimeRange/2.0,
		     nfbins,fbins);
  delete fbins;
  tilemap->GetXaxis()->SetTitle("Time [s]");
  tilemap->GetYaxis()->SetTitle("Frequency [Hz]");
  tilemap->GetZaxis()->SetTitle("SNR");
  tilemap->GetXaxis()->SetLabelSize(0.045);
  tilemap->GetYaxis()->SetLabelSize(0.045);
  tilemap->GetZaxis()->SetLabelSize(0.045);
  tilemap->GetXaxis()->SetTitleSize(0.045);
  tilemap->GetYaxis()->SetTitleSize(0.045);
  tilemap->GetZaxis()->SetTitleSize(0.045);
  tilemap->GetZaxis()->SetRangeUser(1,50);
}

////////////////////////////////////////////////////////////////////////////////////
Otile::~Otile(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"Otile::~Otile"<<endl;
  for(int p=0; p<nq; p++) delete qplanes[p];
  delete qplanes;
  delete tilemap;
}
 
////////////////////////////////////////////////////////////////////////////////////
bool Otile::SetPower(Spectrum *aSpec){
////////////////////////////////////////////////////////////////////////////////////
  if(!aSpec->GetStatus()){
    cerr<<"Otile::SetPower: the Spectrum object is corrupted"<<endl;
    return false;
  }
  for(int p=0; p<nq; p++){
    if(!qplanes[p]->SetPower(aSpec)){
      cerr<<"Otile::SetPower: cannot set power for plane #"<<p<<endl;
      return false;
    }
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::ProjectData(double *aDataRe, double *aDataIm){
////////////////////////////////////////////////////////////////////////////////////
  // project onto q planes
  for(int p=0; p<nq; p++){
    if(!qplanes[p]->ProjectData(aDataRe, aDataIm)){
      cerr<<"Otile::ProjectData: cannot project data onto plane #"<<p<<endl;
      return false;
    }
  }
  //MakeTiling();
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::SaveTriggers(MakeTriggers *aTriggers, const double aSNRThr, const int aLeftTimePad, const int aRightTimePad, const int aT0){
////////////////////////////////////////////////////////////////////////////////////
  if(!aTriggers->Segments::GetStatus()){
    cerr<<"Otile::SaveTriggers: the trigger Segments object is corrupted"<<endl;
    return false;
  }
  if(aLeftTimePad+aRightTimePad>=(double)GetTimeRange()){
    cerr<<"Otile::SaveTriggers: the padding is larger than the time range"<<endl;
    return false;
  }

  // save triggers for each Q plane
  for(int p=0; p<nq; p++){
    if(!qplanes[p]->SaveTriggers(aTriggers, aSNRThr,aLeftTimePad,aRightTimePad,aT0)) return false;// max triggers
  }

  // save segments
  aTriggers->AddSegment((double)aT0-GetTimeRange()/2.0+(double)aLeftTimePad,(double)aT0+GetTimeRange()/2.0-(double)aRightTimePad);
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::SaveMaps(const string aOutdir, const string aName, const int aT0, const string aFormat, vector <int> aWindows){
////////////////////////////////////////////////////////////////////////////////////
  if(!IsDirectory(aOutdir)){
    cerr<<"Otile::SaveMaps: the directory "<<aOutdir<<" is missing"<<endl;
    return false;
  }
  if(!aWindows.size()) aWindows.push_back(GetTimeRange());
   
  if(fVerbosity) cout<<"Otile::SaveMaps: Saving maps for "<<aName<<"..."<<endl;
  ostringstream tmpstream;

  // open root file
  TFile *froot;
  if(aFormat.find("root")!=string::npos){
    tmpstream<<aOutdir<<"/"<<aName<<"_"<<aT0<<"_maps.root";
    froot=TFile::Open(tmpstream.str().c_str(), "recreate");
    tmpstream.clear(); tmpstream.str("");
  }

  // graphix
  vector <string> form;
  if(aFormat.find("gif")!=string::npos) form.push_back("gif");
  if(aFormat.find("png")!=string::npos) form.push_back("png");
  if(aFormat.find("pdf")!=string::npos) form.push_back("pdf");
  if(aFormat.find("ps")!=string::npos)  form.push_back("ps");
  if(aFormat.find("xml")!=string::npos) form.push_back("xml");
  if(aFormat.find("eps")!=string::npos) form.push_back("eps"); 
  if(aFormat.find("jpg")!=string::npos) form.push_back("jpg"); 
  if(aFormat.find("svg")!=string::npos) form.push_back("svg"); 
  
  // save maps for each Q plane
  int xmax, ymax, zmax;
  for(int q=0; q<nq; q++){
    if(fVerbosity>1) cout<<"\t- map Q"<<q<<endl;
    
    // write root
    if(aFormat.find("root")!=string::npos){
      froot->cd();
      qplanes[q]->Write();
    }
      
    if(form.size()){
      // draw map
      if(fVerbosity>1) cout<<"\t\t- Draw map"<<endl;
      qplanes[q]->GetXaxis()->SetRange(-(double)aWindows[(int)aWindows.size()-1]/2.0,(double)aWindows[(int)aWindows.size()-1]/2.0);
      Draw(qplanes[q],"COLZ");
      
      // title
      tmpstream<<aName<<": GPS="<<aT0<<", Q="<<fixed<<setprecision(3)<<qplanes[q]->GetQ();
      qplanes[q]->SetTitle(tmpstream.str().c_str());
      tmpstream.clear(); tmpstream.str("");
      
      // loop over time windows
      if(fVerbosity>1) cout<<"\t\t- Save windowed maps"<<endl;
      for(int w=0; w<(int)aWindows.size(); w++){
	
	// zoom in
	qplanes[q]->GetXaxis()->SetRangeUser(-(double)aWindows[w]/2.0,(double)aWindows[w]/2.0);
	
	// loudest tile
	qplanes[q]->GetMaximumBin(xmax, ymax, zmax);
	tmpstream<<"Loudest: GPS="<<fixed<<setprecision(3)<<(double)aT0+qplanes[q]->GetXaxis()->GetBinCenter(xmax)<<", f="<<qplanes[q]->GetYaxis()->GetBinCenter(ymax)<<" Hz, SNR="<<qplanes[q]->GetBinContent(xmax,ymax);
	AddText(tmpstream.str(), 0.01,0.01,0.03);
	tmpstream.clear(); tmpstream.str("");
	
	// save plot
	for(int f=0; f<(int)form.size(); f++){
	  tmpstream<<aOutdir<<"/"<<aName<<"_"<<aT0<<"_mapQ"<<q<<"dt"<<aWindows[w]<<"."<<form[f];
	  Print(tmpstream.str());
	  tmpstream.clear(); tmpstream.str("");
	}
      }
      
      // unzoom
      qplanes[q]->GetXaxis()->UnZoom();
    }

  }
    
  // close root file
  if(aFormat.find("root")!=string::npos) froot->Close();

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
double Otile::GetQ(const int aQindex){
////////////////////////////////////////////////////////////////////////////////////  
  if(aQindex<0||aQindex>=nq){
    cerr<<"Otile::GetQ: the Q-plane #"<<aQindex<<" does not exist"<<endl;
    return -1.0;
  }
  return qplanes[aQindex]->GetQ();
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::DisplayTiling(const int aQindex){
////////////////////////////////////////////////////////////////////////////////////  
  if(aQindex<0||aQindex>=nq){
    cerr<<"Otile::DisplayTiling: the Q-plane #"<<aQindex<<" does not exist"<<endl;
    return false;
  }
  qplanes[aQindex]->SetTileDisplay();
  Draw(qplanes[aQindex],"COL");
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::DrawPlane(const int aQindex){
////////////////////////////////////////////////////////////////////////////////////  
  if(aQindex<0||aQindex>=nq){
    cerr<<"Otile::DrawPlane: the Q-plane #"<<aQindex<<" does not exist"<<endl;
    return false;
  }
  Draw(qplanes[aQindex],"COLZ");
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
vector <double> Otile::ComputeQs(const double aQMin, const double aQMax, const double aMaximumMismatch){
////////////////////////////////////////////////////////////////////////////////////  
  
  // number of planes
  double QCumulativeMismatch = log(aQMax/aQMin)/sqrt(2.0);// cumulative mismatch across Q range
  double mismatchstep = 2.0*sqrt(aMaximumMismatch/3.0);
  int n = (int)ceil(QCumulativeMismatch/mismatchstep);
  if(n<=0) n=1;
  double Qmismatchstep = QCumulativeMismatch/(double)n;
  
  // compute Q values
  vector <double> qs;
  for(int i=0; i<n; i++) qs.push_back(aQMin * exp(sqrt(2.0) * (0.5+i) * Qmismatchstep));
  return qs;
}

/*
////////////////////////////////////////////////////////////////////////////////////
void Otile::MakeTiling(void){
////////////////////////////////////////////////////////////////////////////////////  
  
  // reset
  maps[0]->Reset();

  // loop over q planes
  for(int q=0; q<(int)Qs.size(); q++){
    for(int f=0; f<qplanes[q]->GetNBands(); f++){
      for(int t=0; t<qplanes[q]->GetBandNtiles(f); t++){
	SetTileContent(qplanes[q]->GetTileTimeStart(t,f),qplanes[q]->GetBandStart(f),qplanes[q]->GetTileTimeEnd(t,f),qplanes[q]->GetBandEnd(f),qplanes[q]->GetTileSNR(t,f));
      }
    }
  }

  
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Otile::SetTileContent(const double t1, const double f1, const double t2, const double f2, const double content){
////////////////////////////////////////////////////////////////////////////////////  
  int tstart=(int)floor((t1-maps[0]->GetXaxis()->GetXmin())/GetTimeRange()*maps[0]->GetNbinsX());
  int tend=(int)ceil((t2-maps[0]->GetXaxis()->GetXmin())/GetTimeRange()*maps[0]->GetNbinsX());
  int fstart=(int)floor(log(f1/maps[0]->GetYaxis()->GetXmin())/log(maps[0]->GetYaxis()->GetXmax()/maps[0]->GetYaxis()->GetXmin())*maps[0]->GetNbinsY());
  int fend=(int)ceil(log(f2/maps[0]->GetYaxis()->GetXmin())/log(maps[0]->GetYaxis()->GetXmax()/maps[0]->GetYaxis()->GetXmin())*maps[0]->GetNbinsY());

  for(int t=tstart; t<tend; t++)
    for(int f=fstart; f<fend; f++)
      if(content>maps[0]->GetBinContent(t+1,f+1)) maps[0]->SetBinContent(t+1,f+1,content);
  return;
}
*/
