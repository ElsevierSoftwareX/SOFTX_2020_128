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
	     const string aPlotStyle, const int aVerbosity): GwollumPlot("otile",aPlotStyle){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // Plot default
  SetLogx(0); SetLogy(1); SetLogz(1);
  snrscale=50;

  // save parameters
  TimeRange=(int)fabs(aTimeRange);
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
  TileFracMax=1;// save all triggers

}

////////////////////////////////////////////////////////////////////////////////////
Otile::~Otile(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"Otile::~Otile"<<endl;
  for(int p=0; p<nq; p++) delete qplanes[p];
  delete qplanes;
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
bool Otile::ProjectData(double *aDataRe, double *aDataIm, const bool aTileDown){
////////////////////////////////////////////////////////////////////////////////////
  // project onto q planes
  for(int p=0; p<nq; p++){
    if(!qplanes[p]->ProjectData(aDataRe, aDataIm)){
      cerr<<"Otile::ProjectData: cannot project data onto plane #"<<p<<endl;
      return false;
    }
  }

  // tile-down
  if(aTileDown) TileDown();

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::SaveTriggers(MakeTriggers *aTriggers, const double aLeftTimePad, const double aRightTimePad, const double aT0){
////////////////////////////////////////////////////////////////////////////////////
  if(!aTriggers->Segments::GetStatus()){
    cerr<<"Otile::SaveTriggers: the trigger Segments object is corrupted"<<endl;
    return false;
  }
  if(aLeftTimePad+aRightTimePad>=(double)TimeRange){
    cerr<<"Otile::SaveTriggers: the padding is larger than the time range"<<endl;
    return false;
  }

  // check tile frac selection
  for(int p=0; p<nq; p++)
    if(qplanes[p]->GetTileFrac()>TileFracMax){
      cerr<<"Otile::SaveTriggers: more than "<<TileFracMax*100.0<<"% of the tile are above the SNR threshold - do not save segment "<<aT0-(double)(TimeRange/2)+aLeftTimePad<<"-"<<aT0+(double)(TimeRange/2)-aRightTimePad<<endl;
      return false;
    }

  // save triggers for each Q plane
  for(int p=0; p<nq; p++)
    if(!qplanes[p]->SaveTriggers(aTriggers, aLeftTimePad,aRightTimePad,aT0)) return false;
  
  // save segments
  aTriggers->AddSegment(aT0-(double)(TimeRange/2)+aLeftTimePad,aT0+(double)(TimeRange/2)-aRightTimePad);
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
double Otile::SaveMaps(const string aOutdir, const string aName, const int aT0, const string aFormat, vector <int> aWindows, const double aSNRThr, const bool aThumb){
////////////////////////////////////////////////////////////////////////////////////
  if(!IsDirectory(aOutdir)){
    cerr<<"Otile::SaveMaps: the directory "<<aOutdir<<" is missing"<<endl;
    return -1.0;
  }
  if(!aWindows.size()) aWindows.push_back(TimeRange);
   
  if(fVerbosity) cout<<"Otile::SaveMaps: Saving maps for "<<aName<<" centered on "<<aT0<<"..."<<endl;
  ostringstream tmpstream;

  // make maps and apply SNR threshold
  double snrmax=-1, snr;
  int n=0;
  for(int q=0; q<nq; q++){
    qplanes[q]->MakeMapContent();
    qplanes[q]->GetXaxis()->SetRangeUser(-(double)aWindows[0]/2.0,(double)aWindows[0]/2.0);
    snr=qplanes[q]->GetBinContent(qplanes[q]->GetMaximumBin());
    if(qplanes[q]->GetBinContent(qplanes[q]->GetMaximumBin())<aSNRThr) n++;
    if(snr>snrmax) snrmax=snr;
    qplanes[q]->GetXaxis()->UnZoom();
  }
  if(n==nq){
    if(fVerbosity) cout<<"Otile::SaveMaps: map "<<aName<<" (+-"<<aWindows[0]/2.0<<"s) is below SNR threshold -> do not save"<<endl;
    return snrmax;
  }

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
    if(fVerbosity>1) cout<<"\t- map Q"<<q<<" = "<<fixed<<setprecision(3)<<qplanes[q]->GetQ()<<endl;
    
    // write root
    if(aFormat.find("root")!=string::npos){
      froot->cd();
      ((TH2D*)qplanes[q])->Write();
    }

    if(form.size()){
      // draw map
      if(fVerbosity>2) cout<<"\t\t- Draw map"<<endl;
      ApplyOffset(qplanes[q],(double)aT0);
      qplanes[q]->GetXaxis()->SetRange((double)aT0-(double)aWindows[(int)aWindows.size()-1]/2.0,(double)aT0+(double)aWindows[(int)aWindows.size()-1]/2.0);
      Draw(qplanes[q],"COLZ");

      // title
      tmpstream<<aName<<": Q="<<fixed<<setprecision(3)<<qplanes[q]->GetQ();
      qplanes[q]->SetTitle(tmpstream.str().c_str());
      tmpstream.clear(); tmpstream.str("");

      // loop over time windows
      if(fVerbosity>2) cout<<"\t\t- Save windowed maps"<<endl;
      for(int w=0; w<(int)aWindows.size(); w++){
	
	// zoom in
	qplanes[q]->GetXaxis()->SetRangeUser((double)aT0-(double)aWindows[w]/2.0,(double)aT0+(double)aWindows[w]/2.0);
	
	// loudest tile
	qplanes[q]->GetMaximumBin(xmax, ymax, zmax);
	tmpstream<<"Loudest: GPS="<<fixed<<setprecision(3)<<qplanes[q]->GetXaxis()->GetBinCenter(xmax)<<", f="<<qplanes[q]->GetYaxis()->GetBinCenter(ymax)<<" Hz, SNR="<<qplanes[q]->GetBinContent(xmax,ymax);
	AddText(tmpstream.str(), 0.01,0.01,0.03);
	tmpstream.clear(); tmpstream.str("");
	      
	// set vertical range
	if(snrscale>1) qplanes[q]->GetZaxis()->SetRangeUser(1,snrscale);
	else qplanes[q]->GetZaxis()->SetRangeUser(1,TMath::Max(qplanes[q]->GetBinContent(xmax,ymax),5.0));

	// save plot
	for(int f=0; f<(int)form.size(); f++){
	  tmpstream<<aOutdir<<"/"<<aName<<"_"<<aT0<<"_mapQ"<<q<<"dt"<<aWindows[w]<<"."<<form[f];
	  Print(tmpstream.str());
	  tmpstream.clear(); tmpstream.str("");
	  if(aThumb){ //thumbnail
	    tmpstream<<aOutdir<<"/"<<aName<<"_"<<aT0<<"_mapQ"<<q<<"dt"<<aWindows[w]<<"th."<<form[f];
	    Print(tmpstream.str(),0.5);
	    tmpstream.clear(); tmpstream.str("");
	  }
	}
      }

      // unzoom
      qplanes[q]->GetXaxis()->UnZoom();
      ApplyOffset(qplanes[q],-(double)aT0);
    }

  }

  // close root file
  if(aFormat.find("root")!=string::npos) froot->Close();

  // full map
  if(form.size()){
    if(fVerbosity>1) cout<<"\t- full map"<<endl;
    TH2D* fullmap;
    for(int w=0; w<(int)aWindows.size(); w++){
      fullmap = MakeFullMap(aWindows[w],(double)aT0);
      
      // draw map
      Draw(fullmap,"COLZ");
      
      // title
      fullmap->SetTitle(aName.c_str());
      
      // loudest tile
      fullmap->GetMaximumBin(xmax, ymax, zmax);
      tmpstream<<"Loudest: GPS="<<fixed<<setprecision(3)<<fullmap->GetXaxis()->GetBinCenter(xmax)<<", f="<<fullmap->GetYaxis()->GetBinCenter(ymax)<<" Hz, SNR="<<fullmap->GetBinContent(xmax,ymax);
      AddText(tmpstream.str(), 0.01,0.01,0.03);
      tmpstream.clear(); tmpstream.str("");
      
      // set vertical range
      if(snrscale>1) fullmap->GetZaxis()->SetRangeUser(1,snrscale);
      else fullmap->GetZaxis()->SetRangeUser(1,TMath::Max(fullmap->GetBinContent(xmax,ymax),5.0));

      // save plot
      for(int f=0; f<(int)form.size(); f++){
	tmpstream<<aOutdir<<"/"<<aName<<"_"<<aT0<<"_fullmapdt"<<aWindows[w]<<"."<<form[f];
	Print(tmpstream.str());
	tmpstream.clear(); tmpstream.str("");
	if(aThumb){ //thumbnail
	  tmpstream<<aOutdir<<"/"<<aName<<"_"<<aT0<<"_fullmapdt"<<aWindows[w]<<"th."<<form[f];
	  Print(tmpstream.str(),0.5);
	  tmpstream.clear(); tmpstream.str("");
	}
      }
    
      delete fullmap;
    }
  }

  return snrmax;
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
bool Otile::DrawMapTiling(const int aQindex){
////////////////////////////////////////////////////////////////////////////////////  
  if(aQindex<0||aQindex>=nq){
    cerr<<"Otile::DrawMapTiling: the Q-plane #"<<aQindex<<" does not exist"<<endl;
    return false;
  }
  qplanes[aQindex]->MakeMapDisplay();
  Draw(qplanes[aQindex],"COL");
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::DrawMapContent(const int aQindex){
////////////////////////////////////////////////////////////////////////////////////  
  if(aQindex<0||aQindex>=nq){
    cerr<<"Otile::DrawMapContent: the Q-plane #"<<aQindex<<" does not exist"<<endl;
    return false;
  }
  qplanes[aQindex]->MakeMapContent();
  Draw(qplanes[aQindex],"COLZ");
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::DrawMapPhase(const int aQindex){
////////////////////////////////////////////////////////////////////////////////////  
  if(aQindex<0||aQindex>=nq){
    cerr<<"Otile::DrawMapPhase: the Q-plane #"<<aQindex<<" does not exist"<<endl;
    return false;
  }
  qplanes[aQindex]->MakeMapPhase();
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

////////////////////////////////////////////////////////////////////////////////////
TH2D* Otile::MakeFullMap(const int aTimeRange, const double aT0){
////////////////////////////////////////////////////////////////////////////////////  
  
  // create combined tiling
  int nfbins = 2*qplanes[nq-1]->GetNBands();
  double *fbins = new double [nfbins+1];
  double FrequencyMin=qplanes[0]->GetFrequencyMin();
  double FrequencyMax=qplanes[nq-1]->GetFrequencyMax();
  double FrequencyLogStep = log(FrequencyMax/FrequencyMin) / (double)nfbins;
  for(int f=0; f<=nfbins; f++) fbins[f] = FrequencyMin * exp((double)f*FrequencyLogStep);
  TH2D *fullmap = new TH2D("fullmap","Full map",350,aT0-(double)aTimeRange/2.0,aT0+(double)aTimeRange/2.0,nfbins,fbins);
  delete fbins;
  fullmap->GetXaxis()->SetTitle("Time [s]");
  fullmap->GetYaxis()->SetTitle("Frequency [Hz]");
  fullmap->GetZaxis()->SetTitle("SNR");
  fullmap->GetXaxis()->SetNoExponent();
  fullmap->GetXaxis()->SetNdivisions(4,5,0);
  fullmap->GetXaxis()->SetLabelSize(0.045);
  fullmap->GetYaxis()->SetLabelSize(0.045);
  fullmap->GetZaxis()->SetLabelSize(0.045);
  fullmap->GetXaxis()->SetTitleSize(0.045);
  fullmap->GetYaxis()->SetTitleSize(0.045);
  fullmap->GetZaxis()->SetTitleSize(0.045);
  fullmap->GetZaxis()->SetRangeUser(1,50);

  int ttstart, ttend, ffstart, ffend;
  double content;

  // loop over q planes
  for(int q=0; q<nq; q++){

    for(int f=0; f<qplanes[q]->GetNBands(); f++){
      ffstart=fullmap->GetYaxis()->FindBin(qplanes[q]->GetBandStart(f));
      ffend=fullmap->GetYaxis()->FindBin(qplanes[q]->GetBandEnd(f));

      for(int t=0; t<qplanes[q]->GetBandNtiles(f); t++){
	if(!qplanes[q]->GetTileTag(t,f)) continue;
	ttstart=fullmap->GetXaxis()->FindBin(qplanes[q]->GetTileTimeStart(t,f)+aT0);
	ttend=fullmap->GetXaxis()->FindBin(qplanes[q]->GetTileTimeEnd(t,f)+aT0);
	
	content=qplanes[q]->GetTileContent(t,f);
	
	for(int tt=ttstart; tt<=ttend; tt++)
	  for(int ff=ffstart; ff<=ffend; ff++)
	    if(content>fullmap->GetBinContent(tt,ff)) fullmap->SetBinContent(tt,ff,content);
	
      }
    }
    
  }
  
  return fullmap;
}

////////////////////////////////////////////////////////////////////////////////////
void Otile::TileDown(void){
////////////////////////////////////////////////////////////////////////////////////  
    
  double tstart, tend, fstart, fend; // reference tile
  int ttstart, ttend, ffstart, ffend;// test tile
  double content; // reference tile
  bool winner;

  // loop over reference q planes
  for(int q=0; q<nq; q++){
    
    // loop over reference tiles
    for(int f=0; f<qplanes[q]->GetNBands(); f++){
      fstart=qplanes[q]->GetBandStart(f);
      fend=qplanes[q]->GetBandEnd(f);
      for(int t=0; t<qplanes[q]->GetBandNtiles(f); t++){
	content=qplanes[q]->GetTileContent(t,f);
	tstart=qplanes[q]->GetTileTimeStart(t,f);
	tend=qplanes[q]->GetTileTimeEnd(t,f);
	qplanes[q]->SetTileTag(t,f,false);

	// drill through and see if the reference tile is the winner
	winner=true;
	for(int qq=0; qq<nq&&winner; qq++){
	  if(qq==q) continue;
	  ffstart=TMath::Max(qplanes[qq]->GetBandIndex(fstart),0);
	  ffend=TMath::Min(qplanes[qq]->GetBandIndex(fend),qplanes[qq]->GetNBands()-1);
	  for(int ff=ffstart; ff<=ffend; ff++){
	    if(!winner) break;
	    ttstart=qplanes[qq]->GetTimeTileIndex(ff,tstart);
	    ttend=qplanes[qq]->GetTimeTileIndex(ff,tend);
	    for(int tt=ttstart; tt<=ttend; tt++){
	      if(qplanes[qq]->GetTileContent(tt,ff)>content){// looser!
		winner=false;
		break;
	      }
	    }
	  }
	}

	// the reference tile is a winner
	if(winner) qplanes[q]->SetTileTag(t,f,true);

      }
    }
    
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Otile::ApplyOffset(TH2D *aMap, const double aOffset){
////////////////////////////////////////////////////////////////////////////////////
  TArrayD X(*(aMap->GetXaxis()->GetXbins()));
  for(int i = 0; i < X.GetSize(); i++) X[i] += aOffset;
  aMap->GetXaxis()->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
  return;
}

