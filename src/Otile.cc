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
  SetLogx(0); SetLogy(1);

  // save parameters
  TimeRange=(int)fabs(aTimeRange);
  double QMin=fabs(aQMin);
  double QMax=fabs(aQMax);
  double FrequencyMin=fabs(aFrequencyMin);
  double FrequencyMax=fabs(aFrequencyMax);
  int SampleFrequency=(int)fabs(aSampleFrequency);
  MaximumMismatch=TMath::Max(0.01,fabs(aMaximumMismatch));
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
  t_snrmax=new int* [nq];
  f_snrmax=new int* [nq];
  for(int q=0; q<nq; q++){
    qplanes[q]=new Oqplane(Qs[q],SampleFrequency,TimeRange,FrequencyMin,FrequencyMax,MismatchStep);
    if(aVerbosity>1) qplanes[q]->PrintParameters();
  }
  Qs.clear();

  // update parameters  
  TimeRange=qplanes[0]->GetTimeRange();

  // set default SNR threshold
  SetSNRThr();

  // set default fill type
  SetMapFill();
  SetRangez();

  // Sequence
  SeqInSegments  = new Segments();
  SeqOutSegments = new Segments();
  SeqOverlap=0;
  SeqOverlapCurrent=SeqOverlap;
  SeqT0=0;
  SeqSeg=0;
}

////////////////////////////////////////////////////////////////////////////////////
Otile::~Otile(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"Otile::~Otile"<<endl;
  for(int q=0; q<nq; q++){
    delete qplanes[q];
  }
  delete qplanes;
  delete t_snrmax;
  delete f_snrmax;
  delete SeqInSegments;
  delete SeqOutSegments;
}
 

////////////////////////////////////////////////////////////////////////////////////
bool Otile::SetSegments(Segments *aInSeg, Segments *aOutSeg){
////////////////////////////////////////////////////////////////////////////////////
  if(aInSeg==NULL || !aInSeg->GetStatus()){
    cerr<<"Otile::SetSegments: input segments are corrupted"<<endl;
    return false;
  }

  // reset sequence
  SeqInSegments->Reset(); SeqOutSegments->Reset();
  SeqOverlapCurrent=SeqOverlap;
  SeqT0=0;// for initialization in NewChunk()
  SeqSeg=0;

  // update sequence
  SeqInSegments->Append(aInSeg);
  if(aOutSeg==NULL) SeqOutSegments->Append(SeqInSegments);// no selection
  else SeqOutSegments->Append(aOutSeg);

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::NewChunk(bool &aNewSegFlag){
////////////////////////////////////////////////////////////////////////////////////
  if(SeqSeg>=SeqInSegments->GetNsegments()){
    cerr<<"Otile::NewChunk: end of segments"<<endl;
    return false;
  }
  
  // current segment is too short
  if((int)SeqInSegments->GetEnd(SeqSeg)-(int)SeqInSegments->GetStart(SeqSeg)<TimeRange){
    SeqSeg++; //  --> move to next segment
    SeqT0=0;
    return NewChunk(aNewSegFlag);
  }

  // end of current segment
  if(SeqT0+TimeRange/2==(int)SeqInSegments->GetEnd(SeqSeg)){
    SeqSeg++; //  --> move to next segment
    SeqT0=0;
    return NewChunk(aNewSegFlag);
  }

  // initialization = start of current segment
  if(!SeqT0){
    SeqT0=(int)SeqInSegments->GetStart(SeqSeg)-TimeRange/2+SeqOverlap;
    aNewSegFlag=true;
  }
  else aNewSegFlag=false;

  // new test chunk
  SeqOverlapCurrent = SeqOverlap;// reset current overlap
  int start_test    = SeqT0+TimeRange/2-SeqOverlap;
  int stop_test     = start_test+TimeRange;

  // chunk ends after current segment end --> adjust overlap
  if(stop_test>(int)SeqInSegments->GetEnd(SeqSeg)){
    SeqT0=(int)SeqInSegments->GetEnd(SeqSeg)-TimeRange/2;
    SeqOverlapCurrent=start_test+SeqOverlap-SeqT0+TimeRange/2;// --> adjust overlap
    return true;
  }

  // OK  
  SeqT0=start_test+TimeRange/2;
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::SetPower(Spectrum *aSpec1, Spectrum *aSpec2){
////////////////////////////////////////////////////////////////////////////////////
  if(!aSpec1->GetStatus()||!aSpec2->GetStatus()){
    cerr<<"Otile::SetPower: the Spectrum object is corrupted"<<endl;
    return false;
  }
  for(int p=0; p<nq; p++){
    if(!qplanes[p]->SetPower(aSpec1, aSpec2)){
      cerr<<"Otile::SetPower: cannot set power for plane #"<<p<<endl;
      return false;
    }
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
int Otile::ProjectData(fft *aDataFft){
////////////////////////////////////////////////////////////////////////////////////

  int nt=0;// number of tiles above threshold
  
  // project onto q planes
  for(int p=0; p<nq; p++){
    nt+=qplanes[p]->ProjectData(aDataFft,(double)(SeqOverlap/2));
  }

  return nt;
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::SaveTriggers(TriggerBuffer *aTriggers){
////////////////////////////////////////////////////////////////////////////////////
  if(!aTriggers->Segments::GetStatus()){
    cerr<<"Otile::SaveTriggers: the trigger Segments object is corrupted"<<endl;
    return false;
  }
  
  if(fVerbosity) cout<<"Otile::SaveTriggers: Saving triggers for "<<aTriggers->GetName()<<" chunk centered on "<<SeqT0<<"..."<<endl;

  // output segments
  Segments *seg = new Segments((double)(SeqT0-TimeRange/2+SeqOverlapCurrent-SeqOverlap/2),(double)(SeqT0+TimeRange/2-SeqOverlap/2));// remove overlaps
  if(!seg->GetStatus()){
    cerr<<"Otile::SaveTriggers: the output segment list is corrupted"<<endl;
    return false;
  }
  seg->Intersect(SeqOutSegments);// apply user-defined output selection
  if(!seg->GetLiveTime()) {delete seg; return true; } // nothing to do

  // if buffer, segments must be provided first
  if(aTriggers->GetBufferSize()) aTriggers->SetBufferSegments(seg);
  
  // save triggers for each Q plane
  for(int p=0; p<nq; p++){
    if(!qplanes[p]->SaveTriggers(aTriggers,(double)SeqT0, seg)){ delete seg; return false; }
  }
    
  // save trigger segments (if no buffer)
  if(!aTriggers->GetBufferSize()){
    if(aTriggers->GetNsegments()&&seg->GetStart(0)>=aTriggers->GetStart(aTriggers->GetNsegments()-1))// after last segments
      aTriggers->Append(seg);
    else
      for(int s=0; s<seg->GetNsegments(); s++) aTriggers->AddSegment(seg->GetStart(s),seg->GetEnd(s));
  }
  
  delete seg;
  return true;
}


////////////////////////////////////////////////////////////////////////////////////
double Otile::SaveMaps(const string aOutdir, const string aName, const string aFormat, vector <int> aWindows, const double aTimeOffset, const bool aThumb){
////////////////////////////////////////////////////////////////////////////////////
  if(!IsDirectory(aOutdir)){
    cerr<<"Otile::SaveMaps: the directory "<<aOutdir<<" is missing"<<endl;
    return -1.0;
  }
  if(!aWindows.size()) aWindows.push_back(TimeRange);
   
  if(fVerbosity) cout<<"Otile::SaveMaps: Saving maps for "<<aName<<" centered on "<<SeqT0+aTimeOffset<<"..."<<endl;
  ostringstream tmpstream;

  // get SNR-max bins /qplane and /window
  int tstart, tend;
  double snr2max, snr2;
  for(int q=0; q<nq; q++){
    t_snrmax[q] = new int [(int)aWindows.size()];
    f_snrmax[q] = new int [(int)aWindows.size()];

    for(int w=0; w<(int)aWindows.size(); w++){
      snr2max =-1.0;
      t_snrmax[q][w] = 0;
      f_snrmax[q][w] = 0;

      for(int f=0; f<qplanes[q]->GetNBands(); f++){
	tstart=qplanes[q]->GetTimeTileIndex(f,aTimeOffset-(double)aWindows[w]/2.0);
	tend=qplanes[q]->GetTimeTileIndex(f,aTimeOffset+(double)aWindows[w]/2.0);

	for(int t=tstart; t<tend; t++){
	  snr2=qplanes[q]->GetTileSNR2(t,f);
	  if(snr2>snr2max){
	    snr2max=snr2;
	    t_snrmax[q][w] = t;
	    f_snrmax[q][w] = f;
	  }
	}
      }
    }
  }

  // get qplane with SNR-max /window (for full map)
  int *q_snrmax = new int [(int)aWindows.size()];
  for(int w=0; w<(int)aWindows.size(); w++){
    snr2max=-1.0;
    q_snrmax[w]=0;
    for(int q=0; q<nq; q++){
      snr2=qplanes[q]->GetTileSNR2(t_snrmax[q][w],f_snrmax[q][w]);
      if(snr2>snr2max){
	snr2max=snr2;
	q_snrmax[w]=q;
      }
    }
  }


  // apply map SNR threshold (first window only!)
  if(qplanes[q_snrmax[0]]->GetTileSNR2(t_snrmax[q_snrmax[0]][0],f_snrmax[q_snrmax[0]][0])<SNRThr_map*SNRThr_map){
    if(fVerbosity) cout<<"Otile::SaveMaps: map "<<aName<<" (+-"<<aWindows[0]/2.0<<"s) is below SNR threshold -> do not save"<<endl;
    double ss = qplanes[q_snrmax[0]]->GetTileSNR2(t_snrmax[q_snrmax[0]][0],f_snrmax[q_snrmax[0]][0]);
    delete q_snrmax;
    for(int q=0; q<nq; q++){
      delete t_snrmax[q];
      delete f_snrmax[q];
    }
    return sqrt(ss);
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
  if(aFormat.find("root")!=string::npos) form.push_back("root");
  
  // no graphical format --> stop here
  if(!form.size()){
    double ss = qplanes[q_snrmax[0]]->GetTileSNR2(t_snrmax[q_snrmax[0]][0],f_snrmax[q_snrmax[0]][0]);
    delete q_snrmax;
    for(int q=0; q<nq; q++){
      delete t_snrmax[q];
      delete f_snrmax[q];
    }
    return sqrt(ss);
  }
  
  // fill maps
  for(int q=0; q<nq; q++) qplanes[q]->FillMap(mapfill,aTimeOffset-(double)aWindows[(int)aWindows.size()-1]/2.0,aTimeOffset+(double)aWindows[(int)aWindows.size()-1]/2.0);
  
  // save maps for each Q plane
  for(int q=0; q<nq; q++){
    if(fVerbosity>1) cout<<"\t- map Q"<<q<<" = "<<fixed<<setprecision(3)<<qplanes[q]->GetQ()<<endl;
    
    // draw map
    if(fVerbosity>2) cout<<"\t\t- Draw map"<<endl;
    ApplyOffset(qplanes[q],(double)SeqT0);
    qplanes[q]->GetXaxis()->SetRange((double)SeqT0+aTimeOffset-(double)aWindows[(int)aWindows.size()-1]/2.0,(double)SeqT0+aTimeOffset+(double)aWindows[(int)aWindows.size()-1]/2.0);
    Draw(qplanes[q],"COLZ");
    
    // title
    tmpstream<<aName<<": Q="<<fixed<<setprecision(3)<<qplanes[q]->GetQ();
    qplanes[q]->SetTitle(tmpstream.str().c_str());
    tmpstream.clear(); tmpstream.str("");
    
    // loop over time windows
    if(fVerbosity>2) cout<<"\t\t- Save windowed maps"<<endl;
    for(int w=0; w<(int)aWindows.size(); w++){
      
      // zoom in
      qplanes[q]->GetXaxis()->SetRangeUser((double)SeqT0+aTimeOffset-(double)aWindows[w]/2.0,(double)SeqT0+aTimeOffset+(double)aWindows[w]/2.0);
      
      // loudest tile
      if(!mapfill.compare("amplitude")) tmpstream<<"Loudest: GPS="<<fixed<<setprecision(3)<<qplanes[q]->GetTileTime(t_snrmax[q][w],f_snrmax[q][w])<<", f="<<qplanes[q]->GetBandFrequency(f_snrmax[q][w])<<" Hz, "<<mapfill<<"="<<scientific<<qplanes[q]->GetTileContent(t_snrmax[q][w],f_snrmax[q][w]);
      else tmpstream<<"Loudest: GPS="<<fixed<<setprecision(3)<<qplanes[q]->GetTileTime(t_snrmax[q][w],f_snrmax[q][w])<<", f="<<qplanes[q]->GetBandFrequency(f_snrmax[q][w])<<" Hz, "<<mapfill<<"="<<qplanes[q]->GetTileContent(t_snrmax[q][w],f_snrmax[q][w]);
      AddText(tmpstream.str(), 0.01,0.01,0.03);
      tmpstream.clear(); tmpstream.str("");
      
      // set vertical range
      if(vrange[0]<vrange[1]) qplanes[q]->GetZaxis()->SetRangeUser(vrange[0],vrange[1]);
      else qplanes[q]->GetZaxis()->UnZoom();
      
      // save plot
      for(int f=0; f<(int)form.size(); f++){
	tmpstream<<aOutdir<<"/"<<aName<<"MAPQ"<<q<<"-"<<SeqT0<<"-"<<aWindows[w]<<"."<<form[f];
	Print(tmpstream.str());
	tmpstream.clear(); tmpstream.str("");
	if(aThumb){ //thumbnail
	  tmpstream<<aOutdir<<"/th"<<aName<<"MAPQ"<<q<<"-"<<SeqT0<<"-"<<aWindows[w]<<"."<<form[f];
	  Print(tmpstream.str(),0.5);
	  tmpstream.clear(); tmpstream.str("");
	}
      }
    }
    
    // unzoom
    qplanes[q]->GetXaxis()->UnZoom();
    ApplyOffset(qplanes[q],-(double)SeqT0);
  }
  
  // full map
  if(fVerbosity>1) cout<<"\t- full map"<<endl;
  TH2D* fullmap;
  for(int w=0; w<(int)aWindows.size(); w++){
    fullmap = MakeFullMap(aWindows[w],aTimeOffset);
    
    // draw map
    Draw(fullmap,"COLZ");
    
    // title
    fullmap->SetTitle(aName.c_str());
    
    // loudest tile
    if(!mapfill.compare("amplitude")) tmpstream<<"Loudest: GPS="<<fixed<<setprecision(3)<<(double)SeqT0+qplanes[q_snrmax[w]]->GetTileTime(t_snrmax[q_snrmax[w]][w],f_snrmax[q_snrmax[w]][w])<<", f="<<qplanes[q_snrmax[w]]->GetBandFrequency(f_snrmax[q_snrmax[w]][w])<<" Hz, "<<mapfill<<"="<<scientific<<qplanes[q_snrmax[w]]->GetTileContent(t_snrmax[q_snrmax[w]][w],f_snrmax[q_snrmax[w]][w]);
    else tmpstream<<"Loudest: GPS="<<fixed<<setprecision(3)<<(double)SeqT0+qplanes[q_snrmax[w]]->GetTileTime(t_snrmax[q_snrmax[w]][w],f_snrmax[q_snrmax[w]][w])<<", f="<<qplanes[q_snrmax[w]]->GetBandFrequency(f_snrmax[q_snrmax[w]][w])<<" Hz, "<<mapfill<<"="<<qplanes[q_snrmax[w]]->GetTileContent(t_snrmax[q_snrmax[w]][w],f_snrmax[q_snrmax[w]][w]);
    AddText(tmpstream.str(), 0.01,0.01,0.03);
    tmpstream.clear(); tmpstream.str("");
    
    // set vertical range
    if(vrange[0]<vrange[1]) fullmap->GetZaxis()->SetRangeUser(vrange[0],vrange[1]);
    else fullmap->GetZaxis()->UnZoom();
        
    // save plot
    for(int f=0; f<(int)form.size(); f++){
      tmpstream<<aOutdir<<"/"<<aName<<"MAP"<<"-"<<SeqT0<<"-"<<aWindows[w]<<"."<<form[f];
      Print(tmpstream.str());
      tmpstream.clear(); tmpstream.str("");
      if(aThumb){ //thumbnail
	tmpstream<<aOutdir<<"/th"<<aName<<"MAP"<<"-"<<SeqT0<<"-"<<aWindows[w]<<"."<<form[f];
	Print(tmpstream.str(),0.5);
	tmpstream.clear(); tmpstream.str("");
      }
    }
    
    delete fullmap;
  }

  double ss = qplanes[q_snrmax[0]]->GetTileSNR2(t_snrmax[q_snrmax[0]][0],f_snrmax[q_snrmax[0]][0]);
  for(int q=0; q<nq; q++){
    delete t_snrmax[q];
    delete f_snrmax[q];
  }
  delete q_snrmax;
  return sqrt(ss);
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::DrawMapTiling(const int aQindex){
////////////////////////////////////////////////////////////////////////////////////  
  if(aQindex<0||aQindex>=nq){
    cerr<<"Otile::DrawMapTiling: the Q-plane #"<<aQindex<<" does not exist"<<endl;
    return false;
  }
  qplanes[aQindex]->FillMap("display",-TimeRange/2,TimeRange/2);
  qplanes[aQindex]->GetZaxis()->SetRangeUser(0,100);
  Draw(qplanes[aQindex],"COL");
  /*
  qplanes[aQindex]->GetXaxis()->SetRangeUser(-0.1,0.1);
  TLine *lf;
  for(int f=0; f<qplanes[aQindex]->GetNBands(); f++){
    lf = new TLine(-0.1,qplanes[aQindex]->GetBandStart(f),0.1,qplanes[aQindex]->GetBandStart(f));
    lf->SetLineStyle(3);
    lf->SetLineColor(2);
    Draw(lf,"same");
  }
  */
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
TH2D* Otile::MakeFullMap(const int aTimeRange, const double aTimeOffset){
////////////////////////////////////////////////////////////////////////////////////  
  
  // create combined tiling
  int nfbins = qplanes[nq-1]->GetNBands();
  double *fbins = new double [nfbins+1];
  double FrequencyMin=qplanes[0]->GetFrequencyMin();
  double FrequencyMax=qplanes[nq-1]->GetFrequencyMax();
  double FrequencyLogStep = log(FrequencyMax/FrequencyMin) / (double)nfbins;
  for(int f=0; f<=nfbins; f++) fbins[f] = FrequencyMin * exp((double)f*FrequencyLogStep);
  //int ntbins = 1000;
  int ntbins = 301;

  TH2D *fullmap = new TH2D("fullmap","Full map",ntbins,SeqT0+aTimeOffset-(double)aTimeRange/2.0,SeqT0+aTimeOffset+(double)aTimeRange/2.0,nfbins,fbins);
  delete fbins;
  fullmap->GetXaxis()->SetTitle("Time [s]");
  fullmap->GetYaxis()->SetTitle("Frequency [Hz]");
  fullmap->GetZaxis()->SetTitle(StringToUpper(mapfill).c_str());
  fullmap->GetXaxis()->SetNoExponent();
  fullmap->GetXaxis()->SetNdivisions(4,5,0);
  fullmap->GetXaxis()->SetLabelSize(0.045);
  fullmap->GetYaxis()->SetLabelSize(0.045);
  fullmap->GetZaxis()->SetLabelSize(0.045);
  fullmap->GetXaxis()->SetTitleSize(0.045);
  fullmap->GetYaxis()->SetTitleSize(0.045);
  fullmap->GetZaxis()->SetTitleSize(0.045);
  
  int tstart, tend, ttstart, ttend, ffstart, ffend;
  double content, snr2;

  // loop over q planes
  for(int q=0; q<nq; q++){

    // loop over frequency band of that q plane
    for(int f=0; f<qplanes[q]->GetNBands(); f++){

      // first and last frequency bin of full map to consider
      ffstart=fullmap->GetYaxis()->FindBin(qplanes[q]->GetBandStart(f));
      ffend=fullmap->GetYaxis()->FindBin(qplanes[q]->GetBandEnd(f));

      // qplane time bin indexes to sweep
      tstart=qplanes[q]->GetTimeTileIndex(f,fullmap->GetXaxis()->GetBinLowEdge(1)-SeqT0);
      tend=qplanes[q]->GetTimeTileIndex(f,fullmap->GetXaxis()->GetBinUpEdge(fullmap->GetNbinsX()-1)-SeqT0);

      // loop over relevant qplane time bins
      for(int t=tstart; t<=tend; t++){

	// first and last time bin of full map to consider
	ttstart=fullmap->GetXaxis()->FindBin(qplanes[q]->GetTileTimeStart(t,f)+SeqT0);
	ttend=fullmap->GetXaxis()->FindBin(qplanes[q]->GetTileTimeEnd(t,f)+SeqT0);
	
	snr2=qplanes[q]->GetTileSNR2(t,f);
	content=qplanes[q]->GetTileContent(t,f);
	
	// loop over full map bins over lapping with qplane bins
	for(int tt=ttstart; tt<=ttend; tt++)
	  for(int ff=ffstart; ff<=ffend; ff++)
	    if(snr2>fullmap->GetBinError(tt,ff)){// update content if higher SNR
	      fullmap->SetBinContent(tt,ff,content);
	      fullmap->SetBinError(tt,ff,snr2);
	    }
      }
    }
    
  }
  
  return fullmap;
}

////////////////////////////////////////////////////////////////////////////////////
void Otile::ApplyOffset(TH2D *aMap, const double aOffset){
////////////////////////////////////////////////////////////////////////////////////
  TArrayD X(*(aMap->GetXaxis()->GetXbins()));
  for(int i = 0; i < X.GetSize(); i++) X[i] += aOffset;
  aMap->GetXaxis()->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
  return;
}


