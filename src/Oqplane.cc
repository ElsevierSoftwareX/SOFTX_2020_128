//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Oqplane.h"

ClassImp(Oqplane)

////////////////////////////////////////////////////////////////////////////////////
Oqplane::Oqplane(const double aQ, const int aSampleFrequency, const int aTimeRange, 
		 const double aFrequencyMin, const double aFrequencyMax, 
		 const double aMismatchStep): Omap(){ 
////////////////////////////////////////////////////////////////////////////////////
  
  // save parameters
  Q=aQ;
  int TimeRange          =aTimeRange;
  int SampleFrequency    =aSampleFrequency;
  double FrequencyMin    =aFrequencyMin;
  double FrequencyMax    =aFrequencyMax;
  double MismatchStep    =aMismatchStep;
  
  // adjust time range
  if(TimeRange<4){
    cerr<<"Oqplane::Oqplane: the time range must be at least 4s --> set to 4s"<<endl;
    TimeRange=4;
  }
  if(!IsPowerOfTwo(TimeRange)){
    TimeRange=NextPowerOfTwo(TimeRange);
    cerr<<"Oqplane::Oqplane: the time range must be a power of 2 --> set to "<<TimeRange<<"s"<<endl;
  }

  // derived parameters
  QPrime                           = Q / sqrt(11.0);
  int NyquistFrequency             = SampleFrequency/2;
  double MinimumAllowableFrequency = 4.0 * Q / (2.0 * TMath::Pi() * (double)TimeRange);
  double MaximumAllowableFrequency = (double)NyquistFrequency/(1.0 + sqrt(11.0) / Q);
  
  // adjust frequency range
  if(FrequencyMin<MinimumAllowableFrequency) FrequencyMin = MinimumAllowableFrequency;
  if(FrequencyMax>MaximumAllowableFrequency) FrequencyMax = MaximumAllowableFrequency;
  if(FrequencyMax<=FrequencyMin){
    FrequencyMin=MinimumAllowableFrequency;
    FrequencyMax=MaximumAllowableFrequency;
  }

  // Q plane definition
  ostringstream titlestream;
  titlestream<<"qplane_"<<setprecision(5)<<fixed<<Q;
  SetName(titlestream.str().c_str());
  SetTitle(titlestream.str().c_str());
  titlestream.str(""); titlestream.clear();
  
  // set binning
  Omap::SetBins(Q,FrequencyMin,FrequencyMax,TimeRange,MismatchStep);
    
  // default snr threshold
  SetSNRThr(2.0);

  // band variables
  bandNoiseAmplitude = new double  [GetNBands()];
  bandFFT            = new fft*    [GetNBands()];
  bandWindow         = new double* [GetNBands()];
  bandWindowSize     = new int     [GetNBands()];
  
  double windowargument;
  double winnormalization;
  double ifftnormalization;
  double delta_f;// Connes window 1/2-width
  int k, end;
  //double A1 = GetA1(); // not used
 
  for(int f=0; f<GetNBands(); f++){
    
    // no power
    bandNoiseAmplitude[f]=0.0;
            
    // band fft
    ifftnormalization = 1.0 / (double)TimeRange;
    bandFFT[f] = new fft(GetBandNtiles(f),"FFTW_ESTIMATE", "c2c");

    // Prepare window stuff
    delta_f=GetBandFrequency(f)/QPrime;// from eq. 5.18
    bandWindowSize[f] = 2 * (int)floor(delta_f*(double)TimeRange) + 1;
    bandWindow[f]     = new double [bandWindowSize[f]];
    winnormalization  = sqrt(315.0*QPrime/128.0/GetBandFrequency(f));// eq. 5.26 Localized bursts only!!!

    // bisquare window
    end=(bandWindowSize[f]+1)/2;
    for(k=0; k<end; k++){
      windowargument=2.0*(double)k/(double)(bandWindowSize[f] - 1);
      bandWindow[f][k] = winnormalization*ifftnormalization*(1.0-windowargument*windowargument)*(1.0-windowargument*windowargument);// bisquare window (1-x^2)^2
    }
    // do not save 0s in the center
    end=bandWindowSize[f];
    for(; k<end; k++){
      windowargument=2.0*(double)(k-end)/(double)(bandWindowSize[f] - 1);
      bandWindow[f][k] = winnormalization*ifftnormalization*(1-windowargument*windowargument)*(1-windowargument*windowargument);// bisquare window (1-x^2)^2
    }
  }
  
}

////////////////////////////////////////////////////////////////////////////////////
Oqplane::~Oqplane(void){
////////////////////////////////////////////////////////////////////////////////////
  for(int f=0; f<GetNBands(); f++){
    delete bandFFT[f];
    delete bandWindow[f];
  }
  delete bandWindowSize;
  delete bandFFT;
  delete bandNoiseAmplitude;
  delete bandWindow;
}

////////////////////////////////////////////////////////////////////////////////////
void Oqplane::FillMap(const string aContentType){
////////////////////////////////////////////////////////////////////////////////////

  if(!aContentType.compare("snr")){
    for(int f=0; f<GetNBands(); f++)
      for(int t=0; t<GetBandNtiles(f); t++)
	SetTileContent(t,f,sqrt(GetTileSNR2(t,f)));
  }
  else if(!aContentType.compare("amplitude")){
    for(int f=0; f<GetNBands(); f++)
      for(int t=0; t<GetBandNtiles(f); t++)
	SetTileContent(t,f,GetTileAmplitude(t,f));
  }
  else if(!aContentType.compare("phase")){
    for(int f=0; f<GetNBands(); f++)
      for(int t=0; t<GetBandNtiles(f); t++)
	SetTileContent(t,f,bandFFT[f]->GetPhase_t(t));
  }
  else{
    for(int f=0; f<GetNBands(); f++)
      for(int t=0; t<GetBandNtiles(f); t++)
	SetTileContent(t,f,50.0*(t%2));
  }
  
  return;
}


////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::SaveTriggers(TriggerBuffer *aTriggers, // trigger structure
			   const double aT0,        // time for the center of the map
			   Segments* aSeg           // segments to select the triggers
			   ){
////////////////////////////////////////////////////////////////////////////////////
  int tstart, tend;
  double snr2, snr;
  double SNRThr2 = SNRThr*SNRThr;
  
  for(int f=0; f<GetNBands(); f++){
    
    // remove padding
    tstart=GetTimeTileIndex(f,aSeg->GetStart(0)-aT0);
    tend=GetTimeTileIndex(f,aSeg->GetEnd(aSeg->GetNsegments()-1)-aT0);

    // enforce GWOLLUM convention
    if(GetTileTime(tend,f)<aSeg->GetEnd(aSeg->GetNsegments()-1)-aT0) tend++;

    // fill triggers
    for(int t=tstart; t<tend; t++){

      // get tile snr
      snr2=GetTileSNR2(t,f);

      // apply SNR threshold
      if(snr2<SNRThr2) continue;

      // time select
      if(!aSeg->IsInsideSegment(GetTileTime(t,f)+aT0)) continue;

      // amplitude SNR
      snr=sqrt(snr2);
      
      // save trigger
      if(!aTriggers->AddTrigger(GetTileTime(t,f)+aT0,
				GetBandFrequency(f),
				snr,
				Q,
				GetTileTimeStart(t,f)+aT0,
				GetTileTimeEnd(t,f)+aT0,
				GetBandStart(f),
				GetBandEnd(f),
				snr*bandNoiseAmplitude[f],
				//bandFFT[f]->GetNorm2_t(t),
				bandFFT[f]->GetPhase_t(t)))
	return false;
    }
  }

  // for trigger segments see Otile::SaveTriggers()
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
int Oqplane::ProjectData(fft *aDataFft, const double aPadding){
////////////////////////////////////////////////////////////////////////////////////

  // locals
  int k, Pql, Nt, end, tstart, tend;
  double snrthr2=SNRThr*SNRThr;
  int nt=0; //number of tiles above threshold
   
  // loop over frequency bands
  for(int f=0; f<GetNBands(); f++){

    // number of tiles in this row
    Nt = GetBandNtiles(f);

    // frequency index shift for this row
    Pql=(int)floor(GetBandFrequency(f)*GetTimeRange());

    // populate \tilde{v}
    end=(bandWindowSize[f]+1)/2;
    for(k=0; k<end; k++){
      bandFFT[f]->SetRe_f(k,bandWindow[f][k]*aDataFft->GetRe_f(k+Pql));
      bandFFT[f]->SetIm_f(k,bandWindow[f][k]*aDataFft->GetIm_f(k+Pql));
    }
    end=Nt-(bandWindowSize[f]-1)/2;
    for(; k<end; k++){
      bandFFT[f]->SetRe_f(k,0.0);
      bandFFT[f]->SetIm_f(k,0.0);
    }
    end=Nt;
    for(; k<end; k++){
      bandFFT[f]->SetRe_f(k,bandWindow[f][k-Nt+bandWindowSize[f]]*aDataFft->GetRe_f(abs(Nt-k-Pql)));
      bandFFT[f]->SetIm_f(k,bandWindow[f][k-Nt+bandWindowSize[f]]*aDataFft->GetIm_f(abs(Nt-k-Pql)));
    }

    // fft-backward
    bandFFT[f]->Backward();
    // note the FFT normalization was already included in the window definition

    // from now on, the final Q coefficients are stored in the time vector of bandFFT[f]

    // count tiles above threshold
    tstart=GetTimeTileIndex(f, GetTimeMin()+aPadding);
    tend=GetTimeTileIndex(f, GetTimeMax()-aPadding);
    for(int t=tstart; t<tend; t++){
      if(GetTileSNR2(t,f)>=snrthr2) nt++;
    }
  }
 
  return nt;
}

////////////////////////////////////////////////////////////////////////////////////
double Oqplane::GetMeanEnergy(const int aBandIndex, const double aPadding){
////////////////////////////////////////////////////////////////////////////////////

  // remove segment edges (2 x 25%)
  int tstart=GetTimeTileIndex(aBandIndex, GetTimeMin()+aPadding);
  int tend=GetTimeTileIndex(aBandIndex, GetTimeMax()-aPadding);
  vector <double> v;
  for(int t=tstart; t<tend; t++) v.push_back(bandFFT[aBandIndex]->GetNorm2_t(t));    

  // outlier threshold eq. 5.91 with alpha=2.0
  size_t n25 = 1*v.size() / 4;
  size_t n75 = 3*v.size() / 4;
  nth_element(v.begin(), v.begin()+n25, v.end());
  nth_element(v.begin(), v.begin()+n75, v.end());
  double thr = v[n75]+2.0*(v[n75]-v[n25]);
  v.clear();

  // get mean
  double MeanEnergy=0;
  int n=0;
  for(int t=tstart; t<tend; t++){
    if(bandFFT[aBandIndex]->GetNorm2_t(t)>thr) continue;
    MeanEnergy+=bandFFT[aBandIndex]->GetNorm2_t(t);    
    n++;
  }
  if(n){
    MeanEnergy/=(double)n;
    MeanEnergy*=BIASFACT2;// bias factor for alpha=2.0
    return MeanEnergy;
  }

  return -1.0;
}

////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::SetPower(Spectrum *aSpec1, Spectrum *aSpec2){
////////////////////////////////////////////////////////////////////////////////////

  //double sumofweight;
  double sum, Wb, deltaf, freq, win;
  int i, end;
  
  double dfreq=1.0 / GetTimeRange();

  // set power for each f-row
  for(int f=0; f<GetNBands(); f++){
    sum=0;
    Wb = sqrt(315.0/128.0*QPrime/GetBandFrequency(f));
    deltaf = GetBandFrequency(f) / QPrime;
    i=0;

    // f > phi
    end=(bandWindowSize[f]+1)/2;
    for(; i<end; i++){
      freq=GetBandFrequency(f)+ (double)i / GetTimeRange();
      win = Wb * (1.0-(freq-GetBandFrequency(f))*(freq-GetBandFrequency(f))/deltaf/deltaf) * (1.0-(freq-GetBandFrequency(f))*(freq-GetBandFrequency(f))/deltaf/deltaf);
      sum += win*win /2.0 /sqrt(aSpec1->GetPower(freq)*aSpec2->GetPower(freq)/4.0) * dfreq;
    }

    // f < phi
    end=bandWindowSize[f];
    for(; i<end; i++){
      freq=GetBandFrequency(f)-(double)(end-i) / GetTimeRange();
      win = Wb * (1.0-(freq-GetBandFrequency(f))*(freq-GetBandFrequency(f))/deltaf/deltaf) * (1.0-(freq-GetBandFrequency(f))*(freq-GetBandFrequency(f))/deltaf/deltaf);
      sum += win*win /2.0 /sqrt(aSpec1->GetPower(freq)*aSpec2->GetPower(freq)/4.0) * dfreq;
    }
        
    bandNoiseAmplitude[f]=1.0/sum;
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
double Oqplane::GetA1(void){
////////////////////////////////////////////////////////////////////////////////////  

  if(QPrime > 10.0) return 1.0; // use asymptotic value of planeNormalization (Fig. 5.3)

  double logfact = log((QPrime + 1) / (QPrime - 1));
  double a =
    logfact*pow(QPrime,8)
    -2.0*pow(QPrime,7)
    -4.0*logfact*pow(QPrime,6)
    +22.0/3.0*pow(QPrime,5)
    +6.0*logfact*pow(QPrime,4)
    -146.0/15.0*pow(QPrime,3)
    -4.0*logfact*pow(QPrime,2)
    +186.0/35.0*QPrime+
    +logfact;// eq. 5.55

  return sqrt(256.0 / 315.0 / QPrime / a);
}

////////////////////////////////////////////////////////////////////////////////////
void Oqplane::PrintParameters(void){
////////////////////////////////////////////////////////////////////////////////////
  cout<<"Oqplane::PrintParameters: Q plane properties"<<endl;
  cout<<"\t- Q                         = "<<Q<<endl;
  cout<<"\t- Time range                = "<<GetTimeRange()<<" s"<<endl;
  cout<<"\t- Time resolution           = "<<GetTileDuration(GetNBands()-1)<<" s - "<<GetTileDuration(0)<<" s"<<endl;
  cout<<"\t- Frequency range           = "<<GetFrequencyMin()<<"-"<<GetFrequencyMax()<<" Hz"<<endl;
  cout<<"\t- Number of frequency rows  = "<<GetNBands()<<endl;
  cout<<"\t- Frequency resolution      = "<<GetBandWidth(0)<<" Hz - "<<GetBandWidth(GetNBands()-1)<<" Hz"<<endl;
  cout<<"\t- Number of tiles           = "<<Ntiles<<endl;
  cout<<"\t- Number of bins (internal) = "<<GetNbinsX()*GetNBands()<<endl;
  return;
}

