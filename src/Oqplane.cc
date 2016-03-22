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
  double MinimumAllowableFrequency = 50.0 * Q / (2.0 * TMath::Pi() * (double)TimeRange);
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

  // set binning
  Omap::SetBins(Q,FrequencyMin,FrequencyMax,TimeRange,MismatchStep);
    
  // default snr threshold
  SetSNRThr(2.0);
  nTriggers=0;
  
  // band variables
  bandPower      = new double  [GetNBands()];
  bandFFT        = new fft*    [GetNBands()];
  bandWindow     = new double* [GetNBands()];
  bandWindowSize = new int     [GetNBands()];

  double windowargument;
  double winnormalization;
  double ifftnormalization;
  double delta_f;// Connes window 1/2-width
  int k, end;
  //double A1 = GetA1(); // not used
 
  for(int f=0; f<GetNBands(); f++){
    
    // no power
    bandPower[f]=0.0;
            
    // band fft
    //ifftnormalization = (double)GetBandNtiles(f) / (double)TimeRange;
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
      bandWindow[f][k] = winnormalization*ifftnormalization*(1-windowargument*windowargument)*(1-windowargument*windowargument);// bisquare window (1-x^2)^2
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
  delete bandPower;
  delete bandWindow;
}


////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::SaveTriggers(MakeTriggers *aTriggers,
			   const double aLeftTimePad, const double aRightTimePad,
			   const double aT0){
////////////////////////////////////////////////////////////////////////////////////
  int tstart, tend;
  for(int f=0; f<GetNBands(); f++){
    
    // remove padding
    tstart=GetTimeTileIndex(f,-GetTimeRange()/2.0+aLeftTimePad)+1;
    tend=GetTimeTileIndex(f,GetTimeRange()/2.0-aRightTimePad);
    if(GetTileTime(tend,f)<GetTimeRange()/2.0-aRightTimePad) tend++;// GWOLLUM convention

    // fill triggers
    for(int t=tstart; t<tend; t++){
      if(GetTileContent(t,f)<SNRThr) continue;// apply SNR threshold
      if(!GetTileTag(t,f)) continue; // apply down-tiling
      
      if(!aTriggers->AddTrigger(GetTileTime(t,f)+aT0,
				GetBandFrequency(f),
				GetTileContent(t,f),
				Q,
				GetTileTimeStart(t,f)+aT0,
				GetTileTimeEnd(t,f)+aT0,
				GetBandStart(f),
				GetBandEnd(f),
				GetTileAmplitude(t,f),
				GetTilePhase(t,f)))
	return false;
    }
  }

  // for trigger segments see Otile::SaveTriggers()
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::ProjectData(fft *aFft){
////////////////////////////////////////////////////////////////////////////////////

  // locals
  int k, kshift, index, dataindex, end; // indexes
  int ZeroPadSize, LeftZeroPadSize, RightZeroPadSize;// number of zeros to pad
  double *energies;             // vector of energies
  double *phases;               // vector of phases
  double meanenergy;            // mean Gaussian energy
  int Nt;                       // number of time tiles in a row
  double snr;                   // tile SNR

  // reset
  nTriggers=0;
 
  // loop over frequency bands
  for(int f=0; f<GetNBands(); f++){

    // number of tiles in this row
    Nt = GetBandNtiles(f);

    // frequency shift
    kshift=(int)floor(GetBandFrequency(f)*GetTimeRange());
 
    end=(bandWindowSize[f]+1)/2;
    for(k=0; k<end; k++){
      bandFFT[f]->SetRe_f(k,bandWindow[f][k]*aFft->GetRe_f(k+kshift));
      bandFFT[f]->SetIm_f(k,bandWindow[f][k]*aFft->GetIm_f(k+kshift));
    }
    end=Nt-(bandWindowSize[f]-1)/2;
    for(; k<end; k++){
      bandFFT[f]->SetRe_f(k,0.0);
      bandFFT[f]->SetIm_f(k,0.0);
    }
    end=Nt;
    for(; k<end; k++){
      bandFFT[f]->SetRe_f(k,bandWindow[f][k-end+bandWindowSize[f]]*aFft->GetRe_f(k+kshift));
      bandFFT[f]->SetIm_f(k,bandWindow[f][k-end+bandWindowSize[f]]*aFft->GetIm_f(k+kshift));
    }

    /*
    // padding sizes
    ZeroPadSize = Nt-bandWindowSize[f];
    LeftZeroPadSize = (ZeroPadSize - 1) / 2;
    RightZeroPadSize = (ZeroPadSize + 1) / 2;

    // fill windowed vector with zero padding
    //
    // |----| (bandWindowSize+1)/2 (right side)
    //      |-------------------------------| ZeroPadSize
    //      |---------------| RightZeroPadSize
    //                      |---------------| LeftZeroPadSize = RightZeroPadSize + 1
    //     (bandWindowSize-1)/2 (left side) |----|
    //
    // |----|---------------|---------------|----| Nt
    
    i=0, index=0, end=0;
    end=Nt/2-RightZeroPadSize;
    for(; i<end; i++){
      index=i+Nt/2-LeftZeroPadSize;
      dataindex=(int)floor((double)(-(bandWindowSize[f]-1)/2 + index)+ 1.5 + GetBandFrequency(f) * GetTimeRange());
      bandFFT[f]->SetRe_f(i,bandWindow[f][index]*aFft->GetRe_f(dataindex));
      bandFFT[f]->SetIm_f(i,bandWindow[f][index]*aFft->GetIm_f(dataindex));
    }
    end+=ZeroPadSize;
    for(; i<end; i++){
      bandFFT[f]->SetRe_f(i,0.0);
      bandFFT[f]->SetIm_f(i,0.0);
    }
    end=Nt;
    for(; i<end; i++){
      index=i-(LeftZeroPadSize+Nt/2);
      dataindex=(int)floor((double)(-(bandWindowSize[f]-1)/2 + index)+ 1.5 + GetBandFrequency(f) * GetTimeRange());
      bandFFT[f]->SetRe_f(i,bandWindow[f][index]*aFft->GetRe_f(dataindex));
      bandFFT[f]->SetIm_f(i,bandWindow[f][index]*aFft->GetIm_f(dataindex));
    }
    */
    // fft-backward
    bandFFT[f]->Backward();
 
    // get energies/phases
    energies = bandFFT[f]->GetNorm2_t();
    phases   = bandFFT[f]->GetPhase_t();

    // get mean energy
    meanenergy=GetMeanEnergy(f,energies);
    if(meanenergy<=0.0){
      cerr<<"Oqplane::ProjectData: cannot normalize energies"<<endl;
      delete energies;
      delete phases;
      return false;
    }

    // fill tile content
    for(int t=0; t<GetBandNtiles(f); t++){
      if(2.0*energies[t]>meanenergy){
        snr=sqrt(2.0*energies[t]/meanenergy-1);
        SetTileContent(t,f,snr,phases[t],true);// eq. 5.79
        if(snr>=SNRThr) nTriggers++;
      }
      else SetTileContent(t,f,0.0);
    }

    delete energies;
    delete phases;
  }
 
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
double Oqplane::GetMeanEnergy(const int aBandIndex, double *aEnergies){
////////////////////////////////////////////////////////////////////////////////////

  // remove segment edges (2 x 25%)
  int tstart=GetTimeTileIndex(aBandIndex, -GetTimeRange()/4.0);
  int tend=GetTimeTileIndex(aBandIndex, GetTimeRange()/4.0);
  vector <double> v;
  for(int t=tstart; t<tend; t++) v.push_back(aEnergies[t]);    

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
    if(aEnergies[t]>thr) continue;
    MeanEnergy+=aEnergies[t];    
    n++;
  }
  if(n){
    MeanEnergy/=(double)n;
    MeanEnergy/=BIASFACT2;// bias factor for alpha=2.0
    return MeanEnergy;
  }

  return -1.0;
}
/*
////////////////////////////////////////////////////////////////////////////////////
double Oqplane::GetMeanEnergy(const int aBandIndex, double *aEnergies){
////////////////////////////////////////////////////////////////////////////////////

  // remove segment edges (2 x 25%)
  int tstart=GetTimeTileIndex(aBandIndex, -GetTimeRange()/4.0);
  int tend=GetTimeTileIndex(aBandIndex, GetTimeRange()/4.0);
  
  // get mean
  double MeanEnergy=0, VarEnergy=0;
  int n=0;
  for(int t=tstart; t<tend; t++){
    MeanEnergy+=aEnergies[t];
    VarEnergy+=aEnergies[t]*aEnergies[t];
    n++;
  }

  double thr=MeanEnergy/(double)n + 2*sqrt(VarEnergy/(double)n-MeanEnergy*MeanEnergy/(double)n/(double)n);
  MeanEnergy=0;
  VarEnergy=0;
  n=0;
  for(int t=tstart; t<tend; t++){
    if(aEnergies[t]>thr) continue;
    MeanEnergy+=aEnergies[t];
    n++;
  }

  if(n){
    MeanEnergy/=(double)n;
    return MeanEnergy*1.00270710469359381;
  }

  return -1.0;
}
*/
////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::SetPower(Spectrum *aSpec){
////////////////////////////////////////////////////////////////////////////////////

  double sumofweight;
  double power, psdval;
  int end;

  // set power for each f-row
  for(int f=0; f<GetNBands(); f++){
    power=0;
    sumofweight=0;

    // weighted average over the window
    end=(bandWindowSize[f]+1)/2;
    for(int k=0; k<end; k++){
      psdval=aSpec->GetPower(GetBandFrequency(f)+(double)k / GetTimeRange());
      if(psdval<0) return false;
      sumofweight+=(bandWindow[f][k]*bandWindow[f][k]);
      power+=psdval*(bandWindow[f][k]*bandWindow[f][k]);
    }
    end=GetBandNtiles(f);
    for(int k=end-(bandWindowSize[f]-1)/2; k<end; k++){
      psdval=aSpec->GetPower(GetBandFrequency(f)-(double)(k-end) / GetTimeRange());
      if(psdval<0) return false;
      sumofweight+=(bandWindow[f][k]*bandWindow[f][k]);
      power+=psdval*(bandWindow[f][k]*bandWindow[f][k]);
    }
    bandPower[f]=power/sumofweight;
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

