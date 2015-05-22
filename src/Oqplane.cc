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
  double MinimumFrequencyStep      = 1.0 / (double)TimeRange;
 
  // adjust frequency range
  if(FrequencyMin<MinimumAllowableFrequency) FrequencyMin = MinimumAllowableFrequency;
  if(FrequencyMax>MaximumAllowableFrequency) FrequencyMax = MaximumAllowableFrequency;
  if(FrequencyMax<=FrequencyMin){
    FrequencyMin=MinimumAllowableFrequency;
    FrequencyMax=MaximumAllowableFrequency;
  }

  // get plane normalization
  GetPlaneNormalization();

  // Q plane definition
  ostringstream titlestream;
  titlestream<<"qplane_"<<setprecision(5)<<fixed<<Q;
  SetName(titlestream.str().c_str());
  SetTitle(titlestream.str().c_str());

  // set binning
  Omap::SetBins(Q,FrequencyMin,FrequencyMax,TimeRange,MismatchStep);

  // band variables
  bandPower      = new double  [GetNBands()];
  bandFFT        = new fft*    [GetNBands()];
  bandWindow     = new double* [GetNBands()];
  bandWindowFreq = new double* [GetNBands()];
  bandWindowSize = new int     [GetNBands()];

  double windowargument;
  double rownormalization;
  double ifftnormalization;

  for(int f=0; f<GetNBands(); f++){
    
    // no power
    bandPower[f]=0.0;
           
    // band fft
    bandFFT[f] = new fft(GetBandNtiles(f),"FFTW_ESTIMATE");

    // Gaussian window
    bandWindowSize[f] = 2 * (int)floor(GetBandFrequency(f)/QPrime*TimeRange) + 1;
    bandWindow[f]     = new double [bandWindowSize[f]];
    bandWindowFreq[f] = new double [bandWindowSize[f]];
    rownormalization  = sqrt(315.0*QPrime/128.0/GetBandFrequency(f));
    ifftnormalization = (double)GetBandNtiles(f) / ((double)SampleFrequency*(double)TimeRange);
    for(int i=0; i<bandWindowSize[f]; i++){
      bandWindowFreq[f][i]=(double)(-(bandWindowSize[f]-1)/2 + i) * MinimumFrequencyStep;
      windowargument=bandWindowFreq[f][i]*QPrime/GetBandFrequency(f);
      bandWindow[f][i] = rownormalization*ifftnormalization*(1-windowargument*windowargument)*(1-windowargument*windowargument);
      bandWindowFreq[f][i]+=GetBandFrequency(f);
    }
  }
  
}

////////////////////////////////////////////////////////////////////////////////////
Oqplane::~Oqplane(void){
////////////////////////////////////////////////////////////////////////////////////
  for(int f=0; f<GetNBands(); f++){
    delete bandFFT[f];
    delete bandWindow[f];
    delete bandWindowFreq[f];
  }
  delete bandWindowSize;
  delete bandFFT;
  delete bandPower;
  delete bandWindow;
  delete bandWindowFreq;
}


////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::SaveTriggers(MakeTriggers *aTriggers, const double aSNRThr, 
			   const int aLeftTimePad, const int aRightTimePad, const int aT0){
////////////////////////////////////////////////////////////////////////////////////
  double tiletime;
  int tstart, tend;
  for(int f=0; f<GetNBands(); f++){
    
    // remove padding
    tstart=GetTimeTileIndex(f,(double)(aLeftTimePad-GetTimeRange()/2))+1;
    tend=GetTimeTileIndex(f,(double)(GetTimeRange()/2.0-aRightTimePad));
    if(GetTileTime(tend,f)<(double)(GetTimeRange()/2.0-aRightTimePad)) tend++;// GWOLLUM convention

    // fill triggers
    for(int t=tstart; t<tend; t++){
      if(!GetTileTag(t,f)) continue; // apply down-tiling
      if(GetTileContent(t,f)<aSNRThr) continue;// apply SNR threshold
      tiletime=GetTileTime(t,f)+(double)aT0;
      if(!aTriggers->AddTrigger(tiletime,
				GetBandFrequency(f),
				GetTileContent(t,f),
				Q,
				tiletime-GetTileDuration(f)/2.0,
				tiletime+GetTileDuration(f)/2.0,
				GetBandFrequency(f)-GetBandWidth(f)/2.0,
				GetBandFrequency(f)+GetBandWidth(f)/2.0,
				GetTileAmplitude(t,f),
				GetTilePhase(t,f))
	 ) return false; // max triggers
    }
  }

  // for trigger segments see Otile::SaveTriggers()
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::ProjectData(double *aDataRe, double *aDataIm){
////////////////////////////////////////////////////////////////////////////////////

  // locals
  double *working_vector[2];    // working vector RE&IM
  int i, index, dataindex, end; // indexes
  int ZeroPadSize, LeftZeroPadSize, RightZeroPadSize;// number of zeros to pad
  double *energies;             // vector of energies
  double *phases;               // vector of phases
  double Thr;                   // Gaussian threshold
  double meanenergy;            // mean Gaussian energy

  // loop over frequency bands
  for(int f=0; f<GetNBands(); f++){
 
    // make working vector
    working_vector[0] = new double [GetBandNtiles(f)];
    working_vector[1] = new double [GetBandNtiles(f)];
    
    // padding sizes
    ZeroPadSize = GetBandNtiles(f)-bandWindowSize[f];
    LeftZeroPadSize = (ZeroPadSize - 1) / 2;
    RightZeroPadSize = (ZeroPadSize + 1) / 2;

    // fill windowed vector with zero padding
    i=0, index=0, end=0;
    end=GetBandNtiles(f)/2-RightZeroPadSize;
    for(; i<end; i++){
      index=i+GetBandNtiles(f)/2-LeftZeroPadSize;
      dataindex=(int)floor((double)(-(bandWindowSize[f]-1)/2 + index)+ 1.5 + GetBandFrequency(f) * GetTimeRange());
      working_vector[0][i]=bandWindow[f][index]*aDataRe[dataindex];// window data
      working_vector[1][i]=bandWindow[f][index]*aDataIm[dataindex];// window data
    }
    end+=LeftZeroPadSize+RightZeroPadSize;
    for(; i<end; i++){
      working_vector[0][i]=0.0;
      working_vector[1][i]=0.0;
    }
    end=GetBandNtiles(f);
    for(; i<end; i++){
      index=i-(LeftZeroPadSize+end/2);
      dataindex=(int)floor((double)(-(bandWindowSize[f]-1)/2 + index)+ 1.5 + GetBandFrequency(f) * GetTimeRange());
      working_vector[0][i]=bandWindow[f][index]*aDataRe[dataindex];// window data
      working_vector[1][i]=bandWindow[f][index]*aDataIm[dataindex];// window data
    }

    // fft-backward
    if(!bandFFT[f]->Backward(working_vector[0],working_vector[1])){
      delete working_vector[0];
      delete working_vector[1];
      return false;
    }
    delete working_vector[0];
    delete working_vector[1];
  
    // get energies/phases
    energies = bandFFT[f]->GetNorm2();
    phases   = bandFFT[f]->GetPhase();

    // make Gaussian threshold
    Thr=1e20;
    UpdateThreshold(f,energies,Thr);
    //UpdateThreshold(f,energies,Thr);
    meanenergy=UpdateThreshold(f,energies,Thr);
    
    // fill tile content
    for(int t=0; t<GetBandNtiles(f); t++){
      SetTileContent(t,f,sqrt(2.0*energies[t]/meanenergy),phases[t]);
      SetTileTag(t,f,1.0);
    }

    delete energies;
    delete phases;
  }
 
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
double Oqplane::UpdateThreshold(const int aBandIndex, double *aEnergies, double &aThreshold){
////////////////////////////////////////////////////////////////////////////////////
  double MeanEnergy=0, RMSEnergy=0;
  int n=0;
  int tstart=GetTimeTileIndex(aBandIndex, -GetTimeRange()/4.0);
  int tend=GetTimeTileIndex(aBandIndex, GetTimeRange()/4.0);
  for(int t=tstart; t<tend; t++){
    if(aEnergies[t]>aThreshold) continue;
    MeanEnergy+=aEnergies[t];    
    RMSEnergy+=aEnergies[t]*aEnergies[t];    
    n++;
  }
  if(n){
    MeanEnergy/=(double)n;
    RMSEnergy/=(double)n;
    RMSEnergy=sqrt(fabs(RMSEnergy-MeanEnergy*MeanEnergy));
    aThreshold = MeanEnergy + 2.0 * RMSEnergy;// 2-sigma threshold
    MeanEnergy/=0.954499736104;//correct gaussian bias 2sigma
  }

  return MeanEnergy;
}

////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::SetPower(Spectrum *aSpec){
////////////////////////////////////////////////////////////////////////////////////

  double sumofweight;
  double power, psdval;

  // set power for each f-row
  for(int f=0; f<GetNBands(); f++){
    power=0;
    sumofweight=0;

    // weighted average over the window
    for(int i=0; i<bandWindowSize[f]; i++){
      psdval=aSpec->GetPower(bandWindowFreq[f][i]);
      if(psdval<0) return false;
      sumofweight+=(bandWindow[f][i]*bandWindow[f][i]);
      power+=psdval*(bandWindow[f][i]*bandWindow[f][i]);
    }
    bandPower[f]=power/sumofweight;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
void Oqplane::GetPlaneNormalization(void){
////////////////////////////////////////////////////////////////////////////////////  
  double coefficients[9] = {0,-2,0,22.0/3.0,0,-146.0/15.0,0,186.0/35.0,0};

  // for large qPrime
  if(QPrime > 10) PlaneNormalization = 1.0; // use asymptotic value of planeNormalization
  else{
    //polynomial coefficients for plane normalization factor
    coefficients[0] = log((QPrime + 1) / (QPrime - 1));
    coefficients[2] = - 4*log((QPrime + 1) / (QPrime - 1));
    coefficients[4] = 6*log((QPrime + 1) / (QPrime - 1));
    coefficients[6] = - 4*log((QPrime + 1) / (QPrime - 1));
    coefficients[8] = log((QPrime + 1) / (QPrime - 1));

    // plane normalization factor
    PlaneNormalization = sqrt(256 / (315 * QPrime * (
						      coefficients[0]*pow(QPrime,8)+
						      coefficients[1]*pow(QPrime,7)+
						      coefficients[2]*pow(QPrime,6)+
						      coefficients[3]*pow(QPrime,5)+
						      coefficients[4]*pow(QPrime,4)+
						      coefficients[5]*pow(QPrime,3)+
						      coefficients[6]*pow(QPrime,2)+
						      coefficients[7]*QPrime+
						      coefficients[8])));
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Oqplane::PrintParameters(void){
////////////////////////////////////////////////////////////////////////////////////
  cout<<"Oqplane::PrintParameters: Q plane properties"<<endl;
  cout<<"\t- Q                         = "<<Q<<endl;
  cout<<"\t- Time range                = "<<GetTimeRange()<<" s"<<endl;
  cout<<"\t- Frequency range           = "<<GetFrequencyMin()<<"-"<<GetFrequencyMax()<<" Hz"<<endl;
  cout<<"\t- Number of frequency rows  = "<<GetNBands()<<endl;
  cout<<"\t- Number of tiles           = "<<Ntiles<<endl;
  cout<<"\t- Number of bins (internal) = "<<GetNbinsX()*GetNBands()<<endl;
  return;
}

