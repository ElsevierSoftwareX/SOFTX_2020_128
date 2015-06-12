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
  double df                        = 1.0 / (double)TimeRange;
 
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
    
  // default snr threshold
  SetSNRThr(2.0);
  nTriggers=0;
  
  // band variables
  bandPower      = new double  [GetNBands()];
  bandFFT        = new fft*    [GetNBands()];
  bandWindow     = new double* [GetNBands()];
  bandWindowFreq = new double* [GetNBands()];
  bandWindowSize = new int     [GetNBands()];

  double windowargument;
  double winnormalization;
  double ifftnormalization;
  double delta_f;

  for(int f=0; f<GetNBands(); f++){
    
    // no power
    bandPower[f]=0.0;

            
    // band fft
    ifftnormalization = (double)GetBandNtiles(f) / ((double)SampleFrequency*(double)TimeRange);
    bandFFT[f] = new fft(GetBandNtiles(f),"FFTW_ESTIMATE");

    // Prepare window stuff
    delta_f=GetBandFrequency(f)/QPrime;// from eq. 5.18
    bandWindowSize[f] = 2 * (int)floor(delta_f/df) + 1;
    bandWindow[f]     = new double [bandWindowSize[f]];
    bandWindowFreq[f] = new double [bandWindowSize[f]];
    winnormalization  = sqrt(315.0*QPrime/128.0/GetBandFrequency(f));// eq. 5.26
 
    // Connes window = A * ( 1 - (f/delta_f)^2 )^2 for |f| < delta_f
    for(int i=0; i<bandWindowSize[f]; i++){
      bandWindowFreq[f][i]=(double)(-(bandWindowSize[f]-1)/2 + i) * df;
      windowargument=bandWindowFreq[f][i]/delta_f;// f/delta_f
      bandWindow[f][i] = winnormalization*ifftnormalization*(1-windowargument*windowargument)*(1-windowargument*windowargument);// connes window (1-x^2)^2
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
bool Oqplane::ProjectData(double *aDataRe, double *aDataIm){
////////////////////////////////////////////////////////////////////////////////////

  // locals
  double *working_vector[2];    // working vector RE&IM
  int i, index, dataindex, end; // indexes
  int ZeroPadSize, LeftZeroPadSize, RightZeroPadSize;// number of zeros to pad
  double *energies;             // vector of energies
  double *phases;               // vector of phases
  double meanenergy;            // mean Gaussian energy
  int Nt;                       // number of time tiles
  double snr;                   // tile SNR

  // reset
  nTriggers=0;
  
  // loop over frequency bands
  for(int f=0; f<GetNBands(); f++){

    // number of tiles in this row
    Nt = GetBandNtiles(f);
    
    // make working vector
    working_vector[0] = new double [Nt];
    working_vector[1] = new double [Nt];
    
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
      working_vector[0][i]=bandWindow[f][index]*aDataRe[dataindex];// window data
      working_vector[1][i]=bandWindow[f][index]*aDataIm[dataindex];// window data
    }
    end+=ZeroPadSize;
    for(; i<end; i++){
      working_vector[0][i]=0.0;
      working_vector[1][i]=0.0;
    }
    end=Nt;
    for(; i<end; i++){
      index=i-(LeftZeroPadSize+Nt/2);
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

  //polynomial coefficients for plane normalization factor
  double logfact = log((QPrime + 1) / (QPrime - 1));
  double coefficients[9] = {logfact,
			    -2.0,
			    -4.0*logfact,
			    22.0/3.0,
			    6.0*logfact,
			    -146.0/15.0,
			    -4.0*logfact,
			    186.0/35.0,
			    logfact};// see eq. 5.55

  // for large qPrime
  if(QPrime > 10.0) PlaneNormalization = 1.0; // use asymptotic value of planeNormalization (Fig. 5.3)

  else{
    PlaneNormalization = sqrt(256.0 / (315.0 * QPrime * (
						      coefficients[0]*pow(QPrime,8)+
						      coefficients[1]*pow(QPrime,7)+
						      coefficients[2]*pow(QPrime,6)+
						      coefficients[3]*pow(QPrime,5)+
						      coefficients[4]*pow(QPrime,4)+
						      coefficients[5]*pow(QPrime,3)+
						      coefficients[6]*pow(QPrime,2)+
						      coefficients[7]*QPrime+
						      coefficients[8])));// eq. 5.55
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

