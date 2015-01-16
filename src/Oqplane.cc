//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Oqplane.h"

ClassImp(Oqplane)

////////////////////////////////////////////////////////////////////////////////////
Oqplane::Oqplane(const double aQ, const int aSampleFrequency, const int aTimeRange, 
		 const double aFrequencyMin, const double aFrequencyMax, 
		 const double aMismatchStep){ 
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

  // frequency binning calculation
  double FrequencyCumulativeMismatch = log(FrequencyMax/FrequencyMin) * sqrt(2.0 + Q*Q) / 2.0;
  int Nf = (int)ceil(FrequencyCumulativeMismatch / MismatchStep);
  if(Nf<=0.0) Nf = 1.0;
  double FrequencyLogStep = log(FrequencyMax/FrequencyMin) / (double)Nf;
  double *fbins = new double [Nf+1];
  for(int f=0; f<=Nf; f++) fbins[f] = FrequencyMin * exp((double)f*FrequencyLogStep);

  // time binning calculation (based on the highest frequency band)
  double TimeCumulativeMismatch = (double)TimeRange * 2.0*TMath::Pi() * sqrt(fbins[Nf-1]*fbins[Nf]) / Q;
  int Nt = NextPowerOfTwo(TimeCumulativeMismatch / MismatchStep);

  // Q plane definition
  ostringstream titlestream;
  titlestream<<"qplane_"<<setprecision(5)<<fixed<<Q;
  qplane = new TH2D(titlestream.str().c_str(),titlestream.str().c_str(),
		    Nt, -(double)TimeRange/2.0,+(double)TimeRange/2.0,
		    Nf, fbins);
  qplane->GetXaxis()->SetTitle("Time [s]");
  qplane->GetYaxis()->SetTitle("Frequency [Hz]");
  qplane->GetZaxis()->SetTitle("SNR");
  qplane->GetXaxis()->SetLabelSize(0.045);
  qplane->GetYaxis()->SetLabelSize(0.045);
  qplane->GetZaxis()->SetLabelSize(0.045);
  qplane->GetXaxis()->SetTitleSize(0.045);
  qplane->GetYaxis()->SetTitleSize(0.045);
  qplane->GetZaxis()->SetTitleSize(0.045);
  qplane->GetZaxis()->SetRangeUser(1,50);
  delete fbins;

  // band variables
  bandMultiple   = new int     [qplane->GetNbinsY()];
  bandPower      = new double  [qplane->GetNbinsY()];
  bandFFT        = new fft*    [qplane->GetNbinsY()];
  bandWindow     = new double* [qplane->GetNbinsY()];
  bandWindowFreq = new double* [qplane->GetNbinsY()];
  bandWindowSize = new int     [qplane->GetNbinsY()];

  double windowargument;
  double rownormalization;
  double ifftnormalization;
  Ntiles=0;

  for(int f=0; f<GetNBands(); f++){
    
    // no power
    bandPower[f]=0.0;
    
    // find multiple
    TimeCumulativeMismatch = (double)TimeRange * 2.0*TMath::Pi() * GetBandFrequency(f) / Q;
    Nt = NextPowerOfTwo(TimeCumulativeMismatch / MismatchStep);
    bandMultiple[f] = qplane->GetNbinsX() / Nt;
    Ntiles+=Nt;
        
    // band fft
    bandFFT[f] = new fft(Nt,"FFTW_ESTIMATE");

    // Gaussian window
    bandWindowSize[f] = 2 * (int)floor(GetBandFrequency(f)/QPrime*TimeRange) + 1;
    bandWindow[f]     = new double [bandWindowSize[f]];
    bandWindowFreq[f] = new double [bandWindowSize[f]];
    rownormalization  = sqrt(315.0*QPrime/128.0/GetBandFrequency(f));
    ifftnormalization = (double)Nt / ((double)SampleFrequency*(double)TimeRange);
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
  for(int f=0; f<qplane->GetNbinsY(); f++){
    delete bandFFT[f];
    delete bandWindow[f];
    delete bandWindowFreq[f];
  }
  delete bandMultiple;
  delete bandWindowSize;
  delete bandFFT;
  delete bandPower;
  delete bandWindow;
  delete bandWindowFreq;
  delete qplane;
}

////////////////////////////////////////////////////////////////////////////////////
void Oqplane::SetTileSNR(const int aTimeTileIndex, const int aBandIndex, const double aSNR){
////////////////////////////////////////////////////////////////////////////////////
  int tstart = aTimeTileIndex * bandMultiple[aBandIndex];
  int t=tstart;
  while(t/bandMultiple[aBandIndex]==tstart/bandMultiple[aBandIndex]){
    qplane->SetBinContent(t+1,aBandIndex+1,aSNR);
    t++;
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Oqplane::SetTileSNR(const double aTime, const double aFrequency, const double aSNR){
////////////////////////////////////////////////////////////////////////////////////
  int t = (int)(aTime-(qplane->GetXaxis()->GetXmax())/qplane->GetXaxis()->GetBinWidth(1));
  int f = qplane->GetYaxis()->FindFixBin(aFrequency)-1;// FIXME: find a direct computation!
  int tstart = (t/bandMultiple[f]) * bandMultiple[f];
  t=tstart;
  while(t/bandMultiple[f] == tstart/bandMultiple[f]){
    qplane->SetBinContent(t+1,f+1,aSNR);
    t++;
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Oqplane::PrintParameters(void){
////////////////////////////////////////////////////////////////////////////////////
  cout<<"\t- Q                         = "<<Q<<endl;
  cout<<"\t- Time range                = "<<qplane->GetXaxis()->GetXmax()-qplane->GetXaxis()->GetXmin()<<" s"<<endl;
  cout<<"\t- Frequency range           = "<<qplane->GetYaxis()->GetXmin()<<"-"<<qplane->GetYaxis()->GetXmax()<<" Hz"<<endl;
  cout<<"\t- Number of frequency rows  = "<<qplane->GetNbinsY()<<endl;
  cout<<"\t- Number of tiles           = "<<Ntiles<<endl;
  cout<<"\t- Number of bins (internal) = "<<qplane->GetNbinsX()*qplane->GetNbinsY()<<endl;
  return;
}

////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::SaveTriggers(MakeTriggers *aTriggers, const double aSNRThr, const int aLeftTimePad, const int aRightTimePad, const int aT0){
////////////////////////////////////////////////////////////////////////////////////
  double tiletime;
  for(int f=0; f<GetNBands(); f++){
    for(int t=0; t<GetBandNtiles(f); t++){
      tiletime=GetTileTime(t,f);
      if(tiletime-GetTileWidth(t,f)/2.0<(double)(aLeftTimePad-GetTimeRange()/2)) continue;
      if(tiletime>(double)(GetTimeRange()/2.0-aRightTimePad)) break;
      if(GetTileSNR(t,f)<aSNRThr) continue;
      if(!aTriggers->AddTrigger(tiletime+(double)aT0,
				GetBandFrequency(f),
				GetTileSNR(t,f),
				Q,
				tiletime+(double)aT0-GetTileWidth(t,f)/2.0,
				tiletime+(double)aT0+GetTileWidth(t,f)/2.0,
				GetBandFrequency(f)-GetBandWidth(f)/2.0,
				GetBandFrequency(f)+GetBandWidth(f)/2.0,
				GetTileAmplitude(t,f))
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
  double *working_vector[2];// working vector RE&IM
  int i, index, dataindex, end;        // indexes
  int ZeroPadSize, LeftZeroPadSize, RightZeroPadSize;// number of zeros to pad at negative frequencies
  double *energies;
  double Thr;
  double meanenergy;

  // loop over frequency bands
  for(int f=0; f<GetNBands(); f++){

    // make working vector
    working_vector[0] = new double [GetBandNtiles(f)];
    working_vector[1] = new double [GetBandNtiles(f)];

    // sizes
    ZeroPadSize = GetBandNtiles(f)-bandWindowSize[f];
    LeftZeroPadSize = (ZeroPadSize - 1) / 2;
    RightZeroPadSize = (ZeroPadSize + 1) / 2;

    // fill windowed vector with zero padding
    i=0, index=0, end=0;
    end=GetBandNtiles(f)/2-RightZeroPadSize;
    for(; i<end; i++){
      index=i+GetBandNtiles(f)/2-LeftZeroPadSize;
      dataindex=(int)floor((double)(-(bandWindowSize[f]-1)/2 + index)+ 1.5 + GetBandFrequency(f) * GetTimeRange());
      working_vector[0][i]=bandWindow[f][index]*aDataRe[dataindex];
      working_vector[1][i]=bandWindow[f][index]*aDataIm[dataindex];
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
      working_vector[0][i]=bandWindow[f][index]*aDataRe[dataindex];
      working_vector[1][i]=bandWindow[f][index]*aDataIm[dataindex];
    }


    // fft-backward
    if(!bandFFT[f]->Backward(working_vector[0],working_vector[1])){
      delete working_vector[0];
      delete working_vector[1];
      return false;
    }
    delete working_vector[0];
    delete working_vector[1];
  
    // get energies
    energies = bandFFT[f]->GetNorm2();

    // make threshold to exclude outliers
    Thr=1e20;
    UpdateThreshold(f,energies,Thr);
    meanenergy=UpdateThreshold(f,energies,Thr);
    
    // fill tile content
    for(int t=0; t<GetBandNtiles(f); t++){
      SetTileSNR(t,f,sqrt(2.0*energies[t]/meanenergy));// SNR
    }
    delete energies;
    
  }
 
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
double Oqplane::UpdateThreshold(const int aBandIndex, double *aEnergies, double &aThreshold){
////////////////////////////////////////////////////////////////////////////////////
  //cout<<aThreshold<<endl;
  double MeanEnergy=0, RMSEnergy=0;
  int numberofvalidtiles=0;
  for(int t=0; t<GetBandNtiles(aBandIndex); t++){
      
    if(GetTileTime(t,aBandIndex)<-GetTimeRange()/4.0) continue;
    if(GetTileTime(t,aBandIndex)>GetTimeRange()/4.0) break;
    if(aEnergies[t]>aThreshold) continue;
    MeanEnergy+=aEnergies[t];    
    RMSEnergy+=aEnergies[t]*aEnergies[t];    
    numberofvalidtiles++;
  }
  if(numberofvalidtiles){
    MeanEnergy/=(double)numberofvalidtiles;
    RMSEnergy/=(double)numberofvalidtiles;
    RMSEnergy-=(MeanEnergy*MeanEnergy);
    RMSEnergy=sqrt(fabs(RMSEnergy));
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

/*
////////////////////////////////////////////////////////////////////////////////////
TH2D* Oqplane::GetMap(double *aDataRe, double *aDataIm, const double time_offset, const bool printamplitude){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Oqplane::GetMap: the Oqplane object is corrupted"<<endl;
    return NULL;
  }

  hplane->Reset();

  // get frequency row snrs
  double *snrs;
  int tfirst, tlast, ffirst, flast;
  for(int f=0; f<fNumberOfRows; f++){
    snrs=freqrow[f]->GetSNRs(aDataRe,aDataIm);
    if(snrs==NULL){
      cerr<<"Oqplane::GetMap: the trigger extraction failed"<<endl;
      return NULL;
    }

    // first and last frequency bin
    ffirst=hplane->GetYaxis()->FindFixBin(freqrow[f]->fF-freqrow[f]->fBandWidth/2.0);
    flast=hplane->GetYaxis()->FindFixBin(freqrow[f]->fF+freqrow[f]->fBandWidth/2.0);

    // loop over tiles
    for(int t=0; t<freqrow[f]->fNumberOfTiles; t++){
      if(!freqrow[f]->ValidIndices[t]) continue;
      if(printamplitude) snrs[t]*=sqrt(freqrow[f]->fPower);

      // first and last time bin
      tfirst=hplane->GetXaxis()->FindFixBin(freqrow[f]->fTime[t]-freqrow[f]->fDuration/2.0-(double)fTimeRange/2.0+time_offset);
      tlast=hplane->GetXaxis()->FindFixBin(freqrow[f]->fTime[t]+freqrow[f]->fDuration/2.0-(double)fTimeRange/2.0+time_offset);

      // fill map
      for(int bt=tfirst; bt<=tlast; bt++)
	for(int bf=ffirst; bf<=flast; bf++)
	  if(snrs[t]>hplane->GetBinContent(bt,bf)) hplane->SetBinContent(bt,bf,snrs[t]);
    }
    delete snrs;
  }

  if(printamplitude) hplane->GetZaxis()->SetTitle("Amplitude");
  else hplane->GetZaxis()->SetTitle("SNR");
   
  return hplane;
}
*/
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

/*
////////////////////////////////////////////////////////////////////////////////////
FreqRow::FreqRow(const int aTimeRange, const int aTimePad, const int aSampleFrequency, 
		 const double aF, const double aQ, const double aMismatchStep,
		 const double aSNRThreshold){ 
////////////////////////////////////////////////////////////////////////////////////

  // save parameters
  fTimeRange=aTimeRange;
  fTimePad=aTimePad;
  fSampleFrequency=aSampleFrequency;
  fF=aF;
  fQ=aQ;
  fMismatchStep=aMismatchStep;
  fSNRThreshold=aSNRThreshold;
  fPower=1;

  // derived parameters
  fBandWidth=2*sqrt(TMath::Pi())*fF/fQ;
  fDuration = 1/(double)fBandWidth;
  double fMinimumFrequencyStep=1/(double)fTimeRange;

  // mismatch between tiles
  double QPrime = fQ / sqrt(11);
  double TimeCumulativeMismatch = (double)fTimeRange * 2.0 * TMath::Pi() * fF / fQ;
  fNumberOfTiles = NextPowerOfTwo(TimeCumulativeMismatch / fMismatchStep);
  double TimeMismatchStep = TimeCumulativeMismatch / (double)fNumberOfTiles;
  double fTimeStep = fQ*TimeMismatchStep/2.0/TMath::Pi()/fF;

  // build times
  fTime.clear();
  for(int i=0; i<fNumberOfTiles; i++) fTime.push_back(i*fTimeStep);

  // build windows
  int HalfWindowLength= (int)floor(fF/QPrime/fMinimumFrequencyStep);
  int fWindowLength = 2 * HalfWindowLength + 1;
  fWindow.clear(); fDataIndices.clear();
  double windowargument;
  double WindowFrequency;
  double RowNormalization = sqrt(315.0*QPrime/128.0/fF);
  double IfftNormalization = (double)fNumberOfTiles / ((double)fSampleFrequency*(double)fTimeRange);
  for(int i=0; i<fWindowLength; i++){
    WindowFrequency=(double)(-HalfWindowLength + i) * fMinimumFrequencyStep;
    fWindowFrequency.push_back(fF+WindowFrequency);
    windowargument=WindowFrequency*QPrime/fF;
    fWindow.push_back(RowNormalization*IfftNormalization*(1-windowargument*windowargument)*(1-windowargument*windowargument));
    fDataIndices.push_back((int)floor(1.5 + fF / fMinimumFrequencyStep - (double)HalfWindowLength + (double)i));
  }
  fZeroPadLength = fNumberOfTiles - fWindowLength;

  // allocate memory
  ValidIndices=new bool [fNumberOfTiles];
  offt = new fft(fNumberOfTiles,"FFTW_ESTIMATE");
}

////////////////////////////////////////////////////////////////////////////////////
bool FreqRow::GetTriggers(MakeTriggers *aTriggers, double *aDataRe, double *aDataIm, const int aStartTime, const int aExtraTimePadMin){
////////////////////////////////////////////////////////////////////////////////////
  
  // get SNRs
  double *snrs = GetSNRs(aDataRe,aDataIm);
  if(snrs==NULL){
    cerr<<"FreqRow::GetTriggers: the trigger extraction failed for row f="<<fF<<endl;
    return false;
  }
  
  // save triggers above SNR thresholds
  double amplitude;
  double peak_time;
  for(int t=0; t<fNumberOfTiles; t++){
    if(!ValidIndices[t]) continue;
    if(fTime[t]-fDuration/2.0<(double)(fTimePad+aExtraTimePadMin)) continue;
    if(snrs[t]<fSNRThreshold) continue;
    amplitude=snrs[t]*sqrt(fPower);
    peak_time=aStartTime+fTime[t];
    if(!aTriggers->AddTrigger(peak_time,
			      fF,
			      snrs[t],
			      fQ,
			      peak_time-fDuration/2.0,
			      peak_time+fDuration/2.0,
			      fF-fBandWidth/2.0,
			      fF+fBandWidth/2.0,
			      amplitude))// maximum number of trigger limit reached
      break;
  }
 
  delete snrs;
 
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
double* FreqRow::GetSNRs(double *aDataRe, double *aDataIm){
////////////////////////////////////////////////////////////////////////////////////
  
  // number of zeros to pad at negative frequencies
  int leftZeroPadLength = (fZeroPadLength - 1) / 2;
  int rightZeroPadLength = (fZeroPadLength + 1) / 2;
    
  if(leftZeroPadLength>fNumberOfTiles/2){
    cerr<<"FreqRow::GetSNRs: the padding is too large???"<<endl;
    return NULL;
  }
  
  double *working_vector[2];
  working_vector[0] = new double [fNumberOfTiles];
  working_vector[1] = new double [fNumberOfTiles];
  
  // fill windowed vector with zero padding
  int i=0, index=0, end=0;
  end=fNumberOfTiles/2-rightZeroPadLength;
  for(; i<end; i++){
    index=i+fNumberOfTiles/2-leftZeroPadLength;
    working_vector[0][i]=fWindow[index]*aDataRe[fDataIndices[index]];
    working_vector[1][i]=fWindow[index]*aDataIm[fDataIndices[index]];
  }
  end+=leftZeroPadLength+rightZeroPadLength;
  for(; i<end; i++){
    working_vector[0][i]=0.0;
    working_vector[1][i]=0.0;
  }
  end=fNumberOfTiles;
  for(; i<end; i++){
    index=i-(leftZeroPadLength+fNumberOfTiles/2);
    working_vector[0][i]=fWindow[index]*aDataRe[fDataIndices[index]];
    working_vector[1][i]=fWindow[index]*aDataIm[fDataIndices[index]];
  }

  // fft-backward
  if(!offt->Backward(working_vector[0],working_vector[1])){
    cerr<<"FreqRow::GetSNRs: FFT failed"<<endl;
    delete working_vector[0];
    delete working_vector[1];
    return NULL;
  }
  delete working_vector[0];
  delete working_vector[1];
  
  // get energies
  double *energies = offt->GetNorm2();
 
  // get threshold to exclude outliers
  double MeanEnergy=0; double RMSEnergy=0; 
  int numberofvalidtiles=0;
  for(int t=0; t<fNumberOfTiles; t++){
    
    // edge selection (GWOLLUM conventions)
    if(fTime[t]-fDuration/2.0<(double)fTimePad) continue;
    if(fTime[t]>(double)(fTimeRange-fTimePad)) continue;
    
    MeanEnergy+=energies[t];    
    RMSEnergy+=energies[t]*energies[t];    
    numberofvalidtiles++;
  }
  MeanEnergy/=(double)numberofvalidtiles;
  RMSEnergy/=(double)numberofvalidtiles;
  RMSEnergy-=(MeanEnergy*MeanEnergy);
  RMSEnergy=sqrt(fabs(RMSEnergy));
  double Thr = MeanEnergy + 2.0 * RMSEnergy;// 2-sigma
  
  // reset and do it again with threshold!!
  MeanEnergy=0; RMSEnergy=0; 
  numberofvalidtiles=0;
  for(int t=0; t<fNumberOfTiles; t++){
    
    // edge selection (GWOLLUM conventions)
    if(fTime[t]-fDuration/2.0<(double)fTimePad) continue;
    if(fTime[t]>(double)(fTimeRange-fTimePad)) continue;
    if(energies[t]>Thr) continue;
    
    MeanEnergy+=energies[t];    
    RMSEnergy+=energies[t]*energies[t];    
    numberofvalidtiles++;
  }
  MeanEnergy/=(double)numberofvalidtiles;
  RMSEnergy/=(double)numberofvalidtiles;
  RMSEnergy-=(MeanEnergy*MeanEnergy);
  RMSEnergy=sqrt(fabs(RMSEnergy));
  Thr = MeanEnergy + 2.0 * RMSEnergy;// 2-sigma
 
  // mean energy without outliers
  MeanEnergy=0;
  numberofvalidtiles=0;
  for(int t=0; t<fNumberOfTiles; t++){
    ValidIndices[t]=true;
        
    // edge selection (GWOLLUM conventions)
    if(fTime[t]-fDuration/2.0<(double)fTimePad){ ValidIndices[t]=false; continue; }
    if(fTime[t]>(double)(fTimeRange-fTimePad)) { ValidIndices[t]=false; continue; }
    if(energies[t]>Thr){ continue; }

    MeanEnergy+=energies[t];    
    numberofvalidtiles++;
  }
  //MeanEnergy/=(numberofvalidtiles*1.00270710469359381);//correct gaussian bias 3sigma
  MeanEnergy/=((double)numberofvalidtiles/0.954499736104);//correct gaussian bias 2sigma
  
  // get energies
  for(int t=0; t<fNumberOfTiles; t++){
    if(!ValidIndices[t]) continue;
    energies[t]=sqrt(2.0*energies[t]/MeanEnergy);// energies is now SNRs
  }


  return energies;
}
*/
