//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Oqplane.h"

ClassImp(Oqplane)

////////////////////////////////////////////////////////////////////////////////////
Oqplane::Oqplane(const double aQ, const int aSampleFrequency, const int aTimeRange, 
		 const int aTimePad, const double aFrequencyMin, const double aFrequencyMax, 
		 const double aMismatchStep, const double aSNRThreshold){ 
////////////////////////////////////////////////////////////////////////////////////
  
  // save parameters
  fQ=aQ;
  fSampleFrequency=aSampleFrequency;
  fTimeRange=aTimeRange;
  fTimePad=aTimePad;
  fFrequencyMin=aFrequencyMin;
  fFrequencyMax=aFrequencyMax;
  fMismatchStep=aMismatchStep;
  fSNRThreshold=aSNRThreshold;

  // init
  status_OK=true;
  fNumberOfTiles=0;

  // derived parameters
  fQPrime = fQ / sqrt(11);
  double NyquistFrequency = fSampleFrequency / 2;
  double MinimumAllowableFrequency = 50 * fQ / (2 * TMath::Pi() * fTimeRange);
  double MaximumAllowableFrequency = NyquistFrequency/(1 + sqrt(11) / fQ);
  double MinimumFrequencyStep = 1 / (double)fTimeRange;

  // use plane specific minimum/maximum allowable frequency if requested
  if(fFrequencyMin<MinimumAllowableFrequency) fFrequencyMin = MinimumAllowableFrequency;
  if(fFrequencyMax>MaximumAllowableFrequency) fFrequencyMax = MaximumAllowableFrequency;
  
  // mismatch between f-rows
  double FrequencyCumulativeMismatch = log(fFrequencyMax/fFrequencyMin) * sqrt(2.0 + fQ*fQ) / 2.0;
  fNumberOfRows = (int)ceil(FrequencyCumulativeMismatch / fMismatchStep);
  if(!fNumberOfRows) fNumberOfRows = 1;
  double fFrequencyMismatchStep = FrequencyCumulativeMismatch / fNumberOfRows;

  // get plan normalization
  GetPlaneNormalization();

  // check parameters
  status_OK*=CheckParameters();

  // nulling
  for(int f=0; f<NFROWMAX; f++) freqrow[f]=NULL;
  hplane=NULL;
  ostringstream titlestream;

  // build q-plane
  if(status_OK){
    fFreq.clear();
    for(int i=0; i<fNumberOfRows; i++)
      fFreq.push_back(floor((fFrequencyMin * exp(2.0/sqrt(2+fQ*fQ) * (0.5+i) * fFrequencyMismatchStep)/MinimumFrequencyStep+0.5))*MinimumFrequencyStep);
    
    // build map
    double *fbin = new double [fNumberOfRows+1];// frequency bins
    double fl=fFreq[0]-sqrt(TMath::Pi())*fFreq[0]/fQ;
    double fh=fFreq[fNumberOfRows-1]+sqrt(TMath::Pi())*fFreq[fNumberOfRows-1]/fQ;
    for(int b=0; b<fNumberOfRows+1; b++) // log
      fbin[b]=fl*pow(10.0,b*log10(fh/fl)/fNumberOfRows);
    titlestream<<"qplane_"<<fQ;
    hplane = new TH2D(titlestream.str().c_str(),titlestream.str().c_str(),(double)(fTimeRange-2*fTimePad)*2.0*sqrt(TMath::Pi())*fFreq[fNumberOfRows-1]/fQ,-(double)fTimeRange/2.0+fTimePad,+(double)fTimeRange/2.0-fTimePad,fNumberOfRows,fbin);
    delete fbin;
    titlestream.str(""); titlestream.clear();
    hplane->GetXaxis()->SetTitle("Time [s]");
    hplane->GetYaxis()->SetTitle("Frequency [Hz]");
    BuildTiles();
  }

  if(!status_OK)
    cerr<<"Oqplane::Oqplane: initialization failed!"<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
Oqplane::~Oqplane(void){
////////////////////////////////////////////////////////////////////////////////////
  fFreq.clear();
  for(int f=0; f<fNumberOfRows; f++) if(freqrow[f]!=NULL) delete freqrow[f];
  if(hplane!=NULL) delete hplane;
}

////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::GetTriggers(Triggers *aTriggers, double *aDataRe, double *aDataIm, const int aTimeStart){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Oqplane::GetTriggers: the Oqplane object is corrupted"<<endl;
    return false;
  }
  for(int f=0; f<fNumberOfRows; f++){
    if(!freqrow[f]->GetTriggers(aTriggers,aDataRe,aDataIm,aTimeStart)){
      cerr<<"Oqplane::GetTriggers: the frequency row "<<f<<" has failed to produce triggers"<<endl;
      return false;
    }
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::SetPowerSpectrum(Spectrum *aSpec){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Oqplane::SetPowerSpectrum: the Oqplane object is corrupted"<<endl;
    return false;
  }

  double sumofweight;
  double power, psdval;

  // set power for each f-row
  for(int f=0; f<fNumberOfRows; f++){
    power=0;
    sumofweight=0;

    // weighted average over the window
    for(int i=0; i<(int)((freqrow[f]->fWindow).size()); i++){
      psdval=aSpec->GetPower(freqrow[f]->fWindowFrequency[i]);
      if(psdval<0){
	cerr<<"Oqplane::SetPowerSpectrum: the power cannot be computed for f="<<freqrow[f]->fWindowFrequency[i]<<endl;
	return false;
      }
      sumofweight+=(freqrow[f]->fWindow[i]*freqrow[f]->fWindow[i]);
      power+=psdval*(freqrow[f]->fWindow[i]*freqrow[f]->fWindow[i]);
    }
    freqrow[f]->SetPower(power/sumofweight);
  }

  return true;
}

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
  for(int f=0; f<fNumberOfRows; f++){
    snrs=freqrow[f]->GetSNRs(aDataRe,aDataIm);
    if(snrs==NULL){
      cerr<<"Oqplane::GetMap: the trigger extraction failed"<<endl;
      return NULL;
    }

    int tfirst, tlast, ffirst, flast;

    // loop over tiles
    for(int t=0; t<freqrow[f]->fNumberOfTiles; t++){
      if(!freqrow[f]->ValidIndices[t]) continue;
      if(printamplitude) snrs[t]*=sqrt(freqrow[f]->fPower);

      // first and last time bin
      tfirst=hplane->GetXaxis()->FindFixBin(freqrow[f]->fTime[t]-freqrow[f]->fDuration/2.0-(double)fTimeRange/2.0+time_offset);
      tlast=hplane->GetXaxis()->FindFixBin(freqrow[f]->fTime[t]+freqrow[f]->fDuration/2.0-(double)fTimeRange/2.0+time_offset);

      // first and last frequency bin
      ffirst=hplane->GetYaxis()->FindFixBin(freqrow[f]->fF-freqrow[f]->fBandWidth/2.0);
      flast=hplane->GetYaxis()->FindFixBin(freqrow[f]->fF+freqrow[f]->fBandWidth/2.0);
      
      // fill map
      for(int bt=tfirst; bt<=tlast; bt++)
	for(int bf=ffirst; bf<=flast; bf++)
	  if(snrs[t]>=hplane->GetBinContent(bt,bf)) hplane->SetBinContent(bt,bf,snrs[t]);
    }

    delete snrs;
  }

  if(printamplitude) hplane->GetZaxis()->SetTitle("Amplitude");
  else hplane->GetZaxis()->SetTitle("SNR");
   
  return hplane;
}

////////////////////////////////////////////////////////////////////////////////////
void Oqplane::GetPlaneNormalization(void){
////////////////////////////////////////////////////////////////////////////////////  
  double coefficients[9] = {0,-2,0,22.0/3.0,0,-146.0/15.0,0,186.0/35.0,0};

  // for large qPrime
  if(fQPrime > 10) fPlaneNormalization = 1.0; // use asymptotic value of planeNormalization
  else{
    //polynomial coefficients for plane normalization factor
    coefficients[0] = log((fQPrime + 1) / (fQPrime - 1));
    coefficients[2] = - 4*log((fQPrime + 1) / (fQPrime - 1));
    coefficients[4] = 6*log((fQPrime + 1) / (fQPrime - 1));
    coefficients[6] = - 4*log((fQPrime + 1) / (fQPrime - 1));
    coefficients[8] = log((fQPrime + 1) / (fQPrime - 1));

    // plane normalization factor
    fPlaneNormalization = sqrt(256 / (315 * fQPrime * (
						      coefficients[0]*pow(fQPrime,8)+
						      coefficients[1]*pow(fQPrime,7)+
						      coefficients[2]*pow(fQPrime,6)+
						      coefficients[3]*pow(fQPrime,5)+
						      coefficients[4]*pow(fQPrime,4)+
						      coefficients[5]*pow(fQPrime,3)+
						      coefficients[6]*pow(fQPrime,2)+
						      coefficients[7]*fQPrime+
						      coefficients[8])));
  }


  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Oqplane::BuildTiles(void){
////////////////////////////////////////////////////////////////////////////////////  

  // for each frequency row
  for(int f=0; f<fNumberOfRows; f++){
    freqrow[f] = new FreqRow(fTimeRange,fTimePad,fSampleFrequency,fFreq[f],fQ,fMismatchStep,fSNRThreshold);
    fNumberOfTiles+=freqrow[f]->fNumberOfTiles;
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////////
bool Oqplane::CheckParameters(void){
////////////////////////////////////////////////////////////////////////////////////  
  if(!status_OK) return false;

  if(fTimeRange<=0){
    cerr<<"Oqplane::CheckParameters: incorrect time range"<<endl;
    return false;
  }

  if(fFrequencyMin >= fFrequencyMax){
    cerr<<"Oqplane::CheckParameters: The time range is too small to cover this frequency range"<<endl;
    return false;
  }

  int x = fSampleFrequency * fTimeRange;
  if(!IsPowerOfTwo(x)){
    cerr<<"Oqplane::CheckParameters: data length is not an integer power of two"<<endl;
    return false;
  }
 
  if(fNumberOfRows > NFROWMAX){
    cerr<<"Oqplane::CheckParameters: the number of frequency rows is too large: "<<fNumberOfRows<<">"<<NFROWMAX<<endl;
    return false;
  }

  return true;
}

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
  double nextpowerof2=ceil(log(TimeCumulativeMismatch / fMismatchStep)/log(2));
  fNumberOfTiles = (int)pow(2.0,nextpowerof2);
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
  double RowNormalization = sqrt(315*QPrime/128/fF);
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
}

////////////////////////////////////////////////////////////////////////////////////
FreqRow::~FreqRow(void){ 
////////////////////////////////////////////////////////////////////////////////////
  fWindow.clear();
  fWindowFrequency.clear();
  fDataIndices.clear();
  fTime.clear();
  delete ValidIndices;
}

////////////////////////////////////////////////////////////////////////////////////
bool FreqRow::GetTriggers(Triggers *aTriggers, double *aDataRe, double *aDataIm, const int aStartTime){
////////////////////////////////////////////////////////////////////////////////////
  
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
  // fill windowed vector with zero padding
  double *working_vector[2];
  working_vector[0] = new double [fNumberOfTiles];
  working_vector[1] = new double [fNumberOfTiles];
  int i=0, index=0, end=0;

  end=fNumberOfTiles/2-rightZeroPadLength;
  for(; i<end; i++){
    index=i+fNumberOfTiles/2-leftZeroPadLength;
    working_vector[0][i]=fWindow[index]*aDataRe[fDataIndices[index]];
    working_vector[1][i]=fWindow[index]*aDataIm[fDataIndices[index]];
  }
  end+=leftZeroPadLength+rightZeroPadLength;
  for(; i<end; i++){
    working_vector[0][i]=0;
    working_vector[1][i]=0;
  }
  end=fNumberOfTiles;
  for(; i<end; i++){
    index=i-(leftZeroPadLength+fNumberOfTiles/2);
    working_vector[0][i]=fWindow[index]*aDataRe[fDataIndices[index]];
    working_vector[1][i]=fWindow[index]*aDataIm[fDataIndices[index]];
  }

  // fft-backward
  fft *offt = new fft(fNumberOfTiles);
  if(!offt->Backward(working_vector[0],working_vector[1])){
    cerr<<"FreqRow::GetSNRs: FFT failed"<<endl;
    delete working_vector[0];
    delete working_vector[1];
    return NULL;
  }
  
  // get energies
  double *energies = offt->GetNorm2();
  delete working_vector[0];
  delete working_vector[1];
  delete offt;
 
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
  
  // save triggers
  for(int t=0; t<fNumberOfTiles; t++){
    if(!ValidIndices[t]) continue;
    energies[t]=sqrt(2.0*energies[t]/MeanEnergy);// energies is now SNRs
  }


  return energies;
}
