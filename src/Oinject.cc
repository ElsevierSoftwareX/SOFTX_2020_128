//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Oinject.h"

ClassImp(Oinject)

////////////////////////////////////////////////////////////////////////////////////
Oinject::Oinject(const double aDuration){ 
////////////////////////////////////////////////////////////////////////////////////

  duration = aDuration;
  
  taumin  = 0;
  taumax  = 0;
  phimin  = 32;
  phimax  = 2048;
  Qmin    = 4;
  Qmax    = 100;
  ampmin  = 1e-21;
  ampmax  = 1e-21;

  // random generator
  randgen = new TRandom3();
  randgen->SetSeed(0); // random seed
  MakeWaveform();
  
}

////////////////////////////////////////////////////////////////////////////////////
Oinject::~Oinject(void){
////////////////////////////////////////////////////////////////////////////////////
  delete randgen;
}

////////////////////////////////////////////////////////////////////////////////////
void Oinject::MakeWaveform(void){
////////////////////////////////////////////////////////////////////////////////////

  // waveform parameters
  GenerateParameters();

  // derived parameters
  Wg       = sqrt(sqrt(2.0/TMath::Pi())*Q/phi);
  sigma_t = Q/phi/TMath::Pi()/sqrt(8.0);
  sigma_f = sqrt(2.0)*phi/Q;

  return;
}

////////////////////////////////////////////////////////////////////////////////////
double Oinject::GetTrueSNR(Spectrum *aSpec){
////////////////////////////////////////////////////////////////////////////////////

  // sum = <|X_n|^2>
  // integrating over positive frequencies and over all frequencies gives the same result,
  // because the window is only non zero over positive frequencies (anti-aliasing)
  double freq, win, sum=0;
  double dfreq=aSpec->GetSpectrumResolution();
 
  for(int i=1; i<aSpec->GetSpectrumSize(); i++){
    freq=aSpec->GetSpectrumFrequency(i);
    win = Wg * exp(-(phi-freq)*(phi-freq)/4.0*Q*Q/phi/phi);
    sum += win*win /sqrt(aSpec->GetPower(freq)/2.0) * dfreq;
  }

  return amp*sum;
}

////////////////////////////////////////////////////////////////////////////////////
void Oinject::GenerateParameters(void){
////////////////////////////////////////////////////////////////////////////////////

  tau   = randgen->Uniform(taumin,taumax);
  phase = randgen->Uniform(0,2*TMath::Pi());
  amp   = ampmin*pow(ampmax/ampmin,randgen->Uniform());
  phi   = randgen->Uniform(phimin, phimax);
  Q     = randgen->Uniform(Qmin, Qmax);


}
