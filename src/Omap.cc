//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omap.h"

ClassImp(Omap)

////////////////////////////////////////////////////////////////////////////////////
Omap::Omap(): TH2D(){ 
////////////////////////////////////////////////////////////////////////////////////
  GetXaxis()->SetTitle("Time [s]");
  GetYaxis()->SetTitle("Frequency [Hz]");
  GetZaxis()->SetTitle("SNR");
  GetXaxis()->SetNoExponent();
  GetXaxis()->SetNdivisions(4,5,0);
  GetXaxis()->SetLabelSize(0.045);
  GetYaxis()->SetLabelSize(0.045);
  GetZaxis()->SetLabelSize(0.045);
  GetXaxis()->SetTitleSize(0.045);
  GetYaxis()->SetTitleSize(0.045);
  GetZaxis()->SetTitleSize(0.045);
  Ntiles=0;
  bandMultiple   = new int[0];
  bandCenter     = new double[0];
}

////////////////////////////////////////////////////////////////////////////////////
Omap::~Omap(void){
////////////////////////////////////////////////////////////////////////////////////
  delete bandCenter;
  delete bandMultiple;
}

////////////////////////////////////////////////////////////////////////////////////
void Omap::SetBins(const double aQ, const double aFrequencyMin, const double aFrequencyMax,
		   const int aTimeRange, const double aMismatchStep){
////////////////////////////////////////////////////////////////////////////////////

  // number of frequency bands
  double FrequencyCumulativeMismatch = log(aFrequencyMax/aFrequencyMin) * sqrt(2.0 + aQ*aQ) / 2.0;
  int Nf = (int)ceil(FrequencyCumulativeMismatch / aMismatchStep);
  if(Nf<=0.0) Nf = 1.0;
  
  // frequency bands (log)
  double FrequencyLogStep = log(aFrequencyMax/aFrequencyMin) / (double)Nf;
  double *fbins = new double [Nf+1];
  for(int f=0; f<=Nf; f++) fbins[f] = aFrequencyMin * exp((double)f*FrequencyLogStep);
  
  // number of time bins
  double TimeCumulativeMismatch = (double)aTimeRange * 2.0*TMath::Pi() * sqrt(fbins[Nf-1]*fbins[Nf]) / aQ;
  int Nt = NextPowerOfTwo(TimeCumulativeMismatch / aMismatchStep);

  // time bins (linear)
  double *tbins = new double [Nt+1];
  for(int t=0; t<=Nt; t++) tbins[t] = -(double)aTimeRange/2.0 + (double)t/(double)Nt*(double)aTimeRange;

  // set binning
  TH1::SetBins(Nt,tbins, Nf,fbins);
  delete fbins;
  delete tbins;
 
  // band parameters
  delete bandCenter;
  bandCenter = new double [GetNBands()];
  delete bandMultiple;
  bandMultiple = new int [GetNBands()];
  Ntiles=0;

  //
  for(int f=0; f<GetNBands(); f++){
    bandCenter[f] = GetYaxis()->GetBinCenterLog(f+1);
    TimeCumulativeMismatch = (double)aTimeRange * 2.0*TMath::Pi() * GetBandFrequency(f) / aQ;
    Nt = NextPowerOfTwo(TimeCumulativeMismatch / aMismatchStep);
    bandMultiple[f] = GetNbinsX() / Nt;
    
    Ntiles+=(long int)Nt;
  }

  return;
}



