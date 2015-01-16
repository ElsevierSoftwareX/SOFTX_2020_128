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
  GetXaxis()->SetLabelSize(0.045);
  GetYaxis()->SetLabelSize(0.045);
  GetZaxis()->SetLabelSize(0.045);
  GetXaxis()->SetTitleSize(0.045);
  GetYaxis()->SetTitleSize(0.045);
  GetZaxis()->SetTitleSize(0.045);
  GetZaxis()->SetRangeUser(1,50);
  Ntiles=0;
  bandMultiple = new int[0];
}

////////////////////////////////////////////////////////////////////////////////////
Omap::~Omap(void){
////////////////////////////////////////////////////////////////////////////////////
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
  
  // frequency bands
  double FrequencyLogStep = log(aFrequencyMax/aFrequencyMin) / (double)Nf;
  double *fbins = new double [Nf+1];
  for(int f=0; f<=Nf; f++) fbins[f] = aFrequencyMin * exp((double)f*FrequencyLogStep);

  // number of time bins
  double TimeCumulativeMismatch = (double)aTimeRange * 2.0*TMath::Pi() * sqrt(fbins[Nf-1]*fbins[Nf]) / aQ;
  int Nt = NextPowerOfTwo(TimeCumulativeMismatch / aMismatchStep);
  double *tbins = new double [Nt+1];
  for(int t=0; t<=Nt; t++) tbins[t] = -(double)aTimeRange/2.0 + (double)t/(double)Nt*(double)aTimeRange;

  // set binning
  TH1::SetBins(Nt,tbins, Nf,fbins);
  delete fbins;
  delete tbins;

  // band variables
  delete bandMultiple;
  bandMultiple   = new int [GetNBands()];
  Ntiles=0;

  for(int f=0; f<GetNBands(); f++){
    TimeCumulativeMismatch = (double)aTimeRange * 2.0*TMath::Pi() * GetBandFrequency(f) / aQ;
    Nt = NextPowerOfTwo(TimeCumulativeMismatch / aMismatchStep);
    bandMultiple[f] = GetNbinsX() / Nt;
    Ntiles+=Nt;
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omap::SetTileContent(const int aTimeTileIndex, const int aBandIndex, const double aContent){
////////////////////////////////////////////////////////////////////////////////////
  int tstart = aTimeTileIndex * bandMultiple[aBandIndex];
  int tend=tstart+bandMultiple[aBandIndex];
  for(int t=tstart; t<tend; t++) SetBinContent(t+1,aBandIndex+1,aContent);
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omap::SetTileDisplay(){
  ////////////////////////////////////////////////////////////////////////////////////
  for(int f=0; f<GetNbinsY(); f++)
    for(int t=0; t<GetNbinsX(); t++)
      SetBinContent(t+1,f+1,50*((t/bandMultiple[f])%2));
  return;
}
