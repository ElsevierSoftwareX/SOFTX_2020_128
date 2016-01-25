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
  bandMultiple = new int[0];
  bandCenter = new double[0];
  tilephase = NULL;
  tilecontent = NULL;
  tiletag = NULL;
}

////////////////////////////////////////////////////////////////////////////////////
Omap::~Omap(void){
////////////////////////////////////////////////////////////////////////////////////
  delete bandCenter;
  delete bandMultiple;
  if(tilecontent!=NULL){
    for(int f=0; f<GetNBands(); f++){
      delete tilecontent[f];
      delete tilephase[f];
      delete tiletag[f];
    }
    delete tilecontent;
    delete tilephase;
    delete tiletag;
  }

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

  // delete content / phase / tag
  if(tilecontent!=NULL){
    for(int f=0; f<GetNBands(); f++) delete tilecontent[f];
    delete tilecontent;
  }
  if(tilephase!=NULL){
    for(int f=0; f<GetNBands(); f++) delete tilephase[f];
    delete tilephase;
  }
  if(tiletag!=NULL){
    for(int f=0; f<GetNBands(); f++) delete tiletag[f];
    delete tiletag;
  }

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

  // create content / phase / tag array
  tilecontent = new double* [GetNBands()];
  tilephase   = new double* [GetNBands()];
  tiletag     = new bool*   [GetNBands()];

  // 
  for(int f=0; f<GetNBands(); f++){
    TimeCumulativeMismatch = (double)aTimeRange * 2.0*TMath::Pi() * GetBandFrequency(f) / aQ;
    Nt = NextPowerOfTwo(TimeCumulativeMismatch / aMismatchStep);
    bandMultiple[f] = GetNbinsX() / Nt;
    bandCenter[f] = GetYaxis()->GetBinCenterLog(aBandIndex+1);
  
    tilecontent[f] = new double [Nt];
    for(int t=0; t<Nt; t++) tilecontent[f][t]=0.0;

    tilephase[f] = new double [Nt];
    for(int t=0; t<Nt; t++) tilephase[f][t]=-100.0;

    tiletag[f] = new bool [Nt];
    for(int t=0; t<Nt; t++) tiletag[f][t]=true;

    Ntiles+=(long int)Nt;
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omap::MakeMapContent(void){
////////////////////////////////////////////////////////////////////////////////////
  for(int bf=0; bf<GetNbinsY(); bf++)
    for(int bt=1; bt<=GetNbinsX(); bt++)
      TH2::SetBinContent(bt,bf+1,(double)tiletag[bf][(bt-1)/bandMultiple[bf]]*tilecontent[bf][(bt-1)/bandMultiple[bf]]);
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omap::MakeMapPhase(void){
////////////////////////////////////////////////////////////////////////////////////
  for(int bf=0; bf<GetNbinsY(); bf++)
    for(int bt=1; bt<=GetNbinsX(); bt++)
      TH2::SetBinContent(bt,bf+1,(double)tiletag[bf][(bt-1)/bandMultiple[bf]]*tilephase[bf][(bt-1)/bandMultiple[bf]]);
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omap::MakeMapDisplay(void){
////////////////////////////////////////////////////////////////////////////////////
  for(int f=0; f<GetNbinsY(); f++)
    for(int t=0; t<GetNbinsX(); t++)
      SetBinContent(t+1,f+1,50*((t/bandMultiple[f])%2));
  return;
}

