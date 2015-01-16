//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Oqplane__
#define __Oqplane__

#include "MakeTriggers.h"
#include "Spectrum.h"
#include "FFT.h"
#include "TH2D.h"

using namespace std;

/**
 * Create a time-frequency Q-plane.
 * This class was designed to create and use a time-frequency Q-plane defined by a Q value. This class is entirely private and can only be used through the Otile class.
 *
 * \author    Florent Robinet
 */
class Oqplane {

 public:
  friend class Otile;  ///< Friendly class

 private:
  
  // constructor - destructor
  Oqplane(const double aQ, 
	  const int aSampleFrequency, 
	  const int aTimeRange, 
	  const double aFrequencyMin, 
	  const double aFrequencyMax, 
	  const double aMismatchStep);
  virtual ~Oqplane(void);

  bool ProjectData(double *aDataRe, double *aDataIm);
  bool SaveTriggers(MakeTriggers *aTriggers, const double aSNRThr, const int aLeftTimePad=0, const int aRightTimePad=0, const int aT0=0);
  inline TH2D* GetMap(void){ return qplane; };

  // GETS
  inline double GetTimeRange(void) { return qplane->GetXaxis()->GetXmax()-qplane->GetXaxis()->GetXmin(); };
  //inline const Double_t* GetBands(void) { return qplane->GetYaxis()->GetXbins()->GetArray(); };
  inline int GetNBands(void) { return qplane->GetNbinsY(); };
  inline int GetBandNtiles(const int aBandIndex) { return qplane->GetNbinsX()/bandMultiple[aBandIndex]; };

  inline double GetTileSNR(const int aTimeTileIndex, const int aFrequencyTileIndex){    
    return qplane->GetBinContent(aTimeTileIndex*bandMultiple[aFrequencyTileIndex]+1,aFrequencyTileIndex+1);
  };
  inline double GetTileAmplitude(const int aTimeTileIndex, const int aFrequencyTileIndex){    
    return GetTileSNR(aTimeTileIndex,aFrequencyTileIndex)*sqrt(bandPower[aFrequencyTileIndex]);
  };
  inline double GetTileTime(const int aTimeTileIndex, const int aFrequencyTileIndex){
    return qplane->GetXaxis()->GetBinCenter(aTimeTileIndex*bandMultiple[aFrequencyTileIndex]+1);
  };
  inline double GetTileTimeStart(const int aTimeTileIndex, const int aFrequencyTileIndex){
    return qplane->GetXaxis()->GetBinLowEdge(aTimeTileIndex*bandMultiple[aFrequencyTileIndex]+1);
  };
  inline double GetTileTimeEnd(const int aTimeTileIndex, const int aFrequencyTileIndex){
    return qplane->GetXaxis()->GetBinUpEdge(aTimeTileIndex*bandMultiple[aFrequencyTileIndex]+1);
  };
  inline double GetTileWidth(const int aTimeTileIndex, const int aFrequencyTileIndex){
    return qplane->GetXaxis()->GetBinWidth(aTimeTileIndex*bandMultiple[aFrequencyTileIndex]+1)*bandMultiple[aFrequencyTileIndex];
  };
  inline double GetBandFrequency(const int aBandIndex){ 
    return qplane->GetYaxis()->GetBinCenterLog(aBandIndex+1);
  };
  inline double GetBandStart(const int aBandIndex){ 
    return qplane->GetYaxis()->GetBinLowEdge(aBandIndex+1);
  };
  inline double GetBandEnd(const int aBandIndex){ 
    return qplane->GetYaxis()->GetBinUpEdge(aBandIndex+1);
  };
  inline double GetBandWidth(const int aBandIndex){ 
    return qplane->GetYaxis()->GetBinWidth(aBandIndex+1);
  };
  
  // SETS
  void SetTileSNR(const int aTimeTileIndex, const int aBandIndex, const double aSNR);
  void SetTileSNR(const double aTime, const double aFrequency, const double aSNR);

  void PrintParameters(void); ///< print plane parameters
  double UpdateThreshold(const int aBandIndex, double *aEnergies, double &aThreshold);

  inline void PresentTile(void){
    for(int f=0; f<qplane->GetNbinsY(); f++)
      for(int t=0; t<qplane->GetNbinsX(); t++)
	qplane->SetBinContent(t+1,f+1,(t/bandMultiple[f])%2);
    return;
  }


    
  //bool GetTriggers(MakeTriggers *aTriggers, double *aDataRe, double *aDataIm, const int aTimeStart, const int aExtraTimePadMin=0);
  bool SetPower(Spectrum *aSpec);
    
  // PARAMETERS
  double Q;                         ///< Q value
  double QPrime;                    ///< Q prime = Q / sqrt(11)
 
  // Q-PLANE
  TH2D *qplane;                     ///< Q plane
  int Ntiles;                       ///< number of tiles in the plane

  void GetPlaneNormalization(void); ///< get plane normalization
  double PlaneNormalization;        ///< plane normalization
  
  // FREQUENCY BANDS
  int *bandMultiple;                ///< band multiple (time resolution)
  int *bandWindowSize;              ///< band Gaussian window size
  double **bandWindow;              ///< band Gaussian window
  double **bandWindowFreq;          ///< band Gaussian window frequency
  double *bandPower;                ///< band power
  fft **bandFFT;                    ///< band fft
  
  ClassDef(Oqplane,0)  
};


#endif


