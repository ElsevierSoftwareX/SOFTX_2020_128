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

  void PrintParameters(void);
  bool ProjectData(double *aDataRe, double *aDataIm);
  bool SaveTriggers(MakeTriggers *aTriggers, 
		    const double aSNRThr, 
		    const int aLeftTimePad=0, 
		    const int aRightTimePad=0, 
		    const int aT0=0);

  // GETS
  inline TH2D* GetMap(void){ 
    return qplane;
  };
  inline double GetTimeRange(void){ 
    return qplane->GetXaxis()->GetXmax()-qplane->GetXaxis()->GetXmin(); 
  };
  inline double GetTimeMin(void){ 
    return qplane->GetXaxis()->GetXmin(); 
  };
  inline double GetTimeMax(void){ 
    return qplane->GetXaxis()->GetXmax(); 
  };
  inline double GetFrequencyMin(void){ 
    return qplane->GetYaxis()->GetXmin(); 
  };
  inline double GetFrequencyMax(void){ 
    return qplane->GetYaxis()->GetXmax(); 
  };
  inline int GetNBands(void){ 
    return qplane->GetNbinsY();
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
  inline double GetTileDuration(const int aBandIndex){
    return qplane->GetXaxis()->GetBinWidth(1)*bandMultiple[aBandIndex];
  };
  inline int GetBandNtiles(const int aBandIndex){ 
    return qplane->GetNbinsX()/bandMultiple[aBandIndex];
  };
  inline double GetTileSNR(const int aTimeTileIndex, const int aBandIndex){    
    return qplane->GetBinContent(aTimeTileIndex*bandMultiple[aBandIndex]+1,aBandIndex+1);
  };
  inline double GetTileAmplitude(const int aTimeTileIndex, const int aBandIndex){    
    return GetTileSNR(aTimeTileIndex,aBandIndex)*sqrt(bandPower[aBandIndex]);
  };
  inline double GetTileTime(const int aTimeTileIndex, const int aBandIndex){
    return qplane->GetXaxis()->GetBinCenter(aTimeTileIndex*bandMultiple[aBandIndex]+1);
  };
  inline double GetTileTimeStart(const int aTimeTileIndex, const int aBandIndex){
    return qplane->GetXaxis()->GetBinLowEdge(aTimeTileIndex*bandMultiple[aBandIndex]+1);
  };
  inline double GetTileTimeEnd(const int aTimeTileIndex, const int aBandIndex){
    return qplane->GetXaxis()->GetBinUpEdge(aTimeTileIndex*bandMultiple[aBandIndex]+1);
  };
  inline int GetTimeTileIndex(const int aBandIndex, const double aTime){
    return (int)floor((aTime-GetTimeMin())/GetTileDuration(aBandIndex));
  };

  // SETS
  bool SetPower(Spectrum *aSpec);
  void SetTileSNR(const int aTimeTileIndex, const int aBandIndex, const double aSNR);
  void SetTileDisplay(void);

  // INTERNAL
  double UpdateThreshold(const int aBandIndex, double *aEnergies, double &aThreshold);
  void GetPlaneNormalization(void);
 
  // Q-PLANE
  double Q;                         ///< Q value
  double QPrime;                    ///< Q prime = Q / sqrt(11)
  TH2D *qplane;                     ///< Q plane
  int Ntiles;                       ///< number of tiles in the plane
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


