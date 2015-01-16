//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Oqplane__
#define __Oqplane__

#include "MakeTriggers.h"
#include "Spectrum.h"
#include "FFT.h"
#include "Omap.h"

using namespace std;

/**
 * Create a time-frequency Q-plane.
 * This class was designed to create and use a time-frequency Q-plane defined by a Q value. This class is entirely private and can only be used through the Otile class.
 *
 * \author    Florent Robinet
 */
class Oqplane: public Omap {

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
  inline double GetQ(void){ 
    return Q;
  };
  inline double GetTileAmplitude(const int aTimeTileIndex, const int aBandIndex){    
    return GetTileContent(aTimeTileIndex,aBandIndex)*sqrt(bandPower[aBandIndex]);
  };

  // SETS
  bool SetPower(Spectrum *aSpec);

  // INTERNAL
  double UpdateThreshold(const int aBandIndex, double *aEnergies, double &aThreshold);
  void GetPlaneNormalization(void);
 
  // Q-PLANE
  double Q;                         ///< Q value
  double QPrime;                    ///< Q prime = Q / sqrt(11)
  double PlaneNormalization;        ///< plane normalization
  
  // FREQUENCY BANDS
  int *bandWindowSize;              ///< band Gaussian window size
  double **bandWindow;              ///< band Gaussian window
  double **bandWindowFreq;          ///< band Gaussian window frequency
  double *bandPower;                ///< band power
  fft **bandFFT;                    ///< band fft
  
  ClassDef(Oqplane,0)  
};


#endif


