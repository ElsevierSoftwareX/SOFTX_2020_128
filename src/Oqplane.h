//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Oqplane__
#define __Oqplane__

#include "MakeTriggers.h"
#include "Spectrum.h"
#include "Omap.h"

// eq 5.95 with alpha=2
#define BIASFACT2 (1.0-log(4.0*3.0*3.0)/(4.0*3.0*3.0-1.0))

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
		    const double aLeftTimePad=0.0, 
		    const double aRightTimePad=0.0, 
		    const double aT0=0.0);

  // GETS
  inline double GetQ(void){ return Q; };
  inline double GetTileAmplitude(const int aTimeTileIndex, const int aBandIndex){    
    return GetTileContent(aTimeTileIndex,aBandIndex)*sqrt(bandPower[aBandIndex]);
  };

  inline double GetTriggerFrac(void){    
    return (double)nTriggers/(double)Ntiles;
  }
  
  // SETS
  bool SetPower(Spectrum *aSpec);
  inline void SetSNRThr(const double aSNRThr){ SNRThr=aSNRThr; };


  // INTERNAL
  double GetMeanEnergy(const int aBandIndex, double *aEnergies);
  void GetPlaneNormalization(void);
 
  // Q-PLANE
  double Q;                         ///< Q value
  double QPrime;                    ///< Q prime = Q / sqrt(11)
  double PlaneNormalization;        ///< plane normalization

  // TRIGGER SELECTION
  double SNRThr;                    ///< SNR threshold to save triggers
  long int nTriggers;               ///< number of tiles above trigger SNR thr
  
  // FREQUENCY BANDS
  int *bandWindowSize;              ///< band 'Gaussian' window size
  double **bandWindow;              ///< band 'Gaussian' window
  double **bandWindowFreq;          ///< band 'Gaussian' window frequency
  double *bandPower;                ///< band power
  fft **bandFFT;                    ///< band fft
  
  ClassDef(Oqplane,0)  
};


#endif


