//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Oinject__
#define __Oinject__
#include "CUtils.h"
#include "TRandom3.h"
#include "TMath.h"

#define sqrt2pi 2.5066282746310002416123552393401042

using namespace std;

/**
 * Inject sinusoidal Gaussian waveforms
 *
 * \author    Florent Robinet
 */
class Oinject{

 public:
  
  Oinject(const double aDuration);
  virtual ~Oinject(void);

  void MakeWaveform(void);
  inline double GetWaveform(const int aIndex, const int aSamplingFrequency){
    return
      amp*
      Wg/sigma_t/sqrt2pi*exp(-(-duration/2.0+(double)aIndex/(double)aSamplingFrequency-tau)*(-duration/2.0+(double)aIndex/(double)aSamplingFrequency-tau)/2.0/sigma_t/sigma_t)* // Gaussian window
      TMath::Cos(2*TMath::Pi()*phi*(-duration/2.0+(double)aIndex/(double)aSamplingFrequency)+phase);// sine
  };

  inline void SetTimeRange(const double aTimeMin, const double aTimeMax){
    taumin=aTimeMin; taumax=aTimeMax;
  };
  inline void SetFrequencyRange(const double aFreqMin, const double aFreqMax){
    phimin=aFreqMin; phimax=aFreqMax;
  };
  inline void SetAmplitudeRange(const double aAmpMin, const double aAmpMax){
    ampmin=aAmpMin; ampmax=aAmpMax;
  };
  inline void SetQRange(const double aQMin, const double aQMax){
    Qmin=aQMin; Qmax=aQMax;
  };

  inline double GetTimeMin(void){ return taumin; };
  inline double GetTimeMax(void){ return taumax; };
  inline double GetFrequencyMin(void){ return phimin; };
  inline double GetFrequencyMax(void){ return phimax; };
  inline double GetQMin(void){ return Qmin; };
  inline double GetQMax(void){ return Qmax; };
  inline double GetAmplitudeMin(void){ return ampmin; };
  inline double GetAmplitudeMax(void){ return ampmax; };


 private:
  
  double duration;                        ///< vector duration
  TRandom3 *randgen;                   ///< random generator

  double Wg;                           ///< Gaussian window normalization
  double sigma_t;                      ///< Gaussian window width

  void GenerateParameters(void);       ///< generate random set of parameters
  double tau;                          ///< injection time
  double taumin;                       ///< minimum injection time
  double taumax;                       ///< maximum injection time
  double phi;                          ///< injection frequency
  double phimin;                       ///< minimum injection frequency
  double phimax;                       ///< maximum injection frequency
  double Q;                            ///< injection Q
  double Qmin;                         ///< minimum injection Q
  double Qmax;                         ///< maximum injection Q
  double amp;                          ///< injection amplitude
  double ampmin;                       ///< minimum injection amplitude
  double ampmax;                       ///< maximum injection amplitude
  double snr;                          ///< injection SNR
  double phase;                        ///< injection phase

  ClassDef(Oinject,0)  
};


#endif


