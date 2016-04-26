//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Oinject__
#define __Oinject__
#include "Spectrum.h"
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
  
  /**
   * @name Constructors and destructors
   @{
  */
  /**
   * Constructor of the Oinject class.
   * The ranges for the parameters are set to default values:
   * - time: fixed at 0
   * - frequency: 32-2048 Hz
   * - Q: 4-100
   * - amplitude: fixed at 1e-21
   *
   * A set of random parameters is generated with MakeWaveform().
   * The user must specify the duration of the injection waveform.
   * @param aDuration waveform duration [s]
   */
  Oinject(const double aDuration);

  /**
   * Destructor of the Oinject class.
   */
  virtual ~Oinject(void);
  /**
     @}
  */

  
  /**
   * Generates a new set of waveform parameters.
   * The waveform parameters are randomly generated:
   * - time: uniform
   * - frequency: uniform
   * - Q: uniform
   * - amplitude: uniform in log
   * - phase: uniform
   */
  void MakeWaveform(void);

  /**
   * Returns the waveform value for a given index.
   * The waveform value is computed according to the set of parameters previously generated with MakeWaveform(). The user must specify the frequency at which the data vector is sampled. This way, combined with the duration set in the constructor, the index can be converted to a time value.
   * @param aIndex waveform vector index
   * @param aSamplingFrequency sampling frequency [Hz]
   */  
  inline double GetWaveform(const int aIndex, const int aSamplingFrequency){
    return
      amp*
      Wg/sigma_t/sqrt2pi*exp(-(-duration/2.0+(double)aIndex/(double)aSamplingFrequency-tau)*(-duration/2.0+(double)aIndex/(double)aSamplingFrequency-tau)/2.0/sigma_t/sigma_t)* // Gaussian window
      TMath::Cos(2*TMath::Pi()*phi*(-duration/2.0+(double)aIndex/(double)aSamplingFrequency)+phase);// sine
  };

  /**
   * Returns the true value of SNR.
   * @param aSpec noise spectrum
   */
  double GetTrueSNR(Spectrum *aSpec);

  /**
   * Sets a new time range.
   * @param aTimeMin minimum time [s]
   * @param aTimeMax maximum time [s]
   */
  inline void SetTimeRange(const double aTimeMin, const double aTimeMax){
    taumin=aTimeMin; taumax=aTimeMax;
  };

  /**
   * Sets a new frequency range.
   * @param aFreqMin minimum frequency [Hz]
   * @param aFreqMax maximum frequency [Hz]
   */
  inline void SetFrequencyRange(const double aFreqMin, const double aFreqMax){
    phimin=aFreqMin; phimax=aFreqMax;
  };

  /**
   * Sets a new amplitude range.
   * @param aAmpMin minimum amplitude
   * @param aAmpMax maximum amplitude
   */
  inline void SetAmplitudeRange(const double aAmpMin, const double aAmpMax){
    ampmin=aAmpMin; ampmax=aAmpMax;
  };

  /**
   * Sets a new Q range.
   * @param aQMin minimum Q
   * @param aQMax maximum Q
   */
  inline void SetQRange(const double aQMin, const double aQMax){
    Qmin=aQMin; Qmax=aQMax;
  };

  /**
   * Returns the minimum time.
   */  
  inline double GetTimeMin(void){ return taumin; };

  /**
   * Returns the maximum time.
   */  
  inline double GetTimeMax(void){ return taumax; };

  /**
   * Returns the minimum frequency.
   */  
  inline double GetFrequencyMin(void){ return phimin; };

  /**
   * Returns the maximum frequency.
   */  
  inline double GetFrequencyMax(void){ return phimax; };

  /**
   * Returns the minimum Q.
   */  
  inline double GetQMin(void){ return Qmin; };

  /**
   * Returns the maximum Q.
   */  
  inline double GetQMax(void){ return Qmax; };

  /**
   * Returns the minimum amplitude.
   */  
  inline double GetAmplitudeMin(void){ return ampmin; };

  /**
   * Returns the maximum amplitude.
   */  
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


