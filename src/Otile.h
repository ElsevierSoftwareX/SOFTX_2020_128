//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Otile__
#define __Otile__

#include "Oqplane.h"
//#include "MakeTriggers.h"
#include "GwollumPlot.h"

using namespace std;

/**
 * Construct a time-frequency-Q tiling.
 * This class was designed to tile the 3-dimensional space in time, frequency and Q. The tiling consists of logarithmically spaced Q-planes. Each of these planes is divided in logarithmically spaced frequency bands. Each of these bands are then linearly divided in time bins.
 * \author    Florent Robinet
 */
class Otile: public GwollumPlot {

 public:

  /**
   * @name Constructors and destructors
   @{
  */
  /**
   * Constructor of the Otile class.
   * The tiling is constructed given the user parameters. The parameter space is defined by a time range, a frequency range and a Q range. The user must specify a maximum mismatch value corresponding to a maximal fractional energy loss from one tile to the next.
   *
   * Some conditions are to be met:
   * - The time range must be a power of 2
   * - The sampling frequency must be a power of 2
   * @param aTimeRange time range [s]
   * @param aQMin minimal Q value
   * @param aQMax maximal Q value
   * @param aFrequencyMin minimal frequency [Hz]
   * @param aFrequencyMax maximal frequency [Hz]
   * @param aSampleFrequency sampling frequency [Hz]
   * @param aMaximumMismatch maximum mismatch between tiles
   * @param aVerbosity verbosity level
   */
  Otile(const int aTimeRange, 
	const double aQMin, 
	const double aQMax, 
	const double aFrequencyMin, 
	const double aFrequencyMax, 
	const int aSampleFrequency, 
	const double aMaximumMismatch, 
	const int aVerbosity=0);

  /**
   * Destructor of the Otile class.
   */
  virtual ~Otile(void);// destructor
  /**
     @}
  */
  
  /**
   * Sets the data power spectrum.
   * This function must be called before the GetTriggers() or GetMap() functions. The power spectrum must be set to compute the trigger amplitude defined as SNR*sqrt(power). power is the weighted average of the PSD over the tile (weighted by the tile Gaussian window).
   * If this function is not called, power is set to 1 and the amplitude is just the SNR.
   *
   * The PSD must be given as a valid Spectrum structure, i.e, the PSD was previously computed.
   * @param aSpec Spectrum structure where the PSD has been computed
   */
  bool SetPower(Spectrum *aSpec);

  bool ProjectData(double *aDataRe, double *aDataIm);
  bool SaveTriggers(MakeTriggers *aTriggers, const double aSNRThr, const int aLeftTimePad=0, const int aRightTimePad=0, const int aT0=0);
  bool SaveMaps(const string aOutdir, const string aName, const int aT0, const string aFormat, vector <int> aWindows);

  /**
   * Saves triggers above SNR threshold.
   * A complex data vector is projected onto all Q-planes. The tiles are populated with the SNR values. A trigger is defined as a tile with a SNR value above the threshold defined in Otile(). This trigger is added to the Triggers object with the <a href="../../Main/convention.html#triggers">GWOLLUM convention</a> with the following parameters:
   * - time = tile central time
   * - frequency = tile central frequency
   * - tstart = tile starting time
   * - tend = tile ending time
   * - fstart = tile starting frequency
   * - fend = tile ending frequency
   * - snr = tile SNR value
   * - amplitude = SNR*sqrt(power) where power is the noise power average in the tile set with SetPowerSpectrum(). 
   *
   * A time offset must be given corresponding to the starting time of the time range.
   * @param aTriggers MakeTriggers object
   * @param aDataRe real part of the data vector (frequency domain)
   * @param aDataIm imaginary part of the data vector (frequency domain)
   * @param aTimeStart time offset
   */
  //bool GetTriggers(MakeTriggers *aTriggers, double *aDataRe, double *aDataIm, const int aTimeStart, const int aExtraTimePadMin=0);

  /**
   * Returns a SNR map of data projected onto a given Q plane. 
   * A complex data vector is projected onto the Q-plane with index qindex. The tiles are populated with the SNR values and the resulting map is returned as a TH2 histogram.
   *
   * The input complex vector must be given in the frequency domain and have the right size, i.e. half of the sampling frequency times the time range.
   *
   * If given, the time offset is added to the time of the tiles.
   * 
   * If printamplitude=true, the amplitude is used instead of the SNR.
   *
   * Do not delete the returned histogram. It will be deleted by the class.
   * @param qindex Q-plane index
   * @param aDataRe real part of the data vector (frequency domain)
   * @param aDataIm imaginary part of the data vector (frequency domain)
   * @param time_offset time offset [s]
   * @param printamplitude switch to amplitude
   */
  //TH2D* GetMap(const int qindex, double *aDataRe, double *aDataIm, const double time_offset=0.0, const bool printamplitude=false);

  /**
   * Returns the Q value of plane with index 'qindex'.
   * -1.0 is returned if this function fails.
   * @param qindex Q-plane index
   */
  double GetQ(const int qindex);

  /**
   * Returns the number of Q planes
   */
  inline int GetNQPlanes(void){ return (int)Qs.size(); };

  /**
   * Returns a given Q plane.
   */
  inline TH2D* GetQPlane(const int qindex){ return qplanes[qindex]->qplane; };

  /**
   * Draws a given Q plane.
   */
  inline void DrawQPlane(const int qindex){ return Draw(qplanes[qindex]->qplane,"COLZ"); };
  
  /**
   * Computes a set of Q values.
   * This function returns a vector of Q values corresponding to a set of parameters:
   * @param aQMin minimal Q value
   * @param aQMax maximal Q value
   * @param aMaximumMismatch maximum mismatch between Q planes
   */
  static vector <double> ComputeQs(const double aQMin, const double aQMax, const double aMaximumMismatch);

  inline double GetTimeRange(void) { return maps[0]->GetXaxis()->GetXmax()-maps[0]->GetXaxis()->GetXmin(); };


 private:

  int fVerbosity;              ///< verbosity level
  TH2D **maps;                 ///< maps
  vector <double> Qs;          ///< vector of Qs
  Oqplane **qplanes;           ///< Q planes

  void MakeTiling(void);       /// fill tiling from Q planes
  void SetTileContent(const double t1, const double f1, const double t2, const double f2, const double content);
    
  ClassDef(Otile,0)  
};

#endif


