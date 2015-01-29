//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Otile__
#define __Otile__

#include "Oqplane.h"
#include "GwollumPlot.h"

using namespace std;

/**
 * Construct a time-frequency-Q tiling.
 * This class was designed to tile the 3-dimensional space in time, frequency and Q. The tiling consists of logarithmically spaced Q-planes. Each of these planes is divided in logarithmically spaced frequency bands. Each of these bands are then linearly divided in time bins. Once constructed, the planes can be used to apply a Q-transform to data. This class offers a graphical interface (GwollumPlot inheritance) and plotting functions to display the tiles and the data.
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
   * Some conditions are to be met to use this class:
   * - The time range must be a power of 2 and at least 4s long
   * - The sampling frequency must be a power of 2
   * - The Q value cannot be smaller than sqrt(11)
   * - The maximum mismatch cannot be larger than 0.5
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
   * Returns the number of Q planes.
   */
  inline int GetNQ(void){ return nq; };

  /**
   * Returns the Q value of a given plane.
   * -1.0 is returned if this function fails.
   * @param aQindex Q-plane index
   */
  double GetQ(const int aQindex);

  /**
   * Displays a canonical representation of a given Q-plane.
   * @param aQindex Q-plane index
   */
  bool DisplayTiling(const int aQindex);

  /**
   * Draws a given Q-plane.
   * @param aQindex Q-plane index
   */
  bool DrawPlane(const int aQindex);

  /**
   * Sets the data power spectrum.
   * This function must be called to compute the amplitude in a given tile: amplitude = SNR * sqrt(power).
   * The power of a tile is given by the input power spectrum weighted by the Gaussian window.
   * If this function is never called, the power is set to 0. 
   *
   * The PSD must be given as a valid Spectrum structure, i.e, the PSD was previously computed.
   * @param aSpec Spectrum structure where the PSD has been computed
   */
  bool SetPower(Spectrum *aSpec);

  /**
   * Projects a data vector onto the Q planes.
   * A complex data vector is projected onto all the Q-planes. The tiles are populated with the resulting SNR values.
   *
   * IMPORTANT: the input data vector must the right size, i.e. SampleFrequency/2 as defined in the constructor. No check will be performed!
   * @param aDataRe real part of the data vector (frequency domain)
   * @param aDataIm imaginary part of the data vector (frequency domain)
   */
  bool ProjectData(double *aDataRe, double *aDataIm);

  /**
   * Saves tiles above a SNR threshold in a MakeTriggers structure.
   * The triggers Segments are also saved follwing the GWOLLUM convention for triggers. A padding can be provided to NOT saved triggers on the plane edges. The planes are always centered on 0. A T0 must therefore be provided.
   * @param aTriggers MakeTriggers object
   * @param aSNRThr SNR threshold
   * @param aLeftTimePad duration of the left padding
   * @param aRightTimePad duration of the right padding
   * @param aT0 plane central time
   */
  bool SaveTriggers(MakeTriggers *aTriggers, const double aSNRThr, const int aLeftTimePad=0, const int aRightTimePad=0, const int aT0=0);

  /**
   * Saves the maps for each Q-planes in output files.
   * The maps are saved in output files.
   * @param aOutdir output directory path
   * @param aName name identifier
   * @param aT0 plane central time
   * @param aFormat output format string
   * @param aWindows list of time windows
   */
  bool SaveMaps(const string aOutdir, const string aName, const int aT0, const string aFormat, vector <int> aWindows);

  /**
   * Computes a set of Q values.
   * This function returns a vector of Q values corresponding to a set of parameters:
   * @param aQMin minimal Q value
   * @param aQMax maximal Q value
   * @param aMaximumMismatch maximum mismatch between Q planes
   */
  static vector <double> ComputeQs(const double aQMin, const double aQMax, const double aMaximumMismatch);

  
  /**
   * Returns the tiling time range.
   */
  inline double GetTimeRange(void){ return tilemap->GetXaxis()->GetXmax()-tilemap->GetXaxis()->GetXmin(); };


 private:

  int fVerbosity;              ///< verbosity level
  Omap *tilemap;               ///< maps: 0=combined, then individual Qplanes
  Oqplane **qplanes;           ///< Q planes
  int nq;                      ///< number of q planes

  //void MakeTiling(void);       /// fill tiling from Q planes
  //void SetTileContent(const double t1, const double f1, const double t2, const double f2, const double content);
    
  ClassDef(Otile,0)  
};

#endif


