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
   * @param aPlotStyle plotting style
   * @param aVerbosity verbosity level
   */
  Otile(const int aTimeRange, 
	const double aQMin, 
	const double aQMax, 
	const double aFrequencyMin, 
	const double aFrequencyMax, 
	const int aSampleFrequency, 
	const double aMaximumMismatch, 
	const string aPlotStyle="GWOLLUM", 
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
  inline double GetQ(const int aQindex){
    if(aQindex<0) return -1.0;
    if(aQindex>=nq) return -1.0;
    return qplanes[aQindex]->GetQ();
  };

  /**
   * Displays a canonical representation of a given Q-plane.
   * @param aQindex Q-plane index
   */
  bool DrawMapTiling(const int aQindex);

  /**
   * Draws a given Q-plane.
   * @param aQindex Q-plane index
   */
  bool DrawMapContent(const int aQindex);

  /**
   * Draws phase for a given Q-plane.
   * @param aQindex Q-plane index
   */
  bool DrawMapPhase(const int aQindex);

  /**
   * Sets the data power spectrum.
   * This function must be called to compute the amplitude in a given tile: amplitude = SNR * sqrt(power).
   * The power of a tile is given by the input power spectrum weighted by the Cone window.
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
   * If requested, tiles are de-activated if they overlap (in time or frequency) another tile. This process is call down-tiling.
   *
   * IMPORTANT: the input data vector must the right size, i.e. SampleFrequency/2 as defined in the constructor. No check will be performed!
   * @param aDataRe real part of the data vector (frequency domain)
   * @param aDataIm imaginary part of the data vector (frequency domain)
   * @param aTileDown apply down-tiling if set to true
   */
  bool ProjectData(double *aDataRe, double *aDataIm, const bool aTileDown=true);

  /**
   * Saves active tiles in a MakeTriggers structure.
   * The triggers Segments are also saved follwing the GWOLLUM convention for triggers. A padding can be provided to NOT saved triggers on the plane edges. The planes are always centered on 0. A T0 must therefore be provided.
   *
   * See also SetSaveSelection().
   * @param aTriggers MakeTriggers object
   * @param aLeftTimePad duration of the left padding
   * @param aRightTimePad duration of the right padding
   * @param aT0 plane central time
   */
  bool SaveTriggers(MakeTriggers *aTriggers, const double aLeftTimePad=0.0, const double aRightTimePad=0.0, const double aT0=0.0);

  /**
   * Saves the maps for each Q-planes in output files.
   * The maps are saved in output files.
   * An additionnal map called 'fullmap' is also saved. It represents active tiles projected in the time-frequency plane.
   * IMPORTANT NOTE: Maps are not saved if the maximum SNR within the first window time range is below the SNR threshold, see SetSaveSelection().
   * The returned value is the maximum SNR value within the first window time range. -1.0 is returned if this function fails.
   *
   * See also SetSaveSelection().
   * @param aOutdir output directory path
   * @param aName name identifier
   * @param aT0 plane central time
   * @param aFormat output format string
   * @param aWindows list of time windows
   * @param aThumb also produce thumbnails if set to true
   */
  double SaveMaps(const string aOutdir, const string aName, const int aT0, const string aFormat, vector <int> aWindows, const bool aThumb=false);

  /**
   * Computes a set of Q values.
   * This function returns a vector of Q values corresponding to a set of parameters:
   * @param aQMin minimal Q value
   * @param aQMax maximal Q value
   * @param aMaximumMismatch maximum mismatch between Q planes
   */
  static vector <double> ComputeQs(const double aQMin, const double aQMax, const double aMaximumMismatch);

  /**
   * Sets the maximum for the map SNR vertical scale.
   * If the SNR value is smaller than 1, an automatic scale is used.
   * @param aSNRScale maximal SNR value
   */
  inline void SetSNRScale(const int aSNRScale){ snrscale=aSNRScale; };

  /**
   * Sets a selection when saving maps and triggers.
   * This selection is applied when calling the SaveMaps() or SaveTriggers() functions.
   * @param aSNRThr_map when calling SaveMaps(), a map is not saved if the loudest tile is below that threshold
   * @param aSNRThr_trigger tiles with a SNR value below that threshold are not saved when calling SaveTriggers()
   * @param aNTriggerMax if, for a given Qplane, the number of tiles with a SNR above 'aSNRThr_trigger' is larger than this value, not a single tile of all planes can be saved. In other words, SaveTriggers() will return false
   */
  inline void SetSaveSelection(const double aSNRThr_map=-1.0, const double aSNRThr_trigger=2.0, const int aNTriggerMax=1000000){
    NTriggerMax=aNTriggerMax;
    SNRThr_map=aSNRThr_map;
    for(int q=0; q<nq; q++) qplanes[q]->SetSNRThr(aSNRThr_trigger);
  };

  /**
   * Returns the current SNR threshold for maps.
   * See SetSaveSelection().
   */
  inline double GetSNRMapThr(void){ return SNRThr_map; };

  /**
   * Returns the current SNR threshold for triggers.
   * See SetSaveSelection().
   */
  inline double GetSNRTriggerThr(void){ return qplanes[0]->GetSNRThr(); };

  /**
   * Returns the current limit on the number of triggers.
   * See SetSaveSelection().
   */
  inline int GetNTriggerMax(void){ return NTriggerMax; };

  /**
   * Returns the lowest frequency of this tiling.
   * The minimum frequency of the lowest Q plane is returned.
   */
  inline double GetFrequencyMin(void){ return qplanes[0]->GetFrequencyMin(); };

  /**
   * Returns the highest frequency of this tiling.
   * The maximum frequency of the highest Q plane is returned.
   */
  inline double GetFrequencyMax(void){ return qplanes[nq-1]->GetFrequencyMax(); };

  /**
   * Returns the maximum mismatch between tiles.
   */
  inline double GetMismatchMax(void){ return MaximumMismatch; };

 private:

  int fVerbosity;           ///< verbosity level
  double MaximumMismatch;   ///< maximum mismatch
  Oqplane **qplanes;        ///< Q planes
  int nq;                   ///< number of q planes
  int TimeRange;            ///< map time range
  int snrscale;             ///< map snr scale
  double SNRThr_map;        ///< map SNR threshold
  int NTriggerMax;          ///< max. number of tiles to save (triggers only)
  
  TH2D* MakeFullMap(const int aTimeRange, const double aT0=0.0); ///< make full map
  void TileDown(void);         ///< tile-down
  void ApplyOffset(TH2D *aMap, const double aOffset);

  ClassDef(Otile,0)  
};

#endif


