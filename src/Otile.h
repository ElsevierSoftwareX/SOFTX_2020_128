//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Otile__
#define __Otile__

#include "Oqplane.h"
#include "GwollumPlot.h"

using namespace std;

/**
 * Construct and apply a time-frequency-Q analysis.
 * This class was designed to tile the 3-dimensional space in time, frequency and Q. The tiling consists of logarithmically spaced Q-planes. Each of these planes is divided in logarithmically spaced frequency bands. Each of these bands are then linearly divided in time bins. Once constructed, the planes can be used to apply a Q-transform data segments.
 *
 * This class offers an algorithm, called Sequence, to read an input segment list sequentially. The Segments object is divided into overlapping time chunks matching the tiling duration. The chunks are loaded sequentially any time the NewChunk() function is called. The chunk sequence can be represented in the following way:
 * \verbatim
------------------------------------------------------------ current segment
 |------------------| chunk i-1
                |------------------| chunk i
                               |------------------| chunk i+1
 
                |---| overlap
 \endverbatim
 *
 * In general, the Segments object contain multiple time segments. The sequence described above does not necessarily match the size of the input segments. The sequence algorithm is designed to deal with such edge effects. Firstly, segments shorter than the tiling duration are skipped. When calling NewChunk() for the last chunk of a segment, the overlap duration is adjusted to fit the leftover:
 * \verbatim
 -----------------------------------------|   <--- input segment under processing

    |--------------------------|              <--- penultimate chunk 
  
 ###### call NextChunk() to cover the left-over

               |--------------------------|   <--- last chunk
	       |---------------|              <--- adjusted overlap
 * \endverbatim 
 * Obviously, the user must be careful about this special case as the overlap duration is modified (the chunk duration is never changed). Some functions are available to monitor the overlap size.
 *
 * When moving to a new segment, the overlap duration is set back to nominal values.
 *
 * This class offers a graphical interface (GwollumPlot inheritance) and plotting functions to display the tiles and the data.
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
   *
   * The tiling parameters (time, Q and frequency range and the mismatch max) can be internally modified to match these conditions.
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
   * The data are provided through a fft object. The fft:Forward() must be done before calling this function.
   *
   * IMPORTANT: the input data vector must the right size, i.e. SampleFrequency/2 as defined in the constructor. No check will be performed!
   * @param aFft fft structure containing the data to project
   */
  bool ProjectData(fft *aDataFft);

  /**
   * Saves active tiles in a MakeTriggers structure.
   * The triggers Segments are also saved following the GWOLLUM convention for triggers. If the Sequence algorithm is in use, the current timing is applied to the tiling.
   *
   * A time selection is performed if specific output segments were previously set with SetSegments(): triggers starting outside the output segment list are not saved.
   *
   * See also SetSNRThr() and SetSegments().
   * @param aTriggers MakeTriggers object
   */
  bool SaveTriggers(MakeTriggers *aTriggers);

  /**
   * Saves the maps for each Q-planes in output files.
   * The maps are saved in output files.
   * An additionnal map called 'fullmap' is also saved. It represents active tiles projected in the time-frequency plane.
   * IMPORTANT NOTE: Maps are not saved if the maximum SNR within the first window time range is below the SNR threshold, see SetSNRThr().
   * The returned value is the maximum SNR value within the first window time range. -1.0 is returned if this function fails.
   *
   * See also SetSNRThr().
   * @param aOutdir output directory path
   * @param aName name identifier
   * @param aFormat output format string
   * @param aWindows list of time windows
   * @param aThumb also produce thumbnails if set to true
   */
  double SaveMaps(const string aOutdir, const string aName, const string aFormat, vector <int> aWindows, const bool aThumb=false);

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
   * Sets a SNR threshold when saving maps and triggers.
   * The thresholds are applied when calling the SaveMaps() or SaveTriggers() functions.
   * @param aSNRThr_map when calling SaveMaps(), a map is not saved if the loudest tile is below that threshold
   * @param aSNRThr_trigger tiles with a SNR value below that threshold are not saved when calling SaveTriggers()
   */
  inline void SetSNRThr(const double aSNRThr_map=0.0, const double aSNRThr_trigger=2.0){
    SNRThr_map=aSNRThr_map;
    for(int q=0; q<nq; q++) qplanes[q]->SetSNRThr(aSNRThr_trigger);
  };

  /**
   * Returns the current SNR threshold for maps.
   * See SetSNRThr().
   */
  inline double GetSNRMapThr(void){ return SNRThr_map; };

  /**
   * Returns the current SNR threshold for triggers.
   * See SetSNRThr().
   */
  inline double GetSNRTriggerThr(void){ return qplanes[0]->GetSNRThr(); };

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

  /**
   * Returns the time range.
   */
  inline int GetTimeRange(void){ return TimeRange; };

  /**
   * Sets new input/output segments.
   * The input list of segments will be read sequencially using Sequence.
   * The input segment times must be integer numbers. They will be considered as such!
   *
   * Optionally, an ouput segment list can be provided to select triggers when calling SaveTriggers(). If pointing to NULL, no time selection is performed.
   * @param aInSeg input segment list
   * @param aOutSeg output segment list
   */
  bool SetSegments(Segments *aInSeg, Segments *aOutSeg=NULL);

  /**
   * Sets a new sequence overlap duration.
   * The input parameter can be modified to be an even number.
   * @param aOverlapDuration new overlap duration [s]
   */
  inline void SetOverlapDuration(const int aOverlapDuration){
    SeqOverlap=aOverlapDuration+aOverlapDuration%2;
  };

  /**
   * Loads a new sequence chunk.
   * The chunks are loaded following the definition presented in the description of this class. This function should be called iteratively to cover the full data set defined with SetSegments(). The returned value indicates the status of this operation:
   * - true : a new chunk has been loaded
   * - false : no more chunk to load or an error occured
   * @param aNewSegFlag set to true if a new segment is started 
   */
  bool NewChunk(bool &aNewSegFlag);

  /**
   * Returns the GPS center time of current chunk.
   */
  inline int GetChunkTimeCenter(void){ return SeqT0; };
  
  /**
   * Returns the GPS starting time of current chunk.
   */
  inline int GetChunkTimeStart(void){ return SeqT0-TimeRange/2; };
  
  /**
   * Returns the GPS ending time of current chunk.
   */
  inline int GetChunkTimeEnd(void){ return SeqT0+TimeRange/2; };
  
  /**
   * Returns the current overlap duration.
   * In most cases the overlap duration is nominal unless the special case of the end of an input segment is hit.
   */
  inline int GetCurrentOverlapDuration(void){ return SeqOverlapCurrent; };

  /**
   * Returns the nominal overlap duration.
   */
  inline int GetOverlapDuration(void){ return SeqOverlap; };


 private:

  int fVerbosity;               ///< verbosity level
  double MaximumMismatch;       ///< maximum mismatch
  Oqplane **qplanes;            ///< Q planes
  int nq;                       ///< number of q planes
  int TimeRange;                ///< map time range
  int snrscale;                 ///< map snr scale
  double SNRThr_map;            ///< map SNR threshold
    
  TH2D* MakeFullMap(const int aTimeRange); ///< make full map
  void ApplyOffset(TH2D *aMap, const double aOffset);

  // SEQUENCE
  Segments *SeqOutSegments;     ///< output trigger segments (current - request)
  Segments *SeqInSegments;      ///< input segments (current - request)
  int SeqOverlap;               ///< nominal overlap duration
  int SeqOverlapCurrent;        ///< current overlap duration
  int SeqT0;                    ///< current chunk center
  int SeqSeg;                   ///< current segment index

  ClassDef(Otile,0)  
};

#endif


