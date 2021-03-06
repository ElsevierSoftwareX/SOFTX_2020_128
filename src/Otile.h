//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Otile__
#define __Otile__

#include "Oqplane.h"
#include <GwollumPlot.h>

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
   * This function must be called to compute the amplitude in a given tile: amplitude = SNR * (band noise amplitude).
   * The noise power of a tile is given by the input noise power spectrum averaged across the bisquare window.
   * If this function is never called, the power is set to 0. 
   * Two spectra must be given as a double whitening is into play.
   *
   * The PSDs must be given as valid Spectrum structures, i.e, the PSDs were previously computed.
   * @param aSpec1 Spectrum structure where the PSD has been computed (1st)
   * @param aSpec2 Spectrum structure where the PSD has been computed (2nd)
   */
  bool SetPower(Spectrum *aSpec1, Spectrum *aSpec2);

  /**
   * Projects a data vector onto the Q planes.
   * The complex data vector is projected onto all the Q-planes.
   * The data are provided through a fft object. The fft:Forward() must be done before calling this function.
   * The number of tiles (excluding overlaps/2) above the SNR threshold is returned.
   *
   * IMPORTANT: the input data vector must the right size, i.e. SampleFrequency/2 as defined in the constructor. No check will be performed!
   * @param aDataFft fft structure containing the data to project
   */
  int ProjectData(fft *aDataFft);

  /**
   * Saves tiles in a MakeTriggers structure.
   * Tiles with a SNR value above the SNR threshold are saved in the input trigger structure.
   * The corresponding triggers Segments are also saved following the GWOLLUM convention for triggers. If the Sequence algorithm is in use, the current timing is applied to the tiling.
   *
   * A time selection is performed if specific output segments were previously set with SetSegments(): triggers the time of which is outside the output segment list are not saved.
   *
   * See also SetSNRThr() and SetSegments().
   * @param aTriggers MakeTriggers object
   */
  bool SaveTriggers(TriggerBuffer *aTriggers);

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
   * @param aTimeOffset Time offset.
   * @param aThumb also produce thumbnails if set to true
   */
  double SaveMaps(const string aOutdir, const string aName, const string aFormat, vector <int> aWindows, const double aTimeOffset=0.0, const bool aThumb=false);

  /**
   * Computes a set of Q values.
   * This function returns a vector of Q values corresponding to a set of parameters:
   * @param aQMin minimal Q value
   * @param aQMax maximal Q value
   * @param aMaximumMismatch maximum mismatch between Q planes
   */
  static vector <double> ComputeQs(const double aQMin, const double aQMax, const double aMaximumMismatch);

  /**
   * Defines how to fill the maps.
   * - "snr": fill with SNR values
   * - "amplitude": fill with amplitude values
   * - "phase": fill with phase values
   * @param aMapFill fill type
   */
  inline void SetMapFill(const string aMapFill="snr"){
    if(!aMapFill.compare("amplitude")) mapfill="amplitude";
    else if(!aMapFill.compare("phase")) mapfill="phase";
    else mapfill="snr";
    for(int q=0; q<nq; q++) qplanes[q]->GetZaxis()->SetTitle(StringToUpper(mapfill).c_str());
    return;
  };

  /**
   * Sets the map vertical range.
   * If aZmin>=aZmax, the map is automatically ranged.
   * @param aZmin minimum Z value
   * @param aZmax maximum Z value
   */
  inline void SetRangez(const double aZmin=-1.0, const double aZmax=-1.0){
    vrange[0]=aZmin;
    vrange[1]=aZmax;
  };

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
   * Returns the current map fill type.
   * See SetMapFill().
   */
  inline string GetMapFill(void){ return mapfill; };

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

  
  /**
   * Sets the chirp mass [solar mass].
   * Use a negative value, not to draw the chirp track.
   * If the GPS time is negative, the merger time is taken at the center of the  timing window.
   * @param[in] aMchirp Chirp mass in solar masses.
   * @param[in] aMergerTime merger GPS time.
   */
  inline void SetChirp(const double aMchirp=-1.0, const double aMergerTime=-1.0){
    mchirp=aMchirp;
    tchirp=aMergerTime;
    chirp->SetParameter(0, -8.0*96.0/3.0/5.0*TMath::Power(TMath::Pi(),8.0/3.0)*TMath::Power(TMath::G()*mchirp*1.989e30/TMath::C()/TMath::C()/TMath::C(), 5.0/3.0));
    if(tchirp>0) chirp->SetParameter(1, tchirp);
  }

  /**
   * Returns the chirp mass [solar mass].
   */
  inline double GetChirpMass(void){ return mchirp; };

 private:

  int fVerbosity;               ///< verbosity level
  double MaximumMismatch;       ///< maximum mismatch
  Oqplane **qplanes;            ///< Q planes
  int nq;                       ///< number of q planes
  int TimeRange;                ///< map time range
  double vrange[2];             ///< map vertical range
  double SNRThr_map;            ///< map SNR threshold
  string mapfill;               ///< map fill type
  int **t_snrmax;               ///< loudest time tile (SNR)
  int **f_snrmax;               ///< loudest frequency tile (SNR)
  TF1 *chirp;                   ///< chirp track
  double mchirp;                ///< chirp mass in solar masses
  double tchirp;                ///< chirp merger GPS time.

  TH2D* MakeFullMap(const int aTimeRange, const double aTimeOffset); ///< make full map
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


