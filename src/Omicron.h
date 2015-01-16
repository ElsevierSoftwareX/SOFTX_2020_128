//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Omicron__
#define __Omicron__

#include "GwollumPlot.h"
#include "IO.h"
#include "ffl.h"
#include "Streams.h"
#include "Spectrum.h"
#include "Sample.h"
#include "MakeTriggers.h"
#include "EventMap.h"
#include "TMath.h"
#include "Otile.h"
#include "Odata.h"
#include "Date.h"

using namespace std;


/**
 * Process data with the Omicron algorithm.
 * An introduction to Omicron is available: <a href="../../Friends/omicron.html">Omicron introduction</a>
 *
 * This class was designed to offer various ways to run Omicron methods; either step-by step or in a all-in-one way. If the algorithm remains the same, Omicron provides two different outputs: triggers or maps. The triggers corresponds to tiles with a SNR value above a given threshold . The maps are a time-frequency decomposition of the input signal. Triggers are produced with the Process() method while maps are obtained with the Scan() method.
 *
 * \author Florent Robinet
 */
class Omicron {

 public:
  
  /**
   * @name Constructors and destructors
   @{
  */
  /**
   * Constructor of the Omicron class.
   * This constructor initializes all the components to run Omicron: data structures, data streams, tiling, triggers, maps, injections, monitoring, etc.
   * An option file is required to define all the parameters to run Omicron. For more details about Omicron configuration, see <a href="../../Friends/omicron.html">this page</a>.
   *
   * The Omicron class offers wrapping functions (see Process(), Scan(), ScanTriggers()) where every analysis steps are included. Omicron can also be used step-by-step for tailored applications (like a low-latency searches where data are provided sequentially when they are available).
   * @param aOptionFile path to the option file
   */
  Omicron(const string aOptionFile);

  /**
   * Destructor of the Omicron class.
   */
  virtual ~Omicron(void);
  /**
     @}
  */


  /**
   * Runs the full trigger analysis of data segments.
   * This function runs the Omicron algorithm over the data defined by the input segments. The data are segmented, conditioned, projected on the tiling structure and resulting triggers are saved on disk. This function calls the following sequence of Omicron functions:
   * - InitSegments() to load the segments to process
   * - MakeDirectories() to create the output directory structure
   * - NewChunk() to walk through the chunks
   * - NewChannel() to walk through the channels
   * - LoadData() to load a data vector
   * - ConditionVector() to condition the data vector
   * - ExtractTriggers() to project data onto the tiles and fill the Triggers structure
   * - WriteTriggers() to save triggers on disk
   *
   * This function is only available if a FFL structure has been previously declared.
   * @param aSeg Segments to process.
   */
  //bool Process(Segments *aSeg);

  /**
   * Runs the full scan analysis of a GPS time.
   * This function runs the Omicron algorithm over the data defined by a central time. The data are conditioned, projected on the tiling structure and resulting maps are saved on disk. This function calls the following sequence of Omicron functions:
   * - InitSegments() to load the chunk to process
   * - MakeDirectories() to create the output directory structure
   * - NewChunk() to load the chunk
   * - NewChannel() to walk through the channels
   * - LoadData() to load a data vector
   * - ConditionVector() to condition the data vector
   * - MakeMaps() to project data onto the tiles and fill the maps
   * - WriteMaps() to save maps on disk
   *
   * This function is only available if a FFL structure has been previously declared.
   * @param aTimeCenter central time of the maps
   */
  //bool Scan(const double aTimeCenter);

  /**
   * Initializes the segments to process.
   * This function should always be called before any type of processing. This function can also be used to introduce a time offset (for the plots only!).
   * 
   * WARNING: the input Segments object is not copied, only the pointer is used. This means that the Segments structure pointed by aSeg should not be modified or deleted before the end of the processing.
   * @param aSeg pointer to the input Segments structure
   * @param aTimeOffset time offset to define a new time origin for the plots
   */
  bool InitSegments(Segments *aSeg, const double aTimeOffset=0.0);

  /**
   * Modifies the output directory structure.
   * Two directory structures are possible:
   * - [path_to_outdir]/aGPS/[channel_name] if aGPS is not 0
   * - [path_to_outdir]/[channel_name] if aGPS is 0
   *
   * where [path_to_outdir] is the output directory specified by the user in the option file and [channel_name] is the channel name being processed.
   * 
   * The GPS value is rounded to the ms digit.
   *
   * If this function is never called, all the output is dumped in the output directory specified by the user in the option file.
   * @param aGPS GPS time
   */
  bool MakeDirectories(const double aGPS = 0);

  /**
   * Loads a new chunk.
   * The chunks are loaded following the time structure defined in the option file and the Segments object defined with InitSegments(). When there is not enough data to fill one chunk (end of a segment), the chunk duration is shortened as explained in Odata. This function should be called iteratively to cover the full data set. The segmentation procedure is detailed in Odata. The returned value indicates the status of this operation:
   * - true : a new chunk has been successfully loaded
   * - false : no more chunk to load
   */
  bool NewChunk(void);

  /**
   * Loads a new channel.
   * The channels defined in the option file are loaded incrementally. If this function is called after the last channel, false is returned and the channel sequence is reset: the next call will load the first channel.
   * This function should be called before any data processing.
   */
  bool NewChannel(void);

  /**
   * Loads a data vector.
   * The data vector of the current channel and the current chunk is loaded. If requested, the injection channel is added to the vector. This function loads the data from the frames listed in the FFL. The FFL option is therefore mandatory.
   * It is the user's responsibility to delete the returned data vector.
   *
   * If this function fails, a pointer to NULL is returned.
   * @param aDataVector pointer to the returned data vector
   * @param aDataVector sample size of the returned data vector
   */
  bool LoadData(double **aDataVector, int *aSize);

  /**
   * Conditions a data vector.
   * Before projecting the data onto the tiles, the data are conditioned with this funtion. The input data chunk is first re-sampled and highpassed. Then the data is used to estimate the noise (PSD). Finally, the data subsegments in the chunk are Tukey-windowed, Fourier-transformed and normalized by the ASD.
   *
   * IMPORTANT: The input vector size MUST MATCH the current chunk size loaded with NewChunk(). NO check is performed against that!
   *
   * If the returned integer value is negative, it means that a fatal error occured and the Omicron object got corrupted. If it is positive, it means that the conditioning ran into some errors but the Omicron object is still valid. If it is 0, the conditioning ended correctly. The error code is the following:
   * - -1 = the Omicron object is corrupted.
   * - -2 = the frequency parameters update of the Sample object failed.
   * -  0 = OK
   * -  1 = the input vector is null
   * -  2 = the input vector is empty
   * -  3 = the input vector is flat
   * -  4 = the vector transformation failed (resampling+highpassing)
   * -  5 = the spectrum could not be computed
   * -  6 = the tiling could not be normalized
   * @param aInVectSize input vector size
   * @param aInVect input data vector (time domain)
   */
  int Condition(const int aInVectSize, double *aInVect);
  
  /**
   * Projects whitened data onto the tiles and fills output structures.
   * Instead of calling specific output functions, like ExtractTriggers() or MakeMaps(), this function automatically produce the data products requested to be saved to disk (see WriteOutput())
   */
  bool ProjectAndSave(void);
  
  /**
   * Writes output to disk.
   * The output data products are selected by the user in the option file.
   * - the chunk triggers: triggers are saved if the the number of triggers is below the limit specified in the option file.
   * - the chunk maps
   * - the chunk time-series
   * - the chunk ASD
   * - the chunk PSD
   *
   * Note that this function writes the data product no matter if they were previously built or not (see  ExtractTriggers(), MakeMaps() or WriteOutput())
   */
  //bool WriteOutput(void);
  
  /**
   * Projects conditioned data onto the tiles and fills the Triggers structure.
   * The Triggers stucture is then sorted (by tstart) and the clustering algorithm is applied if requested.
   * This function returns the current number of triggers in memory. A negative value is returned if this function fails.
   */
  //int ExtractTriggers(void);
  
  
  /**
   * Writes current maps to disk.
   * The maps built with MakeMaps() are saved to disk with the formats defined in Omicron(). When graphical formats are selected, an additional file is saved on disk: '[channel]_mapsummary.root'. This file contains a TTree summarizing the properties of the maps.
   *
   * In addition to maps, some complementary plots are also saved (graphical formats only):
   * - SNR vs. time for the frequency where the SNR is maximal
   * - SNR vs. frequency for the time where the SNR is maximal
   * 
   * This function can also be used with an external set of Q-maps. In this case, additional information must be given in arguments:
   * @param aChNumber channel number of the current Q-maps (always required)
   * @param aTimeCenter GPS time where to set the time origin
   * @param aQ vector of Q values corresponding to each map
   * @param aQmap pointer to an array of Q-maps vor each Q values
   */
  //bool WriteMaps(const vector <double>& aQ = vector<double>(), TH2D **aQmap=NULL);

  /**
   * To be documented.
   * @param aTimeCenter central time of the maps
   */
  //bool ScanTriggers(const double aTimeCenter);

  /**
   * Create a html report for a scan.
   * @param aScanDir path to scan directory
   */
  //bool ScanReport(const string aScanDir="");


  /**
   * Returns the segments associated to the trigger time coverage.
   * See Triggers::GetTriggerSegments().
   *
   * NOTE: This function should be called somewhere after ExtractTriggers() and before WriteTriggers() so the triggers are present in memory.
   * @param aThr threshold object
   * @param aInfValue value above which the threshold is considered infinite
   */
  Segments* GetTriggerSegments(TH1D *aThr=NULL, const double aInfValue=1e20);
  
  /**
   * Returns list of channels to process.
   */
  inline vector <string> GetChannelList(void){return fChannels;};

  /**
   * Returns chunk duration [s].
   */
  inline int GetChunkDuration(void){return fChunkDuration;};

  /**
   * Returns segment/block duration [s].
   */
  inline int GetSegmentDuration(void){return fSegmentDuration;};

  /**
   * Returns overlap duration [s].
   */
  inline int GetOverlapDuration(void){return fOverlapDuration;};

  /**
   * Returns working sampling frequency [Hz]
   */
  inline int GetSampleFrequency(void){return fSampleFrequency;};
 
  /**
   * Prints a formatted message with a timer.
   * @param aMessage message to print
   */
  void PrintMessage(const string aMessage);

  /**
   * Prints a progress report of the processing.
   */
  void PrintStatusInfo(void);
  
  /**
   * Returns class status.
   */
  inline bool GetStatus(void){ return status_OK; };

  /**
   * Returns verbosity level.
   */
  inline int GetVerbosity(void){ return fVerbosity; };

 private:

  // PROCESS STATUS
  bool status_OK;               ///< general status
  time_t timer;                 ///< timer
  time_t timer_start;           ///< timer start
  struct tm * ptm;              ///< gmt time
  int chanindex;                ///< current channel index
  double timeoffset;            ///< current time offset

  // PARAMETERS
  bool ReadOptions(void);       ///< to parse option card
  void AdjustParameters(void);  ///< adjust default parameters
  string fOptionFile;           ///< option file name
  vector <string> fOptionName;  ///< option name (metadata)
  vector <string> fOptionType;  ///< option type (metadata)
  string fMaindir;              ///< main output directory
  vector <string> fOutdir;      ///< output directories per channel
  vector <string> fChannels;    ///< list of channel names
  vector <string> fInjChan;     ///< injection channel names
  vector <double> fInjFact;     ///< injection factors
  string fFflFile;              ///< path to FFL file
  string fTrigDir;              ///< path to trigger directory
  int fSampleFrequency;         ///< sampling frequency of input data
  vector <double> fFreqRange;   ///< frequency range
  vector <double> fQRange;      ///< Q range
  int fChunkDuration;           ///< chunk duration (varies!)
  int fSegmentDuration;         ///< segment duration
  int fOverlapDuration;         ///< overlap duration
  double fMismatchMax;          ///< maximum mismatch
  vector <int> fWindows;        ///< scan windows
  int fWindowMax;               ///< maximal window value
  string ffftplan;              ///< fft plan
  double fSNRThreshold;         ///< SNR Threshold
  int fNtriggerMax;             ///< trigger limit
  vector <string> fClusterAlgo; ///< clustering modes
  double fcldt;                 ///< clustering dt
  int fVerbosity;               ///< verbosity level
  string fOutFormat;            ///< output format
  string fStyle;                ///< plotting style
  string fOutProducts;          ///< output products
  string fWriteMode;            ///< write mode

  // PROCESS MONITORING
  Segments *inSegments;         ///< cumulative requested segments
  Segments **outSegments;       ///< segments currently written on disk
  int chunk_ctr;                ///< number of loaded chunks
  int *chan_ctr;                ///< number of times a channel was loaded
  int *chan_data_ctr;           ///< number of times a channel was found data
  int *chan_cond_ctr;           ///< number of times a channel was conditioned
  int *chan_wmaps_ctr;          ///< number of times a channel's maps were saved
  int *chan_wtrig_ctr;          ///< number of times a channel's triggers were saved

  // COMPONENTS
  Odata *dataseq;               ///< data sequence
  Sample **sample;              ///< sampling structures
  Streams **streams;            ///< streams
  Spectrum *spectrum;           ///< spectrum structure
  ffl *FFL;                     ///< ffl
  Otile *tile;                  ///< tiling structure
  MakeTriggers **triggers;      ///< output triggers

  // RAW DATA
  double *ChunkVect;            ///< chunk raw data (time domain)
  double *SegVect;              ///< subsegment raw data (time domain)
  
  // CONDITIONING
  bool Whiten(double **aDataRe, double **aDataIm); ///< whiten data vector
  double* GetTukeyWindow(const int aSize, const int aFractionSize); ///< create tukey window
  double *TukeyWindow;          ///< tukey window
  fft *offt;                    ///< FFT plan to FFT the input data
  double **dataRe;              ///< conditioned data (Re)
  double **dataIm;              ///< conditioned data (Im)

  // OUTPUT
  void SaveTriggers(void);          ///< save triggers

  // SCANS
  //TH2D **Qmap_full;             ///< combined Q-maps
  //int *loudest_qmap;            ///< Q-map containing the loudest tile
  string fScandir;              ///< latest scan directory
  //bool MakeMaps(const int aNseg);
  //void MakeFullMaps(const int aNq, TH2D **aQmap);

  // MISC
  GwollumPlot *GPlot;           ///< Gwollum plots
  void SaveAPSD(const string type="PSD");    ///< Save current PSD
  void SaveTS(void);///< Save current chunk time series
  bool *first_save;             ///< flags the first save
  void PrintASCIIlogo(void);    ///< print ascii logo

  ClassDef(Omicron,0)  
};

#endif


