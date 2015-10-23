//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Omicron__
#define __Omicron__

#include "IO.h"
#include "Otile.h"
#include "Odata.h"
#include "Date.h"
#include "InjEct.h"
#include "ffl.h"

using namespace std;


/**
 * Process data with the Omicron algorithm.
 * An introduction to Omicron is available: <a href="../../Friends/omicron.html">Omicron introduction</a>
 *
 * This class was designed to offer various methods to conduct an Omicron analysis.
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
   * This constructor initializes all the components to run Omicron: data structures, data streams, tiling, triggers, injections, monitoring, etc.
   * An option file is required to define all the parameters to run Omicron. For more details about Omicron configuration, see <a href="../../Friends/omicron.html">this page</a>.
   *
   * After initialization, the Omicron methods should be called sequentially to perform the analysis. Here is a typical sequence:
   * - InitSegments() defines the data segments to process.
   * - MakeDirectories() creates a specific directory tree for the output 
(optional).
   * - NewChunk() loads a new chunk of data (loop #1).
   * - NewChannel() loads a new channel (loop #2).
   * - LoadData() loads the data vector for this chunk and this channel from FFL file (in loop #1/2)
   * - Condition() conditions data vector (in loop #1/2)
   * - Project() projects data onto the tiles (in loop #1/2)
   * - WriteOutput() writes output data productes to disk (in loop #1/2)
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
   * Initializes the segments to process.
   * This function should always be called before any type of processing. This function can also be used to introduce a time offset (for graphical plots only!).
   * 
   * WARNING: the input Segments object is not copied, only the pointer is used. This means that the Segments structure pointed by aSeg should not be modified or deleted before the end of the processing.
   * @param aSeg pointer to the input Segments structure
   * @param aTimeOffset time offset to define a new time origin for the graphical plots
   */
  bool InitSegments(Segments *aSeg, const double aTimeOffset=0.0);

  /**
   * Creates a specific directory tree for the output.
   * Two directory structures are possible:
   * - [path_to_outdir]/aId/[channel_name] if aId is not 0
   * - [path_to_outdir]/[channel_name] if aId is 0
   *
   * where [path_to_outdir] is the output directory specified by the user in the option file and [channel_name] is the channel name being processed.
   * 
   * The GPS value is rounded to the ms digit.
   *
   * If this function is never called, all the output is dumped in the output directory specified by the user in the option file.
   * @param aId directory id
   */
  bool MakeDirectories(const int aId = 0);

  /**
   * Loads a new chunk.
   * The chunks are loaded following the time structure defined in the option file and the Segments object defined with InitSegments(). When there is not enough data to fill one chunk (end of a segment), the chunk duration is shortened as explained in Odata. This function should be called iteratively to cover the full data set. The segmentation procedure is detailed in Odata. The returned value indicates the status of this operation:
   * - true : a new chunk has been successfully loaded
   * - false : no more chunk to load
   */
  bool NewChunk(void);

  /**
   * Loads a new channel.
   * The channels defined in the option file are loaded incrementally. If this function is called after the last channel, false is returned and the channel sequence is reset: the next call will load the first channel again.
   */
  bool NewChannel(void);

  /**
   * Loads a data vector.
   * The data vector of the current channel and the current chunk is loaded. If requested in the option file, the injection channel is added to the vector. This function loads the data from the frames listed in the FFL. The FFL option is therefore mandatory to use this function.
   * It is the user's responsibility to delete the returned data vector.
   *
   * If this function fails, a pointer to NULL is returned.
   * @param aDataVector pointer to the data vector
   * @param aSize sample size of the data vector
   */
  bool LoadData(double **aDataVector, int *aSize);

  /**
   * Conditions a data vector.
   * Before projecting the data onto the tiles, the data are conditioned with this function. The input data chunk is first resampled and highpassed. If requested in the option file, software injection waveforms are added to the data vector. Then the data is used to estimate the noise (PSD). Finally, the data subsegments in the chunk are Tukey-windowed, Fourier-transformed and normalized by the ASD.
   *
   * IMPORTANT: The input vector size MUST MATCH the current chunk size loaded with NewChunk(). NO check is performed against that!
   *
   * If the returned value is negative, it means that a fatal error occured and the Omicron object got corrupted. If it is positive, it means that the conditioning failed but the Omicron object is still valid for further use. If it is 0, the conditioning ended correctly. The error code is the following:
   * - -1 = the Omicron object is corrupted.
   * -  0 = OK
   * -  1 = the input vector is null
   * -  2 = the input vector is empty
   * -  3 = the input vector is flat
   * -  4 = the native frequency is not compatible with frequency settings.
   * -  5 = the vector transformation failed (resampling+highpassing)
   * -  6 = simulated signals could not be injected in the data vector
   * -  7 = the spectrum could not be computed
   * -  8 = the tiling could not be normalized
   * @param aInVectSize input vector size
   * @param aInVect input data vector (time domain)
   */
  int Condition(const int aInVectSize, double *aInVect);
  
  /**
   * Projects whitened data onto the tiles and fills output structures.
   * The data vector is subdivided into subsegments. Subsegments are Tukey-windowed, Fourier-transformed and normalized by the ASD. Finally, sub-segments are projected onto the tiling structure.
   *
   * In this function, the trigger structure is also filled with tiles above SNR threshold. If requested, the maps are written to disk.
   */
  bool Project(void);
  
  /**
   * Writes output to disk.
   * The output data products selected by the user in the option file are written to disk.
   *
   * Note that maps are written in the Project() function as they are built at the segment scale, not at the chunk scale.
   */
  bool WriteOutput(void);
  
  /**
   * Prints a progress report of the processing.
   */
  void PrintStatusInfo(void);
  

  /**
   * Returns the segments associated to the trigger time coverage.
   * See Triggers::GetTriggerSegments().
   *
   * NOTE: This function should be called somewhere after Project() and before WriteOutput() while the triggers are present in memory.
   * @param aThr threshold object
   * @param aInfValue value above which the threshold is considered infinite
   */
  Segments* GetTriggerSegments(TH1D *aThr=NULL, const double aInfValue=1e20);
  
  /**
   * Returns list of channels.
   */
  inline vector <string> GetChannels(void){return fChannels;};

  /**
   * Returns chunk duration [s].
   */
  inline int GetChunkDuration(void){return dataseq->GetChunkDuration();};

  /**
   * Returns segment/block duration [s].
   */
  inline int GetSegmentDuration(void){return dataseq->GetSegmentDuration();};

  /**
   * Returns overlap duration [s].
   */
  inline int GetOverlapDuration(void){return dataseq->GetOverlapDuration();};

  /**
   * Returns working sampling frequency [Hz]
   */
  inline int GetSampleFrequency(void){return triggers[0]->GetWorkingFrequency();};
 
  /**
   * Prints a formatted message with a timer.
   * @param aMessage message to print
   */
  void PrintMessage(const string aMessage);

  /**
   * Returns class status.
   */
  inline bool GetStatus(void){ return status_OK; };

  /**
   * Returns verbosity level.
   */
  inline int GetVerbosity(void){ return fVerbosity; };

  /**
   * Returns Omicron version.
   */
  inline string GetVersion(void){ return "v2r2"; };

 private:

  // PROCESS STATUS
  bool status_OK;               ///< general status
  time_t timer;                 ///< timer
  time_t timer_start;           ///< timer start
  struct tm * ptm;              ///< gmt time
  int chanindex;                ///< current channel index
  double timeoffset;            ///< current time offset

  // OPTIONS
  void ReadOptions(void);       ///< to parse option card
  string fOptionFile;           ///< option file name
  int fVerbosity;               ///< verbosity level
  string fMaindir;              ///< main output directory
  string fOutFormat;            ///< output format string
  string fOutProducts;          ///< output product string
  vector <string> fChannels;    ///< list of channel names
  vector <int> fWindows;        ///< plot time windows
  string fClusterAlgo;          ///< clustering algorithm
  int fTileDown;                ///< tile-down flag
  vector <string> fInjChan;     ///< injection channel names
  vector <double> fInjFact;     ///< injection factors

  // PROCESS MONITORING
  Segments *inSegments;         ///< requested segments
  Segments **outSegments;       ///< segments currently processed
  int chunk_ctr;                ///< number of loaded chunks
  int *chan_ctr;                ///< number of times a channel was loaded
  int *chan_data_ctr;           ///< number of times a channel was found data
  int *chan_cond_ctr;           ///< number of times a channel was conditioned
  int *chan_proj_ctr;           ///< number of times a channel was projected
  int *chan_write_ctr;          ///< number of times a channel's triggers were saved
  double *chan_mapsnrmax;       ///< channel SNR max in maps (only for html)
  vector <int> mapcenter;       ///< save map centers (only for html)
  vector <int> chunkstart;      ///< chunk start for html plots (only for html)
  vector <int> chunkstop;       ///< chunk stop for html plots (only for html)

  // COMPONENTS
  GwollumPlot *GPlot;           ///< Gwollum plots
  Odata *dataseq;               ///< data sequence
  Spectrum *spectrum;           ///< spectrum structure
  ffl *FFL;                     ///< ffl
  Otile *tile;                  ///< tiling structure
  MakeTriggers **triggers;      ///< output triggers
  InjEct **inject;              ///< software injections

  // DATA VECTORS
  double *ChunkVect;            ///< chunk raw data (time domain)
  double *CondChunkVect;        ///< chunk conditioned data (time domain)
  double *SegVect;              ///< subsegment raw data (time domain)
  
  // CONDITIONING
  bool Whiten(double **aDataRe, 
	      double **aDataIm);///< whiten data vector
  double* GetTukeyWindow(const int aSize, const int aFractionSize); ///< create tukey window
  double *TukeyWindow;          ///< tukey window
  fft *offt;                    ///< FFT plan to FFT the input data
  double **dataRe;              ///< conditioned data (Re)
  double **dataIm;              ///< conditioned data (Im)

  // OUTPUT
  string maindir;               ///< output main directory
  vector <string> outdir;       ///< output directories per channel
  void SaveAPSD(const string type="PSD", const bool aCond=false);///< Save current PSD/ASD
  void SaveTS(const bool aCond=false); ///< Save current chunk time series
  void MakeHtml(void);          ///< make html report

  // MISC
  void PrintASCIIlogo(void);    ///< print ascii logo
  static const string colorcode[17];
  string GetColorCode(const double aSNRratio);


  ClassDef(Omicron,0)  
};

#endif


