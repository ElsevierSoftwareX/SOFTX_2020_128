//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Omicron__
#define __Omicron__

#include "IO.h"
#include "Otile.h"
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
   * This constructor initializes all the components to run Omicron: data structures, data streams, tiling, maps, triggers, injections, monitoring, etc.
   * An option file is required to define all the parameters to run Omicron. For more details about Omicron configuration, see <a href="../../Friends/omicron.html">this page</a>.
   *
   * After initialization, the Omicron methods should be called sequentially to perform the analysis. Here is a typical sequence:
   * - InitSegments() defines the data segments to process.
   * - MakeDirectories() creates a specific directory structure for the output 
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
   * This function should always be called before any type of processing. Use NewChunk() to sequence the Omicron analysis.
   * @param aSeg pointer to the input Segments structure
   */
  bool InitSegments(Segments *aSeg);

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
   * Calls a new time chunk.
   * The time chunks are called following the time sequence defined by the Otile class. The returned value indicates the status of this operation:
   * - true : a new time chunk has been successfully called
   * - false : no more chunk to load
   */
  bool NewChunk(void);

  /**
   * Calls a new channel.
   * The channels defined in the option file are called incrementally. If this function is called after the last channel, false is returned and the channel sequence is reset: the next call will call the first channel again.
   */
  bool NewChannel(void);

  /**
   * Loads a data vector.
   * The data vector of the current channel and the current chunk is loaded. If requested in the option file, the injection stream and the software injections are added to the data vector. This function loads the data from the frames listed in the FFL. The FFL option is therefore mandatory to use this function.
   * It is the user's responsibility to delete the returned data vector.
   *
   * If this function fails, a pointer to NULL is returned.
   * @param aDataVector pointer to the data vector
   * @param aSize sample size of the data vector
   */
  bool LoadData(double **aDataVector, int *aSize);

  /**
   * Conditions a data vector.
   * Before projecting the data onto the tiles, the data is conditioned with this function. The input data chunk is first resampled, highpassed and Tukey-windowed. Then the data is used to estimate the noise (PSD).
   *
   * IMPORTANT: The input vector size MUST MATCH the chunk size loaded with NewChunk(). NO check is performed against that!
   *
   * If the returned value is negative, it means that a fatal error occured and the Omicron object got corrupted. If it is positive, it means that the conditioning failed but the Omicron object is still valid for further use. If it is 0, the conditioning ended correctly. The error code is the following:
   * - -1 = the Omicron object is corrupted.
   * -  0 = OK
   * -  1 = the input vector is NULL
   * -  2 = the input vector size is 0
   * -  3 = the input vector appears to be flat
   * -  4 = the vector transformation failed (resampling+highpassing)
   * -  5 = the spectrum could not be updated
   * -  6 = the tiling power could not be computed
   * @param aInVectSize input vector size
   * @param aInVect input data vector (time domain)
   */
  int Condition(const int aInVectSize, double *aInVect);
  
  /**
   * Projects whitened data onto the tiles and fills output structures.
   * The data vector Fourier-transformed and normalized by the ASD. The data are then projected onto the tiling structure.
   *
   * In this function, the trigger structure is also filled with tiles above SNR threshold.
   */
  bool Project(void);
  
  /**
   * Writes output products to disk.
   * The output data products selected by the user in the option file and for the current chunk/channel are written to disk.
   */
  bool WriteOutput(void);

  /**
   * Writes output triggers to disk.
   * The triggers accumulated so far are written to disk (for the current channel). The trigger object in then flushed.
   */
  bool WriteTriggers(void);

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
  inline int GetChunkDuration(void){return tile->GetTimeRange();};

  /**
   * Returns overlap duration [s].
   */
  inline int GetOverlapDuration(void){return tile->GetOverlapDuration();};

  /**
   * Returns working sampling frequency [Hz]
   */
  inline int GetSampleFrequency(void){return triggers[0]->GetWorkingFrequency();};
 
  /**
   * Resets PSD buffer.
   * (for the current channel)
   */
  inline void ResetPSDBuffer(void) { spectrum[chanindex]->Reset(); };

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

  // OPTIONS
  void ReadOptions(void);       ///< to parse option card
  string fOptionFile;           ///< option file name
  int fVerbosity;               ///< verbosity level
  string fMaindir;              ///< main output directory
  string fOutFormat;            ///< output format string
  string fOutProducts;          ///< output product string
  bool fNoLogo;                 ///< no logo flag
  vector <string> fChannels;    ///< list of channel names
  vector <int> fWindows;        ///< plot time windows. FIXME: to move in Otile
  string fClusterAlgo;          ///< clustering algorithm
  int fTileDown;                ///< tile-down flag FIXME: to move in Otile
  vector <string> fInjChan;     ///< injection channel names
  vector <double> fInjFact;     ///< injection factors

  // PROCESS MONITORING
  Segments *inSegments;         ///< requested segments
  Segments **outSegments;       ///< segments currently processed
  int chunk_ctr;                ///< number of called chunks
  int *chan_ctr;                ///< number of called chunks /channel
  int *chan_data_ctr;           ///< number of LoadData() calls /channel
  int *chan_cond_ctr;           ///< number of Condition() calls /channel
  int *chan_proj_ctr;           ///< number of Project() calls /channel
  int *chan_write_ctr;          ///< number of WriteOutput() calls /channel
  double *chan_mapsnrmax;       ///< channel SNR max in maps (only for html)
  vector <int> chunkcenter;     ///< save chunk centers (only for html)
  
  // COMPONENTS
  GwollumPlot *GPlot;           ///< Gwollum plots
  Spectrum **spectrum;          ///< spectrum structure
  ffl *FFL;                     ///< ffl
  Otile *tile;                  ///< tiling structure
  MakeTriggers **triggers;      ///< output triggers
  InjEct **inject;              ///< software injections

  // DATA VECTORS
  double *ChunkVect;            ///< chunk raw data (time domain)
  double *WhiteChunkVect;       ///< chunk for whitened data (time domain)
  double *DataRe;               ///< whitened data vector (Re)
  double *DataIm;               ///< whitened data vector (Im)
   
  // CONDITIONING & WHITENING
  bool Whiten(void);            ///< whiten data vector
  double* GetTukeyWindow(const int aSize, const int aFractionSize); ///< create tukey window
  double *TukeyWindow;          ///< tukey window
  fft *offt;                    ///< FFT plan to FFT the input data

  // OUTPUT
  string maindir;               ///< output main directory
  vector <string> outdir;       ///< output directories per channel
  void SaveAPSD(const string aType);///< Save current PSD/ASD
  void SaveTS(const bool aWhite=false); ///< Save current chunk time series
  void SaveSpectral(void);      ///< Save current spectral plots
  void MakeHtml(void);          ///< make html report

  // MISC
  void PrintASCIIlogo(void);    ///< print ascii logo
  static const string colorcode[17];
  string GetColorCode(const double aSNRratio);


  ClassDef(Omicron,0)  
};

#endif


