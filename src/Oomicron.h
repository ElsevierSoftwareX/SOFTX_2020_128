//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Oomicron__
#define __Oomicron__

#include <IO.h>
#include "Otile.h"
#include "Oinject.h"
#include "Oconfig.h"
#include <Date.h>
#include <InjEct.h>
#include <ffl.h>

using namespace std;


/**
 * Process data with the Omicron algorithm.
 * - An introduction to Omicron is available: <a href="../../Friends/omicron.html">Omicron introduction</a>
 * - <a href="https://tds.ego-gw.it/ql/?c=10651">Omicron technical note</a>
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
   * This constructor initializes all the components to run Omicron: data structures, data streams, spectra, tiling, maps, triggers, injections, monitoring, etc.
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
   * - WriteOutput() writes output data products to disk (in loop #1/2)
   *
   * @param aOptionFile path to the option file
   * @param aGpsRef Reference time to initiate structures (default).
   * @param aStrict strict mode when set to true: the status of the Omicron object is set to false if options are incorrectly provided.
   */
  Omicron(const string aOptionFile, const int aGpsRef=-1, const bool aStrict=false);
  
  /**
   * Destructor of the Omicron class.
   */
  virtual ~Omicron(void);
  /**
     @}
  */
  

  /**
   * Initializes the segments to process and to output.
   * This function should always be called before any type of processing. Use NewChunk() to sequence the Omicron analysis.
   *
   * Optionally, output segments (for triggers only!) can be specified. If so, triggers outside the output segments will not be saved. Use a pointer to NULL to not use this option.
   * @param aInSeg pointer to the input Segments structure
   * @param aOutSeg pointer to the output Segments structure
   */
  bool InitSegments(Segments *aInSeg, Segments *aOutSeg=NULL);

  /**
   * Creates a specific directory tree for the output.
   * Two directory structures are possible:
   * - [path_to_outdir]/aId/[channel_name] if aId is not 0
   * - [path_to_outdir]/[channel_name] if aId is 0
   *
   * where [path_to_outdir] is the output directory specified by the user in the option file and [channel_name] is the channel name being processed.
   * 
   * The aId value is rounded to the third digit.
   *
   * If this function is never called, all the output is dumped in the output directory specified by the user in the option file.
   * @param aId directory id
   */
  bool MakeDirectories(const double aId = 0.0);

  /**
   * Calls a new time chunk.
   * The time chunks are called following the time sequence defined by the Otile class. The returned value indicates the status of this operation:
   * - true : a new time chunk has been successfully called
   * - false : no more chunk to load
   */
  bool NewChunk(void);

  /**
   * Defines a new time chunk.
   * Instead of defining a list of input segments (see InitSegments()) and processing sequentially the data (see NewChunk()), it is possible to define any time chunk.
   *
   * Optionally, it is possible to reset the PSD buffer (for all channels).
   *
   * The chunk duration must match the one defined in the option file.
   * @param aTimeStart GPS start time of the chunk
   * @param aTimeEnd GPS end time of the chunk
   * @param aResetPSDBuffer flag to reset the PSD buffers
   */
  bool DefineNewChunk(const int aTimeStart, const int aTimeEnd, const bool aResetPSDBuffer=false);

  /**
   * Calls a new channel.
   * The channels defined in the option file are called incrementally. If this function is called after the last channel, false is returned and the channel sequence is reset: the next call will call the first channel again.
   */
  bool NewChannel(void);

  /**
   * Returns the name of the current channel.
   * Returns "none" if no channel is defined.
   */
  inline string GetChannelName(void){ if(chanindex>=0) return triggers[chanindex]->GetName(); return "none"; };

  /**
   * Loads a data vector.
   * The data vector of the current channel and the current chunk is loaded. If requested in the option file, the injection stream and the software injections are added to the data vector. This function loads the data from the frames listed in the FFL. The FFL option is therefore mandatory to use this function.
   * It is the user's responsibility to delete the returned data vector.
   *
   * If this function fails, a pointer to NULL is returned.
   * @param aDataVector pointer to the returned data vector
   * @param aSize sample size of the return ed data vector
   */
  bool LoadData(double **aDataVector, int *aSize);

  /**
   * Conditions a data vector.
   * Before projecting the data onto the tiles, the data is conditioned and whitened with this function. The input data chunk is first removed its DC component, then resampled, optionally highpassed, Tukey-windowed, Fourier-transformed and whitened twice. In this process, the conditioned data vector is used to update the estimate of the noise power density (PSD).
   *
   * IMPORTANT: The input vector size MUST MATCH the chunk size loaded with NewChunk(). NO check is performed against that!
   *
   * If the returned value is negative, it means that a fatal error occured and the Omicron object got corrupted. If it is positive, it means that the conditioning failed but the Omicron object is still valid for further use. If it is 0, the conditioning ended correctly. The error code is the following:
   * - -1 = the Omicron object is corrupted.
   * -  0 = OK
   * -  1 = the input vector is pointing to NULL
   * -  2 = the input vector size is 0
   * -  3 = the native sampling frequency cannot be updated
   * -  4 = the vector transformation failed (DC removal+resampling+highpassing)
   * -  5 = the spectrum (1) is not available
   * -  6 = the chunk data could not be FFTed (forward)
   * -  7 = the spectrum (2) is not available
   * -  8 = the tiling power could not be computed
   * @param aInVectSize input vector size
   * @param aInVect input data vector (time domain)
   */
  int Condition(const int aInVectSize, double *aInVect);
  
  /**
   * Projects whitened data onto the tiles and fills output structures.
   * The data are projected onto the tiling structure. The number of tiles above threshold is returned. -1 is returned if this function fails.
   */
  int Project(void);
  
  /**
   * Extract triggers above threshold.
   * Triggers are saved one Q plane at a time. After each Q plane, the number of triggers is checked. If it exceeds a given maximal rate, defined in the user parameter file, false is returned and the next Q planes are not searched for (see also Otile::SaveTriggers()). All triggers are reset.
   */
  bool ExtractTriggers(void);

  /**
   * Write triggers to disk.
   * All triggers collected until now are saved to disk.
   * This function returns the trigger file path. "none" is returned if this function fails.
   *
   * Optionally, the trigger files can be saved in a directory structure defined by a LIGO-virgo convention: [main dir]/[IFO prefix]/[channel]_OMICRON/
   * @param aLVDirConvention set to true to apply the LIGO trigger path convention
   */
  string WriteTriggers(const bool aLVDirConvention = false);

  /**
   * Writes output products to disk.
   * The output data products selected by the user in the option file and for the current chunk/channel are written to disk.
   */
  bool WriteOutput(void);

  /**
   * Prints a progress report of the processing.
   */
  void PrintStatusInfo(void);
  
  /**
   * When time plots are requested in output, a time offset [s] can be added.
   * @param aTimeOffset time offset [s]
   */
  inline void SetPlotTimeOffset(const double aTimeOffset=0.0){ toffset=aTimeOffset; };
  
  /**
   * Returns the segments associated to a selection of triggers.
   * See Triggers::GetTriggerSegments().
   *
   * NOTE: This function should be called somewhere after Project() and before WriteOutput() while the triggers are present in memory.
   * @param aThr threshold object
   * @param aInfValue value above which the threshold is considered to be infinite
   */
  Segments* GetTriggerSegments(TH1D *aThr=NULL, const double aInfValue=1e20);
  
  /**
   * Returns the segments associated to a selection of triggers.
   * See Triggers::GetLastTriggerSegments().
   *
   * NOTE: This function should be called somewhere after Project() and before WriteOutput() while the triggers are present in memory.
   * @param aThr threshold object
   * @param aInfValue value above which the threshold is considered to be infinite
   */
  Segments* GetLastTriggerSegments(TH1D *aThr=NULL, const double aInfValue=1e20);
  
  /**
   * Returns the number of channels.
   */
  inline int GetNChannels(void){ return nchannels; };

  /**
   * Returns list of channels.
   */
  vector <string> GetChannels(void);

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
   * Returns the current number of triggers in memory for the current channel.
   */
  inline int GetNTriggers(void){ return triggers[chanindex]->GetNtriggers(); };

  /**
   * Returns the SNR threshold used for the triggers
   */
  inline int GetSNRThreshold(void){return tile->GetSNRTriggerThr(); };
  
  /**
   * Returns the buffer livetime for the current channel.
   */
  inline double GetTriggerBufferLiveTime(void){ return triggers[chanindex]->GetBufferLiveTime(); };

  /**
   * Reset the trigger buffer.
   */
  inline void ResetTriggerBuffer(void){ return triggers[chanindex]->ResetBuffer(); };

  /**
   * Reset the trigger buffer after a given time.
   * The time value must be aligned with the chunk starts: only chunks can be purged.
   * @param aTime time after which triggers are removed

   */
  inline void ResetTriggerBufferAfter(const double aTime){ return triggers[chanindex]->ResetAfter(aTime); };

  /**
   * Resets PSD buffer.
   * (for the current channel)
   */
  inline void ResetPSDBuffer(void) { spectrum1[chanindex]->Reset(); };

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

 private:

  // PROCESS STATUS
  bool status_OK;               ///< general status
  time_t timer;                 ///< timer
  time_t timer_start;           ///< timer start
  struct tm * ptm;              ///< gmt time
  int chanindex;                ///< current channel index

  // OPTIONS
  void ReadOptions(const int aGpsRef=-1, const bool aStrict=false);
  string fOptionFile;           ///< option file name
  int fVerbosity;               ///< verbosity level
  string fMaindir;              ///< main output directory
  string fOutFormat;            ///< output format string
  string fOutProducts;          ///< output product string
  bool fNoLogo;                 ///< no logo flag
  double toffset;               ///< time offset for plots
  vector <int> fWindows;        ///< plot time windows. FIXME: to move in Otile
  string fClusterAlgo;          ///< clustering algorithm
  string fftplan;               ///< fft plan
  double fratemax;              ///< maximum trigger rate /chunk
  vector <string> fInjChan;     ///< injection channel names
  vector <double> fInjFact;     ///< injection factors
  int fsginj;                   ///< perform SG injections
   
  // PROCESS MONITORING
  Segments *inSegments;         ///< requested segments
  Segments **outSegments;       ///< segments currently processed
  int chunk_ctr;                ///< number of called chunks
  int *chan_ctr;                ///< number of called chunks /channel
  int *chan_data_ctr;           ///< number of LoadData() calls /channel
  int *chan_cond_ctr;           ///< number of Condition() calls /channel
  int *chan_proj_ctr;           ///< number of Project() calls /channel
  int *chan_write_ctr;          ///< number of WriteOutput() calls /channel
  int *trig_ctr;                ///< number of tiles above snr thr /channel
  double *chan_mapsnrmax;       ///< channel SNR max in maps (only for html)
  vector <int> chunkcenter;     ///< save chunk centers (only for html)
  vector <string> chunktfile;   ///< save chunk file (only for html)
  ofstream oinjfile;            ///< file with sg injection parameters
  
  // COMPONENTS
  int nchannels;                ///< number of channels
  vector <string> outdir;       ///< output directories / channel
  Spectrum **spectrum1;         ///< 1st spectrum structure / channel
  Spectrum **spectrum2;         ///< 2nd spectrum structure / channel
  Spectrum *spectrumw;          ///< spectrum structure to test whitening
  ffl *FFL;                     ///< ffl
  ffl *FFL_inject;              ///< ffl for injection signals
  fft *offt;                    ///< FFT plan to FFT the input data
  Otile *tile;                  ///< tiling structure
  TriggerBuffer **triggers;     ///< output triggers / channel
  Oinject *oinj;                ///< software sg injections
  InjEct **inject;              ///< software injections (in frame data) / channel

  // DATA VECTORS
  double *ChunkVect;            ///< chunk raw data (time domain)
    
  // CONDITIONING & WHITENING
  void Whiten(Spectrum *aSpec); ///< whiten data vector
  double* GetTukeyWindow(const int aSize, const int aFractionSize); ///< create tukey window
  double *TukeyWindow;          ///< tukey window
 
  // OUTPUT
  string maindir;               ///< output main directory
  void SaveAPSD(const string aType);///< Save current PSD/ASD
  void SaveWPSD(void);          ///< Save whitened PSD
  void SaveTS(const bool aWhite=false); ///< Save current chunk time series
  void SaveSG(void);            ///< Save current sg injection parameters
  void SaveSpectral(void);      ///< Save current spectral plots
  void MakeHtml(void);          ///< make html report

  // MISC
  void PrintASCIIlogo(void);    ///< print ascii logo
  static const string colorcode[17];
  string GetColorCode(const double aSNRratio);
  bool IsFlat(const int aInVectSize, double *aInVect);

  ClassDef(Omicron,0)  
};

#endif


