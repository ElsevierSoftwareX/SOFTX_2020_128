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
#include "Triggers.h"
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
 * \author    Florent Robinet
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
   * The Omicron class offers a wrapping functions (see Process() and Scan()) where every analysis steps are included. Omicron can also be used piece-by-piece for tailored applications (like a low-latency searches where data are provided sequentially when they are available). In this case, the option file can be minimal.
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
   * - InitSegments(): to load the segments to process
   * - MakeDirectories(): to create the output directory structure
   * - NewChunk(): to walk through the input segments
   * - ConditionVector(): to load and condition the input data vector
   * - MakeTriggers(): to project data onto the tiles
   * - WriteTriggers(): to save triggers on disk
   *
   * This function is only available if a FFL structure has been previously declared.
   * @param aSeg Segments to process.
   */
  bool Process(Segments *aSeg);

  /**
   * Runs the full scan analysis of a GPS time.
   * This function runs the Omicron algorithm over the data defined by a central time. The data are conditioned, projected on the tiling structure and resulting maps are saved on disk. This function calls the following sequence of Omicron functions:
   * - InitSegments(): to load the chunk to process
   * - MakeDirectories(): to create the output directory structure
   * - NewChunk(): to load a chunk centered on aTimeCenter
   * - ConditionVector(): to load and condition the input data vector
   * - MakeMaps(): to project data onto the tiles and fill the maps
   * - WriteMaps(): to save maps on disk
   *
   * This function is only available if a FFL structure has been previously declared.
   * @param aTimeCenter central time of the maps
   */
  bool Scan(const double aTimeCenter);

  /**
   * Initialize the segments to process.
   * This function should be called before any type of processing.
   * 
   * WARNING: the input segment object is not copied, only the pointer is saved. This means that the Segments structure pointed by aSeg should not be modified or deleted before the end of the processing.
   * @param *aSeg pointer to the input Segments structure
   */
  bool InitSegments(Segments *aSeg);

  /**
   * Create the output directory structure.
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
   * Load a new chunk.
   * The chunks are loaded following the structure defined in the option file and the Segments object defined with InitSegments(). When there is not enough data to fill one chunk (end of a segment), the chunk duration is shortened. This function should be called iteratively to cover the full data set. The segmentation is detailed in Odata. The returned value indicates the status of this operation:
   * - true : a new chunk has been loaded
   * - false : no more chunk to load
   */
  bool NewChunk(void);

  /**
   * Conditions a data vector.
   * @param aInVect input data vector (time domain)
   *
   * Before projecting the data onto the tiles, the data are conditioned with this funtion. The input data chunk is first re-sampled. Then the data is used to estimate the noise (PSD). Finally, the data segments are normalized by the ASD.
   *
   * IMPORTANT: The input vector MUST MATCH the current chunk segment loaded with NewChunk(). NO check will be performed against that!
   *
   * The user must provide information about the data vector given in argument:
   * @param aChNumber channel number as previously declared (indexing starts at 0)
   * @param aInVectSize number of samples in the input vector
   *
   * If the returned integer value is negative, it means that a fatal error occured and the Omicron object got corrupted. If it is positive, it means that the conditioning ran into some errors but the Omicron object is still valid. If it is 0, the conditioning ended correctly.
   */
  int ConditionVector(const int aChNumber, const int aInVectSize, double *aInVect);
  
  /**
   * Projects conditioned data onto the tiles and fills the Triggers structure.
   * It returns the current number of triggers in memory. -1 is returned if this function fails.
   * @param aChNumber channel number as previously declared (indexing starts at 0)
   */
  int MakeTriggers(const int aChNumber);
  
  /**
   * Writes triggers to disk.
   * After being saved, triggers are flushed out of memory. Triggers cannot be saved if the maximum number of trigger limit has been reached.
   * @param aChNumber channel number as previously declared (indexing starts at 0)
   */
  bool WriteTriggers(const int aChNumber);
  
  /**
   * Projects conditioned data onto the tiles and fills the maps structure.
   * This function only projects the first data segment of the current chunk onto the tiles. Maps are saved in memory. A time offset is applied to set the origin to aTimeCenter. The following maps are made:
   * - one map for each Q-plane spanning the first segment of the chunk.
   * - one full map combining all the Q-planes, centered on aTimeCenter and spanning each window value.
   * 
   * @param aChNumber channel number as previously declared (indexing starts at 0)
   * @param aTimeCenter GPS time where to set the time origin
   */
  bool MakeMaps(const int aChNumber, const double aTimeCenter);
  
  /**
   * Writes current maps to disk.
   * The maps built with MakeMaps() are saved to disk with the formats defined in Omicron(). When graphical formats are selected, an additional file is saved on disk: '[channel]_mapsummary.root'. This file contains a TTree summarizing the properties of the maps.
   *
   * In addition to maps, some complementary plots are also saved (graphical formats only):
   * - SNR vs. time for the frequency where the SNR is maximal
   * - SNR vs. frequency for the time where the SNR is maximal
   */
  bool WriteMaps(void);

  /**
   * Create a html report for a scan.
   * @param aScanDir path to scan directory
   */
  bool ScanReport(const string aScanDir="");


  //Segments* GetOnlineSegments(const int aChNumber, TH1D *aThr, const double aPadding=0.0, const double aInfValue=1e20);
  
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
   * Prints a progress report of the processing.
   */
  void PrintStatusInfo(void);
  
  /**
   * Returns class status.
   */
  inline bool GetStatus(void){ return status_OK; };

 private:

  // STATUS
  bool status_OK;               ///< general status
  time_t timer;                 ///< timer
  time_t timer_start;           ///< start time
  struct tm * ptm;              ///< gmt time

  // INPUT OPTIONS
  bool ReadOptions(void);       ///< to parse option card
  string fOptionFile;           ///< option file name
  vector <string> fOptionName;  ///< option name (metadata)
  vector <string> fOptionType;  ///< option type (metadata)
  string fMaindir;              ///< main output directory
  vector <string> fOutdir;      ///< output directories per channel
  vector <string> fChannels;    ///< list of channel names
  vector <string> fInjChan;     ///< injection channel names
  vector <double> fInjFact;     ///< injection factors
  string fFflFile;              ///< path to FFL file
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
  bool writepsd;                ///< writing PSD flag
  bool writetimeseries;         ///< writing time series flag
  string fWriteMode;            ///< write mode

  // PROCESS MONITORING
  Segments *inSegments;         ///< cumulative requested segments
  Segments **outSegments;       ///< segments currently written on disk
  int *chunk_ctr;               ///< number of chunks
  int *cor_chunk_ctr;           ///< number of corrupted chunks
  int *cor_data_ctr;            ///< number of corrupted data chunks
  int *max_chunk_ctr;           ///< number of maxed-out chunks

  // COMPONENTS
  Odata *dataseq;               ///< data sequence
  Sample **sample;              ///< sampling structures
  Streams **streams;            ///< streams
  Spectrum *spectrum;           ///< spectrum structure
  ffl *FFL;                     ///< ffl
  Otile *tile;                  ///< tiling structure
  Triggers **triggers;          ///< output triggers

  // DATA
  int ChunkSize;                ///< chunk sample size (varies)
  int OverlapSize;              ///< overlap sample size
  int SegmentSize;              ///< segment sample size
  double *ChunkVect;            ///< chunk data container (time domain)
  double *SegVect;              ///< segment data container (time domain)
  
  // CONDITIONING
  bool Condition(double **aDataRe, double **aDataIm); ///< condition data vector
  double* GetTukeyWindow(const int aSize, const int aFractionSize); ///< create tukey window
  double *TukeyWindow;          ///< tukey window
  fft *offt;                    ///< FFT plan to FFT the input data
  double **dataRe;              ///< conditioned data container (Re)
  double **dataIm;              ///< conditioned data container (Im)

  // SCANS
  int MapChNumber;              ///< channel currently mapped
  TH2D **Qmap;                  ///< set of Q-maps
  TH2D **Qmap_full;             ///< combined Q-maps
  double Qmap_center;           ///< current GPS time of Q maps
  int *loudest_qmap;            ///< Q-map conatining the loudest tile
  GwollumPlot *GPlot;           ///< Gwollum plots
  string fScandir;              ///< latest scan directory

  // MISC
  void SaveAPSD(const int c, const string type="PSD");    ///< Save current PSD
  void SaveTS(const int c, double tcenter=0);///< Save current chunk time series
  bool *first_save;             ///< flags the first save
  void PrintASCIIlogo(void);    ///< print ascii logo

  ClassDef(Omicron,0)  
};

#endif


