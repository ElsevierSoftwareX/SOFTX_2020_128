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

using namespace std;


/**
 * Process data with the Omicron algorithm.
 * This class was designed...
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
   * This constructor initializes all the components to run Omicron: data structures, data streams, tiling, triggers, injections, monitoring, etc.
   * An option file is required to define all the parameters to run Omicron. For more details about Omicron configuration, see <a href="../../Friends/omicron.html">this page</a>.
   *
   * The Omicron class offers a wrapping function (see Process()) where every analysis steps are included. Omicron can also be used piece-by-piece for tailored applications (like a low-latency search where data are provided sequentially when they are available). In this case, the option file can be minimal (or even empty).
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
   * Sets a new timing structures for the analysis segmentation.
   */
  //bool SetTiming(const int aChunkDuration, const int aSegmentDuration, const int aOverlapDuration);


  /**
   * Runs the full trigger analysis of data segments.
   * This function runs the Omicron algorithm over the data defined by the input segments. The data are segmented, conditioned, projected on the tiling structure and resulting triggers are saved on disk.
   *
   * This function is only available if a FFL structure has been previously declared.
   * @param aSeg Segments to process.
   */
  bool Process(Segments *aSeg);

  /**
   * Runs the full acan analysis of a GPS time.
   * This function runs the Omicron algorithm over the data defined by a central time.
   *
   * This function is only available if a FFL structure has been previously declared.
   * @param aTimeCenter
   */
  bool Scan(const double aTimeCenter);

  /**
   * Conditions a data vector.
   * @param aInVect input data vector (time domain)
   *
   * Before projecting the data onto the tiles, the data are conditioned with this funtion.
   *
   * This vector MUST have the duration of a chunk as previously declared. NO check will be performed against that!
   *
   * The user must provide information about the data vector given in argument:
   * @param aChNumber channel number as previously declared (indexing starts at 0)
   * @param aInVectSize number of samples in the input vector
   *
   * If the returned integer value is negative, it means that a fatal error occured and the Omicron object got corrupted. If it is positive, it means that the processing ran into some errors but the Omicron object is still valid. If it is 0, the conditioning ended correctly.
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
   * After being saved, triggers are flushed out of memory.
   * @param aChNumber channel number as previously declared (indexing starts at 0)
   */
  bool WriteTriggers(const int aChNumber);
  
  /**
   * Projects conditioned data onto the tiles and fills the maps structure.
   * This function only projects the first data segment of the current chunk onto the tiles. Maps are saved in memory. See WriteMaps() to save the maps on disk.
   * By default, the map
   * @param aChNumber channel number as previously declared (indexing starts at 0)
   * @param aTimeCenter GPS time where to center the map
   */
  bool MakeMaps(const int aChNumber, const double aTimeCenter);
  
  /**
   * Writes mapss to disk.
   * 
   * @param aChNumber channel number as previously declared (indexing starts at 0)
   */
  bool WriteMaps(const int aChNumber);

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
   * Returns working sampling frequency
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
  int fChunkDuration;           ///< chunk duration
  int fSegmentDuration;         ///< segment duration
  int fOverlapDuration;         ///< overlap duration
  double fMismatchMax;          ///< maximum mismatch
  vector <int> fWindows;        ///< scan windows
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
  Segments **outSegments;       ///< segments currently processed
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
  int ChunkSize;                ///< chunk sample size (not fixed)
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
  TH2D **Qmap;                  ///< set of Q-maps
  double Qmap_center;           ///< current GPS time of Q maps
  GwollumPlot *GPlot;           ///< Gwollum plots

  // MISC
  void SaveAPSD(const int c, const string type="PSD");    ///< Save current PSD in a ROOT file
  void SaveData(const int c, double *aData, const int s, const int e);///< Save time series in a ROOT file
  bool *first_save;             ///< flags the first save

  ClassDef(Omicron,0)  
};

#endif


