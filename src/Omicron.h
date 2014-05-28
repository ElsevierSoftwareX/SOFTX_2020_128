//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Omicron__
#define __Omicron__

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
   * This constructors initializes all the components to run Omicron: data structures, data streams, tiling, triggers, injections, etc.
   *
   * An option file is required to define all the parameters to run Omicron. For more details about Omicron configuration, see <a href="../../Friends/omicron.html">this page</a>.
   *
   * Omicron can be used as a low-latency search and data are provided sequentially when they are available. To activate this mode, the user should remove the DATA FFL/LCF parameters from the option file.
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
   * Runs the analysis of data segments.
   * This function runs the Omicron algorithm over the data defined by the input segments.
   * @param aSeg Segments to process.
   */
  bool Process(Segments *aSeg);

  /**
   * Conditions a data vector.
   * This function is typically used for an online analysis. An online process provides a data vector to be processed by Omicron:
   * @param aInVect input data vector (time domain)
   *
   * Before projecting the data onto the tiles, the data are conditioned with this funtion.
   *
   * This vector MUST have the duration of a chunk as declared in the option file. NO check will be performed against that!
   *
   * The user must provide information about the data vector passed in arguments:
   * @param aChNumber channel number as declared in the option file (indexing starts at 0)
   * @param aInVectSize number of samples in the input vector
   *
   * If the returned integer value is negative, it means that a fatal error occured and the Omicron object is corrupted. If it is positive, it means that the processing ran into some errors but the Omicron object is still valid. If it is 0, the processing ended correctly.
   */
  int ConditionVector(const int aChNumber, const int aInVectSize, double *aInVect);
  
  /**
   * Writes triggers on disk.
   * When called, this function writes the triggers in the current Triggers structure on disk. Typically, this function should be called when processing data with ProcessVector().
   * @param aChNumber channel number as declared in the option file (indexing starts at 0)
   */
  bool WriteTriggers(const int aChNumber);
  
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
  //bool online;                  ///< online running if true

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
  string fFflFormat;            ///< FFL format
  int fSampleFrequency;         ///< sampling frequency of input data
  vector <double> fFreqRange;   ///< frequency range
  vector <double> fQRange;      ///< Q range
  int fChunkDuration;           ///< segment duration
  int fSegmentDuration;         ///< segment duration
  int fOverlapDuration;         ///< overlap duration
  double fMismatchMax;          ///< maximum mismatch
  double fSNRThreshold;         ///< SNR Threshold
  int fNtriggerMax;             ///< trigger limit
  string fClusterAlgo;          ///< clustering mode
  double fcldt;                 ///< clustering dt
  int fVerbosity;               ///< verbosity level
  string fOutFormat;            ///< output format
  bool writepsd;                ///< writing PSD flag
  bool writetimeseries;         ///< writing time series flag

  // MONITORING
  Segments *inSegments;         ///< requested segments
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
  double **dataRe;              ///< conditioned data (Re)
  double **dataIm;              ///< conditioned data (Im)

  // MISC
  void SavePSD(const int c, const int s, const int e);///< Save PSD in a ROOT file
  bool first_PSD;               ///< flags the PSD to write
  void SaveData(const int c, double *aData, const int s, const int e);///< Save time series in a ROOT file
  bool first_Data;              ///< flags the Data to write

  ClassDef(Omicron,0)  
};

#endif


