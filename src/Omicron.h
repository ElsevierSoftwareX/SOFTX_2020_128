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
   * This constructors initializes all the components to run Omicron: data structures, tiling, triggers, injections and network of detectors.
   *
   * A Segments object is required to define the data segments to process. An option file is also required to defined all the parameters to run Omicron. For more details about Omicron configuration, see <a href="../../Friends/omicron.html">this page</a>.
   *
   * Omicron can be used as a low-latency search and data are provided sequentially when they are available. The FFL (or LCF) option must be set to "ONLINE". For this online mode, no Segments input is necessary: a pointer to NULL should used.
   * @param aOptionFile option file
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
   * Sets Segments to process.
   * @param aSeg Segments to process
   * @param aChNumber channel number to process
   */
  bool SetSegments(Segments *aSeg, const int aChNumber);

  /**
   * Runs the analysis of data segments.
   */
  bool Process(Segments *aSeg);

  //int ProcessOnline(const int aChNumber, FrVect *aVect);
  //bool WriteOnline(const int aChNumber);
  
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
 
  //int GetNativeSampleFrequency(const int aChNumber);

  /**
   * Prints a progress report of the processing.
   */
  void PrintStatusInfo(void);

 private:

  // STATUS
  bool status_OK;               ///< general status
  bool online;                  ///< online running if true

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
  //vector <int> fNativeFrequency;///< native sampling frequency
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
  Streams **streams;            ///< streams
  Spectrum *spectrum;           ///< spectrum structure
  ffl *FFL;                     ///< ffl
  Otile *tile;                  ///< tiling structure
  Triggers **triggers;          ///< output triggers

  // DATA
   
  ClassDef(Omicron,0)  
};

#endif


