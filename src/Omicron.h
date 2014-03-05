//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Omicron__
#define __Omicron__

#include "IO.h"
#include "Inject.h"
#include "Triggers.h"
#include "TMath.h"
#include "Otile.h"
#include "Odata.h"

#define NDATASTREAMS 50

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
   * @param aSegments Segments to process
   * @param aOptionFile option file
   */
  Omicron(Segments *aSegments, const string aOptionFile);

  /**
   * Destructor of the Omicron class.
   */
  virtual ~Omicron(void);
  /**
     @}
  */

  /**
   * Runs the analysis of data segments.
   */
  bool Process(void);

  int ProcessOnline(const int aChNumber, FrVect *aVect);
  bool WriteOnline(const int aChNumber);
  
  Segments* GetOnlineSegments(const int aChNumber, TH1D *aThr, const double aPadding=0.0, const double aInfValue=1e20);
  
  inline int GetChunkDuration(void){return fChunkDuration;};
  inline int GetSegmentDuration(void){return fSegmentDuration;};
  inline int GetOverlapDuration(void){return fOverlapDuration;};
  int GetNativeSampleFrequency(const int aChNumber);
  inline int GetSampleFrequency(void){return fSampleFrequency;};
  inline vector <string> GetChannelList(void){return fChannels;};
  bool PrintStatusInfo(void);

 private:

  // STATUS
  bool status_OK;               ///< general status
  bool online;                  ///< online running if true

  // INPUT OPTIONS
  bool ReadOptions(void);       ///< to parse option card
  string fOptionFile;           ///< option file name
  vector <string> fOptionName;  ///< option name (metadata)
  vector <string> fOptionType;  ///< option type (metadata)
  vector <string> fChannels;    ///< list of channel names
  vector <string> fInjChan;     ///< injection channel names
  vector <double> fInjFact;     ///< injection factors
  vector <string> fDetectors;   ///< detectors
  string fFflFile;              ///< path to FFL file (Virgo)
  string fLcfFile;              ///< path to LCF file (LIGO)
  vector <int> fNativeFrequency;///< native sampling frequency
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
  string fOutdir[NDATASTREAMS]; ///< output directories per channel
  string fOutFormat;            ///< output format
  bool writepsd;                ///< writing PSD flag
  bool writetimeseries;         ///< writing time series flag

  // MONITORING
  Segments **outSegments;       ///< segments currently processed
  int *chunk_ctr;               ///< number of chunks
  int *cor_chunk_ctr;           ///< number of corrupted chunks
  int *max_chunk_ctr;           ///< number of maxed-out chunks

  // NETWORK (optional)
  string fInjFile;              ///< injection file
  Network *Net;                 ///< network
  Inject *Inj;                  ///< software injections

  // TILING
  Otile *tile;                  ///< tiling structure

  // DATA
  bool LCF2FFL(const string lcf_file, const string ffl_file);
  Segments *fSegments;          ///< segments to process - DO NOT DELETE
  Odata *odata[NDATASTREAMS];   ///< data structures
  double *psd;                  ///< psd vector - DO NOT DELETE

  // OUTPUT
  Triggers *triggers[NDATASTREAMS];///< output triggers
  
  ClassDef(Omicron,0)  
};

#endif


