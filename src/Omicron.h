//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Omicron__
#define __Omicron__

#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <stdio.h>
#include "CUtils.h"
#include "IO.h"
#include "Segments.h"
#include "Triggers.h"
#include "TMath.h"
#include "Otile.h"
#include "Odata.h"

#define NDATASTREAMS 50

using namespace std;

//////////////////////////////////////////////////////////////////////////////
class Omicron {

 public:

  Omicron(Segments *aSegments, const string aOptionFile);// constructor
  virtual ~Omicron(void);// destructor

  // PROCESSES
  bool MakeTiling(void);
  bool Process();
  int ProcessOnline(const int aChNumber, FrVect *aVect);

  //INFO
  inline int GetChunkDuration(void){return fChunkDuration;};
  inline int GetSegmentDuration(void){return fSegmentDuration;};
  inline int GetOverlapDuration(void){return fOverlapDuration;};
  int GetNativeSampleFrequency(const int aChNumber);
  inline int GetSampleFrequency(void){return fSampleFrequency;};
  inline vector <string> GetChannelList(void){return fChannels;};

 protected:

  // STATUS
  bool status_OK;  // general status
  bool tiling_OK;  // tiling status
  bool online;     // online running

  // OUTPUT FLAGS
  bool writepsd;          // writing PSD flag
  bool writetimeseries;   // writing time series flag
  bool writewhiteneddata; // writing whiten data flag

  // OPTIONS
  int fVerbosity;
  string fOptionFile;
  IO *fOptions;
  string fFflFile;              // path to FFL file (Virgo)
  string fCacheFile;            // path to Cache file (LIGO)
  vector <string> fChannels;    // list of channels
  vector <int> fFreqRange;      // Frequency range
  vector <double> fQRange;      // Q range
  vector <int> fNativeFrequency;// native sampling frequency
  int fSampleFrequency;         // sample frequency of input data
  int fChunkDuration;           // segment duration
  int fSegmentDuration;         // segment duration
  int fOverlapDuration;         // overlap duration
  double fMismatchMax;          // maximum mismatch
  string fClusterAlgo;          // clustering mode
  string fOutdir[NDATASTREAMS]; // output directories
  string fOutFormat;            // output format

  // TRIGGER
  double fSNRThreshold; // SNR Threshold
  int fNtriggerMax;     // trigger limit
  double fcldt;         // clustering dt

  // TILING
  Otile *tile;

  // DATA
  Segments *fSegments;
  Odata *odata[NDATASTREAMS];
  double *c_data[2];// conditioned data complex vector

  //OUTPUT
  Triggers *triggers[NDATASTREAMS];// output triggers
  int fFileSize; // file size in seconds
  
  bool ReadOptions(void);
  

  ClassDef(Omicron,0)  
};

#endif


