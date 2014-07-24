//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Oqplane__
#define __Oqplane__

#include "CUtils.h"
#include "Triggers.h"
#include "Spectrum.h"
#include "FFT.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"


#define NFROWMAX 1000


using namespace std;

/**
 * Create a frequency row in a time-frequency Q-plane.
 * This class was designed to create and use a frequency row defined by a central frequency and a Q value to build time-frequency Q-planes. This class is private and can only be used by the Oqplane class.
 *
 * \author    Florent Robinet
 */
class FreqRow {

 public:
  friend class Oqplane;  ///< friendly class
  
 private:
  
  FreqRow(const int aTimeRange, 
	  const int aTimePad, 
	  const int aSampleFrequency, 
	  const double aF, 
	  const double aQ, 
	  const double aMismatchStep, 
	  const double aSNRThreshold);
  virtual ~FreqRow(void);

  double* GetSNRs(double *aDataRe, double *aDataIm);// Get SNRs in FD
  bool GetTriggers(Triggers *aTriggers, double *aDataRe, double *aDataIm, const int aStartTime);
  inline void SetPower(const double aPower){ fPower=aPower; };
  
  // PARAMETERS
  double fF;                ///< central frequency [Hz]
  double fQ;                ///< Q of the plane
  int fTimeRange;           ///< duration of analysis [s]
  int fTimePad;             ///< time pad [s]
  int fSampleFrequency;     ///< sample frequency of input data [Hz]
  double fMismatchStep;     ///< maximum mismatch between neighboring tile
  double fSNRThreshold;     ///< SNR Threshold

  // DERIVED PARAMETERS
  double fBandWidth;        ///< f-row tile bandwidth
  double fDuration;         ///< f-row tile duration

  // F-ROW
  vector <double> fTime;    ///< vector of central times
  int fNumberOfTiles;       ///< number of time tiles
  bool *ValidIndices;       ///< energy indices to keep
  double fPower;            ///< power of the f-row

  // Q-TRANSFORM
  vector <double> fWindow;  ///< bi square window function
  vector <double> fWindowFrequency; ///< window frequencies
  vector <int> fDataIndices;///< vector of data indices to inverse fourier transform
  int fZeroPadLength;       ///< number of zeros to append to windowed data
  //double *  working_vector[2];
  fft *offt;
};


/**
 * Create a time-frequency Q-plane.
 * This class was designed to create and use a time-frequency Q-plane defined by a Q value. This class is private and can only be used by the Otile class.
 *
 * \author    Florent Robinet
 */
class Oqplane {

 public:
  friend class Otile;  ///< friendly class

 private:

  Oqplane(const double aQ, 
	  const int aSampleFrequency, 
	  const int aTimeRange, 
	  const int aTimePad, 
	  const double aFrequencyMin, 
	  const double aFrequencyMax, 
	  const double aMismatchStep, 
	  const double aSNRThreshold);
  virtual ~Oqplane(void);
    
  bool GetTriggers(Triggers *aTriggers, double *aDataRe, double *aDataIm, const int aTimeStart);
  bool SetPowerSpectrum(Spectrum *aSpec);
  TH2D* GetMap(double *aDataRe, double *aDataIm, const double time_offset=0.0, const bool printamplitude=false);
  
  // STATUS
  bool status_OK;                  ///< class status

  // PARAMETERS
  double fQ;
  int fTimeRange;                  ///< duration of analysis [s]
  int fTimePad;                    ///< time pad [s]
  int fSampleFrequency;            ///< sampling frequency [Hz]
  double fFrequencyMin,            ///< frequency min
    fFrequencyMax;                 ///< frequency max
  double fMismatchStep;            ///< maximum mismatch between neighboring tiles
  double fSNRThreshold;            ///< SNR Threshold

  // DERIVED PARAMETERS
  double fQPrime;                  ///< Q prime = Q / sqrt(11)

  // Q-PLANE
  TH2D *hplane;                    ///< map
  int fNumberOfRows;               ///< number of frequency rows
  double fPlaneNormalization;      ///< plane normalization
  vector <double> fFreq;           ///< vector of frequencies
  FreqRow *freqrow[NFROWMAX];      ///< f-row objects
  int fNumberOfTiles;              ///< number of tiles
  
  bool CheckParameters(void);      ///< check the validity of the parameters
  void GetPlaneNormalization(void);///< get plane normalization
  void BuildTiles(void);

  ClassDef(Oqplane,0)  
};


#endif


