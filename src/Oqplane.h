//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Oqplane__
#define __Oqplane__

#include "MakeTriggers.h"
#include "Spectrum.h"
#include "FFT.h"
#include "TH2D.h"

using namespace std;

/**
 * Create a frequency row in a time-frequency Q-plane.
 * This class was designed to create and use a frequency row defined by a central frequency and a Q value to build time-frequency Q-planes. This class is private and can only be used by the Oqplane class.
 *
 * \author    Florent Robinet
 */
/*
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
  bool GetTriggers(MakeTriggers *aTriggers, double *aDataRe, double *aDataIm, const int aStartTime, const int aExtraTimePadMin=0);
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
  fft *offt;                ///< row fft
};

*/
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
	  const double aFrequencyMin, 
	  const double aFrequencyMax, 
	  const double aMismatchStep);
 
  virtual ~Oqplane(void);

  bool ProjectData(double *aDataRe, double *aDataIm);

  bool SetTileContent(const int aTimeTileIndex, const int aFrequencyTileIndex, const double aContent);
  bool SetTileContent(const double aTime, const double aFrequency, const double aContent);
  void PrintParameters(void); ///< print plane parameters
  double UpdateThreshold(const int aBandIndex, double *aEnergies, double &aThreshold);

  inline void PresentTile(void){
    for(int f=0; f<qplane->GetNbinsY(); f++)
      for(int t=0; t<qplane->GetNbinsX(); t++)
	qplane->SetBinContent(t+1,f+1,(t/bandMultiple[f])%2);
    return;
  }

  inline int GetNBands(void) { return qplane->GetNbinsY(); };
  inline int GetBandNtiles(const int aBandIndex) { return qplane->GetNbinsX()/bandMultiple[aBandIndex]; };
  inline double GetBandFrequency(const int aBandIndex) { return qplane->GetYaxis()->GetBinCenter(aBandIndex+1); };
  inline double GetTimeRange(void) { return qplane->GetXaxis()->GetXmax()-qplane->GetXaxis()->GetXmin(); };

    
  //bool GetTriggers(MakeTriggers *aTriggers, double *aDataRe, double *aDataIm, const int aTimeStart, const int aExtraTimePadMin=0);
  bool SetPower(Spectrum *aSpec);
    
  // PARAMETERS
  double Q;                         ///< Q value
  double QPrime;                    ///< Q prime = Q / sqrt(11)
  double SampleFrequency;           ///< sampling frequency
  double MismatchStep;              ///< maximum mismatch between neighboring tiles

  // Q-PLANE
  TH2D *qplane;                     ///< Q plane
  int Ntiles;                       ///< number of tiles in the plane

  void GetPlaneNormalization(void); ///< get plane normalization
  double PlaneNormalization;        ///< plane normalization
  
  // FREQUENCY BANDS
  int *bandMultiple;                ///< band multiple (time resolution)
  int *bandWindowSize;              ///< band Gaussian window size
  double **bandWindow;              ///< band Gaussian window
  double **bandWindowFreq;          ///< band Gaussian window frequency
  double *bandPower;                ///< band power
  fft **bandFFT;                    ///< band fft
  
  ClassDef(Oqplane,0)  
};


#endif


