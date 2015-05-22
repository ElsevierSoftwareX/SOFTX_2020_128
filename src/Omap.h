//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Omap__
#define __Omap__

#include "CUtils.h"
#include "TMath.h"
#include "TH2D.h"

using namespace std;

/**
 * Create a time-frequency map.
 * This class was designed to create and use a multi-resolution time-frequency map.
 *
 * \author    Florent Robinet
 */
class Omap: public TH2D {

 public:
  friend class Otile;   ///< Friendly class
  friend class Oqplane; ///< Friendly class

 private:
  
  // constructor - destructor
  Omap();
  virtual ~Omap(void);

  void SetBins(const double aQ, const double aFrequencyMin, const double aFrequencyMax,
	       const int aTimeRange, const double aMismatchStep);

  void SetBins(const int aNf, const double aFrequencyMin, const double aFrequencyMax,
	       const int aNt, const int aTimeRange);


  inline double GetTimeRange(void){
    return GetXaxis()->GetBinUpEdge(GetNbinsX())-GetXaxis()->GetBinLowEdge(1); 
  };
  inline double GetTimeMin(void){ 
    return GetXaxis()->GetXmin(); 
  };
  inline double GetTimeMax(void){ 
    return GetXaxis()->GetXmax(); 
  };
  inline double GetFrequencyMin(void){ 
    return GetYaxis()->GetXmin(); 
  };
  inline double GetFrequencyMax(void){ 
    return GetYaxis()->GetXmax(); 
  };
  inline int GetNBands(void){ 
    return GetNbinsY();
  };
  inline int GetBandIndex(const double aFrequency){
    return GetYaxis()->FindBin(aFrequency)-1;
  };
  inline double GetBandFrequency(const int aBandIndex){ 
    return GetYaxis()->GetBinCenterLog(aBandIndex+1);
  };
  inline double GetBandStart(const int aBandIndex){ 
    return GetYaxis()->GetBinLowEdge(aBandIndex+1);
  };
  inline double GetBandEnd(const int aBandIndex){ 
    return GetYaxis()->GetBinUpEdge(aBandIndex+1);
  };
  inline double GetBandWidth(const int aBandIndex){ 
    return GetYaxis()->GetBinWidth(aBandIndex+1);
  };
  inline double GetTileDuration(const int aBandIndex){
    return GetXaxis()->GetBinWidth(1)*bandMultiple[aBandIndex];
  };
  inline int GetBandNtiles(const int aBandIndex){ 
    return GetNbinsX()/bandMultiple[aBandIndex];
  };
  inline double GetTileContent(const int aTimeTileIndex, const int aBandIndex){    
    return GetBinContent(aTimeTileIndex*bandMultiple[aBandIndex]+1,aBandIndex+1);
  };
  inline double GetTilePhase(const int aTimeTileIndex, const int aBandIndex){    
    return phase[aBandIndex][aTimeTileIndex];
  };
  inline double GetTileTag(const int aTimeTileIndex, const int aBandIndex){    
    return GetBinError(aTimeTileIndex*bandMultiple[aBandIndex]+1,aBandIndex+1);
  };
  inline double GetTileTimeStart(const int aTimeTileIndex, const int aBandIndex){
    return GetXaxis()->GetBinLowEdge(aTimeTileIndex*bandMultiple[aBandIndex]+1);
  };
  inline double GetTileTimeEnd(const int aTimeTileIndex, const int aBandIndex){
    return GetXaxis()->GetBinUpEdge((aTimeTileIndex+1)*bandMultiple[aBandIndex]);
  };
  inline double GetTileTime(const int aTimeTileIndex, const int aBandIndex){
    return GetXaxis()->GetBinLowEdge(aTimeTileIndex*bandMultiple[aBandIndex]+bandMultiple[aBandIndex]/2+1);
  };
  inline int GetTimeTileIndex(const int aBandIndex, const double aTime){
    return (int)floor((aTime-GetTimeMin())/GetTileDuration(aBandIndex));
  };

  // SETS
  void SetTileContent(const int aTimeTileIndex, const int aBandIndex, const double aContent, const double aPhase=-100.0);
  inline void SetTileTag(const int aTimeTileIndex, const int aBandIndex, const double aTag){
    SetBinError(aTimeTileIndex*bandMultiple[aBandIndex]+1,aBandIndex+1,aTag);
  };

  void SetTileDisplay(void);
 
  int Ntiles;                       ///< number of tiles in the plane
  int *bandMultiple;                ///< band multiple (time resolution)
  double **phase;                   ///< tile phase array

  ClassDef(Omap,0)  
};


#endif


