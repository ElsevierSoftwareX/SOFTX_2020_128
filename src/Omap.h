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

  void SetBins(const double aQ,
	       const double aFrequencyMin, const double aFrequencyMax,
	       const int aTimeRange, const double aMismatchStep);

  inline double GetTimeRange(void){
    return GetXaxis()->GetBinUpEdge(GetNbinsX())-GetXaxis()->GetBinLowEdge(1); 
  };
  inline double GetTimeMin(void){ 
    return GetXaxis()->GetBinLowEdge(1); 
  };
  inline double GetTimeMax(void){ 
    return GetXaxis()->GetBinUpEdge(GetNbinsX()); 
  };
  inline double GetFrequencyMin(void){ 
    return GetYaxis()->GetBinLowEdge(1); 
  };
  inline double GetFrequencyMax(void){ 
    return GetYaxis()->GetBinUpEdge(GetNbinsY()); 
  };
  inline int GetNBands(void){ 
    return GetNbinsY();
  };
  inline long int GetNTiles(void){ 
    return Ntiles;
  };
  inline int GetBandIndex(const double aFrequency){
    return GetYaxis()->FindBin(aFrequency)-1;
  };
  inline double GetBandFrequency(const int aBandIndex){
    return bandCenter[aBandIndex];
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

  inline double GetTileTimeStart(const int aTimeTileIndex, const int aBandIndex){
    return GetXaxis()->GetBinLowEdge(aTimeTileIndex*bandMultiple[aBandIndex]+1);
  };
  inline double GetTileTimeEnd(const int aTimeTileIndex, const int aBandIndex){
    return GetXaxis()->GetBinUpEdge((aTimeTileIndex+1)*bandMultiple[aBandIndex]);
  };
  inline double GetTileTime(const int aTimeTileIndex, const int aBandIndex){
    return (GetTileTimeStart(aTimeTileIndex,aBandIndex) + GetTileTimeEnd(aTimeTileIndex,aBandIndex)) / 2.0;
  };
  inline int GetTimeTileIndex(const int aBandIndex, const double aTime){
    return (int)floor((aTime-GetTimeMin())/GetTileDuration(aBandIndex));
  };

  inline double GetTileContent(const int aTimeTileIndex, const int aBandIndex){
    return TH2::GetBinContent(aTimeTileIndex*bandMultiple[aBandIndex]+1,aBandIndex+1);
  };
  inline void SetTileContent(const int aTimeTileIndex, const int aBandIndex, const double aContent){
    int tstart=aTimeTileIndex*bandMultiple[aBandIndex]+1;
    int tend=tstart+bandMultiple[aBandIndex];
    for(int t=tstart; t<tend; t++) TH2::SetBinContent(t,aBandIndex+1,aContent);
  };

  long int Ntiles;                  ///< number of tiles in the plane
  double *bandCenter;               ///< frequency bin center
  int *bandMultiple;                ///< band multiple (time resolution)

  ClassDef(Omap,0)  
};


#endif


