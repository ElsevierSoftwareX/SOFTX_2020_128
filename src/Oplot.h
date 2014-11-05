//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Oplot__
#define __Oplot__

#include "TriggerPlot.h"
#include "EventMap.h"
#include "TGraph.h"

using namespace std;

/**
 * Plot Omicron triggers.
 * This class was specifically designed to plot Omicron triggers. It inherits from the TriggerPlot class. This class is actually just a wrapper of TriggerPlot as it defines 4 trigger collections differing by a SNR threshold. These collection are then assigned a given plotting style (colors, markers...).
 * \author Florent Robinet
 */
class Oplot: public TriggerPlot {

 public:

  /**
   * @name Constructors and destructors
   @{
  */
  /**
   * Constructor of the Oplot class.
   * The 4 TriggerPlot collections are defined over clusters (not triggers!) with 4 different SNR thresholds and 4 different plotting styles.
   * @param aPattern input file pattern
   * @param aDirectory trigger ROOT directory
   * @param aVerbose verbosity level
   */
  Oplot(const string aPattern, const string aDirectory="", const int aVerbose=0);

  /**
   * Destructor of the Oplot class.
   */
  virtual ~Oplot(void);
  /**
     @}
  */
  
  void SetTimeRange(const int aTimeMin, const int aTimeMax);

  void PrintLoudestEventMap(const string aFileName="");

 private:
  double snrthr[4];
  EventMap *Eloud;

  ClassDef(Oplot,0)  
};

#endif


