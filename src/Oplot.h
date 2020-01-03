//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Oplot__
#define __Oplot__

#include <TriggerPlot.h>
#include <EventMap.h>
#include <TGraph.h>

using namespace std;

/**
 * Plot Omicron triggers.
 * This class was specifically designed to plot Omicron triggers. It inherits from the TriggerPlot class. This class is actually just a wrapper of TriggerPlot as it defines 4 trigger collections differing by a SNR threshold. These collections are then assigned a given plotting style (colors, markers...).
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
   * The 4 TriggerPlot collections are defined over TIME clusters (not triggers!) with 4 different SNR thresholds (atomatically set) and 4 different plotting styles.
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
  
  /**
   * Defines new SNR thrsholds.
   * 
   * @param SNR0 first SNR threshold
   * @param SNR1 second SNR threshold
   * @param SNR2 third SNR threshold
   * @param SNR3 fourth SNR threshold
   */ 
  void SetSNRThresholds(const double aSNR0, const double aSNR1, const double aSNR2, const double aSNR3);

  /**
   * Returns the SNR thrshold of given collection.
   * @param aCollIndex collection index
   */ 
  double GetSNRThreshold(const int aCollIndex);

  /**
   * Defines a new time range for the plots.
   * By default, the time range is given by the input trigger files.
   * @param aTimeMin starting time
   * @param aTimeMax stopping time
   */ 
  void SetTimeRange(const int aTimeMin, const int aTimeMax);

  /**
   * Prints the loudest event map.
   * After the plots have been made (see TriggerPlot::MakeCollections()), the loudest event of the plot can be mapped. The resulting map can be saved in a file if aFileName is used.
   * @param aFileName output file name
   */ 
  void PrintLoudestEventMap(const string aFileName="");

 private:

  double snrthr[5];  ///< snr thresholds
  EventMap *Eloud;   ///< loudest event map

  ClassDef(Oplot,0)  
};

#endif


