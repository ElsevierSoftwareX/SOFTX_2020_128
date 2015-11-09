//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#ifndef __Odata__
#define __Odata__

#include "TMath.h"
#include "CUtils.h"
#include "Segments.h"
#include "FFT.h"

using namespace std;

/**
 * Read a segment list sequentially (for Omicron).
 * An input Segments object is divided into chunks to perform a structured analysis of data segments. The input Segments object (set with SetSegments()) is divided into overlapping chunks. The chunks are loaded sequentially any time the NewChunk() function is called. The chunk sequence can be represented in the following way:
 * \verbatim
------------------------------------------------------------ input segment
 |------------------| chunk i-1
                |------------------| chunk i
                               |------------------| chunk i+1
 
                |---| overlap
 \endverbatim
 *
 * In general, the input Segments object contain multiple time segments. The sequence described above does not necessarily match the size of the input segments. The Odata class is designed to deal with such edge effects. Firstly, segments shorter than a chunk duration are skipped. When calling NewChunk() for the last chunk of a segment, the overlap duration is adjusted to fit the leftover:
 * \verbatim
 -----------------------------------------|   <--- input segment under processing

    |--------------------------|              <--- penultimate chunk 
  
 ###### call NextChunk() to cover the left-over

               |--------------------------|   <--- last chunk
	       |---------------|              <--- adjusted overlap
 * \endverbatim 
 * Obviously, the user must be careful about this special case as the overlap duration is modified (the chunk duration is never changed). Some functions are available to monitor the overlap size.
 *
 * When moving to a new segment, the overlap duration is set back to nominal values.
 * \author    Florent Robinet
 */
class Odata{

 public:

  /**
   * @name Constructors and destructors
   @{
  */
  /**
   * Constructor of the Odata class.
   * Nominal durations are defined.
   * @param aChunkDuration nominal chunk duration [s]
   * @param aOverlapDuration nominal overlap duration [s]
   * @param averbose verbosity level
   */
  Odata(const int aChunkDuration, 
	const int aOverlapDuration,
	const int aVerbosity=0);

  /**
   * Destructor of the Odata class.
   */
  virtual ~Odata(void);
  /**
     @}
  */
  
  /**
   * Sets new segments to read.
   * The input Segments structure is loaded and chunks can be called with the NewChunk function.
   * This function returns false if the input segments are not valid.
   * @param aSegments input Segments
   */
  bool SetSegments(Segments *aSegments);

  /**
   * Loads a new (next) chunk.
   * The chunks are loaded following the definition presented in the description of this class. This function should be called iteratively to cover the full data set defined with SetSegments(). The returned value indicates the status of this operation:
   * - true : a new chunk has been loaded
   * - false : no more chunk to load
   */
  bool NewChunk(void);

  /**
   * Returns the GPS starting time of current chunk.
   * Returns -1 if no chunk has been loaded.
   */
  inline int GetChunkTimeStart(void){ return ChunkStart; };
  
  /**
   * Returns the chunk duration.
   */
  inline int GetChunkDuration(void){ return ChunkDuration; };

  /**
   * Returns the current overlap duration.
   * In most cases the overlap duration is nominal unless the special case of the end of an input segment is hit.
   */
  inline int GetCurrentOverlapDuration(void){ return OverlapDurationCurrent; };

  /**
   * Returns the nominal overlap duration.
   */
  inline int GetOverlapDuration(void){ return OverlapDuration; };

 private:
  
  int fVerbosity;             ///< verbosity level

  Segments *fSegments;        ///< input segments
  int ChunkDuration;          ///< nominal chunk duration
  int OverlapDuration;        ///< nominal overlap duration
  int OverlapDurationCurrent; ///< current overlap duration
  int ChunkStart;             ///< current chunk start
  int seg;                    ///< current segment index
  
  ClassDef(Odata,0)  
};

#endif


