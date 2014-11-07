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
 * An input Segments object is segmented into several levels to perform an organized analysis of data segments. The input Segments object (set with SetSegments()) is divided into overlapping chunks. The chunks are loaded sequentially any time the NewChunk() function is called. The chunk sequence can be represented in the following way:
 * \verbatim
------------------------------------------------------------ input segment
 |------------------| chunk i-1
                |------------------| chunk i
                               |------------------| chunk i+1
 
                |---| overlap
 \endverbatim
 *
 * Chunks are subdivided in subsegments overlaping the same way as chunks:
 * \verbatim
 |------------------------------------------| chunk i
 |----------| subsegment 0
         |----------| subsegment 1
                 |----------| subsegment 2
                         |----------| subsegment 3
                                 |----------| subsegment 4
				 
                                         |------------------------------------------| chunk i+1
 \endverbatim 
 *
 * In general, the input Segments object contain multiple time segments. The sequence described above does not necessarily match the size of the input segments. The Odata class is designed to deal with such edge effects. When calling NewChunk(), if the next chunk ends after the end of the current segment, some actions are taken to adjust the chunk size to the left-over (represented below):
 * \verbatim
 ------------------------------------| input segment under processing

 -------| chunk i-1
 
 ###### call NextChunk() to cover the left-over
 
     |--------------------------| chunk i
     |----------| subsegment 0
             |----------| subsegment 1
                     |----------| subsegment 2

 ###### call NextChunk() to cover the left-over
                              
                          |----------| chunk i+1
                          |----------| subsegment 0
 \endverbatim 
 *
 * - the chunk (i) size is reduced until it is fully contained inside the left-over segment. This reduction is perfomed so that the chunk size is a multiple of at least one sub-segment (in the example, there are 3 subsegments instead of the nominal 5).
 * - if the left-over is smaller than the size of one sub-segment, the chunk (i+1) size is set to the size of one subsegment and the chunk stop is set to the end of the current segment. The overlap duration value is modified.
 *
 * Obviously, the user must be careful about these special cases as the chunck and/or the overlap duration is modified (the subsegment duration is never changed). Some functions are available to monitor the different sizes.
 *
 * When moving to a new segment, the chunk and overlap durations are set back to nominal values.
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
   *
   * IMPORTANT: the nominal durations must verify some conditions. If these conditions are not met, the duration are automatically changed to match the Odata class requirements.
   * @param aChunkDuration nominal chunk duration [s]
   * @param aSegmentDuration nominal subsegment duration [s]
   * @param aOverlapDurationnominal overlap duration [s]
   * @param averbose verbosity level
   */
  Odata(const int aChunkDuration, 
	const int aSegmentDuration,
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
   * Sets the segments to read.
   * The input Segments structure is loaded and chunks can be called with the NewChunk function.
   * @param aSegments input Segments
   */
  bool SetSegments(Segments *aSegments);

  /**
   * Loads a new (next) chunk.
   * The chunks are loaded following the structure defined for the Odata class and the Segments object defined with SetSegments(). When there is not enough data to fill one chunk (end of a segment), the chunk duration is shortened. This function should be called iteratively to cover the full data set. The returned value indicates the status of this operation:
   * - true : a new chunk has been loaded
   * - false : no more chunk to load
   */
  bool NewChunk(void);

  /**
   * Returns number of subsegments in the current chunk.
   * The nominal value is returned if no chunk is loaded.
   */
  inline int GetNSegments(void){ return NSegments; };
  
  /**
   * Returns the GPS starting time of current chunk.
   */
  inline int GetChunkTimeStart(void){ return ChunkStart; };
  
  /**
   * Returns the GPS ending time of current chunk.
   */
  inline int GetChunkTimeEnd(void){ return ChunkStop; };
  
  /**
   * Returns the GPS starting time of segment aNseg.
   * Returns -1 if this function fails.
   * @param aNseg segment number
   */
  int GetSegmentTimeStart(const int aNseg);

  /**
   * Returns the GPS ending time of segment aNseg.
   * Returns -1 if this function fails.
   * @param aNseg segment number
   */
  int GetSegmentTimeEnd(const int aNseg);
  
  /**
   * Returns the current chunk duration.
   * In most cases the chunk duration is nominal unless the special case of the end of an input segment is hit.
   */
  inline int GetCurrentChunkDuration(void){ return ChunkStop-ChunkStart; };

  /**
   * Returns the nominal chunk duration.
   */
  inline int GetChunkDuration(void){ return ChunkDuration; };

  /**
   * Returns the segment duration.
   * This value is never changed.
   */
  inline int GetSegmentDuration(void){ return SegmentDuration; };

  /**
   * Returns the current overlap duration.
   * In most cases the overlap duration is nominal unless the special case of the end of an input segment is hit.
   */
  inline int GetCurrentOverlapDuration(void){ return OverlapDurationCurrent; };

  /**
   * Returns the nominal overlap duration.
   * In most cases the overlap duration is nominal unless the special case of the end of an input segment is hit.
   */
  inline int GetOverlapDuration(void){ return OverlapDuration; };

 private:
  
  bool status_OK;        ///< class status
  int fVerbosity;        ///< verbosity level

  Segments *fSegments;   ///< input segments - DO NOT DELETE
  int ChunkDuration;     ///< nominal chunk duration
  int SegmentDuration;   ///< nominal segment duration
  int OverlapDuration;   ///< nominal overlap duration
  int OverlapDurationCurrent;///< current overlap duration
  int ChunkStart;        ///< current chunk start
  int ChunkStop;         ///< current chunk end
  int NSegments;         ///< number of overlapping subsegments in the current chunk
  int seg;               ///< current segment index
  bool TestChunk(void);  ///< test chunk against segments

  ClassDef(Odata,0)  
};

#endif


