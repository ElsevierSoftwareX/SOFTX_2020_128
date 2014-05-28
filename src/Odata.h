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
 * A Segments object is segmented with several levels to perform an organized analysis of data segments. The input segments are divided in chunks. The chunks can be represented in the following way:
 * \verbatim
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
 * The chunk duration is pre-defined. However, if the input segment is not long enough, the chunk size is adjusted to fit the maximum number of subsegments. The segment size cannot be ajusted.
 *
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
   *
   * @param aChunkDuration chunk duration [s]
   * @param aSegmentDuration subsegment duration [s]
   * @param aOverlapDuration overlap duration [s]
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
   * The chunks are loaded following the structure defined for the Odata class and the Segments object defined in Odata(). When there is not enough data to fill one chunk (end of a segment), the chunk duration is shortened. This function should be called iteratively to cover the full data set. The returned value indicates the status of this operation:
   * - true : a new chunk has been loaded
   * - false : no more chunk to load
   */
  bool NewChunk(void);

  /**
   * Returns number of subsegments per chunk.
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
   * Returns whitened data (frequency domain).
   * The data segment number 'aNseg' (after being downsampled and highpassed) is FFTed and normalized using the PSD. After this transformation, the real and imaginary parts are returned as vectors. The user is responsible to delete these vectors after use.
   *
   * A pointer to the PSD vector is returned as well as its size.
   * @param aNseg segment number
   * @param aDataRe returned vector for the real part
   * @param aDataIm returned vector for the imaginary part
   * @param aPSD pointer to the PSD vector
   * @param aPSDsize size of the PSD vector
   */
  //bool GetConditionedData(const int aNseg, double **aDataRe, double **aDataIm, double **aPSD, int &aPSDsize);


  /**
   * Returns current class status.
   */
  inline bool GetStatus(void){ return status_OK; };

  
 private:
  
  bool status_OK;              ///< class status
  int fVerbosity;              ///< verbosity level

  Segments *fSegments;         ///< input segments - DO NOT DELETE
  int ChunkDuration;           ///< chunk duration
  int SegmentDuration;         ///< segment duration
  int OverlapDuration;         ///< overlap duration
  int ChunkStart;              ///< current chunk start
  int ChunkStop;               ///< current chunk end
  int NSegments;               ///< number of 50% overlapping segments in one chunk
  int seg;                     ///< current segment index
  bool TestChunk(void);        ///< test chunk against segments

  ClassDef(Odata,0)  
};

#endif


