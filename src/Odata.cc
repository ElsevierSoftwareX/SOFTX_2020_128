//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Odata.h"

ClassImp(Odata)

////////////////////////////////////////////////////////////////////////////////////
Odata::Odata(const int aChunkDuration, const int aSegmentDuration,
	     const int aOverlapDuration, const int aVerbosity){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // save parameters
  ChunkDuration      = (int)fabs(aChunkDuration);
  SegmentDuration    = (int)fabs(aSegmentDuration);
  OverlapDuration    = (int)fabs(aOverlapDuration);
  fVerbosity         = aVerbosity;
  status_OK=true;

  //***************** CHECKS *****************
  if(OverlapDuration%2){
    cerr<<"Odata::Odata: the overlap duration is not even"<<endl;
    status_OK*=false;
  }
  if(SegmentDuration>ChunkDuration){
    cerr<<"Odata::Odata: the segment duration cannot be larger than the chunk duration"<<endl;
    status_OK*=false;
  }
  if(OverlapDuration>=SegmentDuration){
    cerr<<"Odata::Odata: the overlap duration cannot be larger than the segment duration"<<endl;
    status_OK*=false;
  }
  if((ChunkDuration-OverlapDuration)%(SegmentDuration-OverlapDuration)){
    cerr<<"Odata::Odata: the overlap/segment/chunk durations do not macth. You could use:"<<endl;
    cerr<<"              chunk duration   = "<<((ChunkDuration-OverlapDuration)/(SegmentDuration-OverlapDuration)+1)*(SegmentDuration-OverlapDuration)+OverlapDuration<<" s"<<endl;
    cerr<<"              segment duration = "<<SegmentDuration<<" s"<<endl;
    cerr<<"              overlap duration = "<<OverlapDuration<<" s"<<endl;
    status_OK*=false;
  }

  //******************************************

  // no timing yet
  fSegments  = NULL;
  seg        = -1;
  ChunkStart = -1;
  ChunkStop  = -1;
  NSegments  = (ChunkDuration-OverlapDuration)/(SegmentDuration-OverlapDuration);
}

////////////////////////////////////////////////////////////////////////////////////
Odata::~Odata(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"Odata::~Odata"<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
bool Odata::SetSegments(Segments *aSegments){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Odata::SetSegments: the Odata object is corrupted"<<endl;
    return false;
  }
  if(!aSegments->GetLiveTime()){
    cerr<<"Odata::SetSegments: no live time in input segments"<<endl;
    return false;
  }

  fSegments = aSegments;

  // init timing
  seg        = 0; // sets on first segment
  ChunkStart = (int)(fSegments->GetStart(seg));
  ChunkStop  = -1;// flag for the first chunk
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Odata::NewChunk(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Odata::NewChunk: the Odata object is corrupted"<<endl;
    return false;
  }
  if(ChunkStart<0){
    cerr<<"Odata::NewChunk: no segment to read"<<endl;
    return false;
  }
  
  // First chunk
  if(ChunkStop<0) ChunkStop=ChunkStart+ChunkDuration;
  
  // New chunk
  else{
    ChunkStart=ChunkStop-OverlapDuration;
    ChunkStop=ChunkStart+ChunkDuration;
  }
  
  // test chunk against segments and update chunk if necessary
  while(!TestChunk());
  
  // end of data segments --> stop
  if(seg>=fSegments->GetNsegments()){
    cout<<"Odata::NewChunk: end of data segments"<<endl;
    ChunkStart=-1;
    return false;
  }
  
  // new values
  NSegments=(ChunkStop-ChunkStart-OverlapDuration)/(SegmentDuration-OverlapDuration);

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Odata::TestChunk(void){
 
  // too short chunk -> move to next segment
  if(ChunkStop-ChunkStart<SegmentDuration){
    seg++;// next segment
    if(seg>=fSegments->GetNsegments()) return true;// end of data segments
    ChunkStart=(int)(fSegments->GetStart(seg));
    ChunkStop=ChunkStart+ChunkDuration;
    return false;// the validity of this case needs to be tested again
  }

  // ideal case: inside the current segment
  if(ChunkStart>=fSegments->GetStart(seg)&&ChunkStop<=fSegments->GetEnd(seg))
    return true;

  // the chunk stops after the segment end -> shorten chunk     
  if(ChunkStart<fSegments->GetEnd(seg)&&ChunkStop>fSegments->GetEnd(seg)){
    ChunkStop-=(SegmentDuration-OverlapDuration);
    return false;// the validity of this case needs to be tested again
  }

  // the chunk starts after the segment end
  if(ChunkStart>=fSegments->GetEnd(seg)){
    seg++;// move to next segment
    if(seg>=fSegments->GetNsegments()) return true;// end of data segments
    ChunkStart=(int)(fSegments->GetStart(seg));
    ChunkStop=ChunkStart+ChunkDuration;
    return false;// the validity of this case needs to be tested again
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
int Odata::GetSegmentTimeStart(const int aNseg){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Odata::GetSegmentTimeStart: the Odata object is corrupted"<<endl;
    return -1;
  }
  if(aNseg<0||aNseg>=NSegments){
    cerr<<"Odata::GetSegmentTimeStart: segment "<<aNseg<<" cannot be found in the data chunk"<<endl;
    return -1;
  }

  return ChunkStart+aNseg*(SegmentDuration-OverlapDuration);
}

////////////////////////////////////////////////////////////////////////////////////
int Odata::GetSegmentTimeEnd(const int aNseg){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Odata::GetSegmentTimeEnd: the Odata object is corrupted"<<endl;
    return -1;
  }
  if(aNseg<0||aNseg>=NSegments){
    cerr<<"Odata::GetSegmentTimeEnd: segment "<<aNseg<<" cannot be found in the data chunk"<<endl;
    return -1;
  }

  return ChunkStart+aNseg*(SegmentDuration-OverlapDuration)+SegmentDuration;
}


