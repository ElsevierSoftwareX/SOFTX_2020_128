//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Odata.h"

ClassImp(Odata)

////////////////////////////////////////////////////////////////////////////////////
Odata::Odata(Segments *aSegments, const int aChunkDuration, 
	     const int aSegmentDuration, const int aOverlapDuration,
	     const int aVerbosity){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // save parameters
  fSegments          = aSegments;
  ChunkDuration      = fabs(aChunkDuration);
  SegmentDuration    = fabs(aSegmentDuration);
  OverlapDuration    = fabs(aOverlapDuration);
  fVerbosity         = aVerbosity;
  status_OK=true;

  //***************** CHECKS *****************
  if(!aSegments->GetLiveTime()){
    cerr<<"Odata::Odata: no live time in input segments"<<endl;
    status_OK*=false;
  }
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

  // init timing
  seg=0; // sets on first segment
  ChunkStart = (int)(fSegments->GetStart(seg));
  ChunkStop  = -1;// flag for the first chunk
  NSegments  = (ChunkDuration-OverlapDuration)/(SegmentDuration-OverlapDuration);
}

////////////////////////////////////////////////////////////////////////////////////
Odata::~Odata(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>2) cout<<"Odata::~Odata"<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
bool Odata::NewChunk(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Odata::NewChunk: the Odata object is corrupted"<<endl;
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

/*
////////////////////////////////////////////////////////////////////////////////////
bool Odata::GetConditionedData(const int aNseg, double **aDataRe, double **aDataIm, double **aPSD, int &aPSDsize){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Odata::GetConditionedData: ("<<fChannelName<<") the Odata object is corrupted"<<endl;
    return false;
  }
  if(aNseg<0||aNseg>=NSegments){
    cerr<<"Odata::GetConditionedData: ("<<fChannelName<<") segment "<<aNseg<<" cannot be found in the data chunk"<<endl;
    return false;
  }

  // fill data vector (time-domain) and apply tukey window
  double *datain_tmp = new double [SegmentSize];
  for(int i=0; i<SegmentSize; i++)
    datain_tmp[i] = chunkvect[aNseg*(SegmentSize-OverlapSize)+i] * TukeyWindow[i];

  // fft-forward
  if(!offt->Forward(datain_tmp)){
    cerr<<"Odata::GetConditionedData: ("<<fChannelName<<") FFT-forward failed"<<endl;
    delete datain_tmp;
    return false;
  }
    
  // get ffted data
  *aDataRe=offt->GetRe();
  *aDataIm=offt->GetIm();
  delete datain_tmp;// no longer needed
  if(*aDataRe==NULL||*aDataIm==NULL){
    cerr<<"Odata::GetConditionedData: ("<<fChannelName<<") cannot retrieve FFT data"<<endl;
    return false;
  }

  // zero-out below high-frequency cutoff
  int icutoff = (int)floor(fHighPassFrequency/(double)fSamplingFrequency*(double)SegmentSize);
  for(int i=0; i<icutoff; i++){
    (*aDataRe)[i]=0.0;
    (*aDataIm)[i]=0.0;
  }

  //interpolate PSD vector
  gsl_interp_accel_reset(acc);
  gsl_spline_init(interp_psd, f_PSD, PSD, fPsdSize);

  // normalize data by the PSD
  double new_psd;
  for(int i=icutoff; i<SegmentSize/2; i++){
    if(f_data[i]<=f_PSD[fPsdSize-1]) // interpolate
      new_psd=gsl_spline_eval(interp_psd, f_data[i], acc);
    else // extrapolate
      new_psd=PSD[fPsdSize-1]+(PSD[fPsdSize-1]-PSD[fPsdSize-2])/(f_PSD[fPsdSize-1]-f_PSD[fPsdSize-2])*(f_data[i]-f_PSD[fPsdSize-1]);
        
    (*aDataRe)[i] /= sqrt(new_psd);
    (*aDataIm)[i] /= sqrt(new_psd);
    // NOTE: no FFT normalization here because we do a FFTback in Oqplane later with no normalization either.
  }

  // point to current PSD
  *aPSD=PSD;
  aPSDsize=fPsdSize;

  return true;
}


////////////////////////////////////////////////////////////////////////////////////
double* Odata::GetTukeyWindow(const int aSize, const int aFractionSize){
////////////////////////////////////////////////////////////////////////////////////

  int FracSize_2=aFractionSize/2;
  double *Window = new double [aSize];
  
  double factor = TMath::Pi()/(double)FracSize_2;

  for (int i=0; i<FracSize_2; i++)
    Window[i] = 0.5*(1+TMath::Cos(factor*(double)(i-FracSize_2)));
  for (int i=FracSize_2; i<aSize-FracSize_2; i++)
    Window[i] = 1.0;
  for (int i=aSize-FracSize_2; i<aSize; i++)
    Window[i] = 0.5*(1+TMath::Cos(factor*(double)(i-aSize+FracSize_2)));
 
  return Window;
}
*/

