//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

ClassImp(Omicron)

////////////////////////////////////////////////////////////////////////////////////
Omicron::Omicron(const string aOptionFile){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // input
  fOptionFile=aOptionFile;

  // default options
  fOptionName.push_back("omicron_OUTPUT_DIRECTORY");          fOptionType.push_back("s");
  fMaindir=".";
  fOutdir.clear();
  fOptionName.push_back("omicron_DATA_CHANNEL");              fOptionType.push_back("s");
  fChannels.clear();
  fOptionName.push_back("omicron_INJECTION_CHANNEL");         fOptionType.push_back("s");
  fInjChan.clear();
  fOptionName.push_back("omicron_INJECTION_FACTOR");          fOptionType.push_back("d");
  fInjFact.clear();
  fOptionName.push_back("omicron_DATA_FFLFILE");              fOptionType.push_back("s");
  fFflFile="none";
  fOptionName.push_back("omicron_DATA_FFLFORMAT");            fOptionType.push_back("s");
  fFflFormat="ffl";
  fOptionName.push_back("omicron_DATA_SAMPLEFREQUENCY");      fOptionType.push_back("i");
  fSampleFrequency=4096;
  fOptionName.push_back("omicron_PARAMETER_FMIN");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_FMAX");            fOptionType.push_back("d");
  fFreqRange.clear(); fFreqRange.push_back(32);  fFreqRange.push_back(2048); 
  fOptionName.push_back("omicron_PARAMETER_QMIN");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_QMAX");            fOptionType.push_back("d");
  fQRange.clear(); fQRange.push_back(4);  fQRange.push_back(100); 
  fOptionName.push_back("omicron_PARAMETER_CHUNKDURATION");   fOptionType.push_back("i");
  fChunkDuration=484;
  fOptionName.push_back("omicron_PARAMETER_BLOCKDURATION");   fOptionType.push_back("i");
  fSegmentDuration=64;
  fOptionName.push_back("omicron_PARAMETER_OVERLAPDURATION"); fOptionType.push_back("i");
  fOverlapDuration=8;
  fOptionName.push_back("omicron_PARAMETER_MISMATCHMAX");     fOptionType.push_back("d");
  fMismatchMax=0.2;
  fOptionName.push_back("omicron_TRIGGER_SNRTHRESHOLD");      fOptionType.push_back("d");
  fSNRThreshold=8.0;
  fOptionName.push_back("omicron_TRIGGER_NMAX");              fOptionType.push_back("i");
  fNtriggerMax=1000000;
  fOptionName.push_back("omicron_TRIGGER_CLUSTERING");        fOptionType.push_back("s");
  fClusterAlgo="none";
  fOptionName.push_back("omicron_TRIGGER_CLUSTERDT");         fOptionType.push_back("d");
  fcldt=0.1;
  fOptionName.push_back("omicron_OUTPUT_VERBOSITY");          fOptionType.push_back("i");
  fVerbosity=2;
  fOptionName.push_back("omicron_OUTPUT_FORMAT");             fOptionType.push_back("s");
  fOutFormat="root";
  fOptionName.push_back("omicron_OUTPUT_WRITEPSD");           fOptionType.push_back("i");
  writepsd=false;
  fOptionName.push_back("omicron_OUTPUT_WRITETIMESERIES");    fOptionType.push_back("i");
  writetimeseries=false;

  // read option file
  status_OK=ReadOptions();

  // init process monitoring
  if(fVerbosity) cout<<"Omicron::Omicron: Allocate memory for monitoring container..."<<endl;
  inSegments     = new Segments();
  outSegments    = new Segments* [(int)fChannels.size()];
  chunk_ctr      = new int       [(int)fChannels.size()];
  cor_chunk_ctr  = new int       [(int)fChannels.size()];
  cor_data_ctr   = new int       [(int)fChannels.size()];
  max_chunk_ctr  = new int       [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){
    outSegments[c]    = new Segments();
    chunk_ctr[c]      = 0;
    cor_chunk_ctr[c]  = 0;
    cor_data_ctr[c]   = 0;
    max_chunk_ctr[c]  = 0;
  }

  // DATA
  if(fVerbosity) cout<<"Omicron::Omicron: Allocate memory for data container and procedure..."<<endl;
  ChunkSize   = fSampleFrequency * fChunkDuration;
  OverlapSize = fSampleFrequency * fOverlapDuration;
  SegmentSize = fSampleFrequency * fSegmentDuration;
  ChunkVect   = new double [ChunkSize];
  SegVect     = new double [SegmentSize];
  TukeyWindow = GetTukeyWindow(SegmentSize,OverlapSize);
  offt = new fft(SegmentSize,"FFTW_MEASURE");
  
  // init FFL
  FFL=NULL;
  if(fFflFile.compare("none")){
    if(fVerbosity) cout<<"Omicron::Omicron: Define FFL object..."<<endl;
    FFL = new ffl(fFflFile, fFflFormat, fVerbosity);
    status_OK*=FFL->DefineTmpDir(fMaindir);
    status_OK*=FFL->LoadFrameFile();
  }
  
  // init Streams
  if(fVerbosity) cout<<"Omicron::Omicron: Define Streams object..."<<endl;
  streams = new Streams* [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){
    streams[c] = new Streams(fChannels[c], fVerbosity);
    status_OK*=streams[c]->GetStatus();
  }
  first_Data=true;

  // init Spectrum
  if(fVerbosity) cout<<"Omicron::Omicron: Define Spectrum object..."<<endl;
  double psdsize;
  if(fFreqRange[0]>=1.0) psdsize = fSampleFrequency; // 0.5Hz binning
  else{
    psdsize = (int)floor((double)fSampleFrequency/fFreqRange[0]);// over-binning (factor 2)
    int nextpowerof2=(int)floor(log(psdsize)/log(2));
    psdsize=(int)pow(2.0,(double)nextpowerof2);
  }
  spectrum = new Spectrum(fSampleFrequency,psdsize,fOverlapDuration,fVerbosity);
  status_OK*=spectrum->GetStatus();
  status_OK*=spectrum->GetStatus();
  first_PSD=true;
    
  // init Triggers
  if(fVerbosity) cout<<"Omicron::Omicron: Define Triggers object..."<<endl;
  triggers    = new Triggers* [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){

    // create trigger directory
    system(("mkdir -p "+fOutdir[c]).c_str());

    // init triggers object
    triggers[c] = new Triggers(fOutdir[c],fChannels[c],fOutFormat,fVerbosity);
    triggers[c]->SetNtriggerMax(fNtriggerMax);// maximum number of triggers per file

    // set clustering if any
    if(fClusterAlgo.compare("none")){// clustering is requested
      status_OK*=triggers[c]->SetClustering(fClusterAlgo);// set clustering
      status_OK*=triggers[c]->SetClusterDeltaT(fcldt);// set dt
    }

    // set metadata
    status_OK*=triggers[c]->InitUserMetaData(fOptionName,fOptionType);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],fOutdir[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[1],fChannels[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[2],fInjChan[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],fInjFact[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[4],fFflFile);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[5],fFflFormat);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[6],fSampleFrequency);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[7],fFreqRange[0]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[8],fFreqRange[1]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[9],fQRange[0]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[10],fQRange[1]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[11],fChunkDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[12],fSegmentDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[13],fOverlapDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[14],fMismatchMax);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[15],fSNRThreshold);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[16],fNtriggerMax);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[17],fClusterAlgo);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[18],fcldt);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[19],fVerbosity);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[20],fOutFormat);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[21],(int)writepsd);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[22],(int)writetimeseries);
    triggers[c]->SetMprocessname("Omicron");
    triggers[c]->SetMstreamname(streams[c]->GetName());
    triggers[c]->SetMdetindex(streams[c]->GetDetIndex());
  }

  // init tiles
  if(fVerbosity) cout<<"Omicron::Omicron: Define Tile object..."<<endl;
  tile = new Otile(fSegmentDuration,fOverlapDuration/2,fQRange[0],fQRange[1],fFreqRange[0],fFreqRange[1],fSampleFrequency,fMismatchMax,fSNRThreshold,fVerbosity);
  status_OK*=tile->GetStatus();

  // online specific
  sample_online = new Sample* [(int)fChannels.size()];
  nativesampling_online = new int [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){
    sample_online[c]=NULL;
    nativesampling_online[c]=0;
  }

  if(!status_OK) cerr<<"Omicron::Omicron: Initialization failed"<<endl;
  cout<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
Omicron::~Omicron(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"delete Omicron"<<endl;
  for(int c=0; c<(int)fChannels.size(); c++){
    delete outSegments[c];
    delete streams[c];
    if(sample_online[c]!=NULL) delete sample_online[c];
  }
  delete sample_online;
  delete nativesampling_online;
  delete outSegments;
  delete inSegments;
  delete streams;
  delete chunk_ctr;
  delete cor_chunk_ctr;
  delete cor_data_ctr;
  delete max_chunk_ctr;
  if(FFL!=NULL) delete FFL;
  delete tile;
  delete ChunkVect;
  delete SegVect;
  delete TukeyWindow;
  delete offt;

  fOptionName.clear();
  fOptionType.clear();
  fOutdir.clear();
  fChannels.clear();
  fInjChan.clear();
  fInjFact.clear();
  fFreqRange.clear();
  fQRange.clear();
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::Process(Segments *aSeg){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::Process: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(!aSeg->GetStatus()){
    cerr<<"Omicron::Process: the input Segments is corrupted"<<endl;
    return false;
  }
  if(!aSeg->GetLiveTime()){
    cerr<<"Omicron::Process: the input Segments has a null live time"<<endl;
    return false;
  }

  // add requested segments (try to optimize the update of inSegments)
  if(inSegments->GetLiveTime()&&aSeg->GetStart(0)>=inSegments->GetStart(inSegments->GetNsegments()-1))
    inSegments->Append(aSeg);
  else{
    for(int s=0; s<aSeg->GetNsegments(); s++) inSegments->AddSegment(aSeg->GetStart(s),aSeg->GetEnd(s));
  }

  // data structure
  if(fVerbosity) cout<<"Omicron::Process: initiate data segments..."<<endl;
  Odata *data = new Odata(aSeg, fChunkDuration, fSegmentDuration, fOverlapDuration, fVerbosity);
  if(!data->GetStatus()){
    cerr<<"Omicron::Process: cannot initiate data."<<endl;
    delete data;
    return false;
  }

  // locals
  int dsize;         // data size
  double *dvector;   // data vector before resampling
  double *c_data[2]; // conditionned data vector (Re & Im)
  int sampling_test; // test sampling frequency
  
  // sample structure
  int nativesampling=32768; // just a guess (this is bad if fFreqRange[0]>32768)
  Sample *sample = new Sample(nativesampling,0);
  if(!sample->GetStatus()||!sample->SetWorkingFrequency(fSampleFrequency)||!sample->SetHighPassFrequency(fFreqRange[0])){
    cerr<<"Omicron::Process: cannot initiate sampling."<<endl;
    delete data;
    delete sample;
    return false;
  }
  
  // loop over chunks
  while(data->NewChunk()){
    if(fVerbosity) cout<<"Omicron::Process: chunk "<<data->GetChunkTimeStart()<<"-"<<data->GetChunkTimeEnd()<<endl;
    
    // loop over channels
    for(int c=0; c<(int)fChannels.size(); c++){
      if(fVerbosity) cout<<"Omicron::Process: channel "<<fChannels[c]<<"..."<<endl;
      
      // one more chunk of data
      chunk_ctr[c]++;

      // get data vector
      dvector=FFL->GetData(dsize, fChannels[c], data->GetChunkTimeStart(), data->GetChunkTimeEnd());
      sampling_test=dsize/(data->GetChunkTimeEnd()-data->GetChunkTimeStart());

      // the data chunk was not filled -> skip channel
      if(dsize<=0){
	cerr<<"Omicron::Process: chunk "<<data->GetChunkTimeStart()<<"-"<<data->GetChunkTimeEnd()<<" for channel "<<fChannels[c]<<" is corrupted --> skip"<<endl;
	cor_chunk_ctr[c]++;
	continue;
      }

      // flat data?  -> skip channel
      if(dvector[0]==dvector[dsize-1]){
	cerr<<"Omicron::Process: chunk "<<data->GetChunkTimeStart()<<"-"<<data->GetChunkTimeEnd()<<" for channel "<<fChannels[c]<<" has flat data"<<endl;
	cor_data_ctr[c]++;
	delete dvector;
	continue;
      }

      // Test native sampling frequency
      if(sampling_test!=nativesampling){
	delete sample; sample = new Sample(sampling_test,0);
	nativesampling=sampling_test;// for next round
	if(!sample->GetStatus()||!sample->SetWorkingFrequency(fSampleFrequency)||!sample->SetHighPassFrequency(fFreqRange[0])){
	  cerr<<"Omicron::Process: the Sample structure could not be re-defined"<<endl;
	  delete data;
	  delete sample;
	  delete dvector;
	  return false;
	}
      }

      // transform data vector
      ChunkSize=(data->GetChunkTimeEnd()-data->GetChunkTimeStart())*fSampleFrequency;
      if(!sample->Transform(dsize, dvector, ChunkSize, ChunkVect)){
	cerr<<"Omicron::Process: the data chunk cannot be transformed for "<<fChannels[c]<<" in chunk "<<data->GetChunkTimeStart()<<"-"<<data->GetChunkTimeEnd()<<" --> skip"<<endl;
	cor_data_ctr[c]++;
	delete dvector;
	continue;
      }
      delete dvector;// not needed anymore

      // Save time series if requested
      if(writetimeseries) SaveData(c,ChunkVect,data->GetChunkTimeStart(),data->GetChunkTimeEnd());

      // Make spectrum
      if(!spectrum->LoadData(ChunkSize, ChunkVect)){
	cerr<<"Omicron::Process: the spectrum cannot be computed for "<<fChannels[c]<<" in chunk "<<data->GetChunkTimeStart()<<"-"<<data->GetChunkTimeEnd()<<" --> skip"<<endl;
	cor_data_ctr[c]++;
	continue;
      }

      // Save PSD if requested
      if(writepsd) SavePSD(c,data->GetChunkTimeStart(),data->GetChunkTimeEnd());

      // set new power for this chunk
      if(fVerbosity) cout<<"set power to the tiles for this chunk..."<<endl;
      if(!tile->SetPowerSpectrum(spectrum)){
	cerr<<"Omicron::Process: the spectrum could not be attached to the tile structure for "<<fChannels[c]<<" in chunk "<<data->GetChunkTimeStart()<<"-"<<data->GetChunkTimeEnd()<<" --> skip"<<endl;
	cor_data_ctr[c]++;
	continue;
      }

      //loop over segments
      for(int s=0; s<data->GetNSegments(); s++){
	if(fVerbosity) cout<<"Omicron::Process: segment "<<data->GetSegmentTimeStart(s)+fOverlapDuration/2<<"-"<<data->GetSegmentTimeEnd(s)-fOverlapDuration/2<<endl;
	
	// fill segment vector (time-domain) and apply tukey window
	for(int i=0; i<SegmentSize; i++)
	  SegVect[i] = ChunkVect[s*(SegmentSize-OverlapSize)+i] * TukeyWindow[i];
	
	// get conditioned data
	if(!Condition(&(c_data[0]), &(c_data[1]))){
	  cerr<<"Omicron::Process: conditionned data are corrupted!"<<endl;
	  delete data;
	  delete sample;
	  return false;
	}

	// get triggers above SNR threshold
	if(!tile->GetTriggers(triggers[c],c_data[0],c_data[1],data->GetSegmentTimeStart(s))){
	  cerr<<"Omicron::Process: could not get triggers for channel "<<fChannels[c]
	      <<" in segment starting at "<<data->GetSegmentTimeStart(s)<<endl;
	  delete c_data[0];
	  delete c_data[1];
	  continue;
	}

	// add processed segment
	else
	  triggers[c]->AddSegment(data->GetSegmentTimeStart(s)+fOverlapDuration/2,data->GetSegmentTimeEnd(s)-fOverlapDuration/2);

	delete c_data[0];
	delete c_data[1];

	//stop getting triggers if max flag
	if(triggers[c]->GetMaxFlag()) break;
      }
      
      // don't save if max flag
      if(triggers[c]->GetMaxFlag()){
	cerr<<"Omicron::Process: channel "<<fChannels[c]<<" is maxed-out. This chunk of triggers is not saved"<<endl;
	triggers[c]->Reset();
	max_chunk_ctr[c]++;
	continue; // move to next channel
      }
      
      // save triggers for this chunk
      if(!triggers[c]->Write("ALL", "default").compare("none")){
	cerr<<"Omicron::Process: writing events failed for channel "<<fChannels[c]<<endl;
	return false;
      }

      // processed chunks
      else
	outSegments[c]->AddSegment(data->GetChunkTimeStart()+fOverlapDuration/2,data->GetChunkTimeEnd()-fOverlapDuration/2);
	
    }
  }

  delete sample;
  delete data;
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::Condition(double **aDataRe, double **aDataIm){
////////////////////////////////////////////////////////////////////////////////////

  // fft-forward
  if(!offt->Forward(SegVect)){
    cerr<<"Omicron::Condition: FFT-forward failed"<<endl;
    return false;
  }
    
  // get ffted data
  *aDataRe=offt->GetRe();
  *aDataIm=offt->GetIm();
  if(*aDataRe==NULL||*aDataIm==NULL){
    cerr<<"Omicron::Condition: cannot retrieve FFT data"<<endl;
    return false;
  }

  // zero-out below high-frequency cutoff
  int icutoff = (int)floor(fFreqRange[0]/(double)fSampleFrequency*(double)SegmentSize);
  for(int i=0; i<icutoff; i++){
    (*aDataRe)[i]=0.0;
    (*aDataIm)[i]=0.0;
  }

  // normalize data by the ASD
  double asdval;
  for(int i=icutoff; i<SegmentSize/2; i++){
    asdval=spectrum->GetPower(i*(double)fSampleFrequency/(double)SegmentSize);
    if(asdval<0){
      cerr<<"Omicron::Condition: could not retrieve power"<<endl;
      delete *aDataRe; delete *aDataIm;
      *aDataRe=NULL; *aDataIm=NULL;
      return false;
    }
    asdval=sqrt(asdval);
    (*aDataRe)[i] /= sqrt(asdval);
    (*aDataIm)[i] /= sqrt(asdval);
    // NOTE: no FFT normalization here because we do a FFTback in Oqplane later with no normalization either.
  }

  return true;
}


////////////////////////////////////////////////////////////////////////////////////
int Omicron::ProcessVector(const int aChNumber, const int aInVectSize, double *aInVect, const int aTimeStart){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::ProcessVector: the Omicron object is corrupted"<<endl;
    return -1;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::ProcessVector: channel number "<<aChNumber<<" does not exist"<<endl;
    return -2;
  }
  if(aInVect==NULL){
    cerr<<"Omicron::ProcessVector: input vector is NULL"<<endl;
    return -3;
  }
  if(fChunkDuration!=fSegmentDuration){
    cerr<<"Omicron::ProcessVector: this function can only be used if the chunks and segments have the same size"<<endl;
    return -4;

  }

  // test sampling frequency and update sample if necessary
  int nativesampling_test = aInVectSize/fChunkDuration;
  if(nativesampling_test!=nativesampling_online[aChNumber]){
    nativesampling_online[aChNumber]=nativesampling_test;
    if(sample_online[aChNumber]!=NULL) delete sample_online;
    sample_online[aChNumber] = new Sample(nativesampling_online[aChNumber],0);
    if(!sample_online[aChNumber]->GetStatus()||!sample_online[aChNumber]->SetWorkingFrequency(fSampleFrequency)||!sample_online[aChNumber]->SetHighPassFrequency(fFreqRange[0])){
      cerr<<"Omicron::ProcessVector: cannot initiate sampling for channel "<<fChannels[aChNumber]<<endl;
      delete sample_online[aChNumber]; sample_online[aChNumber]=NULL;
      nativesampling_online[aChNumber]=0;
      return -5;
    }
  }

  // one more chunk of data
  chunk_ctr[aChNumber]++;

  // flat data?  -> skip channel
  if(aInVect[0]==aInVect[aInVectSize-1]){
    cerr<<"Omicron::ProcessVector: chunk "<<aTimeStart<<"-"<<aTimeStart+fChunkDuration<<" for channel "<<fChannels[aChNumber]<<" has flat data"<<endl;
    cor_data_ctr[aChNumber]++;
    return 1;
  }

  // transform data vector
  ChunkSize=fChunkDuration*fSampleFrequency;
  if(!sample_online[aChNumber]->Transform(aInVectSize, aInVect, ChunkSize, ChunkVect)){
    cerr<<"Omicron::ProcessVector: chunk "<<aTimeStart<<"-"<<aTimeStart+fChunkDuration<<" for channel "<<fChannels[aChNumber]<<" cannot be transformed"<<endl;
    cor_data_ctr[aChNumber]++;
    return 2;
  }

  // Save time series if requested
  if(writetimeseries) SaveData(aChNumber,ChunkVect,aTimeStart,aTimeStart+fChunkDuration);

  // Make spectrum
  if(!spectrum->LoadData(ChunkSize, ChunkVect)){
    cerr<<"Omicron::ProcessVector: PSD for chunk "<<aTimeStart<<"-"<<aTimeStart+fChunkDuration<<" for channel "<<fChannels[aChNumber]<<" cannot be computed"<<endl;
    cor_data_ctr[aChNumber]++;
    return 3;
  }

  // Save PSD if requested
  if(writepsd) SavePSD(aChNumber,aTimeStart,aTimeStart+fChunkDuration);

  // set new power for this chunk
  if(!tile->SetPowerSpectrum(spectrum)){
    cerr<<"Omicron::ProcessVector: PSD for chunk "<<aTimeStart<<"-"<<aTimeStart+fChunkDuration<<" for channel "<<fChannels[aChNumber]<<" cannot be attached to the tile structure"<<endl;
    cor_data_ctr[aChNumber]++;
    return 4;
  }

  // fill segment vector (time-domain) and apply tukey window
  for(int i=0; i<SegmentSize; i++)
    SegVect[i] = ChunkVect[i] * TukeyWindow[i];

  // get conditioned data
  double *c_data[2]; // conditionned data vector (Re & Im)
  if(!Condition(&(c_data[0]), &(c_data[1]))){
    cerr<<"Omicron::ProcessVector: conditionned data for chunk "<<aTimeStart<<"-"<<aTimeStart+fChunkDuration<<" for channel "<<fChannels[aChNumber]<<" are corrupted"<<endl;
    cor_data_ctr[aChNumber]++;
    return 5;
  }

  // get triggers above SNR threshold
  if(!tile->GetTriggers(triggers[aChNumber],c_data[0],c_data[1],aTimeStart)){
    cerr<<"Omicron::ProcessVector: could not get triggers for chunk "<<aTimeStart<<"-"<<aTimeStart+fChunkDuration<<" for channel "<<fChannels[aChNumber]<<endl;
    delete c_data[0];
    delete c_data[1];
    return 6;
  }

  // add processed segment
  else
    triggers[aChNumber]->AddSegment(aTimeStart+fOverlapDuration/2,aTimeStart+fChunkDuration-fOverlapDuration/2);
  
  delete c_data[0];
  delete c_data[1];

  //stop getting triggers if max flag
  if(triggers[aChNumber]->GetMaxFlag()){
    cerr<<"Omicron::ProcessVector: channel "<<fChannels[aChNumber]<<" is maxed-out. This chunk is erased"<<endl;
    triggers[aChNumber]->Reset();
    return 7;
  }
  
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::WriteTriggers(const int aChNumber){
  ////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::WriteOnline: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::WriteOnline: channel number "<<aChNumber<<" does not exist"<<endl;
    return -3;
  }
  if(fVerbosity>0) cout<<"writing triggers for channel "<<fChannels[aChNumber]<<"..."<<endl;
  if(triggers[aChNumber]->Write("PROC","default").compare("none")) return true;
  return false;
}
/*
////////////////////////////////////////////////////////////////////////////////////
Segments* Omicron::GetOnlineSegments(const int aChNumber, TH1D *aThr, const double aPadding, const double aInfValue){
////////////////////////////////////////////////////////////////////////////////////
  Segments* empty = new Segments();
  if(!status_OK){
    cerr<<"Omicron::GetOnlineSegments: the Omicron object is corrupted"<<endl;
    return empty;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::GetOnlineSegments: channel number "<<aChNumber<<" does not exist"<<endl;
    return empty;
  }
  delete empty;

  return triggers[aChNumber]->GetTriggerSegments(aThr,aPadding,aInfValue);
}
*/

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SavePSD(const int c, const int s, const int e){
////////////////////////////////////////////////////////////////////////////////////

  TGraph *GPSD = spectrum->GetPSD();
  if(GPSD==NULL) return;

  stringstream ss;
  ss<<"PSD_"<<s<<"_"<<e;
  GPSD->SetName(ss.str().c_str());
  
  TFile *fpsd;
  if(first_PSD){
    fpsd=new TFile((fOutdir[c]+"/PSD_"+fChannels[c]+".root").c_str(),"RECREATE");
    first_PSD=false;
  }
  else fpsd=new TFile((fOutdir[c]+"/PSD_"+fChannels[c]+".root").c_str(),"UPDATE");
  fpsd->cd();
  GPSD->Write();
  fpsd->Close();
  
  delete GPSD;
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveData(const int c, double *aData, const int s, const int e){
////////////////////////////////////////////////////////////////////////////////////

  int N = (e-s)*fSampleFrequency;
  TGraph *GDATA = new TGraph(N);
  if(GDATA==NULL) return;

  stringstream ss;
  ss<<"data_"<<s<<"_"<<e;
  GDATA->SetName(ss.str().c_str());

  for(int i=0; i<N; i++) GDATA->SetPoint(i,(double)s+(double)i*(double)(e-s)/(double)N,aData[i]);
  
  TFile *fdata;
  if(first_Data){
    fdata=new TFile((fOutdir[c]+"/DATA_"+fChannels[c]+".root").c_str(),"RECREATE");
    first_Data=false;
  }
  else fdata=new TFile((fOutdir[c]+"/DATA_"+fChannels[c]+".root").c_str(),"UPDATE");
  fdata->cd();
  GDATA->Write();
  fdata->Close();
  
  delete GDATA;
  return;
}

////////////////////////////////////////////////////////////////////////////////////
double* Omicron::GetTukeyWindow(const int aSize, const int aFractionSize){
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

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::ReadOptions(void){
////////////////////////////////////////////////////////////////////////////////////

  // check that the option file exists
  if(!IsTextFile(fOptionFile)){
    cerr<<"Omicron::ReadOptions: option file "<<fOptionFile<<" cannot be found"<<endl;
    return false;
  }

  // create parser
  IO *io = new IO(fOptionFile.c_str());
  
  //***** output directory *****
  if(!io->GetOpt("OUTPUT","DIRECTORY", fMaindir)) fMaindir=".";
  if(!IsDirectory(fMaindir)){
    cerr<<"Omicron::ReadOptions: output directory cannot be found: OUTPUT/DIRECTORY"<<endl;
    return false;
  }
  //*****************************

  //***** List of channels *****
  if(!io->GetOpt("DATA","CHANNELS", fChannels)){
    cerr<<"Omicron::ReadOptions: you must provide at least one channel: DATA/CHANNELS"<<endl;
    return false;
  }
  for(int c=0; c<(int)fChannels.size(); c++) // trigger output dir
    fOutdir.push_back(fMaindir+"/"+fChannels[c]);
  //*****************************

  //***** injection channels *****
  if(io->GetOpt("INJECTION","CHANNELS", fInjChan)){
    if(fInjChan.size()!=fChannels.size()){
      cerr<<"Omicron::ReadOptions: INJECTION/CHANNELS is inconsistent with the number of channels"<<endl;
      return false;
    }
    if(io->GetOpt("INJECTION","FACTORS", fInjFact)){
      if(fInjFact.size()!=fInjChan.size()){
	cerr<<"Omicron::ReadOptions: INJECTION/FACTORS is inconsistent with the number of channels"<<endl;
	return false;
      }
    }
    else{
      for(int i=0; i<(int)fChannels.size(); i++) fInjFact.push_back(1.0);
    }
  }
  else{
    for(int i=0; i<(int)fChannels.size(); i++){
      fInjChan.push_back("none");
      fInjFact.push_back(0.0);
    }
  }
  //*****************************

  /*
  //***** network mode *****
  fDetectors.clear();
  fInjFile="none";
  io.GetOpt("NETWORK","DETECTORS", fDetectors);
  if(fDetectors.size()){// network mode activated
    if(fDetectors.size()!=fChannels.size()){
      cerr<<"Omicron::ReadOptions: NETWORK/DETECTORS is inconsistent with the number of channels"<<endl;
      return false;
    }
    io.GetOpt("NETWORK","INJFILE", fInjFile);
  }
  
  fOptionName.push_back("omicron_NETWORK_DETECTORS");
  fOptionType.push_back("s");
  fOptionName.push_back("omicron_NETWORK_INJFILE");
  fOptionType.push_back("s");
  */
  //*****************************
  
  //***** ffl file *****
  if(io->GetOpt("DATA","LCF", fFflFile)) fFflFormat="lcf";
  if(io->GetOpt("DATA","FFL", fFflFile)) fFflFormat="ffl";
  else{
    fFflFile="none";
    fFflFormat="none";
  }
  //*****************************

  //***** Sampling frequency *****
  if(!io->GetOpt("DATA","SAMPLEFREQUENCY", fSampleFrequency)){
    cerr<<"Omicron::ReadOptions: A working sampling frequency (PARAMETER/SAMPLEFREQUENCY) is required"<<endl;
    return false;
  }
  if(fSampleFrequency<16||fSampleFrequency>40000){
    cerr<<"Omicron::ReadOptions: Sampling frequency "<<fSampleFrequency<<"Hz (PARAMETER/SAMPLEFREQUENCY) is not reasonable"<<endl;
    return false;
  }
  //*****************************

  //***** Frequency range *****
  fFreqRange.clear();
  if(!io->GetOpt("PARAMETER","FREQUENCYRANGE", fFreqRange)){
    cerr<<"Omicron::ReadOptions: A search frequency range (PARAMETER/FREQUENCYRANGE) is required"<<endl;
    return false;
  }
  if(fFreqRange.size()!=2){
    cerr<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) is not correct"<<endl;
    return false;
  }
  //*****************************

  //***** Q range *****
  fQRange.clear();
  if(!io->GetOpt("PARAMETER","QRANGE", fQRange)){
    cerr<<"Omicron::ReadOptions: A search Q range (PARAMETER/QRANGE) is required"<<endl;
    return false;
  }
  if(fQRange.size()!=2){
    cerr<<"Omicron::ReadOptions: Q range (PARAMETER/QRANGE) is not correct"<<endl;
    return false;
  }
  //*****************************

  //***** timing *****
  if(!io->GetOpt("PARAMETER","CHUNKDURATION", fChunkDuration)) fChunkDuration=484;
  if(!io->GetOpt("PARAMETER","BLOCKDURATION", fSegmentDuration)) fSegmentDuration=64;
  if(!io->GetOpt("PARAMETER","OVERLAPDURATION", fOverlapDuration)) fOverlapDuration=8;
  //*****************************

  //***** maximum mismatch *****
  if(!io->GetOpt("PARAMETER","MISMATCHMAX", fMismatchMax)) fMismatchMax=0.2;
  //*****************************

  //***** SNR Threshold *****
  if(!io->GetOpt("TRIGGER","SNRTHRESHOLD", fSNRThreshold)) fSNRThreshold=8.0;
  //*****************************
  
  //***** set trigger limit *****
  fNtriggerMax=1000000;
  double ratemax;
  if(io->GetOpt("TRIGGER","RATEMAX", ratemax)) fNtriggerMax=(int)ceil(ratemax*fChunkDuration);
  else io->GetOpt("TRIGGER","NMAX", fNtriggerMax);
  //*****************************

  //***** set clustering *****
  if(!io->GetOpt("TRIGGER","CLUSTERING", fClusterAlgo)) fClusterAlgo="none";
  if(!io->GetOpt("TRIGGER","CLUSTERDT", fcldt))         fcldt=0.1;
  //*****************************

  //***** set verbosity *****
  if(!io->GetOpt("OUTPUT","VERBOSITY", fVerbosity)) fVerbosity=0;
  //*****************************

  //***** set output format ***** 
  if(!io->GetOpt("OUTPUT","FORMAT", fOutFormat)) fOutFormat="root";
  //*****************************

  //***** set writing flags *****
  if(!io->GetOpt("OUTPUT","WRITEPSD", writepsd)) writepsd=0;
  if(!io->GetOpt("OUTPUT","WRITETIMESERIES", writetimeseries)) writetimeseries=0;;
  //*****************************

  // dump options
  if(fVerbosity>1) io->Dump(cout);
  delete io;

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::PrintStatusInfo(void){
////////////////////////////////////////////////////////////////////////////////////

  
  cout<<"\n************* Omicron status info *************"<<endl;
  cout<<"requested start = "<<inSegments->GetStart(0)<<endl;
  cout<<"requested end   = "<<inSegments->GetEnd(inSegments->GetNsegments()-1)<<endl;
  cout<<"requested livetime   = "<<inSegments->GetEnd(inSegments->GetNsegments()-1)<<endl;
  for(int c=0; c<(int)fChannels.size(); c++){
    cout<<"\n*** "<<fChannels[c]<<endl;
    if(outSegments[c]->GetNsegments()){
      cout<<"start_out                  = "<<outSegments[c]->GetStart(0)<<endl;
      cout<<"end_out                    = "<<outSegments[c]->GetEnd(outSegments[c]->GetNsegments()-1)<<endl;
      cout<<"processed livetime         = "<<outSegments[c]->GetLiveTime()<<" ("<<outSegments[c]->GetLiveTime()/inSegments->GetLiveTime()*100<<"%)"<<endl;
    }
    else{
      cout<<"start_out                  = -1"<<endl;
      cout<<"end_out                    = -1"<<endl;
      cout<<"processed livetime         = "<<"0 (0%)"<<endl;
    }
    cout<<"number of chunks           = "<<chunk_ctr[c]<<endl;
    cout<<"number of corrupted chunks = "<<cor_chunk_ctr[c]<<endl;
    cout<<"number of maxed-out chunks = "<<max_chunk_ctr[c]<<endl;
  }
  cout<<"***********************************************\n"<<endl;

  return;
}
