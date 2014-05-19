//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

ClassImp(Omicron)

////////////////////////////////////////////////////////////////////////////////////
Omicron::Omicron(const string aOptionFile){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // inputs
  fOptionFile=aOptionFile;
  
  // read option file
  status_OK=ReadOptions();
 
  // init process monitoring
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

  // init FFL
  if(fFflFile.compare("none")){
    FFL = new ffl(fFflFile, fFflFormat, fVerbosity);
    status_OK*=FFL->DefineTmpDir(fMaindir);
    status_OK*=FFL->LoadFrameFile();
  }

  // init Streams
  streams = new Streams* [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++)
    streams[c] = new Streams(fChannels[c], fVerbosity);
  
  // init Spectrum
  double psdsize;
  if(fFreqRange[0]<=0||fFreqRange[0]>=1.0) psdsize = fSamplingFrequency; // 0.5Hz binning
  else{
    psdsize = (int)floor((double)fSamplingFrequency/fFreqRange[0]);// or over-binning (factor 2)
    int nextpowerof2=(int)floor(log(psdsize)/log(2));
    psdsize=(int)pow(2.0,(double)nextpowerof2);
  }
  spectrum = new Spectrum(fSampleFrequency,psdsize,fOverlapDuration,fVerbosity);
  status_OK*=spectrum->GetStatus();
    
  // init Triggers
  triggers    = new Triggers* [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){
    system(("mkdir -p "+fOutdir[c]).c_str());// output dir
    triggers[c] = new Triggers(fOutdir[c],fChannels[c],fOutFormat,fVerbosity);
    status_OK*=triggers[c]->SetNtriggerMax(fNtriggerMax);// maximum number of triggers per file
    if(fClusterAlgo.compare("none")){// clustering is requested
      status_OK*=triggers[c]->SetClustering(fClusterAlgo);// set clustering
      status_OK*=triggers[c]->SetClusterDeltaT(fcldt);// set dt
    }
    status_OK*=triggers[c]->InitUserMetaData(fOptionName,fOptionType);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],fOutdir[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[1],fChannels[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[2],fInjChan[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],fInjFact[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[4],fFflFile);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[5],fFflFormat);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[9],fSampleFrequency);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[10],fFreqRange[0]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[11],fFreqRange[1]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[12],fQRange[0]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[13],fQRange[1]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[14],fChunkDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[15],fSegmentDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[16],fOverlapDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[17],fMismatchMax);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[18],fSNRThreshold);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[19],fNtriggerMax);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[20],fClusterAlgo);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[21],fcldt);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[22],fVerbosity);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[23],fOutFormat);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[24],(int)writepsd);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[25],(int)writetimeseries);
    status_OK*=triggers[c]->SetMprocessname("Omicron");
    status_OK*=triggers[c]->SetMstreamname(streams[c]->GetName());
    status_OK*=triggers[c]->SetMdetindex(streams[c]->GetDetIndex());
  }

  // init tiles
  tile = new Otile(fSegmentDuration,fOverlapDuration/2,fQRange[0],fQRange[1],fFreqRange[0],fFreqRange[1],fSampleFrequency,fMismatchMax,fSNRThreshold,fVerbosity);
  status_OK=tile->GetStatus();


  if(!status_OK) cerr<<"Omicron::Omicron: Initialization failed"<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
Omicron::~Omicron(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"delete Omicron"<<endl;
  for(int c=0; c<(int)fChannels.size(); c++){
    delete outSegments[c];
    delete streams[c];
  }
  delete outSegments;
  delete inSegments;
  delete streams;
  delete chunk_ctr;
  delete cor_chunk_ctr;
  delete cor_data_ctr;
  delete max_chunk_ctr;
  if(FFL!=NULL) delete FFL;
  delete tile;

  fOptionName.clear();
  fOptionType.clear();
  fChannels.clear();
  fInjChan.clear();
  fInjFact.clear();
  fDetectors.clear();
  fNativeFrequency.clear();
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
  if(!aSeg->GetLivetime()){
    cerr<<"Omicron::Process: the input Segments have a null live time"<<endl;
    return true;
  }

  // add requested segments (try to optimize the update of inSegments)
  if(inSegments->GetLiveTime()&&aSeg->GetStart(0)>=inSegments->GetStart(inSegments->GetNsegments()-1))
    inSegments->Append(aSeg);
  else{
    for(int s=0; s<aSeg->GetNsegments(); s++) inSegments->AddSegment(aSeg->GetStart(s),aSeg->GetEnd(s));
  }

  // locals
  // int psdsize;
  int dsize;        // data size
  double *dvector;  // data vector
  double *c_data[2];// conditionned data vector (Re & Im)
  double *psd;      // psd vector - DO NOT DELETE
  int psdsize;      // psd vector size
  int sampling;     // sampling frequency
  int newsize;      // vector size after conditioning

  // data structure
  Odata *data = new Odata(aSeg, fChunkDuration, fSegmentDuration, fOverlapDuration);

  // sample structure
  int nativesampling=4096; // just a guess (this is bad if fFreqRange[0]>4096)
  Sample *sample = new Sample(nativesampling,0);
  sample->SetWorkingFrequency(fSampleFrequency);
  sample->SetHighPassFrequency(fFreqRange[0]);

  if(!odata[c]->GetStatus()){
    cerr<<"Omicron::Process: cannot initiate the data."<<endl;
    delete data;
    return false;
  }

  // loop over chunks
  while(data->NewChunk()>0){
    if(fVerbosity) cout<<"Omicron::Process: chunk "<<data->GetChunkTimeStart()<<"-"<<data->GetChunkTimeEnd()<<endl;
    
    // loop over channels
    for(int c=0; c<(int)fChannels.size(); c++){
      if(fVerbosity) cout<<"Omicron::Process: channel "<<fChannels[c]<<"..."<<endl;
      
      // one more chunk of data
      chunk_ctr[c]++;

      // get data vector
      dvector=FFL->GetData(dsize, fChannels[c], data->GetChunkTimeStart(), data->GetChunkTimeEnd());
      sampling=dsize/(data->GetChunkTimeEnd()-data->GetChunkTimeStart());

      // the data chunk is corrupted -> skip channel
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
      if(sampling!=nativesampling){
	delete sample; sample = new Sample(nativesampling,0);
	if(!sample->SetWorkingFrequency(fSampleFrequency)||!sample->SetHighPassFrequency(fFreqRange[0])){
	  delete data;
	  delete dvector;
	  delete sample;
	  return false;
	}
      }

      // Sample data vector
      newsize=(data->GetChunkTimeEnd()-data->GetChunkTimeStart())*fSampleFrequency;
      sample->Transform(dsize, dvector, newsize, double *datavector);
      delete dvector;// not needed anymore

      // Make spectrum
      if(!spectrum->LoadData(newsize, datavector)){
	cerr<<"Omicron::Process: the spectrum cannot be computed for "<<fChannels[c]<<" in chunk "<<data->GetChunkTimeStart()<<"-"<<data->GetChunkTimeEnd()<<" --> skip"<<endl;
	cor_data_ctr[c]++;
	continue;
      }


      // write chunk info if requested
      //if(writetimeseries) odata[c]->WriteTimeSeries(fOutdir[c]);
      //if(writepsd)        odata[c]->WritePSD(fOutdir[c]);

      //loop over segments
      for(int s=0; s<data->GetNSegments(); s++){
	if(fVerbosity) cout<<"Omicron::Process: segment "<<data->GetSegmentTimeStart(s)+fOverlapDuration/2<<"-"<<data->GetSegmentTimeEnd(s)-fOverlapDuration/2<<endl;
	
	// get conditioned data
	if(!data->GetConditionedData(&(c_data[0]), &(c_data[1]), s, spectrum)){
	  cerr<<"Omicron::Process: conditionned data are corrupted!"<<endl;
	  return false;
	}

	// set new power for this chunk
	if(!s){
	  if(fVerbosity) cout<<"set power to the tiles for this chunk..."<<endl;
	  tile->SetPowerSpectrum(psd,psdsize);
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
	cerr<<"Omicron::Process: channel "<<fChannels[c]<<" is maxed-out. This chunk is not saved"<<endl;
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
	outSegments[c]->AddSegment(data[c]->GetChunkTimeStart()+fOverlapDuration/2,data[c]->GetChunkTimeEnd()-fOverlapDuration/2);
	
    }
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
int Omicron::ProcessOnline(const int aChNumber, FrVect *aVect){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::ProcessOnline: the Omicron object is corrupted"<<endl;
    return -1;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::ProcessOnline: channel number "<<aChNumber<<" does not exist"<<endl;
    return -3;
  }

  if(fVerbosity>0) cout<<"processing channel "<<fChannels[aChNumber]<<"..."<<endl;
      
  // read and condition data
  if(!odata[aChNumber]->ReadVect(aVect)){
    cerr<<"Omicron::ProcessOnline: channel number "<<aChNumber<<": vector cannot be read"<<endl;
    return 1;
  }
                  
  // write chunk info if requested
  if(writetimeseries) status_OK*=odata[aChNumber]->WriteTimeSeries(fOutdir[aChNumber]);
  if(writepsd)        status_OK*=odata[aChNumber]->WritePSD(fOutdir[aChNumber]);

  // get conditioned data
  int psdsize;
  double *c_data[2];
  if(!odata[aChNumber]->GetConditionedData(0,&(c_data[0]),&(c_data[1]),&psd,psdsize)){
    cerr<<"Omicron::ProcessOnline: conditionned data are corrupted!"<<endl;
    return 2;
  }
  	
  // set new power spectrum
  tile->SetPowerSpectrum(psd,psdsize);

  //get triggers
  cout<<" "<<fChannels[aChNumber]<<" Extracting triggers in "<<odata[aChNumber]->GetSegmentTimeStart(0)+fOverlapDuration/2<<"-"<<odata[aChNumber]->GetSegmentTimeStart(0)+fSegmentDuration-fOverlapDuration/2<<endl;
  if(!tile->GetTriggers(triggers[aChNumber],c_data[0],c_data[1],odata[aChNumber]->GetSegmentTimeStart(0))){
    cerr<<"Omicron::ProcessOnline: could not get trigger for channel "<<fChannels[aChNumber]
	<<" in segment starting at "<<odata[aChNumber]->GetSegmentTimeStart(0)<<endl;
    return 3;
  }
  else
    triggers[aChNumber]->AddSegment(odata[aChNumber]->GetSegmentTimeStart(0)+fOverlapDuration/2,odata[aChNumber]->GetSegmentTimeStart(0)+fSegmentDuration-fOverlapDuration/2);
    
  delete c_data[0];
  delete c_data[1];

  // don't save if max flag
  if(triggers[aChNumber]->GetMaxFlag()){
    cerr<<"Omicron::ProcessOnline: channel "<<fChannels[aChNumber]<<" is maxed-out. This chunk is not saved"<<endl;
    triggers[aChNumber]->Reset();
    return 4;
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::WriteOnline(const int aChNumber){
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

////////////////////////////////////////////////////////////////////////////////////
int Omicron::GetNativeSampleFrequency(const int aChNumber){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::GetNativeSampleFrequency: the Omicron object is corrupted"<<endl;
    return -1;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::GetNativeSampleFrequency: channel number "<<aChNumber<<" does not exist"<<endl;
    return -1;
  }
  return fNativeFrequency[aChNumber];
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
  
  // output meta-data
  fOptionName.clear();
  fOptionType.clear();

  //***** output directory *****
  if(!io->GetOpt("OUTPUT","DIRECTORY", maindir)) fMaindir=".";
  if(!IsDirectory(fMaindir)){
    cerr<<"Omicron::ReadOptions: output directory cannot be found: OUTPUT/DIRECTORY"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_OUTPUT_DIRECTORY");
  fOptionType.push_back("s");
  //*****************************

  //***** List of channels *****
  fChannels.clear();
  if(!io->GetOpt("DATA","CHANNELS", fChannels)){
    cerr<<"Omicron::ReadOptions: you must provide at least one channel: DATA/CHANNELS"<<endl;
    return false;
  }
  if(fChannels.size()>NDATASTREAMS){
    cerr<<"Omicron::ReadOptions: The number of channels cannot exceed "<<NDATASTREAMS<<endl;
    return false;
  }
  for(int c=0; c<(int)fChannels.size(); c++) fOutdir.push_back(fMaindir+"/"+fChannels[c]);
  fOptionName.push_back("omicron_DATA_CHANNEL");
  fOptionType.push_back("s");
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
      for(int i=0; i<(int)fInjChan.size(); i++) fInjFact.push_back(1.0);
    }
  }
  fOptionName.push_back("omicron_INJECTION_CHANNEL");
  fOptionType.push_back("s");
  fOptionName.push_back("omicron_INJECTION_FACTOR");
  fOptionType.push_back("d");
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
  */
  fOptionName.push_back("omicron_NETWORK_DETECTORS");
  fOptionType.push_back("s");
  fOptionName.push_back("omicron_NETWORK_INJFILE");
  fOptionType.push_back("s");
  //*****************************
  
  //***** ffl file *****
  if(io->GetOpt("DATA","LCF", fFflFile)) fFflFormat="lcf";
  if(io->GetOpt("DATA","FFL", fFflFile)) fFflFormat="ffl";
  else{
    fFflFile="none";
    fFflFormat="none";
  }
  fOptionName.push_back("omicron_DATA_FFL");
  fOptionType.push_back("s");
  //*****************************
    
  //***** Native channel frequencies *****
  //fNativeFrequency.clear();
  //io.GetOpt("DATA","NATIVEFREQUENCY", fNativeFrequency);
  //if(fNativeFrequency.size()!=fChannels.size()){
  //cerr<<"Omicron::ReadOptions: Each channel must be given a native frequency: PARAMETER/NATIVEFREQUENCY"<<endl;
  //return false;
  //}
  //fOptionName.push_back("omicron_DATA_NATIVEFREQUENCY");
  //fOptionType.push_back("i");
  //*****************************

  //***** Sampling frequency *****
  if(!io->GetOpt("DATA","SAMPLEFREQUENCY", fSampleFrequency)){
    cerr<<"Omicron::ReadOptions: A working sampling frequency (PARAMETER/SAMPLEFREQUENCY) is required"<<endl;
    return false;
  }
  if(fSampleFrequency<16||fSampleFrequency>40000){
    cerr<<"Omicron::ReadOptions: Sampling frequency "<<fSampleFrequency<<" (PARAMETER/SAMPLEFREQUENCY) is not reasonable"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_DATA_SAMPLEFREQUENCY");
  fOptionType.push_back("i");
  //*****************************

  //***** Frequency range *****
  if(!io->GetOpt("PARAMETER","FREQUENCYRANGE", fFreqRange)){
    cerr<<"Omicron::ReadOptions: A search frequency range (PARAMETER/FREQUENCYRANGE) is required"<<endl;
    return false;
  }
  if(fFreqRange.size()!=2){
    cerr<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) is not correct"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_PARAMETER_FMIN");
  fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_FMAX");
  fOptionType.push_back("d");
  //*****************************

  //***** Q range *****
  if(!io.GetOpt("PARAMETER","QRANGE", fQRange)){
    cerr<<"Omicron::ReadOptions: A search Q range (PARAMETER/QRANGE) is required"<<endl;
    return false;
  }
  if(fQRange.size()!=2){
    cerr<<"Omicron::ReadOptions: Q range (PARAMETER/QRANGE) is not correct"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_PARAMETER_QMIN");
  fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_QMAX");
  fOptionType.push_back("d");
  //*****************************

  //***** timing *****
  if(!io->GetOpt("PARAMETER","CHUNKDURATION", fChunkDuration)) fChunkDuration=484;
  if(!io->GetOpt("PARAMETER","BLOCKDURATION", fSegmentDuration)) fSegmentDuration=64;
  if(!io->GetOpt("PARAMETER","OVERLAPDURATION", fOverlapDuration)) fOverlapDuration=8;
  /*
  if(!IsPowerOfTwo(fSegmentDuration)){
    cerr<<"Omicron::ReadOptions: Block duration (PARAMETER/BLOCKDURATION) should be a power of 2"<<endl;
    return false;
  }
  if(fOverlapDuration%2){
    cerr<<"Omicron::ReadOptions: Overlap duration (PARAMETER/OVERLAPDURATION) should be an even integer value"<<endl;
    return false;
  }
  if((fChunkDuration-fOverlapDuration)%(fSegmentDuration-fOverlapDuration)){
    cerr<<"Omicron::ReadOptions: inconsistency in the segments structure"<<endl;
    cerr<<"For example, you could use:"<<endl;
    cerr<<"PARAMETER/CHUNKDURATION:   "<<(fChunkDuration-fOverlapDuration)/(fSegmentDuration-fOverlapDuration)*(fSegmentDuration-fOverlapDuration)+fOverlapDuration<<endl;
    cerr<<"PARAMETER/BLOCKDURATION:   "<<fSegmentDuration<<endl;
    cerr<<"PARAMETER/OVERLAPDURATION: "<<fOverlapDuration<<endl;
    return false;
  }
  */
  fOptionName.push_back("omicron_PARAMETER_CHUNKDURATION");
  fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_BLOCKDURATION");
  fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_OVERLAPDURATION");
  fOptionType.push_back("i");
  //*****************************

  //***** maximum mismatch *****
  if(!io->GetOpt("PARAMETER","MISMATCHMAX", fMismatchMax)) fMismatchMax=0.2;
  fOptionName.push_back("omicron_PARAMETER_MISMATCHMAX");
  fOptionType.push_back("d");
  //*****************************

  //***** SNR Threshold *****
  if(!io->GetOpt("TRIGGER","SNRTHRESHOLD", fSNRThreshold)) fSNRThreshold=8.0;
  fOptionName.push_back("omicron_TRIGGER_SNRTHRESHOLD");
  fOptionType.push_back("d");
  //*****************************
  
  //***** set trigger limit *****
  fNtriggerMax=1000000;
  double ratemax;
  if(io->GetOpt("TRIGGER","RATEMAX", ratemax)) fNtriggerMax=(int)ceil(ratemax*fChunkDuration);
  else io->GetOpt("TRIGGER","NMAX", fNtriggerMax);
  fOptionName.push_back("omicron_TRIGGER_NMAX");
  fOptionType.push_back("i");
  //*****************************

  //***** set clustering *****
  if(!io->GetOpt("TRIGGER","CLUSTERING", fClusterAlgo)) fClusterAlgo="none";
  if(!io->GetOpt("TRIGGER","CLUSTERDT", fcldt))         fcldt=0.1;
  fOptionName.push_back("omicron_TRIGGER_CLUSTERING");
  fOptionType.push_back("s");
  fOptionName.push_back("omicron_TRIGGER_CLUSTERDT");
  fOptionType.push_back("d");
  //*****************************

  //***** set verbosity *****
  if(!io->GetOpt("OUTPUT","VERBOSITY", fVerbosity)) fVerbosity=0;
  fOptionName.push_back("omicron_OUTPUT_VERBOSITY");
  fOptionType.push_back("i");
  //*****************************

  //***** set verbosity ***** 
  if(!io->GetOpt("OUTPUT","FORMAT", fOutFormat)) fOutFormat="root";
  fOptionName.push_back("omicron_OUTPUT_FORMAT");
  fOptionType.push_back("s");
  //*****************************

  //***** set writing flags *****
  if(!io->GetOpt("OUTPUT","WRITEPSD", writepsd)) writepsd=0;
  if(!io->GetOpt("OUTPUT","WRITETIMESERIES", writetimeseries)) writetimeseries=0;;
  fOptionName.push_back("omicron_OUTPUT_WRITEPSD");
  fOptionType.push_back("i");
  fOptionName.push_back("omicron_OUTPUT_WRITETIMESERIES");
  fOptionType.push_back("i");
  //*****************************

  // dump options
  if(fVerbosity>1) io.Dump(cout);
  delete io;

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::PrintStatusInfo(void){
////////////////////////////////////////////////////////////////////////////////////

  
  cout<<"\n************* Omicron status info *************"<<endl;
  cout<<"requested start = "<<fSegments->GetStart(0)<<endl;
  cout<<"requested end   = "<<fSegments->GetEnd(fSegments->GetNsegments()-1)<<endl;
  cout<<"requested livetime   = "<<fSegments->GetEnd(fSegments->GetNsegments()-1)<<endl;
  for(int c=0; c<(int)fChannels.size(); c++){
    cout<<"\n*** "<<fChannels[c]<<endl;
    if(outSegments[c]->GetNsegments()){
      cout<<"start_out                  = "<<outSegments[c]->GetStart(0)<<endl;
      cout<<"end_out                    = "<<outSegments[c]->GetEnd(outSegments[c]->GetNsegments()-1)<<endl;
      cout<<"processed livetime         = "<<outSegments[c]->GetLiveTime()<<" ("<<outSegments[c]->GetLiveTime()/fSegments->GetLiveTime()*100<<"%)"<<endl;
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
