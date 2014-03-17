//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

ClassImp(Omicron)

////////////////////////////////////////////////////////////////////////////////////
Omicron::Omicron(Segments *aSegments, const string aOptionFile){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // inputs
  fSegments = aSegments;
  fOptionFile=aOptionFile;

  // read option file
  status_OK=ReadOptions();
 
  // init process monitoring
  outSegments    = new Segments* [(int)fChannels.size()];
  chunk_ctr      = new int       [(int)fChannels.size()];
  cor_chunk_ctr  = new int       [(int)fChannels.size()];
  max_chunk_ctr  = new int       [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){
    outSegments[c]    = new Segments();
    chunk_ctr[c]      = 0;
    cor_chunk_ctr[c]  = 0;
    max_chunk_ctr[c]  = 0;
  }

  // init network mode
  Net=NULL;
  string NetName="";
  for(int d=0; d<(int)fDetectors.size(); d++) NetName+=fDetectors[d];
  if(fDetectors.size()){
    if(fVerbosity) cout<<"Omicron::Omicron: define network..."<<endl;
    Net = new Network(NetName, fVerbosity);
    status_OK*=Net->GetNdetectors();
  }
  
  // init inject object at the working sampling frequency
  Inj=NULL;
  if(fDetectors.size()&&fInjFile.compare("none")){
    if(fVerbosity) cout<<"Omicron::Omicron: define injection engine..."<<endl;
    Inj = new Inject(Net,fSampleFrequency,fVerbosity);
    status_OK*=Inj->SetInjectionSet(fInjFile);
  }
  
  // init trigger objects
  for(int c=0; c<(int)fChannels.size(); c++){
    if(fVerbosity) cout<<"Omicron::Omicron: define trigger structures for "<<fChannels[c]<<"..."<<endl;
    triggers[c] = new Triggers(fOutdir[c],fChannels[c],fOutFormat,fVerbosity);
    triggers[c]->SetNtriggerMax(fNtriggerMax);// maximum number of triggers per file
    if(!fClusterAlgo.empty()){// clustering is requested
      status_OK*=triggers[c]->SetClustering(fClusterAlgo);// set clustering
      status_OK*=triggers[c]->SetClusterDeltaT(fcldt);// set dt
    }
  }

  // save Omicron meta-data
  for(int c=0; c<(int)fChannels.size(); c++){
    status_OK*=triggers[c]->InitUserMetaData(fOptionName,fOptionType);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],fOutdir[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[1],fChannels[c]);
    if(fInjChan.size()) status_OK*=triggers[c]->SetUserMetaData(fOptionName[2],fInjChan[c]);
    else status_OK*=triggers[c]->SetUserMetaData(fOptionName[2],"none");
    if(fInjFact.size()) status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],fInjFact[c]);
    else  status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],0.0);
    if(fDetectors.size()) status_OK*=triggers[c]->SetUserMetaData(fOptionName[4],fDetectors[c]);
    else status_OK*=triggers[c]->SetUserMetaData(fOptionName[4],"none");
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[5],fInjFile);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[6],fFflFile);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[7],fNativeFrequency[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[8],fSampleFrequency);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[9],fFreqRange[0]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[10],fFreqRange[1]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[11],fQRange[0]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[12],fQRange[1]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[13],fChunkDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[14],fSegmentDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[15],fOverlapDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[16],fMismatchMax);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[17],fSNRThreshold);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[18],fNtriggerMax);
    if(fClusterAlgo.empty()) status_OK*=triggers[c]->SetUserMetaData(fOptionName[19],"none");
    else status_OK*=triggers[c]->SetUserMetaData(fOptionName[19],fClusterAlgo);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[20],fcldt);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[21],fVerbosity);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[22],fOutFormat);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[23],(int)writepsd);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[24],(int)writetimeseries);
    triggers[c]->SetMprocessname("Omicron");
    triggers[c]->SetMstreamname(fChannels[c]);
    triggers[c]->SetMdetindex(GetDetIndex(fChannels[c].substr(0,2)));
  }

  // init data objects
  for(int c=0; c<(int)fChannels.size(); c++){
    if(fVerbosity) cout<<"Omicron::Omicron: define data structures for "<<fChannels[c]<<"..."<<endl;
    odata[c] = new Odata(fFflFile,fChannels[c],fNativeFrequency[c],fSampleFrequency,fSegments,fChunkDuration,fSegmentDuration,fOverlapDuration,fFreqRange[0],fVerbosity);
    if(fInjChan.size()) 
      status_OK*=odata[c]->SetInjectionChannel(fInjChan[c],fInjFact[c]);// add injection channel
    if(Inj!=NULL) 
      status_OK*=odata[c]->SetInject(Inj,Net->GetDetIndex(fDetectors[c]));// add inject object
    status_OK*=odata[c]->GetStatus();// update status
  }

  // make tiling
  if(fVerbosity) cout<<"Omicron::Omicron: make tiling..."<<endl;
  tile = new Otile(fSegmentDuration,fOverlapDuration/2,fQRange[0],fQRange[1],fFreqRange[0],fFreqRange[1],fSampleFrequency,fMismatchMax,fSNRThreshold,fVerbosity);
  status_OK=tile->GetStatus();
  
  if(!status_OK) cerr<<"Omicron::Omicron: Initialization failed"<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
Omicron::~Omicron(void){
////////////////////////////////////////////////////////////////////////////////////
  cout<<"delete Omicron"<<endl;

  delete tile;
  for(int c=0; c<(int)fChannels.size(); c++){
    delete outSegments[c];
    delete odata[c];
    delete triggers[c];
  }
  delete outSegments;
  delete chunk_ctr;
  delete cor_chunk_ctr;
  delete max_chunk_ctr;
  if(Net!=NULL) delete Net;
  if(Inj!=NULL) delete Inj;
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
bool Omicron::Process(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::Process: the Omicron object is corrupted"<<endl;
    return false;
  }

  // locals
  int psdsize;
  double *c_data[2];
  
  // looping status
  int keep_looping=0; 

  // loop over chunks
  while(keep_looping>=0){
    
    // loop over channels
    for(int c=0; c<(int)fChannels.size(); c++){
      if(fVerbosity) cout<<"processing channel "<<fChannels[c]<<"..."<<endl;
      
      // get new data chunk
      keep_looping=odata[c]->NewChunk();// get fresh data
     
      // check loop status
      if(keep_looping==1) return false; // critical error  --> STOP Omicron
      if(keep_looping==2){              // end of segments --> STOP loop
	keep_looping=-1; break;
      }
      if(keep_looping==3){              // end of FFL file --> STOP loop
	keep_looping=-1; break;
      }
      // increment counters
      chunk_ctr[c]++;         // one more chunk of data
      if(keep_looping==4){    // the data chunk is corrupted -> skip channel
	cor_chunk_ctr[c]++;
	continue;
      }
  
      // write chunk info if requested
      if(writetimeseries) odata[c]->WriteTimeSeries(fOutdir[c]);
      if(writepsd)        odata[c]->WritePSD(fOutdir[c]);

      //loop over segments
      for(int s=0; s<odata[c]->GetNSegments(); s++){

	// get conditioned data
	if(!odata[c]->GetConditionedData(s,&(c_data[0]),&(c_data[1]),&psd,psdsize)){
	  cerr<<"Omicron::Process: conditionned data are corrupted!"<<endl;
	  return false;
	}

	// set new power for this chunk
	if(!s){
	  if(fVerbosity) cout<<"set power to the tiles for this chunk..."<<endl;
	  tile->SetPowerSpectrum(psd,psdsize);
	}

	//get triggers
	cout<<" "<<fChannels[c]<<": Extracting triggers in "<<odata[c]->GetSegmentTimeStart(s)+fOverlapDuration/2<<"-"<<odata[c]->GetSegmentTimeStart(s)+fSegmentDuration-fOverlapDuration/2<<endl;
	if(!tile->GetTriggers(triggers[c],c_data[0],c_data[1],odata[c]->GetSegmentTimeStart(s))){
	  cerr<<"Omicron::Process: could not get triggers for channel "<<fChannels[c]
	      <<" in segment starting at "<<odata[c]->GetSegmentTimeStart(s)<<endl;
	  delete c_data[0];
	  delete c_data[1];
	  continue;
	}
	else
	  triggers[c]->AddSegment(odata[c]->GetSegmentTimeStart(s)+fOverlapDuration/2,odata[c]->GetSegmentTimeStart(s)+fSegmentDuration-fOverlapDuration/2);

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
	continue;
      }
      
      // save triggers
      if(!triggers[c]->Write("ALL", "default").compare("none")){
	cerr<<"Omicron::Process: writing events failed for channel "<<fChannels[c]<<endl;
	return false;
      }
      else
	outSegments[c]->AddSegment(odata[c]->GetChunkTimeStart()+fOverlapDuration/2,odata[c]->GetChunkTimeEnd()-fOverlapDuration/2);
	
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
  IO *fOptions = new IO(fOptionFile.c_str());
  IO& io = *fOptions;

  // output meta-data
  fOptionName.clear();
  fOptionType.clear();

  //***** output directory *****
  string maindir="";
  io.GetOpt("OUTPUT","DIRECTORY", maindir);
  if(maindir.empty()){
    cerr<<"Omicron::ReadOptions: you must provide an output directory: OUTPUT/DIRECTORY"<<endl;
    return false;
  }
  if(!IsDirectory(maindir)){
    cerr<<"Omicron::ReadOptions: output directory cannot be found: OUTPUT/DIRECTORY"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_OUTPUT_DIRECTORY");
  fOptionType.push_back("s");
  //*****************************

  //***** List of channels *****
  fChannels.clear();
  io.GetOpt("DATA","CHANNELS", fChannels);
  if(!fChannels.size()){
    cerr<<"Omicron::ReadOptions: you must provide a list of channels: DATA/CHANNELS"<<endl;
    return false;
  }
  if(fChannels.size()>NDATASTREAMS){
    cerr<<"Omicron::ReadOptions: The number of channels cannot exceed "<<NDATASTREAMS<<endl;
    return false;
  }
  fOptionName.push_back("omicron_DATA_CHANNEL");
  fOptionType.push_back("s");
  //*****************************

  //***** injection channels *****
  fInjChan.clear();
  io.GetOpt("INJECTION","CHANNELS", fInjChan);
  if(fInjChan.size()&&fInjChan.size()!=fChannels.size()){
    cerr<<"Omicron::ReadOptions: INJECTION/CHANNELS is inconsistent with the number of channels"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_INJECTION_CHANNEL");
  fOptionType.push_back("s");
  fInjFact.clear();
  io.GetOpt("INJECTION","FACTORS", fInjFact);
  if(fInjFact.size()!=fInjChan.size()){
    cerr<<"Omicron::ReadOptions: INJECTION/FACTORS is inconsistent with the number of channels"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_INJECTION_FACTOR");
  fOptionType.push_back("d");
  //*****************************

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
  //*****************************

  //***** lcf file *****
  online=false;
  bool lcfset=false;
  io.GetOpt("DATA","LCF", fLcfFile);
  if(!fLcfFile.compare("ONLINE")){// online keyword
    online=true;
    lcfset=true;
  }
  else if(fLcfFile.empty()){
    lcfset=false;
  }
  else if(!IsTextFile(fLcfFile)){
    cerr<<"Omicron::ReadOptions: the LCF file "<<fLcfFile<<" cannot be found"<<endl;
    return false;
  }
  else{// OK
    lcfset=true;

    // convert to FFL
    srand (time(NULL));
    int randint = rand();
    ostringstream tmpstream;
    tmpstream<<maindir<<"/converted_lcf_"<<randint<<".ffl";
    LCF2FFL(fLcfFile,tmpstream.str());
    fFflFile=tmpstream.str();
  }
  //*****************************

  //***** ffl file *****
  if(!lcfset){
    online=false;
    io.GetOpt("DATA","FFL", fFflFile);
    if(!fFflFile.compare("ONLINE"))// online keyword
      online=true;
    else if(fFflFile.empty()){
      cerr<<"Omicron::ReadOptions: you must provide a FFL file: DATA/FFL"<<endl;
      return false;
    }
    else if(!IsTextFile(fFflFile)){
      cerr<<"Omicron::ReadOptions: the FFL file "<<fFflFile<<" cannot be found"<<endl;
      return false;
    }
    else// OK
      online=false;
  }
  fOptionName.push_back("omicron_DATA_FFL");
  fOptionType.push_back("s");
  //*****************************
    
  //***** Native channel frequencies *****
  fNativeFrequency.clear();
  io.GetOpt("DATA","NATIVEFREQUENCY", fNativeFrequency);
  if(fNativeFrequency.size()!=fChannels.size()){
    cerr<<"Omicron::ReadOptions: Each channel must be given a native frequency: PARAMETER/NATIVEFREQUENCY"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_DATA_NATIVEFREQUENCY");
  fOptionType.push_back("i");
  //*****************************

  //***** Sampling frequency *****
  fSampleFrequency=-1;
  io.GetOpt("DATA","SAMPLEFREQUENCY", fSampleFrequency);
  if(fSampleFrequency<16||fSampleFrequency>40000){
    cerr<<"Omicron::ReadOptions: Sampling frequency "<<fSampleFrequency<<" (PARAMETER/SAMPLEFREQUENCY) is not reasonable"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_DATA_SAMPLEFREQUENCY");
  fOptionType.push_back("i");
  //*****************************

  //***** Frequency range *****
  fFreqRange.clear();
  io.GetOpt("PARAMETER","FREQUENCYRANGE", fFreqRange);
  if(fFreqRange.size()!=2||fFreqRange[0]>=fFreqRange[1]){
    cerr<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) is not correct"<<endl;
    return false;
  }
  if(fFreqRange[1]>fSampleFrequency/2){
    cout<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) goes beyond Nyquist frequency: "<<fFreqRange[1]<<">"<<fSampleFrequency/2<<" --> Nyquist frequency will be used"<<endl;
    fFreqRange.pop_back(); fFreqRange.push_back((double)fSampleFrequency/2.0);
  }
  if(fFreqRange[0]>=fFreqRange[1]){
    cerr<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) is not correct"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_PARAMETER_FMIN");
  fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_FMAX");
  fOptionType.push_back("d");
  //*****************************

  //***** Q range *****
  fQRange.clear();
  io.GetOpt("PARAMETER","QRANGE", fQRange);
  if(fQRange.size()!=2||fQRange[0]>fQRange[1]){
    cerr<<"Omicron::ReadOptions: Q range (PARAMETER/QRANGE) is not correct"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_PARAMETER_QMIN");
  fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_QMAX");
  fOptionType.push_back("d");
  //*****************************

  //***** timing *****
  fChunkDuration=484;
  fSegmentDuration=64;
  fOverlapDuration=8;
  io.GetOpt("PARAMETER","CHUNKDURATION", fChunkDuration);
  io.GetOpt("PARAMETER","BLOCKDURATION", fSegmentDuration);
  io.GetOpt("PARAMETER","OVERLAPDURATION", fOverlapDuration);
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
  fOptionName.push_back("omicron_PARAMETER_CHUNKDURATION");
  fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_BLOCKDURATION");
  fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_OVERLAPDURATION");
  fOptionType.push_back("i");
  //*****************************

  //***** maximum mismatch *****
  fMismatchMax=0.2;
  io.GetOpt("PARAMETER","MISMATCHMAX", fMismatchMax);
  if(fMismatchMax<=0||fMismatchMax>0.5){
    cerr<<"Omicron::ReadOptions: maximum mismatch (PARAMETER/MISMATCHMAX) is not reasonable"<<endl;
    return false;
  }
  fOptionName.push_back("omicron_PARAMETER_MISMATCHMAX");
  fOptionType.push_back("d");
  //*****************************

  //***** SNR Threshold *****
  fSNRThreshold=8.0;
  io.GetOpt("TRIGGER","SNRTHRESHOLD", fSNRThreshold);
  fOptionName.push_back("omicron_TRIGGER_SNRTHRESHOLD");
  fOptionType.push_back("d");
  //*****************************

  //***** set trigger limit *****
  fNtriggerMax=1000000;
  double ratemax;
  if(io.GetOpt("TRIGGER","RATEMAX", ratemax)) fNtriggerMax=(int)ceil(ratemax*fChunkDuration);
  else io.GetOpt("TRIGGER","NMAX", fNtriggerMax);
  fOptionName.push_back("omicron_TRIGGER_NMAX");
  fOptionType.push_back("i");
  //*****************************

  //***** set clustering *****
  fcldt=0.1;
  io.GetOpt("TRIGGER","CLUSTERING", fClusterAlgo);
  io.GetOpt("TRIGGER","CLUSTERDT", fcldt);
  fOptionName.push_back("omicron_TRIGGER_CLUSTERING");
  fOptionType.push_back("s");
  fOptionName.push_back("omicron_TRIGGER_CLUSTERDT");
  fOptionType.push_back("d");
  //*****************************

  //***** set verbosity *****
  fVerbosity=0;
  io.GetOpt("OUTPUT","VERBOSITY", fVerbosity);
  fOptionName.push_back("omicron_OUTPUT_VERBOSITY");
  fOptionType.push_back("i");
  //*****************************

  //***** set verbosity ***** 
  io.GetOpt("OUTPUT","FORMAT", fOutFormat);
  if(fOutFormat.compare("root")&&fOutFormat.compare("xml")&&fOutFormat.compare("txt")){
    cerr<<"Omicron::ReadOptions: possible output formats: root, xml or txt"<<endl;
    cerr<<"                      --> root format will be used"<<endl;
    fOutFormat="root";
  }
  fOptionName.push_back("omicron_OUTPUT_FORMAT");
  fOptionType.push_back("s");
  //*****************************

  //***** set writing flags *****
  writepsd=0, writetimeseries=0;
  io.GetOpt("OUTPUT","WRITEPSD", writepsd);
  io.GetOpt("OUTPUT","WRITETIMESERIES", writetimeseries);
  fOptionName.push_back("omicron_OUTPUT_WRITEPSD");
  fOptionType.push_back("i");
  fOptionName.push_back("omicron_OUTPUT_WRITETIMESERIES");
  fOptionType.push_back("i");
  //*****************************

  // create output directories
  for(int c=0; c<(int)fChannels.size(); c++){
    fOutdir[c]=maindir+"/"+fChannels[c];
    system(("mkdir -p "+fOutdir[c]).c_str());
  }

  // dump options
  if(fVerbosity>1) io.Dump(cout);
  delete fOptions;

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::LCF2FFL(const string lcf_file, const string ffl_file){
////////////////////////////////////////////////////////////////////////////////////

  // check that the option file exists
  if(!IsTextFile(lcf_file)){
    cerr<<"Omicron::LCF2FFL: input file "<<lcf_file<<" cannot be found"<<endl;
    return false;
  }

  ReadAscii *RA = new ReadAscii(lcf_file,"s;s;i;i;s");
  if(RA->GetNCol()!=5){
    cerr<<"Omicron::LCF2FFL: input file "<<lcf_file<<" does not look like a LCF file"<<endl;
    delete RA;
    return false;
  }

  ofstream ffile(ffl_file.c_str());
  int start, duration;
  string framefile, framefilename;
  for(int l=0; l<RA->GetNRow(); l++){
    if(!RA->GetElement(start,l,2)){
      cerr<<"Omicron::LCF2FFL: cannot retrieve starting time"<<endl;
      delete RA;
      return false;
    }
    if(!RA->GetElement(duration,l,3)){
      cerr<<"Omicron::LCF2FFL: cannot retrieve duration"<<endl;
      delete RA;
      return false;
    }
    if(!RA->GetElement(framefile,l,4)){
      cerr<<"Omicron::LCF2FFL: cannot retrieve frame file name"<<endl;
      delete RA;
      return false;
    }

    // FIXME: This could be better!!
    if(string::npos != framefile.find("localhost")){
      unsigned pos = framefile.find("localhost");
      framefilename=framefile.substr(pos+9);
    }
    else framefilename=framefile;

    ffile<<framefilename<<" "<<start<<" "<<duration<<" 0 0"<<endl;
  }

  ffile.close();
  delete RA;
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::PrintStatusInfo(void){
////////////////////////////////////////////////////////////////////////////////////

  
  cout<<"\n************* Omicron status info *************"<<endl;
  cout<<"start_in = "<<fSegments->GetStart(0)<<endl;
  cout<<"end_in   = "<<fSegments->GetEnd(fSegments->GetNsegments()-1)<<endl;
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

  return true;
}
