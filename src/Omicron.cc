//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

ClassImp(Omicron)

////////////////////////////////////////////////////////////////////////////////////
Omicron::Omicron(const string aOptionFile){ 
////////////////////////////////////////////////////////////////////////////////////

  cout<<"############################################################################"<<endl;
  cout<<"############################################################################"<<endl;
  cout<<endl;
  cout<<"          ██████╗ ███╗   ███╗██╗ ██████╗██████╗  ██████╗ ███╗   ██╗"<<endl;
  cout<<"         ██╔═══██╗████╗ ████║██║██╔════╝██╔══██╗██╔═══██╗████╗  ██║"<<endl;
  cout<<"         ██║   ██║██╔████╔██║██║██║     ██████╔╝██║   ██║██╔██╗ ██║"<<endl;
  cout<<"         ██║   ██║██║╚██╔╝██║██║██║     ██╔══██╗██║   ██║██║╚██╗██║"<<endl;
  cout<<"         ╚██████╔╝██║ ╚═╝ ██║██║╚██████╗██║  ██║╚██████╔╝██║ ╚████║"<<endl;
  cout<<"          ╚═════╝ ╚═╝     ╚═╝╚═╝ ╚═════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝"<<endl;
  cout<<endl;
  cout<<"############################################################################"<<endl;
  cout<<"############################################################################"<<endl;
 
  // input
  fOptionFile=aOptionFile;

  // metadata
  fOptionName.push_back("omicron_OUTPUT_DIRECTORY");          fOptionType.push_back("s");
  fOptionName.push_back("omicron_DATA_CHANNEL");              fOptionType.push_back("s");
  fOptionName.push_back("omicron_INJECTION_CHANNEL");         fOptionType.push_back("s");
  fOptionName.push_back("omicron_INJECTION_FACTOR");          fOptionType.push_back("d");
  fOptionName.push_back("omicron_DATA_FFLFILE");              fOptionType.push_back("s");
  fOptionName.push_back("omicron_DATA_SAMPLEFREQUENCY");      fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_FMIN");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_FMAX");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_QMIN");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_QMAX");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_CHUNKDURATION");   fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_BLOCKDURATION");   fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_OVERLAPDURATION"); fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_MISMATCHMAX");     fOptionType.push_back("d");
  fOptionName.push_back("omicron_TRIGGER_SNRTHRESHOLD");      fOptionType.push_back("d");
  fOptionName.push_back("omicron_TRIGGER_NMAX");              fOptionType.push_back("i");
  fOptionName.push_back("omicron_TRIGGER_CLUSTERING");        fOptionType.push_back("s");
  fOptionName.push_back("omicron_TRIGGER_CLUSTERDT");         fOptionType.push_back("d");
  fOptionName.push_back("omicron_OUTPUT_VERBOSITY");          fOptionType.push_back("i");
  fOptionName.push_back("omicron_OUTPUT_FORMAT");             fOptionType.push_back("s");
  fOptionName.push_back("omicron_OUTPUT_WRITEPSD");           fOptionType.push_back("i");
  fOptionName.push_back("omicron_OUTPUT_WRITETIMESERIES");    fOptionType.push_back("i");

  // read option file
  status_OK=ReadOptions();

  // init plotting
  GPlot = new GwollumPlot ("Omicron",fStyle);

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
  if(fVerbosity) cout<<"Omicron::Omicron: Allocate memory for data container and procedures..."<<endl;
  ChunkSize   = fSampleFrequency * fChunkDuration;
  OverlapSize = fSampleFrequency * fOverlapDuration;
  SegmentSize = fSampleFrequency * fSegmentDuration;
  ChunkVect   = new double [ChunkSize];
  SegVect     = new double [SegmentSize];
  TukeyWindow = GetTukeyWindow(SegmentSize,OverlapSize);
  offt = new fft(SegmentSize,"FFTW_"+ffftplan);
  dataseq = new Odata(fChunkDuration, fSegmentDuration, fOverlapDuration, fVerbosity);
  status_OK*=dataseq->GetStatus();
  dataRe = new double* [dataseq->GetNSegments()];
  dataIm = new double* [dataseq->GetNSegments()];

  // init FFL
  FFL=NULL;
  if(fFflFile.compare("none")){
    if(fVerbosity) cout<<"Omicron::Omicron: Define FFL object..."<<endl;
    FFL = new ffl(fFflFile, fVerbosity);
    status_OK*=FFL->DefineTmpDir(fMaindir);
    status_OK*=FFL->LoadFrameFile();
  }
  
  // init Streams
  if(fVerbosity) cout<<"Omicron::Omicron: Define Streams object..."<<endl;
  streams = new Streams* [(int)fChannels.size()];
  first_save = new bool [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){
    streams[c] = new Streams(fChannels[c], fVerbosity);
    status_OK*=streams[c]->GetStatus();
    streams[c]->MakeLVDetector();// just an attempt
    first_save[c]=true;
  }
  
  // init Sample
  sample = new Sample* [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++) sample[c] = new Sample(0);
  
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
    
  // init Triggers
  if(fVerbosity) cout<<"Omicron::Omicron: Define Triggers object..."<<endl;
  triggers    = new Triggers* [(int)fChannels.size()];
  string form = "";
  if(fOutFormat.find("xml")!=string::npos)  form+="xml";
  if(fOutFormat.find("txt")!=string::npos)  form+="txt";
  if(fOutFormat.find("root")!=string::npos) form+="root";
  if(!form.compare("")) form="root";
  for(int c=0; c<(int)fChannels.size(); c++){
    
    // create trigger directory
    system(("mkdir -p "+fOutdir[c]).c_str());
    
    // init triggers object
    triggers[c] = new Triggers(fOutdir[c],fChannels[c],form,fVerbosity);
    triggers[c]->SetNtriggerMax(fNtriggerMax);// maximum number of triggers per file
    
    // set clustering if any
    if(fClusterAlgo[0].compare("none")){// clustering is requested
      status_OK*=triggers[c]->SetClustering(fClusterAlgo[0]);// set clustering
      status_OK*=triggers[c]->SetClusterDeltaT(fcldt);// set dt
    }
    
    // set metadata
    status_OK*=triggers[c]->InitUserMetaData(fOptionName,fOptionType);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],fOutdir[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[1],fChannels[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[2],fInjChan[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],fInjFact[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[4],fFflFile);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[5],fSampleFrequency);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[6],fFreqRange[0]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[7],fFreqRange[1]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[8],fQRange[0]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[9],fQRange[1]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[10],fChunkDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[11],fSegmentDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[12],fOverlapDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[13],fMismatchMax);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[14],fSNRThreshold);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[15],fNtriggerMax);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[16],fClusterAlgo[0]+"_"+fClusterAlgo[1]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[17],fcldt);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[18],fVerbosity);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[19],fOutFormat);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[20],(int)writepsd);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[21],(int)writetimeseries);
    triggers[c]->SetMprocessname("Omicron");
    triggers[c]->SetMstreamname(streams[c]->GetName());
    triggers[c]->SetMdetindex(streams[c]->GetDetIndex());
  }
  
  // init tiles
  if(fVerbosity) cout<<"Omicron::Omicron: Define Tile object..."<<endl;
  tile = new Otile(fSegmentDuration,fOverlapDuration/2,fQRange[0],fQRange[1],fFreqRange[0],fFreqRange[1],fSampleFrequency,fMismatchMax,fSNRThreshold,fVerbosity);
  status_OK*=tile->GetStatus();

  // init Qmaps
  if(fVerbosity) cout<<"Omicron::Omicron: Define Qmaps..."<<endl;
  Qmap = new TH2D* [tile->GetNQPlanes()];
  for(int q=0; q<tile->GetNQPlanes(); q++) Qmap[q]=NULL;
  Qmap_center=-1.0;

  if(!status_OK) cerr<<"Omicron::Omicron: Initialization failed"<<endl;
  cout<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
Omicron::~Omicron(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"delete Omicron"<<endl;
  delete inSegments;
  delete chunk_ctr;
  delete cor_chunk_ctr;
  delete cor_data_ctr;
  delete max_chunk_ctr;
  for(int c=0; c<(int)fChannels.size(); c++){
    delete outSegments[c];
    delete sample[c];
    delete streams[c];
    delete triggers[c];
  }
  delete outSegments;
  delete dataRe;
  delete dataIm;
  delete sample;
  delete streams;
  delete spectrum;
  delete triggers;
  if(FFL!=NULL) delete FFL;
  delete dataseq;
  delete tile;
  delete Qmap;
  delete GPlot;
  delete ChunkVect;
  delete SegVect;
  delete TukeyWindow;
  delete offt;
  delete first_save;

  fOptionName.clear();
  fOptionType.clear();
  fOutdir.clear();
  fChannels.clear();
  fInjChan.clear();
  fInjFact.clear();
  fFreqRange.clear();
  fQRange.clear();
  fWindows.clear();
  fClusterAlgo.clear();
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
  if(FFL==NULL){
    cerr<<"Omicron::Process: this function can only be used with a valid FFL object"<<endl;
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
  if(!dataseq->SetSegments(aSeg)){
    cerr<<"Omicron::Process: cannot initiate data segments."<<endl;
    return false;
  }

  // locals
  int dsize;         // native data size
  double *dvector;   // data vector before resampling
  int res;
  
  // loop over chunks
  while(dataseq->NewChunk()){
    if(fVerbosity) cout<<"Omicron::Process: chunk "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<endl;
    
    // loop over channels
    for(int c=0; c<(int)fChannels.size(); c++){
      if(fVerbosity) cout<<"Omicron::Process: channel "<<fChannels[c]<<"..."<<endl;
      
      // one more chunk of data
      chunk_ctr[c]++;

      // get data vector
      dvector=FFL->GetData(dsize, fChannels[c], dataseq->GetChunkTimeStart(), dataseq->GetChunkTimeEnd());

      if(dsize<=0){
	cor_chunk_ctr[c]++;
	cerr<<"Omicron::Process: cannot retrieve data ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
	continue;
      }

      // size update
      fChunkDuration=dataseq->GetChunkTimeEnd()-dataseq->GetChunkTimeStart();
      ChunkSize=fChunkDuration*fSampleFrequency;

      // process data vector
      res=ConditionVector(c, dsize, dvector);
      if(res<0){
	cerr<<"Omicron::Process: conditioning failed."<<endl;
	delete dvector;// not needed anymore
	cor_data_ctr[c]++;
	return false;
      }
      if(res>0){
	cerr<<"Omicron::Process: conditioning failed."<<endl;
	delete dvector;// not needed anymore
	cor_data_ctr[c]++;
	continue;
      }

      delete dvector;// not needed anymore

      // Save info if requested
      if(writetimeseries) SaveData(c,ChunkVect,dataseq->GetChunkTimeStart(),dataseq->GetChunkTimeEnd());
      if(writepsd)        SaveAPSD(c,"PSD");

      // get triggers above SNR threshold
      if(MakeTriggers(c)<0){
	cerr<<"Omicron::Process: cannot make triggers for channel "<<fChannels[c]<<endl;
	cor_data_ctr[c]++;
	continue;
      }
     
      // don't save if max flag
      if(triggers[c]->GetMaxFlag()){
	cerr<<"Omicron::Process: channel "<<fChannels[c]<<" is maxed-out. This chunk of triggers is not saved"<<endl;
	triggers[c]->Reset();
	max_chunk_ctr[c]++;
	continue; // move to next channel
      }
      
      // save triggers for this chunk
      if(!triggers[c]->Write(fWriteMode, "default").compare("none")){
	cerr<<"Omicron::Process: writing events failed for channel "<<fChannels[c]<<endl;
	return false;
      }

      // processed chunks
      else
	outSegments[c]->AddSegment(dataseq->GetChunkTimeStart()+fOverlapDuration/2,dataseq->GetChunkTimeEnd()-fOverlapDuration/2);
	
    }
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::Scan(const double aTimeCenter){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::Scan: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aTimeCenter<700000000){
    cerr<<"Omicron::Scan: the input time is not reasonable"<<endl;
    return false;
  }
  if(fChunkDuration!=fSegmentDuration){
    cerr<<"Omicron::Scan: this function can only be called if the chunk and segment durations are identical"<<endl;
    return false;
  }

  // locals
  int dsize;         // native data size
  double *dvector;   // data vector before resampling
  int res;

  // add requested segments
  if(fVerbosity) cout<<"Omicron::Scan: initiate data segments..."<<endl;
  inSegments->Reset();
  if(!inSegments->AddSegment((int)aTimeCenter-fChunkDuration/2,(int)aTimeCenter+fChunkDuration/2)){
    cerr<<"Omicron::Scan: the input segment is corrupted"<<endl;
    return false;
  }

  // data structure
  if(!dataseq->SetSegments(inSegments)||!dataseq->NewChunk()){
    cerr<<"Omicron::Scan: cannot initiate data segments."<<endl;
    return false;
  }

  // loop over channels
  for(int c=0; c<(int)fChannels.size(); c++){
    if(fVerbosity) cout<<"Omicron::Scan: channel "<<fChannels[c]<<"..."<<endl;
      
    // get data vector
    dvector=FFL->GetData(dsize, fChannels[c], dataseq->GetChunkTimeStart(), dataseq->GetChunkTimeEnd());

    if(dsize<=0){
      cerr<<"Omicron::Scan: cannot retrieve data ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")."<<endl;
      continue;
    }

    // process data vector
    res=ConditionVector(c, dsize, dvector);
    if(res<0){
      cerr<<"Omicron::Scan: conditioning failed."<<endl;
      delete dvector;
      return false;
    }
    if(res>0){
      cerr<<"Omicron::Scan: conditioning failed."<<endl;
      delete dvector;
      continue;
    }

    delete dvector;// not needed anymore

    // get maps
    if(!MakeMaps(c, aTimeCenter)){
      cerr<<"Omicron::Scan: cannot make maps for channel "<<fChannels[c]<<endl;
      continue;
    }

    // save maps on disk
    if(!WriteMaps(c)){
      cerr<<"Omicron::Scan: cannot write maps for channel "<<fChannels[c]<<endl;
      continue;
    }

    // save PSD on disk
    SaveAPSD(c,"ASD");

   
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
int Omicron::ConditionVector(const int aChNumber, const int aInVectSize, double *aInVect){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::ConditionVector: the Omicron object is corrupted"<<endl;
    return -1;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::ConditionVector: channel number "<<aChNumber<<" does not exist"<<endl;
    return 1;
  }

  // Input data checks
  if(aInVect==NULL){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): input vector is NULL"<<endl;
    return 2;
  }
  if(aInVectSize<=0){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): input vector empty"<<endl;
    return 3;
  }
  if(aInVect[0]==aInVect[aInVectSize-1]){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): input vector is flat"<<endl;
    return 4;
  }

  // transform data vector
  int nativesampling = aInVectSize/fChunkDuration;
  if(!sample[aChNumber]->SetFrequencies(nativesampling,fSampleFrequency,fFreqRange[0])){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): the Sample structure could not be defined"<<endl;
    return -2;
  }
  if(!sample[aChNumber]->Transform(aInVectSize, aInVect, ChunkSize, ChunkVect)){
    return 5;
  }

  // Make spectrum
  if(!spectrum->LoadData(ChunkSize, ChunkVect)){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): the spectrum could not be estimated"<<endl;
    return 6;
  }

  // set new power for this chunk
  if(!tile->SetPowerSpectrum(spectrum)){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): the tiling could not be normalized"<<endl;
    return 7;
  }

  //loop over segments
  for(int s=0; s<dataseq->GetNSegments(); s++){
    	
    // fill segment vector (time-domain) and apply tukey window
    for(int i=0; i<SegmentSize; i++)
      SegVect[i] = ChunkVect[s*(SegmentSize-OverlapSize)+i] * TukeyWindow[i];
    
    // get conditioned data
    if(!Condition(&(dataRe[s]), &(dataIm[s]))){
      cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): the conditioning failed"<<endl;
      return 8;
    }
  }
  
  return 0;
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
    if(asdval<=0){
      cerr<<"Omicron::Condition: could not retrieve power"<<endl;
      delete *aDataRe; delete *aDataIm;
      *aDataRe=NULL; *aDataIm=NULL;
      return false;
    }
    asdval=sqrt(asdval);
    (*aDataRe)[i] /= asdval;
    (*aDataIm)[i] /= asdval;
    // NOTE: no FFT normalization here because we do a FFTback in Oqplane later with no normalization either.
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
int Omicron::MakeTriggers(const int aChNumber){
  ////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::MakeTriggers: the Omicron object is corrupted"<<endl;
    return -1;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::MakeTriggers: channel number "<<aChNumber<<" does not exist"<<endl;
    return -1;
  }
  if(fVerbosity) cout<<"Omicron::MakeTriggers: make triggers for channel "<<fChannels[aChNumber]<<" starting at "<<dataseq->GetChunkTimeStart()<<endl;
  
  int s;
  for(s=0; s<dataseq->GetNSegments(); s++){
    if(!tile->GetTriggers(triggers[aChNumber],dataRe[s],dataIm[s],dataseq->GetSegmentTimeStart(s))){
      cerr<<"Omicron::MakeTriggers: could not make triggers for channel "<<fChannels[aChNumber]<<endl;
      return -1;
    }
    else{
      triggers[aChNumber]->AddSegment(dataseq->GetSegmentTimeStart(s)+fOverlapDuration/2,dataseq->GetSegmentTimeEnd(s)-fOverlapDuration/2);
      delete dataRe[s];
      delete dataIm[s];
      if(triggers[aChNumber]->GetMaxFlag()) break;
    }
  }
  for(int ss=s+1; ss<dataseq->GetNSegments(); ss++){
    delete dataRe[s];
    delete dataIm[s];
  } 
   
  return triggers[aChNumber]->GetNTrig();

}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::WriteTriggers(const int aChNumber){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::WriteTriggers: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::WriteTriggers: channel number "<<aChNumber<<" does not exist"<<endl;
    return false;
  }
  if(fVerbosity) cout<<"Omicron::WriteTriggers: write triggers for channel "<<fChannels[aChNumber]<<endl;
  return !(triggers[aChNumber]->Write("PROC","default").compare("none"));
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::MakeMaps(const int aChNumber, const double aTimeCenter){
  ////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::MakeMaps: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::MakeMaps: channel number "<<aChNumber<<" does not exist"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::MakeMaps: make maps for channel "<<fChannels[aChNumber]<<" starting at "<<dataseq->GetChunkTimeStart()<<endl;

  Qmap_center=aTimeCenter;
  double toffset=Qmap_center-dataseq->GetSegmentTimeStart(0)-(double)fSegmentDuration/2.0;
  ostringstream tmpstream;

  // get maps
  for(int q=0; q<tile->GetNQPlanes(); q++){
    Qmap[q] = tile->GetMap(q, dataRe[0], dataIm[0], -toffset);
    if(Qmap[q]==NULL){
      cerr<<"Omicron::MakeMaps: maps for channel number "<<aChNumber<<" are corrupted"<<endl;
      delete dataRe[0];
      delete dataIm[0];
      return false;
    }

    // map title
    tmpstream<<fChannels[aChNumber]<<": GPS="<<setprecision(12)<<Qmap_center<<", Q="<<setprecision(5)<<tile->GetQ(q);
    Qmap[q]->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
    
    // cosmetics
    Qmap[q]->GetXaxis()->SetTitle("Time [s]");
    Qmap[q]->GetYaxis()->SetTitle("Frequency [Hz]");
    Qmap[q]->GetZaxis()->SetTitle("SNR [-]");
    Qmap[q]->GetXaxis()->SetTitleOffset(1.1);
    Qmap[q]->GetXaxis()->SetLabelSize(0.045);
    Qmap[q]->GetYaxis()->SetLabelSize(0.045);
    Qmap[q]->GetXaxis()->SetTitleSize(0.045);
    Qmap[q]->GetYaxis()->SetTitleSize(0.045);
    Qmap[q]->GetZaxis()->SetTitleSize(0.05);
  }
  delete dataRe[0];
  delete dataIm[0];
      
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::WriteMaps(const int aChNumber){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::WriteMaps: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::WriteMaps: channel number "<<aChNumber<<" does not exist"<<endl;
    return false;
  }
  for(int q=0; q<tile->GetNQPlanes(); q++){
    if(Qmap[q]==NULL){
      cerr<<"Omicron::WriteMaps: Q-maps are not been filled"<<endl;
      return false;
    }
  }

  if(fVerbosity) cout<<"Omicron::WriteMaps: write maps for channel "<<fChannels[aChNumber]<<endl;
  
  // canvas style
  GPlot->SetLogx(0);
  GPlot->SetLogy(1);
  GPlot->SetLogz(1);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);

  // locals
  double center;
  ostringstream tmpstream;
  int xmax, ymax, zmax, bin_start, bin_stop, dummy, bin_start_t, bin_stop_t, bin_stop_f, bin_start_f;
  double content;

  // create overall maps
  TH2D **map = new TH2D * [(int)fWindows.size()];
  int nfbins;
  if(fFreqRange[0]>=1) nfbins = fabs(fFreqRange[1]-fFreqRange[0])+1;
  else nfbins = 1000;
  double *f_bin = new double [nfbins+1];
  for(int f=0; f<nfbins+1; f++) f_bin[f]=fFreqRange[0]*pow(10,f*log10(fFreqRange[1]/fFreqRange[0])/nfbins);
  for(int w=0; w<(int)fWindows.size(); w++){
    tmpstream<<"map_"<<fWindows[w];
    map[w] = new TH2D(tmpstream.str().c_str(),tmpstream.str().c_str(),1000,-(double)fWindows[w]/2.0,+(double)fWindows[w]/2.0,nfbins,f_bin);
    tmpstream.str(""); tmpstream.clear();
    tmpstream<<fChannels[aChNumber]<<": GPS="<<setprecision(12)<<Qmap_center;
    map[w]->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
    map[w]->GetXaxis()->SetTitle("Time [s]");
    map[w]->GetYaxis()->SetTitle("Frequency [Hz]");
    map[w]->GetZaxis()->SetTitle("SNR [-]");
    map[w]->GetXaxis()->SetTitleOffset(1.1);
    map[w]->GetXaxis()->SetLabelSize(0.045);
    map[w]->GetYaxis()->SetLabelSize(0.045);
    map[w]->GetXaxis()->SetTitleSize(0.045);
    map[w]->GetYaxis()->SetTitleSize(0.045);
    map[w]->GetZaxis()->SetTitleSize(0.05);
    map[w]->GetZaxis()->SetRangeUser(1,50);
  }
  delete f_bin;
    
  //populate overall maps
  for(int q=0; q<tile->GetNQPlanes(); q++){
    for(int bt=1; bt<=Qmap[q]->GetNbinsX(); bt++){
      for(int bf=1; bf<=Qmap[q]->GetNbinsY(); bf++){
	content=Qmap[q]->GetBinContent(bt,bf);
	for(int w=0; w<(int)fWindows.size(); w++){
	  if(Qmap[q]->GetXaxis()->GetBinLowEdge(bt)<-(double)fWindows[w]/2.0) continue;
	  if(Qmap[q]->GetXaxis()->GetBinUpEdge(bt)>(double)fWindows[w]/2.0) continue;
	  if(Qmap[q]->GetYaxis()->GetBinLowEdge(bf)<map[w]->GetYaxis()->GetBinLowEdge(1)) continue;
	  if(Qmap[q]->GetYaxis()->GetBinUpEdge(bf)>map[w]->GetYaxis()->GetBinUpEdge(nfbins)) continue;
	  bin_start = map[w]->FindBin(Qmap[q]->GetXaxis()->GetBinLowEdge(bt),Qmap[q]->GetYaxis()->GetBinLowEdge(bf));
	  bin_stop = map[w]->FindBin(Qmap[q]->GetXaxis()->GetBinUpEdge(bt),Qmap[q]->GetYaxis()->GetBinUpEdge(bf));
	  map[w]->GetBinXYZ(bin_start,bin_start_t,bin_start_f,dummy);
	  map[w]->GetBinXYZ(bin_stop,bin_stop_t,bin_stop_f,dummy);
	  for(int bbt=bin_start_t; bbt<=bin_stop_t; bbt++){// time-sweep the tile
	    for(int bbf=bin_start_f; bbf<=bin_stop_f; bbf++){// freq-sweep the tile
	      if(content>map[w]->GetBinContent(bbt,bbf)) map[w]->SetBinContent(bbt,bbf,content);
	    }
	  }
	}
      }
    }
  }


  // draw Qmaps
  for(int q=0; q<tile->GetNQPlanes(); q++){

    // SNR range
    Qmap[q]->GetZaxis()->SetRangeUser(1,50);
    
    // plot
    GPlot->Draw(Qmap[q],"COLZ");
    
    // window resize
    center=(Qmap[q]->GetXaxis()->GetXmax()+Qmap[q]->GetXaxis()->GetXmin())/2.0;
    for(int w=0; w<(int)fWindows.size(); w++){
      
      // zoom
      Qmap[q]->GetXaxis()->SetRangeUser(center-(double)fWindows[w]/2.0,center+(double)fWindows[w]/2.0);

      // get max bin
      Qmap[q]->GetMaximumBin(xmax, ymax, zmax);
         
      // loudest tile
      tmpstream<<"Loudest tile: GPS="<<setprecision(12)<<Qmap_center+Qmap[q]->GetXaxis()->GetBinCenter(xmax)<<setprecision(5)<<", f="<<Qmap[q]->GetYaxis()->GetBinCenter(ymax)<<" Hz, SNR="<<Qmap[q]->GetBinContent(xmax,ymax);
      GPlot->AddText(tmpstream.str(), 0.01,0.01,0.03);
      tmpstream.str(""); tmpstream.clear();

      // save qmaps
      tmpstream<<fOutdir[aChNumber]<<"/"<<fChannels[aChNumber]<<"_mapQ"<<q<<"_dt"<<fWindows[w]<<".gif";
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();
    }

  }

  // draw overall maps
  for(int w=0; w<(int)fWindows.size(); w++){
    GPlot->Draw(map[w],"COLZ");

    // get max bin
    map[w]->GetMaximumBin(xmax, ymax, zmax);
         
    // loudest tile
    tmpstream<<"Loudest tile: GPS="<<setprecision(12)<<Qmap_center+map[w]->GetXaxis()->GetBinCenter(xmax)<<setprecision(5)<<", f="<<map[w]->GetYaxis()->GetBinCenter(ymax)<<" Hz, SNR="<<map[w]->GetBinContent(xmax,ymax);
    GPlot->AddText(tmpstream.str(), 0.01,0.01,0.03);
    tmpstream.str(""); tmpstream.clear();

    // save maps
    tmpstream<<fOutdir[aChNumber]<<"/"<<fChannels[aChNumber]<<"_map_dt"<<fWindows[w]<<".gif";
    GPlot->Print(tmpstream.str());
    tmpstream.str(""); tmpstream.clear();
  }

  for(int w=0; w<(int)fWindows.size(); w++) delete map[w];
  delete map;
  return true;
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
void Omicron::SaveAPSD(const int c, const string type){
////////////////////////////////////////////////////////////////////////////////////

  // extract A/PSD 
  TGraph *GAPSD;
  if(!type.compare("ASD")) GAPSD = spectrum->GetASD();
  else GAPSD = spectrum->GetPSD();
  if(GAPSD==NULL) return;

  // cosmetics
  GPlot->SetLogx(1);
  GPlot->SetLogy(1);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);
  GPlot->Draw(GAPSD,"APL");
  GAPSD->GetHistogram()->SetXTitle("Frequency [Hz]");
  if(type.compare("ASD")){
    GAPSD->GetHistogram()->SetYTitle("Power [Amp^{2}/Hz]");
    GAPSD->SetTitle((fChannels[c]+": Power spectrum density").c_str());
  }
  else{
    GAPSD->GetHistogram()->SetYTitle("Amplitude [Amp/#sqrt{Hz}]");
    GAPSD->SetTitle((fChannels[c]+": Amplitude spectrum density").c_str());
  }
  GAPSD->SetLineWidth(2);
  GAPSD->GetXaxis()->SetLimits(fFreqRange[0],fFreqRange[1]);
  GAPSD->GetXaxis()->SetTitleOffset(1.1);
  GAPSD->GetXaxis()->SetLabelSize(0.045);
  GAPSD->GetYaxis()->SetLabelSize(0.045);
  GAPSD->GetXaxis()->SetTitleSize(0.045);
  GAPSD->GetYaxis()->SetTitleSize(0.045);
   
  // set new name
  stringstream ss;
  if(!type.compare("ASD")) ss<<"ASD_"<<fChannels[c]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd();
  else ss<<"PSD_"<<fChannels[c]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd();
  GAPSD->SetName(ss.str().c_str());
  ss.str(""); ss.clear();

  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fpsd;
    if(first_save[c]){
      fpsd=new TFile((fOutdir[c]+"/data_"+fChannels[c]+".root").c_str(),"RECREATE");
      first_save[c]=false;
    }
    else fpsd=new TFile((fOutdir[c]+"/data_"+fChannels[c]+".root").c_str(),"UPDATE");
    fpsd->cd();
    GAPSD->Write();
    fpsd->Close();
  }

  // Graphix
  vector <string> form;
  if(fOutFormat.find("gif")!=string::npos) form.push_back("gif");
  if(fOutFormat.find("png")!=string::npos) form.push_back("png");
  if(fOutFormat.find("pdf")!=string::npos) form.push_back("pdf");
  if(fOutFormat.find("ps")!=string::npos)  form.push_back("ps");
  if(fOutFormat.find("xml")!=string::npos) form.push_back("xml");
  if(fOutFormat.find("eps")!=string::npos) form.push_back("eps"); 
  if(form.size()){
    for(int f=0; f<(int)form.size(); f++){
      if(!type.compare("ASD")) ss<<fOutdir[c]<<"/asd_"<<fChannels[c]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"."<<form[f];
      else ss<<fOutdir[c]<<"/psd_"<<fChannels[c]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"."<<form[f];
      GPlot->Print(ss.str().c_str());
      ss.str(""); ss.clear();
    }
  }
  
  form.clear();
  delete GAPSD;
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
  if(first_save[c]){
    fdata=new TFile((fOutdir[c]+"/data_"+fChannels[c]+".root").c_str(),"RECREATE");
    first_save[c]=false;
  }
  else fdata=new TFile((fOutdir[c]+"/data_"+fChannels[c]+".root").c_str(),"UPDATE");
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
void Omicron::PrintStatusInfo(void){
////////////////////////////////////////////////////////////////////////////////////

  
  cout<<"\n************* Omicron status info *************"<<endl;
  cout<<"requested start = "<<inSegments->GetStart(0)<<endl;
  cout<<"requested end   = "<<inSegments->GetEnd(inSegments->GetNsegments()-1)<<endl;
  cout<<"requested livetime   = "<<inSegments->GetLiveTime()<<endl;
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
