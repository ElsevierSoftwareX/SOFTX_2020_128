//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

ClassImp(Omicron)

////////////////////////////////////////////////////////////////////////////////////
Omicron::Omicron(const string aOptionFile){ 
////////////////////////////////////////////////////////////////////////////////////
  PrintASCIIlogo();
  gErrorIgnoreLevel = 3000;

  // init timer
  time ( &timer );
  timer_start=timer;

  // input
  fOptionFile=aOptionFile;

  // metadata field definition
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
  fOptionName.push_back("omicron_OUTPUT_PRODUCTS");           fOptionType.push_back("s");

  // read option file
  status_OK=ReadOptions();

  // init plotting
  GPlot = new GwollumPlot ("Omicron",fStyle);

  // init process monitoring
  if(fVerbosity) cout<<"Omicron::Omicron: init monitoring..."<<endl;
  chanindex      = -1;
  timeoffset     = 0.0;
  inSegments     = new Segments();
  outSegments    = new Segments* [(int)fChannels.size()];
  chunk_ctr      = 0;
  chan_ctr       = new int       [(int)fChannels.size()];
  chan_data_ctr  = new int       [(int)fChannels.size()];
  chan_cond_ctr  = new int       [(int)fChannels.size()];
  chan_wmaps_ctr = new int       [(int)fChannels.size()];
  chan_wtrig_ctr = new int       [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){
    outSegments[c]    = new Segments();
    chan_ctr[c]       = 0;
    chan_data_ctr[c]  = 0;
    chan_cond_ctr[c]  = 0;
    chan_wmaps_ctr[c] = 0;
    chan_wtrig_ctr[c] = 0;
  }
  
  // init FFL
  FFL=NULL;
  if(fFflFile.compare("none")){
    if(fVerbosity) cout<<"Omicron::Omicron: init FFL..."<<endl;
    FFL = new ffl(fFflFile, fStyle, fVerbosity);
    status_OK*=FFL->DefineTmpDir(fMaindir);
    status_OK*=FFL->LoadFrameFile();
  }
    
  // adjust default parameters using input data
  AdjustParameters();

  // init data containers
  if(fVerbosity) cout<<"Omicron::Omicron: init data container and procedures..."<<endl;
  dataseq = new Odata(fChunkDuration, fSegmentDuration, fOverlapDuration, fVerbosity);
  fChunkDuration=dataseq->GetChunkDuration();
  fSegmentDuration=dataseq->GetSegmentDuration();
  fOverlapDuration=dataseq->GetOverlapDuration();
  ChunkVect   = new double [fChunkDuration*fSampleFrequency];
  SegVect     = new double [fSegmentDuration*fSampleFrequency];
  TukeyWindow = GetTukeyWindow(fSegmentDuration*fSampleFrequency,fOverlapDuration*fSampleFrequency);
  offt = new fft(fSegmentDuration*fSampleFrequency,"FFTW_"+ffftplan);
  dataRe = new double* [dataseq->GetNSegments()];
  dataIm = new double* [dataseq->GetNSegments()];
  
  // init Streams
  if(fVerbosity) cout<<"Omicron::Omicron: init Streams..."<<endl;
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
  if(fVerbosity) cout<<"Omicron::Omicron: init Spectra..."<<endl;
  int psdsize;
  if(fFreqRange[0]>=1.0) psdsize = fSampleFrequency; // 0.5Hz binning
  else{
    psdsize = 10*(int)floor((double)fSampleFrequency/fFreqRange[0]);// over-binning (factor 20)
    int nextpowerof2=(int)floor(log(psdsize)/log(2));
    psdsize=(int)pow(2.0,(double)nextpowerof2);
  }
  spectrum = new Spectrum(fSampleFrequency,psdsize,0,fVerbosity);
  status_OK*=spectrum->GetStatus();

  // default output directory
  for(int c=0; c<(int)fChannels.size(); c++) fOutdir.push_back(fMaindir);
  fScandir="";// used for scans

  // init Triggers
  if(fVerbosity) cout<<"Omicron::Omicron: init Triggers..."<<endl;
  triggers    = new MakeTriggers* [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){
    
    // init triggers object
    triggers[c] = new MakeTriggers(fOutdir[c],fChannels[c],fOutFormat,fVerbosity);
    triggers[c]->SetNtriggerMax(fNtriggerMax);// maximum number of triggers per file
    
    // set clustering parameters
    status_OK*=triggers[c]->SetClusterDeltaT(fcldt);
        
    // set metadata
    status_OK*=triggers[c]->InitUserMetaData(fOptionName,fOptionType);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],fMaindir);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[1],fChannels[c]);
    if(fInjChan.size()) status_OK*=triggers[c]->SetUserMetaData(fOptionName[2],fInjChan[c]);
    else status_OK*=triggers[c]->SetUserMetaData(fOptionName[2],"none");
    if(fInjChan.size()) status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],fInjFact[c]);
    else status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],0.0);
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
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[20],fOutProducts);
    triggers[c]->SetMprocessname("Omicron");
    triggers[c]->SetMstreamname(streams[c]->GetName());
    triggers[c]->SetMdetindex(streams[c]->GetDetIndex());
  }
  
  // init tiles
  if(fVerbosity) cout<<"Omicron::Omicron: init Tile..."<<endl;
  tile = new Otile(fSegmentDuration,fOverlapDuration/2,fQRange[0],fQRange[1],fFreqRange[0],fFreqRange[1],fSampleFrequency,fMismatchMax,fSNRThreshold,fVerbosity);
  status_OK*=tile->GetStatus();

  // init Qmaps
  if(fVerbosity) cout<<"Omicron::Omicron: init Qmaps..."<<endl;
  Qmap = new TH2D* [tile->GetNQPlanes()];
  for(int q=0; q<tile->GetNQPlanes(); q++) Qmap[q]=NULL;
  
  // init full maps
  Qmap_full = new TH2D* [(int)fWindows.size()];
  int nfbins; ostringstream tmpstream;
  if(fFreqRange[0]>=1) nfbins = (int)fabs(fFreqRange[1]-fFreqRange[0])+1;
  else nfbins = 500;
  double *f_bin = new double [nfbins+1];
  for(int f=0; f<nfbins+1; f++) f_bin[f]=fFreqRange[0]*pow(10,f*log10(fFreqRange[1]/fFreqRange[0])/nfbins);
  for(int w=0; w<(int)fWindows.size(); w++){
    tmpstream<<"map_"<<fWindows[w];
    Qmap_full[w] = new TH2D(tmpstream.str().c_str(),tmpstream.str().c_str(),1000,-(double)fWindows[w]/2.0,+(double)fWindows[w]/2.0,nfbins,f_bin);
    tmpstream.str(""); tmpstream.clear();
    Qmap_full[w]->GetXaxis()->SetTitle("Time [s]");
    Qmap_full[w]->GetYaxis()->SetTitle("Frequency [Hz]");
    Qmap_full[w]->GetZaxis()->SetTitle("SNR");
    Qmap_full[w]->GetXaxis()->SetTitleOffset(1.1);
    Qmap_full[w]->GetXaxis()->SetLabelSize(0.045);
    Qmap_full[w]->GetYaxis()->SetLabelSize(0.045);
    Qmap_full[w]->GetXaxis()->SetTitleSize(0.045);
    Qmap_full[w]->GetYaxis()->SetTitleSize(0.045);
    Qmap_full[w]->GetZaxis()->SetTitleSize(0.05);
  }
  delete f_bin;
  loudest_qmap = new int [(int)fWindows.size()];
}

////////////////////////////////////////////////////////////////////////////////////
Omicron::~Omicron(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"Omicron::~Omicron"<<endl;
  delete inSegments;
  delete chan_ctr;
  delete chan_data_ctr;
  delete chan_cond_ctr;
  delete chan_wmaps_ctr;
  delete chan_wtrig_ctr;
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
  for(int w=0; w<(int)fWindows.size(); w++) delete Qmap_full[w];
  delete Qmap_full;
  delete loudest_qmap;
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
bool Omicron::InitSegments(Segments *aSeg, const double aTimeOffset){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::InitSegments: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aSeg==NULL){
    cerr<<"Omicron::InitSegments: the input segment is NULL"<<endl;
    return false;
  }
  if(!aSeg->GetStatus()){
    cerr<<"Omicron::InitSegments: the input segment is corrupted"<<endl;
    return false;
  }
  if(!aSeg->GetNsegments()){
    cerr<<"Omicron::InitSegments: there is no input segment"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::InitSegments: initiate data segments..."<<endl;

  // update time offset
  timeoffset = aTimeOffset;

  // input segment monitor
  // (try to optimize the update of inSegments)
  if(inSegments->GetLiveTime()&&aSeg->GetStart(0)>=inSegments->GetStart(inSegments->GetNsegments()-1))
    inSegments->Append(aSeg);
  else{
    for(int s=0; s<aSeg->GetNsegments(); s++) inSegments->AddSegment(aSeg->GetStart(s),aSeg->GetEnd(s));
  }

  // data structure
  if(!dataseq->SetSegments(aSeg)){
    cerr<<"Omicron::InitSegments: cannot initiate data segments."<<endl;
    return false;
  }

  // update FFL
  if(FFL!=NULL){
    if(!FFL->ExtractChannels(aSeg->GetStart(0))){
      cerr<<"Omicron::InitSegments: cannot update FFL info."<<endl;
      return false;
    }
  }
  
  // update time offset
  timeoffset = aTimeOffset;

  if(fVerbosity>1){
    cout<<"\t- N segments       = "<<aSeg->GetNsegments()<<endl;
    cout<<"\t- Livetime         = "<<aSeg->GetLiveTime()<<endl;
    cout<<"\t- Plot time offset = "<<timeoffset<<endl;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::MakeDirectories(const double aGPS){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::MakeDirectories: the Omicron object is corrupted"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::MakeDirectories: make directory structure..."<<endl;
  fOutdir.clear();
  ostringstream tmpstream;

  if(!aGPS){// basic trigger directories
    for(int c=0; c<(int)fChannels.size(); c++){
      fOutdir.push_back(fMaindir+"/"+fChannels[c]);
      if(fVerbosity>1) cout<<"\t- "<<fOutdir[c]<<endl;
      system(("mkdir -p "+fOutdir[c]).c_str());
      triggers[c]->SetOutputDirectory(fOutdir[c]);
    }
  }
  else{
    tmpstream<<fMaindir<<"/"<<setprecision(3)<<fixed<<aGPS;
    fScandir=tmpstream.str();
    tmpstream.clear(); tmpstream.str("");
    for(int c=0; c<(int)fChannels.size(); c++){
      fOutdir.push_back(fScandir+"/"+fChannels[c]);
      if(fVerbosity>1) cout<<"\t- "<<fOutdir[c]<<endl;
      system(("mkdir -p "+fOutdir[c]).c_str());
      triggers[c]->SetOutputDirectory(fOutdir[c]);
    }
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::NewChunk(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::NewChunk: the Omicron object is corrupted"<<endl;
    return false;
  }

  // load new chunk
  if(fVerbosity) cout<<"Omicron::NewChunk: load a new chunk..."<<endl;
  if(!dataseq->NewChunk()) return false;
  if(fVerbosity>1) cout<<"\tchunk "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<" is loaded"<<endl;

  chunk_ctr++;// one more chunk
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::NewChannel(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::NewChannel: the Omicron object is corrupted"<<endl;
    return false;
  }

  // load new channel
  if(fVerbosity) cout<<"Omicron::NewChannel: load a new channel..."<<endl;
  chanindex++;

  // last channel
  if(chanindex==(int)fChannels.size()){
    chanindex=-1;
    if(fVerbosity>1) cout<<"\tno more channels to load"<<endl;
    return false; 
  }

  if(fVerbosity>1) cout<<"\tchannel "<<fChannels[chanindex]<<" is loaded"<<endl;
  chan_ctr[chanindex]++;
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::LoadData(double **aDataVector, int *aSize){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::LoadData: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(FFL==NULL){
    cerr<<"Omicron::LoadData: this function can only be used with a valid FFL object"<<endl;
    return false;
  }

  if(fVerbosity){
    if(fInjChan.size()) cout<<"Omicron::LoadData: load data vector and add injections..."<<endl;
    else cout<<"Omicron::LoadData: load data vector..."<<endl;
  }

  // get data vector
  *aDataVector = FFL->GetData(*aSize, fChannels[chanindex], dataseq->GetChunkTimeStart(), dataseq->GetChunkTimeEnd());

  // cannot retrieve data
  if(*aSize<=0){
    cerr<<"Omicron::LoadData: cannot retrieve data ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
    return false;
  }
  
  // get injection vector and inject it
  if(fInjChan.size()){
    int dsize_inj;
    double *dvector_inj = FFL->GetData(dsize_inj, fInjChan[chanindex], dataseq->GetChunkTimeStart(), dataseq->GetChunkTimeEnd());

    // cannot retrieve data
    if(dsize_inj<=0){
      cerr<<"Omicron::LoadData: cannot retrieve injection data ("<<fInjChan[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
      delete aDataVector; aDataVector=NULL; *aSize=0;
      return false;
    }

    // size mismatch
    if(dsize_inj!=*aSize){
      cerr<<"Omicron::LoadData: the sampling of the injection channel is not the same as the sampling of the main channel ("<<fInjChan[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
      delete aDataVector; aDataVector=NULL; *aSize=0;
      delete dvector_inj;
      return false;
    }

    // add injections in the data
    for(int d=0; d<*aSize; d++) *aDataVector[d]+=(fInjFact[chanindex]*dvector_inj[d]);
    delete dvector_inj;
  }

  chan_data_ctr[chanindex]++;
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
int Omicron::Condition(const int aInVectSize, double *aInVect){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::Condition: the Omicron object is corrupted"<<endl;
    return -1;
  }

  // Input data checks
  if(aInVect==NULL){
    cerr<<"Omicron::Condition: input vector is NULL ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
    return 1;
  }
  if(aInVectSize<=0){
    cerr<<"Omicron::Condition: input vector is empty ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
    return 2;
  }
  if(aInVect[0]==aInVect[aInVectSize-1]){
    cerr<<"Omicron::Condition: input vector seems to be flat ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
    return 3;
  }
    
  if(fVerbosity) cout<<"Omicron::Condition: condition data vector..."<<endl;

  // test native sampling (and update if necessary)
  int nativesampling = aInVectSize/(dataseq->GetCurrentChunkDuration());
  if(!sample[chanindex]->SetFrequencies(nativesampling,fSampleFrequency,fFreqRange[0])){
    cerr<<"Omicron::Condition: the Sample structure could not be parameterized ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
    status_OK=false;
    return -2;
  }

  // transform data vector
  if(fVerbosity>1) cout<<"\t- transform data vector..."<<endl;
  if(!sample[chanindex]->Transform(aInVectSize, aInVect, dataseq->GetCurrentChunkDuration()*fSampleFrequency, ChunkVect)) return 4;
 
  // Make spectrum
  if(fVerbosity>1) cout<<"\t- make spectrum..."<<endl;
  if(!spectrum->LoadData(dataseq->GetCurrentChunkDuration()*fSampleFrequency, ChunkVect)) return 5;

  // set new power for this chunk
  if(fVerbosity>1) cout<<"\t- normalize tiling..."<<endl;
  if(!tile->SetPowerSpectrum(spectrum)) return 6;

  chan_cond_ctr[chanindex]++;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::Project(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::Project: the Omicron object is corrupted"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::Project: project data onto the tiles..."<<endl;

  // loop over segments
  int segsize=dataseq->GetSegmentDuration()*fSampleFrequency;
  int ovsize=dataseq->GetOverlapDuration()*fSampleFrequency;
  for(int s=0; s<dataseq->GetNSegments(); s++){
    
    // fill segment vector (time-domain) and apply tukey window
    for(int i=0; i<segsize; i++)
      SegVect[i] = ChunkVect[s*(segsize-ovsize)+i] * TukeyWindow[i];
    
    // whiten data
    if(!Whiten(&(dataRe[s]), &(dataIm[s]))) return false;
    


    /////////////////////
    //*** MAPS
    if(fOutProducts.find("maps")!=string::npos) MakeMaps(s);

    //*** TRIGGERS
    if(fOutProducts.find("triggers")!=string::npos&&!triggers[chanindex]->GetMaxFlag()){
      
      if(tile->GetTriggers(triggers[chanindex],dataRe[s],dataIm[s],dataseq->GetSegmentTimeStart(s),dataseq->GetCurrentOverlapDuration()-fOverlapDuration))
	triggers[chanindex]->AddSegment(dataseq->GetSegmentTimeStart(s)+dataseq->GetCurrentOverlapDuration()-fOverlapDuration/2,dataseq->GetSegmentTimeEnd(s)-fOverlapDuration/2);// add processed segment
    }
  
    delete dataRe[s];
    delete dataIm[s];
  } 

  //if(!triggers[chanindex]->SortTriggers()) return -3;

  // clustering if any
  //if(fClusterAlgo[0].compare("none")){
  //if(!triggers[chanindex]->Clusterize(fClusterAlgo[0])) return -4;
  // if(fVerbosity>1) cout<<"\t- "<<triggers[chanindex]->GetNClusters()<<" clusters were found"<<endl;
  //}
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::WriteOutput(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::WriteOutput: the Omicron object is corrupted"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::WriteOutput: write output..."<<endl;

  //*** TRIGGERS
  if(fOutProducts.find("triggers")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write triggers..."<<endl;

    // don't save if max flag
    if(triggers[chanindex]->GetMaxFlag()){
      cerr<<"Omicron::WriteOutput: the number of triggers is maxed-out ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
      triggers[chanindex]->Reset();
    }
    else{      
      // save triggers for this chunk
      if(triggers[chanindex]->Write(fWriteMode).compare("none")){
	chan_wtrig_ctr[chanindex]++;
	outSegments[chanindex]->AddSegment(dataseq->GetChunkTimeStart()+fOverlapDuration/2,dataseq->GetChunkTimeEnd()-fOverlapDuration/2);// FIXME: what about maps?
      }
    }
  }
  
  //*** MAPS
  if(fOutProducts.find("maps")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write ASD..."<<endl;
    if(!WriteMaps()){
      cerr<<"Omicron::WriteOutput: the maps cannot be saved ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
    }
    else chan_wmaps_ctr[chanindex]++;
  }

  //*** ASD
  if(fOutProducts.find("asd")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write ASD..."<<endl;
    SaveAPSD("ASD");
  }
  
  //*** PSD
  if(fOutProducts.find("psd")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write PSD..."<<endl;
    SaveAPSD("PSD");
  }
  
  //*** TS
  if(fOutProducts.find("timeseries")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write time-series..."<<endl;
    SaveTS();
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::PrintStatusInfo(void){
////////////////////////////////////////////////////////////////////////////////////

  cout<<"\n************* Omicron status info *************"<<endl;
  cout<<"requested start         = "<<inSegments->GetStart(0)<<endl;
  cout<<"requested end           = "<<inSegments->GetEnd(inSegments->GetNsegments()-1)<<endl;
  cout<<"requested livetime      = "<<inSegments->GetLiveTime()<<"s"<<endl;
  cout<<"number of loaded chunks = "<<chunk_ctr<<endl;

  for(int c=0; c<(int)fChannels.size(); c++){
    cout<<"\n*** "<<fChannels[c]<<endl;
    if(outSegments[c]->GetNsegments()){
      cout<<"start_out           = "<<outSegments[c]->GetStart(0)<<endl;
      cout<<"end_out             = "<<outSegments[c]->GetEnd(outSegments[c]->GetNsegments()-1)<<endl;
      cout<<"trigger livetime    = "<<outSegments[c]->GetLiveTime()<<"s ("<<outSegments[c]->GetLiveTime()/inSegments->GetLiveTime()*100<<"%)"<<endl;
    }
    else{
      cout<<"start_out             = -1"<<endl;
      cout<<"end_out               = -1"<<endl;
      cout<<"processed livetime    = "<<"0s (0%)"<<endl;
    }
    cout<<"number of calls                = "<<chan_ctr[c]<<endl;
    cout<<"number of data calls           = "<<chan_data_ctr[c]<<endl;
    cout<<"number of conditioning calls   = "<<chan_cond_ctr[c]<<endl;
    cout<<"number of output map sets      = "<<chan_wmaps_ctr[c]<<endl;
    cout<<"number of output trigger files = "<<chan_wtrig_ctr[c]<<endl;
  }
  cout<<"***********************************************\n"<<endl;

  return;
}

////////////////////////////////////////////////////////////////////////////////////
Segments* Omicron::GetTriggerSegments(TH1D *aThr, const double aInfValue){
////////////////////////////////////////////////////////////////////////////////////
  Segments* empty = new Segments();
  if(!status_OK){
    cerr<<"Omicron::GetTriggerSegments: the Omicron object is corrupted"<<endl;
    return empty;
  }
  delete empty;

  return triggers[chanindex]->GetTriggerSegments(aThr,aInfValue);
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::Whiten(double **aDataRe, double **aDataIm){
////////////////////////////////////////////////////////////////////////////////////
  // fft-forward
  if(!offt->Forward(SegVect)){
    cerr<<"Omicron::Whiten: FFT-forward failed"<<endl;
    return false;
  }
    
  // get ffted data
  *aDataRe=offt->GetRe();
  *aDataIm=offt->GetIm();
  if(*aDataRe==NULL||*aDataIm==NULL){
    cerr<<"Omicron::Whiten: cannot retrieve FFT data"<<endl;
    return false;
  }

  // zero-out below high-frequency cutoff
  int icutoff = (int)floor(fFreqRange[0]*(double)(dataseq->GetSegmentDuration()));
  for(int i=0; i<icutoff; i++){
    (*aDataRe)[i]=0.0;
    (*aDataIm)[i]=0.0;
  }

  // normalize data by the ASD
  double asdval;
  int ssize=dataseq->GetSegmentDuration()*fSampleFrequency/2;
  for(int i=icutoff; i<ssize; i++){
    asdval=spectrum->GetPower(i*(double)fSampleFrequency/(double)(2*ssize));
    if(asdval<=0){
      cerr<<"Omicron::Whiten: could not retrieve power"<<endl;
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
void Omicron::SaveAPSD(const string type){
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
    GAPSD->SetTitle((fChannels[chanindex]+": Power spectrum density").c_str());
  }
  else{
    GAPSD->GetHistogram()->SetYTitle("Amplitude [Amp/#sqrt{Hz}]");
    GAPSD->SetTitle((fChannels[chanindex]+": Amplitude spectrum density").c_str());
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
  if(!type.compare("ASD")) ss<<"ASD_"<<fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd();
  else ss<<"PSD_"<<fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd();
  GAPSD->SetName(ss.str().c_str());
  ss.str(""); ss.clear();

  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fpsd;
    if(first_save[chanindex]){
      fpsd=new TFile((fOutdir[chanindex]+"/"+fChannels[chanindex]+"_data.root").c_str(),"RECREATE");
      first_save[chanindex]=false;
    }
    else fpsd=new TFile((fOutdir[chanindex]+"/"+fChannels[chanindex]+"_data.root").c_str(),"UPDATE");
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
  if(fOutFormat.find("jpg")!=string::npos) form.push_back("jpg"); 
  if(fOutFormat.find("svg")!=string::npos) form.push_back("svg"); 
  if(form.size()){
    for(int f=0; f<(int)form.size(); f++){
      if(!type.compare("ASD")) ss<<fOutdir[chanindex]<<"/"<<fChannels[chanindex]<<"_asd_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"."<<form[f];
      else ss<<fOutdir[chanindex]<<"/"<<fChannels[chanindex]<<"_psd_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"."<<form[f];
      GPlot->Print(ss.str().c_str());
      ss.str(""); ss.clear();
      if(!type.compare("ASD")) ss<<fOutdir[chanindex]<<"/"<<fChannels[chanindex]<<"_asd_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"th."<<form[f];
      else ss<<fOutdir[chanindex]<<"/"<<fChannels[chanindex]<<"_psd_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"th."<<form[f];
      GPlot->Print(ss.str(),0.5);
      ss.str(""); ss.clear();
    }
  }
  
  form.clear();
  delete GAPSD;
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveTS(void){
////////////////////////////////////////////////////////////////////////////////////

  TGraph *GDATA = new TGraph(fSampleFrequency*dataseq->GetCurrentChunkDuration());
  if(GDATA==NULL) return;

  stringstream ss;
  ss<<"ts_"<<fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd();
  GDATA->SetName(ss.str().c_str());
  ss.str(""); ss.clear();

  for(int i=0; i<fSampleFrequency*dataseq->GetCurrentChunkDuration(); i++) GDATA->SetPoint(i,(double)dataseq->GetChunkTimeStart()+(double)i*(double)(fSampleFrequency)-timeoffset,ChunkVect[i]);
     
  // cosmetics
  GPlot->SetLogx(0);
  GPlot->SetLogy(0);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);
  GPlot->Draw(GDATA,"APL");
  GDATA->GetHistogram()->SetXTitle("Time [s]");
  GDATA->GetHistogram()->SetYTitle("Amplitude [?]");
  GDATA->SetTitle((fChannels[chanindex]+": amplitude time series").c_str());
  GDATA->SetLineWidth(1);
  GDATA->GetXaxis()->SetNoExponent();
  GDATA->GetXaxis()->SetTitleOffset(1.1);
  GDATA->GetXaxis()->SetLabelSize(0.045);
  GDATA->GetYaxis()->SetLabelSize(0.045);
  GDATA->GetXaxis()->SetTitleSize(0.045);
  GDATA->GetYaxis()->SetTitleSize(0.045);

  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fdata;
    if(first_save[chanindex]){
      fdata=new TFile((fOutdir[chanindex]+"/"+fChannels[chanindex]+"_data.root").c_str(),"RECREATE");
      first_save[chanindex]=false;
    }
    else fdata=new TFile((fOutdir[chanindex]+"/"+fChannels[chanindex]+"_data.root").c_str(),"UPDATE");
    fdata->cd();
    GDATA->Write();
    fdata->Close();
  }

  // Graphix
  vector <string> form;
  if(fOutFormat.find("gif")!=string::npos) form.push_back("gif");
  if(fOutFormat.find("png")!=string::npos) form.push_back("png");
  if(fOutFormat.find("pdf")!=string::npos) form.push_back("pdf");
  if(fOutFormat.find("ps")!=string::npos)  form.push_back("ps");
  if(fOutFormat.find("xml")!=string::npos) form.push_back("xml");
  if(fOutFormat.find("eps")!=string::npos) form.push_back("eps"); 
  if(fOutFormat.find("jpg")!=string::npos) form.push_back("jpg"); 
  if(fOutFormat.find("svg")!=string::npos) form.push_back("svg"); 
  if(form.size()){
    
    // give the time offset info in the title
    if(timeoffset){
      ss<<fChannels[chanindex]+": amplitude time series centered at "<<fixed<<setprecision(3)<<timeoffset;
      //tcenter=0;
      GDATA->SetTitle(ss.str().c_str());
      ss.str(""); ss.clear();
    }

    // zoom
    double tcenter;
    if(timeoffset) tcenter=0.0;
    else tcenter=(double)(dataseq->GetChunkTimeStart()+dataseq->GetChunkTimeEnd())/2.0;
    for(int w=(int)fWindows.size()-1; w>=0; w--){
      GDATA->GetXaxis()->SetLimits(tcenter-(double)fWindows[w]/2.0,tcenter+(double)fWindows[w]/2.0);
      for(int f=0; f<(int)form.size(); f++){
	ss<<fOutdir[chanindex]<<"/"<<fChannels[chanindex]<<"_ts_dt"<<fWindows[w]<<"."<<form[f];
	GPlot->Print(ss.str().c_str());
	ss.str(""); ss.clear();
	ss<<fOutdir[chanindex]<<"/"<<fChannels[chanindex]<<"_tsth_dt"<<fWindows[w]<<"."<<form[f];
	GPlot->Print(ss.str(),0.5);
	ss.str(""); ss.clear();
      }
    }
  }
  
  form.clear();
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
void Omicron::PrintASCIIlogo(void){
////////////////////////////////////////////////////////////////////////////////////

  cout<<endl;
  cout<<endl;
  cout<<"############################################################################"<<endl;
  cout<<"############################################################################"<<endl;
  cout<<endl;
  cout<<"          ██████╗ ███╗   ███╗██╗ ██████╗██████╗  ██████╗ ███╗   ██╗"<<endl;
  cout<<"         ██╔═══██╗████╗ ████║██║██╔════╝██╔══██╗██╔═══██╗████╗  ██║"<<endl;
  cout<<"         ██║   ██║██╔████╔██║██║██║     ██████╔╝██║   ██║██╔██╗ ██║"<<endl;
  cout<<"         ██║   ██║██║╚██╔╝██║██║██║     ██╔══██╗██║   ██║██║╚██╗██║"<<endl;
  cout<<"         ╚██████╔╝██║ ╚═╝ ██║██║╚██████╗██║  ██║╚██████╔╝██║ ╚████║"<<endl;
  cout<<"          ╚═════╝ ╚═╝     ╚═╝╚═╝ ╚═════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝"<<endl;
  cout<<"                                                               V1R4"<<endl;
  cout<<"############################################################################"<<endl;
  cout<<"############################################################################"<<endl;
  cout<<endl;
  cout<<endl;

  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::PrintMessage(const string aMessage){
////////////////////////////////////////////////////////////////////////////////////
  time ( &timer );
  ptm = gmtime ( &timer );
  cout<<"\n("<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) [+"<<timer-timer_start<<"s]: "<<aMessage<<endl;
  return;
}
