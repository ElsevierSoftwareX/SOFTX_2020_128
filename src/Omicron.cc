//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

ClassImp(Omicron)


////////////////////////////////////////////////////////////////////////////////////
Omicron::Omicron(const string aOptionFile, const int aGpsRef){ 
////////////////////////////////////////////////////////////////////////////////////
  PrintASCIIlogo();
  gErrorIgnoreLevel = 3000;
  status_OK=true;

  // init timer
  time ( &timer );
  timer_start=timer;

  // parse option file
  if(fVerbosity) cout<<"Omicron::Omicron: init options..."<<endl;
  fOptionFile=aOptionFile;
  ReadOptions(aGpsRef);

  // sort time windows
  // FIXME: to move in Otile
  std::sort(fWindows.begin(), fWindows.end());

  // spectrum status
  for(int c=0; c<nchannels; c++) status_OK*=spectrum1[c]->GetStatus()*spectrum2[c]->GetStatus()*spectrumw->GetStatus();

  // chunk FFT
  offt = new fft(tile->GetTimeRange()*triggers[0]->GetWorkingFrequency(), fftplan, "r2c");
 
  // data containers
  ChunkVect     = new double    [offt->GetSize_t()];
  TukeyWindow   = GetTukeyWindow(offt->GetSize_t(),
				 tile->GetOverlapDuration()*triggers[0]->GetWorkingFrequency());
 
  // metadata field definition
  vector <string> fOptionName;
  vector <string> fOptionType;
  fOptionName.push_back("omicron_DATA_FFLFILE");                fOptionType.push_back("s");
  fOptionName.push_back("omicron_DATA_CHANNEL");                fOptionType.push_back("s");
  fOptionName.push_back("omicron_DATA_SAMPLEFREQUENCY");        fOptionType.push_back("i");
  fOptionName.push_back("omicron_INJECTION_FFL");               fOptionType.push_back("s");
  fOptionName.push_back("omicron_INJECTION_CHANNEL");           fOptionType.push_back("s");
  fOptionName.push_back("omicron_INJECTION_FACTOR");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_FILENAME");          fOptionType.push_back("s");
  fOptionName.push_back("omicron_INJECTION_SG");                fOptionType.push_back("i");
  fOptionName.push_back("omicron_INJECTION_SGTIMEMIN");         fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGTIMEMAX");         fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGFREQUENCYMIN");    fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGFREQUENCYMAX");    fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGQMIN");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGQMAX");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGAMPMIN");          fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGAMPMAX");          fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_FMIN");              fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_FMAX");              fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_QMIN");              fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_QMAX");              fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_CHUNKDURATION");     fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_OVERLAPDURATION");   fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_MISMATCHMAX");       fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_SNRTHRESHOLD");      fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_PSDLENGTH");         fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_HIGHPASS");          fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_CLUSTERING");        fOptionType.push_back("s");
  fOptionName.push_back("omicron_PARAMETER_CLUSTERDT");         fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_FFTPLAN");           fOptionType.push_back("s");
  fOptionName.push_back("omicron_PARAMETER_TRIGGERRATEMAX");    fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_TRIGGERBUFFERSIZE"); fOptionType.push_back("i");
  fOptionName.push_back("omicron_OUTPUT_DIRECTORY");            fOptionType.push_back("s");
  fOptionName.push_back("omicron_OUTPUT_VERBOSITY");            fOptionType.push_back("i");
  fOptionName.push_back("omicron_OUTPUT_FORMAT");               fOptionType.push_back("s");
  fOptionName.push_back("omicron_OUTPUT_PRODUCTS");             fOptionType.push_back("s");
  fOptionName.push_back("omicron_OUTPUT_STYLE");                fOptionType.push_back("s");

  // triggers metadata
  for(int c=0; c<nchannels; c++){
    triggers[c]->SetMprocessname("OMICRON");
    triggers[c]->SetProcessVersion(GetVersion());
    status_OK*=triggers[c]->InitUserMetaData(fOptionName,fOptionType);
    if(FFL==NULL) status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],"none");
    else          status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],FFL->GetInputFfl());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[1],triggers[c]->GetName());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[2],triggers[c]->GetWorkingFrequency());
    if(FFL_inject==NULL) status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],"none");
    else                 status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],FFL_inject->GetInputFfl());
    if(fInjChan.size()){
      status_OK*=triggers[c]->SetUserMetaData(fOptionName[4],fInjChan[c]);
      status_OK*=triggers[c]->SetUserMetaData(fOptionName[5],fInjFact[c]);
    }
    else{
      status_OK*=triggers[c]->SetUserMetaData(fOptionName[4],"none");
      status_OK*=triggers[c]->SetUserMetaData(fOptionName[5],0.0);
    }
    if(inject==NULL) status_OK*=triggers[c]->SetUserMetaData(fOptionName[6],"none");
    else             status_OK*=triggers[c]->SetUserMetaData(fOptionName[6],inject[c]->GetInputFilePattern());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[7],fsginj);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[8],oinj->GetTimeMin());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[9],oinj->GetTimeMax());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[10],oinj->GetFrequencyMin());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[11],oinj->GetFrequencyMax());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[12],oinj->GetQMin());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[13],oinj->GetQMax());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[14],oinj->GetAmplitudeMin());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[15],oinj->GetAmplitudeMax());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[16],tile->GetFrequencyMin());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[17],tile->GetFrequencyMax());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[18],tile->GetQ(0));
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[19],tile->GetQ(tile->GetNQ()-1));
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[20],tile->GetTimeRange());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[21],tile->GetOverlapDuration());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[22],tile->GetMismatchMax());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[23],tile->GetSNRTriggerThr());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[24],spectrum1[c]->GetDataBufferLength());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[25],triggers[c]->GetHighPassFrequency());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[26],fClusterAlgo);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[27],triggers[c]->GetClusterizeDt());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[28],fftplan);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[29],fratemax);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[30],triggers[c]->GetBufferSize());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[31],fMaindir);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[32],fVerbosity);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[33],fOutFormat);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[34],fOutProducts);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[35],GPlot->GetCurrentStyle());
  }

  // default output directory: main dir
  maindir=fMaindir;
  for(int c=0; c<nchannels; c++){
    outdir.push_back(maindir);
  }
  
  //SG injections file
  if(fsginj){
    oinjfile.open((maindir+"/sginjections.txt").c_str());
    oinjfile<<"# Time # Frequency # Q # Amplitude # Phase # True SNR # sigma_t # sigma_f"<<endl;
    oinjfile.precision(5);
  }

  // default plottime offset
  SetPlotTimeOffset();

  // process monitoring
  chanindex      = -1;
  inSegments     = new Segments();
  outSegments    = new Segments* [nchannels];
  chunk_ctr      = 0;
  chan_ctr       = new int       [nchannels];
  chan_data_ctr  = new int       [nchannels];
  chan_cond_ctr  = new int       [nchannels];
  chan_proj_ctr  = new int       [nchannels];
  chan_write_ctr = new int       [nchannels];
  chan_mapsnrmax = new double    [nchannels];
  trig_ctr       = new int       [nchannels];
  for(int c=0; c<nchannels; c++){
    outSegments[c]    = new Segments();
    chan_ctr[c]       = 0;
    chan_data_ctr[c]  = 0;
    chan_cond_ctr[c]  = 0;
    chan_proj_ctr[c]  = 0;
    chan_write_ctr[c] = 0;
    chan_mapsnrmax[c] = 0.0;
    trig_ctr[c]       = 0;
  }
}

////////////////////////////////////////////////////////////////////////////////////
Omicron::~Omicron(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"Omicron::~Omicron"<<endl;
  if(status_OK&&fOutProducts.find("html")!=string::npos) MakeHtml(); // print final html report
  delete inSegments;
  delete chan_ctr;
  delete chan_data_ctr;
  delete chan_cond_ctr;
  delete chan_proj_ctr;
  delete chan_write_ctr;
  delete chan_mapsnrmax;
  delete trig_ctr;
  for(int c=0; c<nchannels; c++){
    delete outSegments[c];
    delete triggers[c];
    delete spectrum1[c];
    delete spectrum2[c];
  }
  delete outSegments;
  delete spectrum1;
  delete spectrum2;
  delete spectrumw;
  delete triggers;
  if(FFL_inject!=FFL&&FFL_inject!=NULL) delete FFL_inject;
  if(FFL!=NULL) delete FFL;
  if(inject!=NULL){
    for(int c=0; c<nchannels; c++) delete inject[c];
    delete inject;
  }
  delete oinj;
  delete tile;
  delete GPlot;
  delete ChunkVect;
  delete TukeyWindow;
  delete offt;
  
  if(fsginj) oinjfile.close();
  outdir.clear();
  chunkcenter.clear();
  chunktfile.clear();
  fInjChan.clear();
  fInjFact.clear();
  fWindows.clear();
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::InitSegments(Segments *aInSeg, Segments *aOutSeg){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::InitSegments: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aInSeg==NULL){
    cerr<<"Omicron::InitSegments: the input segment is NULL"<<endl;
    return false;
  }
  if(!aInSeg->GetStatus()){
    cerr<<"Omicron::InitSegments: the input segment is corrupted"<<endl;
    return false;
  }
  if(!aInSeg->GetNsegments()){
    cerr<<"Omicron::InitSegments: there is no input segment"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::InitSegments: initiate data segments..."<<endl;

  // cumulative input segment monitor
  // (try to optimize the update of inSegments)
  if(inSegments->GetLiveTime()&&aInSeg->GetStart(0)>=inSegments->GetStart(inSegments->GetNsegments()-1))
    inSegments->Append(aInSeg);
  else
    for(int s=0; s<aInSeg->GetNsegments(); s++) inSegments->AddSegment(aInSeg->GetStart(s),aInSeg->GetEnd(s));
  
  // data structure
  if(!tile->SetSegments(aInSeg, aOutSeg)){
    cerr<<"Omicron::InitSegments: cannot initiate data segments."<<endl;
    return false;
  }

  // update channel list
  if(FFL!=NULL){
    if(!FFL->ExtractChannels(aInSeg->GetStart(0))){
      cerr<<"Omicron::InitSegments: cannot update FFL info."<<endl;
      return false;
    }
  }
  if(FFL_inject!=NULL&&FFL_inject!=FFL){
    if(!FFL_inject->ExtractChannels(aInSeg->GetStart(0))){
      cerr<<"Omicron::InitSegments: cannot update FFL info for injections."<<endl;
      return false;
    }
  }

  if(fVerbosity>1){
    cout<<"\t- N input segments   = "<<aInSeg->GetNsegments()<<endl;
    cout<<"\t- Input livetime     = "<<aInSeg->GetLiveTime()<<endl;
    if(aOutSeg!=NULL){
      cout<<"\t- N output segments  = "<<aOutSeg->GetNsegments()<<endl;
      cout<<"\t- Output livetime    = "<<aOutSeg->GetLiveTime()<<endl;
    }
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::MakeDirectories(const double aId){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::MakeDirectories: the Omicron object is corrupted"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::MakeDirectories: make directory structure..."<<endl;
  ostringstream tmpstream;

  // Main dir
  if(!aId){// single level
    maindir=fMaindir;
  }
  else{// two level
    tmpstream<<fMaindir<<"/"<<setprecision(3)<<fixed<<aId;
    maindir=tmpstream.str();
    tmpstream.clear(); tmpstream.str("");
    if(system(("mkdir -p "+maindir).c_str())){
      cerr<<"Omicron::MakeDirectories: the output directory cannot be created"<<endl;
      return false;
    }
  }

  // channel directories
  outdir.clear();
  for(int c=0; c<nchannels; c++){
    outdir.push_back(maindir+"/"+triggers[c]->GetName());
    if(system(("mkdir -p "+outdir[c]).c_str())){
      cerr<<"Omicron::MakeDirectories: the output directory cannot be created"<<endl;
      return false;
    }
  }


  //SG injections file
  if(fsginj){
    oinjfile.close();
    oinjfile.open((maindir+"/sginjections.txt").c_str());
    oinjfile<<"# Time # Frequency # Q # Amplitude # Phase # True SNR # sigma_t # sigma_f"<<endl;
    oinjfile.precision(5);
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

  bool newseg;
  
  // load new chunk
  if(fVerbosity) cout<<"Omicron::NewChunk: call a new chunk..."<<endl;
  if(!tile->NewChunk(newseg)) return false;
  if(fVerbosity>1){
    if(newseg) cout<<"\t- chunk "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<" is loaded (start a new segment)"<<endl;
    else cout<<"\t- chunk "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<" is loaded"<<endl;
  }

  // save info for html report
  if(fOutProducts.find("html")!=string::npos) chunkcenter.push_back(tile->GetChunkTimeCenter());

  // new segment --> reset PSD buffer
  if(newseg)
    for(int c=0; c<nchannels; c++){ spectrum1[c]->Reset(); spectrum2[c]->Reset(); }  

  // generate SG parameters
  if(fsginj) oinj->MakeWaveform();

  chunk_ctr++;// one more chunk
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::DefineNewChunk(const int aTimeStart, const int aTimeEnd, const bool aResetPSDBuffer){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::DefineNewChunk: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aTimeEnd-aTimeStart!=tile->GetTimeRange()){
    cerr<<"Omicron::DefineNewChunk: the input chunk is not the right duration"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::DefineNewChunk: "<<aTimeStart<<"-"<<aTimeEnd<<"..."<<endl;

  // make segment for this chunk
  Segments *Stmp = new Segments(aTimeStart,aTimeEnd);
  if(!InitSegments(Stmp)){
    cerr<<"Omicron::DefineNewChunk: the input time segment is not correct"<<endl;
    delete Stmp;
    return false;
  }
  delete Stmp;

  // set tiling for this segment
  bool newseg; // not used
  if(!tile->NewChunk(newseg)) return false;
    
  // save info for html report
  if(fOutProducts.find("html")!=string::npos) chunkcenter.push_back(tile->GetChunkTimeCenter());

  // reset PSD buffer
  if(aResetPSDBuffer)
    for(int c=0; c<nchannels; c++){ spectrum1[c]->Reset(); spectrum2[c]->Reset(); }
      
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
  if(chanindex==nchannels){
    chanindex=-1;
    if(fVerbosity>1) cout<<"\t- no more channels to load"<<endl;
    return false; 
  }

  // new channel
  if(fVerbosity>1) cout<<"\t- channel "<<triggers[chanindex]->GetName()<<" is loaded"<<endl;
  chan_ctr[chanindex]++;

  // reset number of tiles above threshold
  trig_ctr[chanindex]=0;
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::LoadData(double **aDataVector, int *aSize){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::LoadData: the Omicron object is corrupted"<<endl;
    aDataVector=NULL; *aSize=0;
    return false;
  }
  if(FFL==NULL){
    cerr<<"Omicron::LoadData: this function can only be used with a valid FFL object"<<endl;
    aDataVector=NULL; *aSize=0;
    return false;
  }
  if(!tile->GetChunkTimeCenter()){
    cerr<<"Omicron::LoadData: no chunk called yet"<<endl;
    aDataVector=NULL; *aSize=0;
    return false;
  }
  if(chanindex<0){
    cerr<<"Omicron::LoadData: no channel called yet"<<endl;
    aDataVector=NULL; *aSize=0;
    return false;
  }

  if(fVerbosity){
    if(fInjChan.size()) cout<<"Omicron::LoadData: load data vector and add injections..."<<endl;
    else cout<<"Omicron::LoadData: load data vector..."<<endl;
  }

  // get data vector
  if(fVerbosity>1) cout<<"\t- get data from frames..."<<endl;
  *aDataVector = FFL->GetData(*aSize, triggers[chanindex]->GetName(), tile->GetChunkTimeStart(), tile->GetChunkTimeEnd());

  // cannot retrieve data
  if(*aSize<=0){
    cerr<<"Omicron::LoadData: cannot retrieve data ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    aDataVector=NULL; *aSize=0;
    return false;
  }

  // test native sampling (and update if necessary)
  int nativesampling = *aSize/(tile->GetTimeRange());
  if(!triggers[chanindex]->SetNativeFrequency(nativesampling)){
    cerr<<"Omicron::LoadData: incompatible native/working frequency ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return false;
  }

  // add software injections if any
  if(inject!=NULL){
    if(fVerbosity>1) cout<<"\t- perform software injections..."<<endl;
    inject[chanindex]->UpdateNativeSamplingFrequency();
    if(!inject[chanindex]->Inject(*aSize, *aDataVector, tile->GetChunkTimeStart())){
      cerr<<"Omicron::LoadData: failed to inject ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
      return false;
    }
  }
  
  // get injection vector and inject it
  if(fInjChan.size()){
    if(fVerbosity>1) cout<<"\t- perform stream injections..."<<endl;
    int dsize_inj;
    double *dvector_inj = FFL_inject->GetData(dsize_inj, fInjChan[chanindex], tile->GetChunkTimeStart(), tile->GetChunkTimeEnd());

    // cannot retrieve data
    if(dsize_inj<=0){
      cerr<<"Omicron::LoadData: cannot retrieve injection data ("<<fInjChan[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
      delete aDataVector; aDataVector=NULL; *aSize=0;
      return false;
    }

    // size mismatch
    if(dsize_inj!=*aSize){
      cerr<<"Omicron::LoadData: the sampling of the injection channel is not the same as the sampling of the main channel ("<<fInjChan[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
      delete aDataVector; aDataVector=NULL; *aSize=0;
      delete dvector_inj;
      return false;
    }

    // add injections in the data
    for(int d=0; d<*aSize; d++) (*aDataVector)[d]+=(fInjFact[chanindex]*dvector_inj[d]);
    delete dvector_inj;
  }

  // add sg injections
  // FIXME: could be moved in Condition()
  if(fsginj){
    for(int d=0; d<*aSize; d++) (*aDataVector)[d]+=oinj->GetWaveform(d,nativesampling);
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
    cerr<<"Omicron::Condition: input vector is NULL ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return 1;
  }
  if(aInVectSize<=0){
    cerr<<"Omicron::Condition: input vector is empty ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return 2;
  }
  if(aInVect[0]==aInVect[aInVectSize-1]&&IsFlat(aInVectSize,aInVect))
    cerr<<"Omicron::Condition: input vector is flat ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
  
  if(fVerbosity) cout<<"Omicron::Condition: condition data vector..."<<endl;

  // test native sampling (and update if necessary)
  int nativesampling = aInVectSize/(tile->GetTimeRange());
  if(!triggers[chanindex]->SetNativeFrequency(nativesampling)){
    cerr<<"Omicron::Condition: incompatible native/working frequency ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return 3;
  }

  // transform data vector
  if(fVerbosity>1) cout<<"\t- transform data vector..."<<endl;
  if(!triggers[chanindex]->Transform(aInVectSize, aInVect, offt->GetSize_t(), ChunkVect)) return 4;

  // apply Tukey Window
  if(fVerbosity>1) cout<<"\t- apply Tukey window..."<<endl;
  for(int i=0; i<offt->GetSize_t(); i++) ChunkVect[i] *= TukeyWindow[i];

  // update first spectrum (if enough data)
  if(fVerbosity>1) cout<<"\t- update spectrum 1..."<<endl;
  int dstart = (tile->GetCurrentOverlapDuration()-tile->GetOverlapDuration()/2)*triggers[chanindex]->GetWorkingFrequency(); // start of 'sane' data
  int dsize = (tile->GetTimeRange()-tile->GetCurrentOverlapDuration())*triggers[chanindex]->GetWorkingFrequency(); // size of 'sane' data
  if(!spectrum1[chanindex]->AddData(dsize, ChunkVect, dstart))
     cerr<<"Omicron::Condition: warning: this chunk is not used for PSD(1) estimation ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
  if(spectrum1[chanindex]->IsBufferEmpty()){
    cerr<<"Omicron::Condition: No PSD is available ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return 5;
  }

  // fft-forward the chunk data
  if(fVerbosity>1) cout<<"\t- move the data in the frequency domain..."<<endl;
  if(!offt->Forward(ChunkVect)) return 6;

  // 1st whitening
  if(fVerbosity>1) cout<<"\t- whiten chunk (1)"<<endl;
  Whiten(spectrum1[chanindex]);

  // back in the time domain
  if(fVerbosity>1) cout<<"\t- move the data back in the time domain..."<<endl;
  offt->Backward();

  // apply FFT (forward and backward) normalization
  for(int i=0; i<offt->GetSize_t(); i++) offt->SetRe_t(i, offt->GetRe_t(i)/(double)offt->GetSize_t());

  // update second spectrum
  if(fVerbosity>1) cout<<"\t- update spectrum 2..."<<endl;
  double *rvec = offt->GetRe_t();
  if(!spectrum2[chanindex]->AddData(dsize, rvec, dstart))
    cerr<<"Omicron::Condition: warning: this chunk is not used for PSD(2) estimation ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
  if(spectrum2[chanindex]->IsBufferEmpty()){// should never happen if it worked for the 1st spectrum
    delete rvec;
    return 7;
  }
  delete rvec;
  
  // fft-forward the chunk data
  if(fVerbosity>1) cout<<"\t- move the data in the frequency domain..."<<endl;
  offt->Forward();// This is mandatory because the frequency-domain vector was messed-up by the previous backward (r2c)

  // apply FFT (forward) normalization
  for(int i=0; i<offt->GetSize_f(); i++){
    offt->SetRe_f(i, offt->GetRe_f(i)/(double)triggers[chanindex]->GetWorkingFrequency());
    offt->SetIm_f(i, offt->GetIm_f(i)/(double)triggers[chanindex]->GetWorkingFrequency());
  }

  // 2nd whitening
  if(fVerbosity>1) cout<<"\t- whiten chunk (2)"<<endl;
  Whiten(spectrum2[chanindex]);

  // compute tiling power
  if(fVerbosity>1) cout<<"\t- compute tiling power..."<<endl;
  if(!tile->SetPower(spectrum1[chanindex],spectrum2[chanindex])) return 8;

  // save sg injection parameters
  // must be done here because the spectra are needed
  if(fsginj) SaveSG();

  chan_cond_ctr[chanindex]++;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////
int Omicron::Project(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::Project: the Omicron object is corrupted"<<endl;
    return -1;
  }

  if(fVerbosity) cout<<"Omicron::Project: project data onto the tiles..."<<endl;
  chan_proj_ctr[chanindex]++;
  trig_ctr[chanindex]=tile->ProjectData(offt);
  return trig_ctr[chanindex];
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::WriteOutput(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::WriteOutput: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(fVerbosity) cout<<"Omicron::WriteOutput: write chunk output..."<<endl;
 
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
    
  //*** CONDITIONNED TS
  if(fOutProducts.find("timeseries")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write conditionned time-series..."<<endl;
    SaveTS(false);
  }

  //*** WHITENED TS
  if(fOutProducts.find("white")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write whitened data..."<<endl;
    offt->Backward();// Back in time domain
    // IMPORTANT: after that, the frequency-domain vector of offt is corrupted (r2c)
    // apply FFT normalization
    for(int i=0; i<offt->GetSize_t(); i++) offt->SetRe_t(i, offt->GetRe_t(i)*triggers[chanindex]->GetWorkingFrequency()/(double)offt->GetSize_t());
    SaveTS(true);

    if(fOutProducts.find("whitepsd")!=string::npos){
      int dstart = (tile->GetCurrentOverlapDuration()-tile->GetOverlapDuration()/2)*triggers[chanindex]->GetWorkingFrequency(); // start of 'sane' data
      int dsize = (tile->GetTimeRange()-tile->GetCurrentOverlapDuration())*triggers[chanindex]->GetWorkingFrequency(); // size of 'sane' data
      if(spectrumw->LoadData(dsize, offt->GetRe_t(), dstart)) SaveWPSD();
    }
  }
  
  //*** MAPS
  if(fOutProducts.find("map")!=string::npos){
    double snr;// snr max of the full map
    if(fVerbosity>1) cout<<"\t- write maps"<<endl;
    snr=tile->SaveMaps(outdir[chanindex],
		       triggers[chanindex]->GetNameConv()+"_OMICRON",
		       fOutFormat,fWindows,toffset,(bool)(fOutProducts.find("html")+1));
    if(snr>chan_mapsnrmax[chanindex]) chan_mapsnrmax[chanindex]=snr;// get snr max over all chunks
  }

  //*** TRIGGERS
  if(fOutProducts.find("triggers")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write triggers "<<endl;
    if(!ExtractTriggers()){// extract triggers
      chunktfile.push_back("none");
      return false;
    }
    chunktfile.push_back(GetFileName(WriteTriggers())); // write triggers to disk
  }

  // update monitoring segments
  // -- do not include output segment selection --
  outSegments[chanindex]->AddSegment((double)(tile->GetChunkTimeStart()+tile->GetCurrentOverlapDuration()-tile->GetOverlapDuration()/2),(double)(tile->GetChunkTimeEnd()-tile->GetOverlapDuration()/2));

  chan_write_ctr[chanindex]++;
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::ExtractTriggers(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::ExtractTriggers: the Omicron object is corrupted"<<endl;
    return false;
  }

  // trigger rate: evaluated over the chunk excluding nominal overlaps/2
  double trate = (double)trig_ctr[chanindex]/(double)(tile->GetTimeRange()-tile->GetOverlapDuration());

  // check against trigger rate limit
  if(trate>fratemax){
    cerr<<"Omicron::ExtractTriggers: the maximum trigger rate ("<<fratemax<<" Hz) is exceeded ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return false;
  }

  // save tiles above SNR threshold
  if(!tile->SaveTriggers(triggers[chanindex])){
    if(triggers[chanindex]->GetBufferSize()) triggers[chanindex]->ResetBuffer(); // reset buffer
    else triggers[chanindex]->Reset(); // reset TTrees
    return false;
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
string Omicron::WriteTriggers(const bool aLVDirConvention){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::WriteTriggers: the Omicron object is corrupted"<<endl;
    return "none";
  }

  // flush buffer if any
  if(triggers[chanindex]->GetBufferSize()) triggers[chanindex]->FlushBuffer();

  // sort triggers
  if(!triggers[chanindex]->SortTriggers()){
    cerr<<"Omicron::WriteTriggers: triggers cannot be sorted ("<<triggers[chanindex]->GetName()<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return "none";
  }
    
  // clustering if any
  if(fClusterAlgo.compare("none")) triggers[chanindex]->Clusterize();

  // LIGO directory convention
  if(aLVDirConvention){
    int gps_5=0;
    if(triggers[chanindex]->GetNsegments()) gps_5 = (int) triggers[chanindex]->GetStart(0) / 100000;
    stringstream lv_dir;
    lv_dir << maindir << "/" << triggers[chanindex]->GetNamePrefix() << "/" << triggers[chanindex]->GetNameSuffixUnderScore() << "_OMICRON/" << gps_5;
    system(("mkdir -p " + lv_dir.str()).c_str());
    return triggers[chanindex]->Write(lv_dir.str(),fOutFormat);
  }
  
  // write triggers to disk
  return triggers[chanindex]->Write(outdir[chanindex],fOutFormat);
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::Whiten(Spectrum *aSpec){
////////////////////////////////////////////////////////////////////////////////////

  int i=0; // frequency index
 
  // zero-out DC
  offt->SetRe_f(i,0.0);
  offt->SetIm_f(i,0.0);
  i++;
  
  // zero-out below highpass frequency
  int n = (int)(triggers[chanindex]->GetHighPassFrequency()*tile->GetTimeRange());
  for(; i<n; i++){
    offt->SetRe_f(i,0.0);
    offt->SetIm_f(i,0.0);
  }
 
  // normalize data by the ASD
  double asdval;
  for(; i<offt->GetSize_f(); i++){
    asdval=aSpec->GetAmplitude((double)i/(double)tile->GetTimeRange())/sqrt(2.0);
    if(!asdval){
      offt->SetRe_f(i,0.0);
      offt->SetIm_f(i,0.0);
      continue;
    }
    offt->SetRe_f(i,offt->GetRe_f(i) / asdval);
    offt->SetIm_f(i,offt->GetIm_f(i) / asdval);
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveAPSD(const string aType){
////////////////////////////////////////////////////////////////////////////////////

  // normalization factor
  double factor;
  if(!aType.compare("ASD")) factor=sqrt(2.0);
  else                      factor=2.0;
  
  // extract A/PSD 1
  TGraph *GAPSD1;
  if(!aType.compare("ASD")) GAPSD1 = spectrum1[chanindex]->GetASD(tile->GetFrequencyMin(),tile->GetFrequencyMax());
  else                      GAPSD1 = spectrum1[chanindex]->GetPSD(tile->GetFrequencyMin(),tile->GetFrequencyMax());
  if(GAPSD1==NULL) return;
 
  // extract A/PSD 2
  TGraph *GAPSD2;
  if(!aType.compare("ASD")) GAPSD2 = spectrum2[chanindex]->GetASD(tile->GetFrequencyMin(),tile->GetFrequencyMax());
  else                     GAPSD2 = spectrum2[chanindex]->GetPSD(tile->GetFrequencyMin(),tile->GetFrequencyMax());
  if(GAPSD2==NULL) return;
  
  // combine the 2 apsd
  for(int i=0; i<GAPSD2->GetN(); i++) GAPSD2->GetY()[i] *= (GAPSD1->GetY()[i]/factor);
   
  // extract sub-segment A/PSD 1
  TGraph **G1;
  G1 = new TGraph* [spectrum1[chanindex]->GetNSubSegmentsMax(0)+spectrum1[chanindex]->GetNSubSegmentsMax(1)];
  if(!aType.compare("ASD")){
    for(int i=0; i<spectrum1[chanindex]->GetNSubSegmentsMax(0); i++) G1[i] = spectrum1[chanindex]->GetSubASD(0,i);
    for(int i=0; i<spectrum1[chanindex]->GetNSubSegmentsMax(1); i++) G1[spectrum1[chanindex]->GetNSubSegmentsMax(0)+i] = spectrum1[chanindex]->GetSubASD(1,i);
  }
  else{
    for(int i=0; i<spectrum1[chanindex]->GetNSubSegmentsMax(0); i++) G1[i] = spectrum1[chanindex]->GetSubPSD(0,i);
    for(int i=0; i<spectrum1[chanindex]->GetNSubSegmentsMax(1); i++) G1[spectrum1[chanindex]->GetNSubSegmentsMax(0)+i] = spectrum1[chanindex]->GetSubPSD(1,i);
  }
 
  // extract sub-segment A/PSD 2
  TGraph **G2;
  G2 = new TGraph* [spectrum2[chanindex]->GetNSubSegmentsMax(0)+spectrum2[chanindex]->GetNSubSegmentsMax(1)];
  if(!aType.compare("ASD")){
    for(int i=0; i<spectrum2[chanindex]->GetNSubSegmentsMax(0); i++) G2[i] = spectrum2[chanindex]->GetSubASD(0,i);
    for(int i=0; i<spectrum2[chanindex]->GetNSubSegmentsMax(1); i++) G2[spectrum2[chanindex]->GetNSubSegmentsMax(0)+i] = spectrum2[chanindex]->GetSubASD(1,i);
  }
  else{
    for(int i=0; i<spectrum2[chanindex]->GetNSubSegmentsMax(0); i++) G2[i] = spectrum2[chanindex]->GetSubPSD(0,i);
    for(int i=0; i<spectrum2[chanindex]->GetNSubSegmentsMax(1); i++) G2[spectrum2[chanindex]->GetNSubSegmentsMax(0)+i] = spectrum2[chanindex]->GetSubPSD(1,i);
  }

  // combine the 2 apsd
  for(int i=0; i<spectrum1[chanindex]->GetNSubSegmentsMax(0)+spectrum1[chanindex]->GetNSubSegmentsMax(1); i++){
    if(G1[i]==NULL||G2[i]==NULL) continue;
    for(int j=0; j<G1[i]->GetN(); j++) G2[i]->GetY()[j] *= (G1[i]->GetY()[j]/factor);
    delete G1[i];
  }
  delete G1;
 
  // cosmetics
  GPlot->SetLogx(1);
  GPlot->SetLogy(1);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);
  GPlot->Draw(GAPSD2,"AP");
  //for(int i=0; i<spectrum1[chanindex]->GetNSubSegmentsMax(0)+spectrum1[chanindex]->GetNSubSegmentsMax(1); i++){
  for(int i=0; i<spectrum1[chanindex]->GetNSubSegmentsMax(1); i++){
    if(G2[i]!=NULL) GPlot->Draw(G2[i],"Lsame");
  }
  GAPSD2->GetHistogram()->SetXTitle("Frequency [Hz]");
  if(aType.compare("ASD")){
    GAPSD2->GetHistogram()->SetYTitle("Power [Amp^{2}/Hz]");
    GAPSD2->SetTitle((triggers[chanindex]->GetName()+": Power spectrum density").c_str());
  }
  else{
    GAPSD2->GetHistogram()->SetYTitle("Amplitude [Amp/#sqrt{Hz}]");
    GAPSD2->SetTitle((triggers[chanindex]->GetName()+": Amplitude spectrum density").c_str());
  }
  GPlot->Draw(GAPSD2,"PLsame");
  GPlot->RedrawAxis();
  GPlot->RedrawAxis("g");
  GAPSD2->SetLineWidth(1);
  GAPSD2->GetXaxis()->SetTitleOffset(0.9);
  GAPSD2->GetYaxis()->SetTitleOffset(1.03);
  GAPSD2->GetXaxis()->SetLabelSize(0.045);
  GAPSD2->GetYaxis()->SetLabelSize(0.045);
  GAPSD2->GetXaxis()->SetTitleSize(0.045);
  GAPSD2->GetYaxis()->SetTitleSize(0.045);
  GAPSD2->GetYaxis()->SetRangeUser(GAPSD2->GetYaxis()->GetXmin()/100.0,GAPSD2->GetYaxis()->GetXmax()*2.0);
  GAPSD1->SetLineWidth(1);
  GPlot->Draw(GAPSD1,"PLsame");

  
  // set new name
  stringstream ss;
  ss<<aType;
  ss<<"_"<<triggers[chanindex]->GetName()<<"_"<<tile->GetChunkTimeCenter();
  GAPSD2->SetName(ss.str().c_str());
  ss.str(""); ss.clear();

  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fpsd;
    ss<<outdir[chanindex]<<"/"+triggers[chanindex]->GetNameConv()<<"_OMICRON"<<aType<<"-"<<tile->GetChunkTimeStart()<<"-"<<tile->GetTimeRange()<<".root";
    fpsd=new TFile((ss.str()).c_str(),"RECREATE");
    ss.str(""); ss.clear();
    fpsd->cd();
    GAPSD2->Write();
    for(int i=0; i<spectrum1[chanindex]->GetNSubSegmentsMax(0)+spectrum1[chanindex]->GetNSubSegmentsMax(1); i++){
      if(G2[i]!=NULL) G2[i]->Write();
    }
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
      ss<<outdir[chanindex]<<"/"+triggers[chanindex]->GetNameConv()<<"_OMICRON"<<aType<<"-"<<tile->GetChunkTimeStart()<<"-"<<tile->GetTimeRange()<<"."<<form[f];
      GPlot->Print(ss.str().c_str());
      ss.str(""); ss.clear();
    }
  }
  
  form.clear();
  delete GAPSD1;
  delete GAPSD2;
  for(int i=0; i<spectrum1[chanindex]->GetNSubSegmentsMax(0)+spectrum1[chanindex]->GetNSubSegmentsMax(1); i++)
    if(G2[i]!=NULL) delete G2[i];
  delete G2;

  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveWPSD(void){
////////////////////////////////////////////////////////////////////////////////////

  // extract PSD
  TGraph *GPSD;
  GPSD = spectrumw->GetPSD(tile->GetFrequencyMin(),tile->GetFrequencyMax());
  if(GPSD==NULL) return;
  
  // cosmetics
  GPlot->SetLogx(1);
  GPlot->SetLogy(1);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);
  GPlot->Draw(GPSD,"APL");
  GPlot->RedrawAxis();
  GPlot->RedrawAxis("g");
  GPSD->SetLineWidth(1);
  GPSD->GetXaxis()->SetTitleOffset(0.9);
  GPSD->GetYaxis()->SetTitleOffset(1.03);
  GPSD->GetXaxis()->SetLabelSize(0.045);
  GPSD->GetYaxis()->SetLabelSize(0.045);
  GPSD->GetXaxis()->SetTitleSize(0.045);
  GPSD->GetYaxis()->SetTitleSize(0.045);
  GPSD->GetHistogram()->SetXTitle("Frequency [Hz]");
  GPSD->GetHistogram()->SetYTitle("Power [Amp^{2}/Hz]");
  GPSD->SetTitle((triggers[chanindex]->GetName()+": Power spectrum density").c_str());
  

  // set new name
  stringstream ss;
  ss<<"PSDW_"<<triggers[chanindex]->GetName()<<"_"<<tile->GetChunkTimeCenter();
  GPSD->SetName(ss.str().c_str());
  ss.str(""); ss.clear();

  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fpsd;
    ss<<outdir[chanindex]<<"/"+triggers[chanindex]->GetNameConv()<<"_OMICRONWPSD-"<<tile->GetChunkTimeStart()<<"-"<<tile->GetTimeRange()<<".root";
    fpsd=new TFile((ss.str()).c_str(),"RECREATE");
    ss.str(""); ss.clear();
    fpsd->cd();
    GPSD->Write();
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
      ss<<outdir[chanindex]<<"/"+triggers[chanindex]->GetNameConv()<<"_OMICRONWPSD-"<<tile->GetChunkTimeStart()<<"-"<<tile->GetTimeRange()<<"."<<form[f];
      GPlot->Print(ss.str().c_str());
      ss.str(""); ss.clear();
    }
  }
  
  form.clear();
  delete GPSD;

  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveSpectral(void){
////////////////////////////////////////////////////////////////////////////////////

  // locals
  stringstream ss;

  // amplitude graph (without DC)
  TGraph *Gamp   = new TGraph(offt->GetSize_f()-1);
  for(int i=1; i<offt->GetSize_f(); i++)
    Gamp->SetPoint(i-1, (double)i*(double)(triggers[chanindex]->GetWorkingFrequency()/2)/(double)offt->GetSize_f(), offt->GetNorm_f(i));
  ss<<"SpecAmp_"<<triggers[chanindex]->GetName()<<"_"<<tile->GetChunkTimeCenter();
  Gamp->SetName(ss.str().c_str());
  ss.str(""); ss.clear();
  Gamp->SetTitle((triggers[chanindex]->GetName()+": spectral density").c_str());
  Gamp->GetHistogram()->SetXTitle("Frequency [Hz]");
  Gamp->GetHistogram()->SetYTitle("Amplitude /#sqrt{Hz}");
  Gamp->SetLineWidth(1);
  Gamp->GetXaxis()->SetTitleOffset(1.03);
  Gamp->GetXaxis()->SetTitleOffset(1.03);
  Gamp->GetXaxis()->SetLabelSize(0.045);
  Gamp->GetYaxis()->SetLabelSize(0.045);
  Gamp->GetXaxis()->SetTitleSize(0.045);
  Gamp->GetYaxis()->SetTitleSize(0.045);
  Gamp->GetXaxis()->SetLimits(tile->GetFrequencyMin(), tile->GetFrequencyMax());

  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fspec;
    ss<<outdir[chanindex]<<"/"+triggers[chanindex]->GetName()<<"_"<<tile->GetChunkTimeCenter()<<"_spec.root";
    fspec=new TFile((ss.str()).c_str(),"RECREATE");
    ss.str(""); ss.clear();
    fspec->cd();
    Gamp->Write();
    fspec->Close();
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

  // draw
  GPlot->Draw(Gamp,"APL");
  GPlot->SetLogx(1);
  GPlot->SetLogy(1);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);

  // save
  if(form.size()){
    for(int f=0; f<(int)form.size(); f++){
      ss<<outdir[chanindex]<<"/"+triggers[chanindex]->GetName()<<"_"<<tile->GetChunkTimeCenter()<<"_spec."<<form[f];
      GPlot->Print(ss.str().c_str());
      ss.str(""); ss.clear();
    }
  }

  // clean and exit
  form.clear();
  delete Gamp;
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveTS(const bool aWhite){
////////////////////////////////////////////////////////////////////////////////////

  // create graph 
  TGraph *GDATA = new TGraph(offt->GetSize_t());
  if(GDATA==NULL) return;

  stringstream ss;
 
  // whitened data
  if(aWhite){
    ss<<"whitets_"<<triggers[chanindex]->GetName()<<"_"<<tile->GetChunkTimeCenter();
    for(int i=0; i<offt->GetSize_t(); i++) GDATA->SetPoint(i,(double)tile->GetChunkTimeStart()+(double)i/(double)(triggers[chanindex]->GetWorkingFrequency()),offt->GetRe_t(i));
  }
  // conditioned data
  else{
    ss<<"ts_"<<triggers[chanindex]->GetName()<<"_"<<tile->GetChunkTimeCenter();
    for(int i=0; i<offt->GetSize_t(); i++) GDATA->SetPoint(i,(double)tile->GetChunkTimeStart()+(double)i/(double)(triggers[chanindex]->GetWorkingFrequency()),ChunkVect[i]);
  }
 
  // plot name
  GDATA->SetName(ss.str().c_str());
  ss.str(""); ss.clear();
     
  // cosmetics
  GPlot->SetLogx(0);
  GPlot->SetLogy(0);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);
  GDATA->GetHistogram()->SetXTitle("Time [s]");
  GDATA->GetHistogram()->SetYTitle("Amplitude");
  if(aWhite) GDATA->SetTitle((triggers[chanindex]->GetName()+": amplitude whitened data time series").c_str());
  else GDATA->SetTitle((triggers[chanindex]->GetName()+": amplitude conditionned time series").c_str());
  GDATA->SetLineWidth(1);
  GDATA->GetXaxis()->SetNoExponent();
  GDATA->GetXaxis()->SetTitleOffset(1.0);
  GDATA->GetYaxis()->SetTitleOffset(1.03);
  GDATA->GetXaxis()->SetLabelSize(0.045);
  GDATA->GetYaxis()->SetLabelSize(0.045);
  GDATA->GetXaxis()->SetTitleSize(0.045);
  GDATA->GetYaxis()->SetTitleSize(0.045);
  GDATA->GetXaxis()->SetNdivisions(4,5,0);
  GPlot->Draw(GDATA,"APL");
 
  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fdata;
    if(aWhite)
      ss<<outdir[chanindex]<<"/"+triggers[chanindex]->GetNameConv()<<"_OMICRONWHITETS-"<<tile->GetChunkTimeStart()<<"-"<<tile->GetTimeRange()<<".root";
    else
      ss<<outdir[chanindex]<<"/"+triggers[chanindex]->GetNameConv()<<"_OMICRONCONDTS-"<<tile->GetChunkTimeStart()<<"-"<<tile->GetTimeRange()<<".root";
    fdata=new TFile((ss.str()).c_str(),"RECREATE");
    ss.str(""); ss.clear();
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
    
    // zoom
    for(int w=(int)fWindows.size()-1; w>=0; w--){
      GDATA->GetXaxis()->SetLimits(tile->GetChunkTimeCenter()+toffset-(double)fWindows[w]/2.0,tile->GetChunkTimeCenter()+toffset+(double)fWindows[w]/2.0);
      for(int f=0; f<(int)form.size(); f++){
	if(aWhite)
	  ss<<outdir[chanindex]<<"/"+triggers[chanindex]->GetNameConv()<<"_OMICRONWHITETS-"<<tile->GetChunkTimeCenter()<<"-"<<fWindows[w]<<"."<<form[f];
	else
	  ss<<outdir[chanindex]<<"/"+triggers[chanindex]->GetNameConv()<<"_OMICRONCONDTS-"<<tile->GetChunkTimeCenter()<<"-"<<fWindows[w]<<"."<<form[f];
	GPlot->Print(ss.str().c_str());
	ss.str(""); ss.clear();

	if(aWhite)
	  ss<<outdir[chanindex]<<"/th"+triggers[chanindex]->GetNameConv()<<"_OMICRONWHITETS-"<<tile->GetChunkTimeCenter()<<"-"<<fWindows[w]<<"."<<form[f];
	else
	  ss<<outdir[chanindex]<<"/th"+triggers[chanindex]->GetNameConv()<<"_OMICRONCONDTS-"<<tile->GetChunkTimeCenter()<<"-"<<fWindows[w]<<"."<<form[f];
	GPlot->Print(ss.str().c_str(),0.5);
	ss.str(""); ss.clear();
      }
    }
  }
      
  form.clear();
  delete GDATA;
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveSG(void){
////////////////////////////////////////////////////////////////////////////////////

  oinjfile<<fixed<<(double)tile->GetChunkTimeCenter()+oinj->GetTime()<<" ";
  oinjfile<<fixed<<oinj->GetFrequency()<<" ";
  oinjfile<<fixed<<oinj->GetQ()<<" ";
  oinjfile<<scientific<<oinj->GetAmplitude()<<" ";
  oinjfile<<fixed<<oinj->GetPhase()<<" ";
  oinjfile<<fixed<<oinj->GetTrueSNR(spectrum1[chanindex], spectrum2[chanindex])<<" ";
  oinjfile<<fixed<<oinj->GetSigmat()<<" ";
  oinjfile<<fixed<<oinj->GetSigmaf()<<" ";
  oinjfile<<endl;
  
  return;
}
