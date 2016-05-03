//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

ClassImp(Omicron)
  
const string Omicron::colorcode[17] = {"#1b02e6","#0220e6","#0257e6","#028ee6","#01c5e6","#01e6cf","#01e698","#01e660","#01e629","#10e601","#7fe601","#b6e601","#e6de00","#e6a600","#e66f00","#e63800","#e60000"};

////////////////////////////////////////////////////////////////////////////////////
Omicron::Omicron(const string aOptionFile){ 
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
  ReadOptions();

  // output directory
  maindir=fMaindir;
  for(int c=0; c<(int)fChannels.size(); c++) outdir.push_back(fMaindir);

  // load ffl if any
  if(FFL!=NULL){
    status_OK*=FFL->DefineTmpDir(fMaindir);
    status_OK*=FFL->LoadFrameFile();
  }
  if(FFL_inject!=NULL&&FFL_inject!=FFL){
    status_OK*=FFL_inject->DefineTmpDir(fMaindir);
    status_OK*=FFL_inject->LoadFrameFile();
  }

  // sort time windows
  // FIXME: to move in Otile
  std::sort(fWindows.begin(), fWindows.end());

  // spectrum status
  for(int c=0; c<(int)fChannels.size(); c++) status_OK*=spectrum[c]->GetStatus();

  // chunk FFT
  offt = new fft(tile->GetTimeRange()*triggers[0]->GetWorkingFrequency(), fftplan, "r2c");
 
  // data containers
  ChunkVect     = new double    [offt->GetSize_t()];
  TukeyWindow   = GetTukeyWindow(offt->GetSize_t(),
				 tile->GetOverlapDuration()*triggers[0]->GetWorkingFrequency());
 
  // metadata field definition
  vector <string> fOptionName;
  vector <string> fOptionType;
  fOptionName.push_back("omicron_DATA_FFLFILE");              fOptionType.push_back("s");
  fOptionName.push_back("omicron_DATA_CHANNEL");              fOptionType.push_back("s");
  fOptionName.push_back("omicron_DATA_SAMPLEFREQUENCY");      fOptionType.push_back("i");
  fOptionName.push_back("omicron_INJECTION_FFL");             fOptionType.push_back("s");
  fOptionName.push_back("omicron_INJECTION_CHANNEL");         fOptionType.push_back("s");
  fOptionName.push_back("omicron_INJECTION_FACTOR");          fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_FILENAME");        fOptionType.push_back("s");
  fOptionName.push_back("omicron_INJECTION_SG");              fOptionType.push_back("i");
  fOptionName.push_back("omicron_INJECTION_SGTIMEMIN");       fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGTIMEMAX");       fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGFREQUENCYMIN");  fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGFREQUENCYMAX");  fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGQMIN");          fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGQMAX");          fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGAMPMIN");        fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_SGAMPMAX");        fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_FMIN");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_FMAX");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_QMIN");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_QMAX");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_CHUNKDURATION");   fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_OVERLAPDURATION"); fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_MISMATCHMAX");     fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_SNRTHRESHOLD");    fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_PSDLENGTH");       fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_CLUSTERING");      fOptionType.push_back("s");
  fOptionName.push_back("omicron_PARAMETER_CLUSTERDT");       fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_FFTPLAN");         fOptionType.push_back("s");
  fOptionName.push_back("omicron_OUTPUT_DIRECTORY");          fOptionType.push_back("s");
  fOptionName.push_back("omicron_OUTPUT_VERBOSITY");          fOptionType.push_back("i");
  fOptionName.push_back("omicron_OUTPUT_FORMAT");             fOptionType.push_back("s");
  fOptionName.push_back("omicron_OUTPUT_PRODUCTS");           fOptionType.push_back("s");
  fOptionName.push_back("omicron_OUTPUT_STYLE");              fOptionType.push_back("s");

  // triggers metadata
  for(int c=0; c<(int)fChannels.size(); c++){
    triggers[c]->SetMprocessname("OMICRON");
    triggers[c]->SetProcessVersion(GetVersion());
    status_OK*=triggers[c]->InitUserMetaData(fOptionName,fOptionType);
    if(FFL==NULL) status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],"none");
    else          status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],FFL->GetInputFfl());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[1],fChannels[c]);
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
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[24],spectrum[c]->GetDataBufferLength());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[25],fClusterAlgo);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[26],triggers[c]->GetClusterizeDt());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[27],fftplan);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[28],fMaindir);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[29],fVerbosity);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[30],fOutFormat);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[31],fOutProducts);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[32],GPlot->GetCurrentStyle());
  }
  
  // process monitoring
  chanindex      = -1;
  inSegments     = new Segments();
  outSegments    = new Segments* [(int)fChannels.size()];
  chunk_ctr      = 0;
  chan_ctr       = new int       [(int)fChannels.size()];
  chan_data_ctr  = new int       [(int)fChannels.size()];
  chan_cond_ctr  = new int       [(int)fChannels.size()];
  chan_proj_ctr  = new int       [(int)fChannels.size()];
  chan_write_ctr = new int       [(int)fChannels.size()];
  chan_mapsnrmax = new double    [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){
    outSegments[c]    = new Segments();
    chan_ctr[c]       = 0;
    chan_data_ctr[c]  = 0;
    chan_cond_ctr[c]  = 0;
    chan_proj_ctr[c]  = 0;
    chan_write_ctr[c] = 0;
    chan_mapsnrmax[c] = 0.0;
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
  for(int c=0; c<(int)fChannels.size(); c++){
    delete outSegments[c];
    delete triggers[c];
    delete spectrum[c];
  }
  delete outSegments;
  delete spectrum;
  delete triggers;
  if(FFL_inject!=FFL&&FFL_inject!=NULL) delete FFL_inject;
  if(FFL!=NULL) delete FFL;
  if(inject!=NULL){
    for(int c=0; c<(int)fChannels.size(); c++) delete inject[c];
    delete inject;
  }
  delete oinj;
  delete tile;
  delete GPlot;
  delete ChunkVect;
  delete TukeyWindow;
  delete offt;

  chunkcenter.clear();
  outdir.clear();
  fChannels.clear();
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
bool Omicron::MakeDirectories(const int aId){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::MakeDirectories: the Omicron object is corrupted"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::MakeDirectories: make specific directory structure..."<<endl;
  outdir.clear();
  ostringstream tmpstream;

  // Main dir
  if(!aId){// single level
    for(int c=0; c<(int)fChannels.size(); c++){
      maindir=fMaindir;
      outdir.push_back(maindir+"/"+fChannels[c]);
    }
  }
  else{// two level
    tmpstream<<fMaindir<<"/"<<setprecision(3)<<fixed<<aId;
    maindir=tmpstream.str();
    tmpstream.clear(); tmpstream.str("");
  }

  // channel dir
  for(int c=0; c<(int)fChannels.size(); c++){
    outdir.push_back(maindir+"/"+fChannels[c]);
    if(fVerbosity>1) cout<<"\t- "<<outdir[c]<<endl;
    system(("mkdir -p "+outdir[c]).c_str());
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
  if(newseg){
    for(int c=0; c<(int)fChannels.size(); c++) spectrum[c]->Reset();
  }

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
    for(int c=0; c<(int)fChannels.size(); c++) spectrum[c]->Reset();
      
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
    if(fVerbosity>1) cout<<"\t- no more channels to load"<<endl;
    return false; 
  }

  if(fVerbosity>1) cout<<"\t- channel "<<fChannels[chanindex]<<" is loaded"<<endl;
  chan_ctr[chanindex]++;
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
  *aDataVector = FFL->GetData(*aSize, fChannels[chanindex], tile->GetChunkTimeStart(), tile->GetChunkTimeEnd());

  // cannot retrieve data
  if(*aSize<=0){
    cerr<<"Omicron::LoadData: cannot retrieve data ("<<fChannels[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    aDataVector=NULL; *aSize=0;
    return false;
  }

  // test native sampling (and update if necessary)
  int nativesampling = *aSize/(tile->GetTimeRange());
  if(!triggers[chanindex]->SetNativeFrequency(nativesampling)){
    cerr<<"Omicron::LoadData: incompatible native/working frequency ("<<fChannels[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return false;
  }

  // add software injections if any
  if(inject!=NULL){
    if(fVerbosity>1) cout<<"\t- perform software injections..."<<endl;
    inject[chanindex]->UpdateNativeSamplingFrequency();
    if(!inject[chanindex]->Inject(*aSize, *aDataVector, tile->GetChunkTimeStart())){
      cerr<<"Omicron::LoadData: failed to inject ("<<fChannels[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
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
    cerr<<"Omicron::Condition: input vector is NULL ("<<fChannels[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return 1;
  }
  if(aInVectSize<=0){
    cerr<<"Omicron::Condition: input vector is empty ("<<fChannels[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return 2;
  }
  if(aInVect[0]==aInVect[aInVectSize-1]){// FIXME: use mean/min/max
    cerr<<"Omicron::Condition: input vector appears to be flat ("<<fChannels[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return 3;
  }
    
  if(fVerbosity) cout<<"Omicron::Condition: condition data vector..."<<endl;

  // test native sampling (and update if necessary)
  int nativesampling = aInVectSize/(tile->GetTimeRange());
  if(!triggers[chanindex]->SetNativeFrequency(nativesampling)){
    cerr<<"Omicron::Condition: incompatible native/working frequency ("<<fChannels[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return 4;
  }

  // transform data vector
  if(fVerbosity>1) cout<<"\t- transform data vector..."<<endl;
  if(!triggers[chanindex]->Transform(aInVectSize, aInVect, offt->GetSize_t(), ChunkVect)) return 5;

  // apply Tukey Window
  if(fVerbosity>1) cout<<"\t- apply Tukey window..."<<endl;
  for(int i=0; i<offt->GetSize_t(); i++) ChunkVect[i] *= TukeyWindow[i];

  // update spectrum
  if(fVerbosity>1) cout<<"\t- update spectrum..."<<endl;
  if(!spectrum[chanindex]->AddData((tile->GetTimeRange()-tile->GetCurrentOverlapDuration())*triggers[chanindex]->GetWorkingFrequency(), ChunkVect, (tile->GetCurrentOverlapDuration()-tile->GetOverlapDuration()/2)*triggers[chanindex]->GetWorkingFrequency())) return 6;

  // compute tiling power
  if(fVerbosity>1) cout<<"\t- compute tiling power..."<<endl;
  if(!tile->SetPower(spectrum[chanindex])) return 7;

  // fft-forward the chunk data
  if(fVerbosity>1) cout<<"\t- move the data in the frequency domain..."<<endl;
  if(!offt->Forward(ChunkVect)) return 8;
  // note: the FFT normalization is performed in Whiten()

  // whitening
  if(fVerbosity>1) cout<<"\t- whiten chunk"<<endl;
  if(!Whiten()) return 9;

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
  if(!tile->ProjectData(offt)) return false;

  chan_proj_ctr[chanindex]++;
  return true;
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

  //*** SPECTRAL DATA
  /*
  if(fOutProducts.find("spectral")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write spectral data..."<<endl;
    SaveSpectral();
  }
  */

  //*** WHITENED TS
  if(fOutProducts.find("white")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write whitened data..."<<endl;
    offt->Backward();// Back in time domain
    // IMPORTANT: after that, the frequency-domain vector of offt is corrupted (r2c)
    // apply FFT normalization
    for(int i=0; i<offt->GetSize_t(); i++) offt->SetRe_t(i, offt->GetRe_t(i)*triggers[chanindex]->GetWorkingFrequency()/(double)offt->GetSize_t());
    SaveTS(true);
  }
  
  //*** INJECTIONS
  if(fsginj==1&&fOutProducts.find("injection")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write injection data..."<<endl;
    SaveSG();
  }

  //*** MAPS
  if(fOutProducts.find("maps")!=string::npos){
    double snr;
    if(fVerbosity>1) cout<<"\t- write maps"<<endl;
    if(fOutProducts.find("html")==string::npos) // need to make thumbnails for html
      snr=tile->SaveMaps(outdir[chanindex],
			 fChannels[chanindex],
			 fOutFormat,fWindows,false);
    else
      snr=tile->SaveMaps(outdir[chanindex],
			 fChannels[chanindex],
			 fOutFormat,fWindows,true);
    if(snr>chan_mapsnrmax[chanindex]) chan_mapsnrmax[chanindex]=snr;
  }

  //*** TRIGGERS
  if(fOutProducts.find("triggers")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write triggers "<<endl;
    if(!tile->SaveTriggers(triggers[chanindex])) return false;
    if(!WriteTriggers()) return false;
  }

  // update monitoring segments
  // -- do not include output segment selection --
  outSegments[chanindex]->AddSegment((double)(tile->GetChunkTimeStart()+tile->GetCurrentOverlapDuration()-tile->GetOverlapDuration()/2),(double)(tile->GetChunkTimeEnd()-tile->GetOverlapDuration()/2));

  chan_write_ctr[chanindex]++;
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::WriteTriggers(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::WriteTriggers: the Omicron object is corrupted"<<endl;
    return false;
  }

  // sort triggers
  if(!triggers[chanindex]->SortTriggers()){
    cerr<<"Omicron::WriteTriggers: triggers cannot be sorted ("<<fChannels[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return false;
  }
    
  // clustering if any
  if(fClusterAlgo.compare("none")) triggers[chanindex]->Clusterize();
  
  // write triggers to disk
  if(!triggers[chanindex]->Write(outdir[chanindex],fOutFormat).compare("none")){
    cerr<<"Omicron::WriteTriggers: triggers cannot be written to disk ("<<fChannels[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return false;
  }

  return true;
}


////////////////////////////////////////////////////////////////////////////////////
void Omicron::PrintStatusInfo(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK) return;
  if(!inSegments->GetNsegments()) return;
  
  cout<<"\n************* Omicron status info *************"<<endl;
  cout<<"requested start         = "<<(int)inSegments->GetStart(0)<<endl;
  cout<<"requested end           = "<<(int)inSegments->GetEnd(inSegments->GetNsegments()-1)<<endl;
  cout<<"requested livetime      = "<<(int)inSegments->GetLiveTime()<<"s"<<endl;
  cout<<"number of loaded chunks = "<<chunk_ctr<<endl;

  for(int c=0; c<(int)fChannels.size(); c++){
    cout<<"\n*** "<<fChannels[c]<<endl;
    cout<<"number of calls                = "<<chan_ctr[c]<<endl;
    cout<<"number of data calls           = "<<chan_data_ctr[c]<<endl;
    cout<<"number of conditioning calls   = "<<chan_cond_ctr[c]<<endl;
    cout<<"number of projection calls     = "<<chan_proj_ctr[c]<<endl;
    cout<<"number of write calls          = "<<chan_write_ctr[c]<<endl;
    if(outSegments[c]->GetNsegments()){
      cout<<"start_out           = "<<(int)outSegments[c]->GetStart(0)<<endl;
      cout<<"end_out             = "<<(int)outSegments[c]->GetEnd(outSegments[c]->GetNsegments()-1)<<endl;
      cout<<"trigger livetime    = "<<(int)outSegments[c]->GetLiveTime()<<"s ("<<outSegments[c]->GetLiveTime()/inSegments->GetLiveTime()*100<<"%)"<<endl;
    }
  }
  cout<<"***********************************************\n"<<endl;

  return;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::Whiten(void){
////////////////////////////////////////////////////////////////////////////////////

  int i=0; // frequency index
 
  // zero-out DC
  offt->SetRe_f(i,0.0);
  offt->SetIm_f(i,0.0);
  i++;
  
  // zero-out below highpass frequency
  int n = (int)(tile->GetFrequencyMin()*tile->GetTimeRange());
  for(; i<n; i++){
    offt->SetRe_f(i,0.0);
    offt->SetIm_f(i,0.0);
  }
 
  // normalize data by the ASD
  double asdval;
  for(; i<offt->GetSize_f(); i++){
    asdval=spectrum[chanindex]->GetPower((double)i/(double)tile->GetTimeRange())/2.0;
    if(asdval<=0){
      cerr<<"Omicron::Whiten: could not retrieve power for f="<<(double)i/(double)tile->GetTimeRange()<<" Hz"<<endl;
      return false;
    }
    asdval=sqrt(asdval);
    offt->SetRe_f(i,offt->GetRe_f(i) / asdval /triggers[chanindex]->GetWorkingFrequency());
    offt->SetIm_f(i,offt->GetIm_f(i) / asdval /triggers[chanindex]->GetWorkingFrequency());
    // /f_w is the FFT normalization, which was not performed in Condition()
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveAPSD(const string aType){
////////////////////////////////////////////////////////////////////////////////////

  // extract A/PSD 
  TGraph *GAPSD;
  if(!aType.compare("ASD")) GAPSD = spectrum[chanindex]->GetASD(tile->GetFrequencyMin(),tile->GetFrequencyMax());
  else                     GAPSD = spectrum[chanindex]->GetPSD(tile->GetFrequencyMin(),tile->GetFrequencyMax());
  if(GAPSD==NULL) return;

  // extract sub-segment A/PSD
  TGraph **G;
  G = new TGraph* [spectrum[chanindex]->GetNSubSegmentsMax(0)+spectrum[chanindex]->GetNSubSegmentsMax(1)];
  if(!aType.compare("ASD")){
    for(int i=0; i<spectrum[chanindex]->GetNSubSegmentsMax(0); i++) G[i] = spectrum[chanindex]->GetSubASD(0,i);
    for(int i=0; i<spectrum[chanindex]->GetNSubSegmentsMax(1); i++) G[spectrum[chanindex]->GetNSubSegmentsMax(0)+i] = spectrum[chanindex]->GetSubASD(1,i);
  }
  else{
    for(int i=0; i<spectrum[chanindex]->GetNSubSegmentsMax(0); i++) G[i] = spectrum[chanindex]->GetSubPSD(0,i);
    for(int i=0; i<spectrum[chanindex]->GetNSubSegmentsMax(1); i++) G[spectrum[chanindex]->GetNSubSegmentsMax(0)+i] = spectrum[chanindex]->GetSubPSD(1,i);
  }

  // cosmetics
  GPlot->SetLogx(1);
  GPlot->SetLogy(1);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);
  GPlot->Draw(GAPSD,"AP");
  for(int i=0; i<spectrum[chanindex]->GetNSubSegmentsMax(0)+spectrum[chanindex]->GetNSubSegmentsMax(1); i++){
    if(G[i]!=NULL) GPlot->Draw(G[i],"Lsame");
  }
  GAPSD->GetHistogram()->SetXTitle("Frequency [Hz]");
  if(aType.compare("ASD")){
    GAPSD->GetHistogram()->SetYTitle("Power [Amp^{2}/Hz]");
    GAPSD->SetTitle((fChannels[chanindex]+": Power spectrum density").c_str());
  }
  else{
    GAPSD->GetHistogram()->SetYTitle("Amplitude [Amp/#sqrt{Hz}]");
    GAPSD->SetTitle((fChannels[chanindex]+": Amplitude spectrum density").c_str());
  }
  GPlot->Draw(GAPSD,"PLsame");
  GPlot->RedrawAxis();
  GPlot->RedrawAxis("g");
  GAPSD->SetLineWidth(1);
  GAPSD->GetXaxis()->SetTitleOffset(1.1);
  GAPSD->GetXaxis()->SetLabelSize(0.045);
  GAPSD->GetYaxis()->SetLabelSize(0.045);
  GAPSD->GetXaxis()->SetTitleSize(0.045);
  GAPSD->GetYaxis()->SetTitleSize(0.045);
  //GAPSD->GetYaxis()->SetRangeUser(1e-25,1e-21);
  
  // set new name
  stringstream ss;
  ss<<aType;
  ss<<"_"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter();
  GAPSD->SetName(ss.str().c_str());
  ss.str(""); ss.clear();

  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fpsd;
    ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_"<<aType;
    ss<<".root";
    fpsd=new TFile((ss.str()).c_str(),"RECREATE");
    ss.str(""); ss.clear();
    fpsd->cd();
    GAPSD->Write();
    for(int i=0; i<spectrum[chanindex]->GetNSubSegmentsMax(0)+spectrum[chanindex]->GetNSubSegmentsMax(1); i++){
      if(G[i]!=NULL) G[i]->Write();
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
      ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_"<<aType;
      ss<<"."<<form[f];
      GPlot->Print(ss.str().c_str());
      ss.str(""); ss.clear();
    }
  }
  
  form.clear();
  delete GAPSD;
  for(int i=0; i<spectrum[chanindex]->GetNSubSegmentsMax(0)+spectrum[chanindex]->GetNSubSegmentsMax(1); i++)
    if(G[i]!=NULL) delete G[i];
  delete G;

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
  ss<<"SpecAmp_"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter();
  Gamp->SetName(ss.str().c_str());
  ss.str(""); ss.clear();
  Gamp->SetTitle((fChannels[chanindex]+": spectral density").c_str());
  Gamp->GetHistogram()->SetXTitle("Frequency [Hz]");
  Gamp->GetHistogram()->SetYTitle("Amplitude /#sqrt{Hz}");
  Gamp->SetLineWidth(1);
  Gamp->GetXaxis()->SetTitleOffset(1.1);
  Gamp->GetXaxis()->SetLabelSize(0.045);
  Gamp->GetYaxis()->SetLabelSize(0.045);
  Gamp->GetXaxis()->SetTitleSize(0.045);
  Gamp->GetYaxis()->SetTitleSize(0.045);
  Gamp->GetXaxis()->SetLimits(tile->GetFrequencyMin(), tile->GetFrequencyMax());

  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fspec;
    ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_spec.root";
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
      ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_spec."<<form[f];
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
    ss<<"whitets_"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter();
    for(int i=0; i<offt->GetSize_t(); i++) GDATA->SetPoint(i,(double)tile->GetChunkTimeStart()+(double)i/(double)(triggers[chanindex]->GetWorkingFrequency()),offt->GetRe_t(i));
  }
  // conditioned data
  else{
    ss<<"ts_"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter();
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
  if(aWhite) GDATA->SetTitle((fChannels[chanindex]+": amplitude whitened data time series").c_str());
  else GDATA->SetTitle((fChannels[chanindex]+": amplitude conditionned time series").c_str());
  GDATA->SetLineWidth(1);
  GDATA->GetXaxis()->SetNoExponent();
  GDATA->GetXaxis()->SetTitleOffset(1.1);
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
      ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_whitets.root";
    else
      ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_ts.root";
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
      GDATA->GetXaxis()->SetLimits(tile->GetChunkTimeCenter()-(double)fWindows[w]/2.0,tile->GetChunkTimeCenter()+(double)fWindows[w]/2.0);
      for(int f=0; f<(int)form.size(); f++){
	if(aWhite)
	  ss<<outdir[chanindex]<<"/"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_whitetsdt"<<fWindows[w]<<"."<<form[f];
	else
	  ss<<outdir[chanindex]<<"/"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_tsdt"<<fWindows[w]<<"."<<form[f];
	GPlot->Print(ss.str().c_str());
	ss.str(""); ss.clear();

	if(aWhite)
	  ss<<outdir[chanindex]<<"/"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_whitetsdt"<<fWindows[w]<<"th."<<form[f];
	else
	  ss<<outdir[chanindex]<<"/"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_tsdt"<<fWindows[w]<<"th."<<form[f];
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

  // text file only
  stringstream ss;
  ss<<outdir[chanindex]<<"/"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_sginjection.txt";

  ofstream outfile(ss.str().c_str());
  outfile.precision(5);
  outfile<<"Oinject: sinusoidal Gaussian waveform:"<<endl;
  outfile<<"         peak time = "<<fixed<<(double)tile->GetChunkTimeCenter()+oinj->GetTime()<<endl;
  outfile<<"         frequency = "<<fixed<<oinj->GetFrequency()<<endl;
  outfile<<"         Q         = "<<fixed<<oinj->GetQ()<<endl;
  outfile<<"         amplitude = "<<scientific<<oinj->GetAmplitude()<<endl;
  outfile<<"         phase     = "<<fixed<<oinj->GetPhase()<<endl;
  outfile<<"         SNR       = "<<fixed<<oinj->GetTrueSNR(spectrum[chanindex])<<endl;
  outfile<<"         sigma_t   = "<<fixed<<oinj->GetSigmat()<<endl;
  outfile<<"         sigma_f   = "<<fixed<<oinj->GetSigmaf()<<endl;
  outfile.close();
  return;
}
