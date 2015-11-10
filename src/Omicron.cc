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

  // load ffl if any (CHECKME: can we load lcf directly?)
  if(FFL!=NULL){
    status_OK*=FFL->DefineTmpDir(fMaindir);// working directory (for LCF conversion)
    status_OK*=FFL->LoadFrameFile();
  }

  // sort time windows
  // FIXME: to move in Otile
  std::sort(fWindows.begin(), fWindows.end());

  // data spectrum
  if(tile->GetFrequencyMin()>1.0) // resolution = 0.5 Hz above 1 Hz
    spectrum = new Spectrum(triggers[0]->GetWorkingFrequency(),tile->GetTimeRange()-tile->GetOverlapDuration(),triggers[0]->GetWorkingFrequency(),fVerbosity);
  else // increase the resolution not to extrapolate the PSD.
    spectrum = new Spectrum(2*(int)floor((double)triggers[0]->GetWorkingFrequency()/tile->GetFrequencyMin()),tile->GetTimeRange()-tile->GetOverlapDuration(),triggers[0]->GetWorkingFrequency(),fVerbosity);
  status_OK*=spectrum->GetStatus();

  // data containers
  ChunkVect     = new double    [tile->GetTimeRange()*triggers[0]->GetWorkingFrequency()];
  WhiteChunkVect= new double    [tile->GetTimeRange()*triggers[0]->GetWorkingFrequency()];
  TukeyWindow   = GetTukeyWindow(tile->GetTimeRange()*triggers[0]->GetWorkingFrequency(),
				 tile->GetOverlapDuration()*triggers[0]->GetWorkingFrequency());
  offt          = new fft(       tile->GetTimeRange()*triggers[0]->GetWorkingFrequency(),
				 "FFTW_MEASURE");
  for(int i=0; i<tile->GetTimeRange()*triggers[0]->GetWorkingFrequency(); i++) WhiteChunkVect[i]=0.0;

  // metadata field definition
  vector <string> fOptionName;  ///< option name (metadata)
  vector <string> fOptionType;  ///< option type (metadata)
  fOptionName.push_back("omicron_DATA_FFLFILE");              fOptionType.push_back("s");
  fOptionName.push_back("omicron_DATA_CHANNEL");              fOptionType.push_back("s");
  fOptionName.push_back("omicron_DATA_SAMPLEFREQUENCY");      fOptionType.push_back("i");
  fOptionName.push_back("omicron_INJECTION_CHANNEL");         fOptionType.push_back("s");
  fOptionName.push_back("omicron_INJECTION_FACTOR");          fOptionType.push_back("d");
  fOptionName.push_back("omicron_INJECTION_FILENAME");        fOptionType.push_back("s");
  fOptionName.push_back("omicron_PARAMETER_FMIN");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_FMAX");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_QMIN");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_QMAX");            fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_CHUNKDURATION");   fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_OVERLAPDURATION"); fOptionType.push_back("i");
  fOptionName.push_back("omicron_PARAMETER_MISMATCHMAX");     fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_SNRTHRESHOLD");    fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_CLUSTERING");      fOptionType.push_back("s");
  fOptionName.push_back("omicron_PARAMETER_CLUSTERDT");       fOptionType.push_back("d");
  fOptionName.push_back("omicron_PARAMETER_TILEDOWN");        fOptionType.push_back("i");
  fOptionName.push_back("omicron_OUTPUT_DIRECTORY");          fOptionType.push_back("s");
  fOptionName.push_back("omicron_OUTPUT_NTRIGGERMAX");        fOptionType.push_back("i");
  fOptionName.push_back("omicron_OUTPUT_VERBOSITY");          fOptionType.push_back("i");
  fOptionName.push_back("omicron_OUTPUT_FORMAT");             fOptionType.push_back("s");
  fOptionName.push_back("omicron_OUTPUT_PRODUCTS");           fOptionType.push_back("s");
  fOptionName.push_back("omicron_OUTPUT_STYLE");              fOptionType.push_back("s");

  // triggers metadata
  for(int c=0; c<(int)fChannels.size(); c++){
    triggers[c]->SetProcessVersion(GetVersion());
    status_OK*=triggers[c]->InitUserMetaData(fOptionName,fOptionType);
    if(FFL==NULL) status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],"none");
    else          status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],FFL->GetInputFfl());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[1],fChannels[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[2],triggers[c]->GetWorkingFrequency());
    if(fInjChan.size()){
      status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],fInjChan[c]);
      status_OK*=triggers[c]->SetUserMetaData(fOptionName[4],fInjFact[c]);
    }
    else{
      status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],"none");
      status_OK*=triggers[c]->SetUserMetaData(fOptionName[4],0.0);
    }
    if(inject==NULL) status_OK*=triggers[c]->SetUserMetaData(fOptionName[5],"none");
    else             status_OK*=triggers[c]->SetUserMetaData(fOptionName[5],inject[c]->GetInputFilePattern());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[6],tile->GetFrequencyMin());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[7],tile->GetFrequencyMax());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[8],tile->GetQ(0));
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[9],tile->GetQ(tile->GetNQ()-1));
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[10],tile->GetTimeRange());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[11],tile->GetOverlapDuration());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[12],tile->GetMismatchMax());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[13],tile->GetSNRTriggerThr());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[14],fClusterAlgo);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[15],triggers[c]->GetClusterizeDt());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[16],fTileDown);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[17],fMaindir);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[18],tile->GetNTriggerMax());
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[19],fVerbosity);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[20],fOutFormat);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[21],fOutProducts);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[22],GPlot->GetCurrentStyle());
    triggers[c]->SetMprocessname("Omicron");
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
  }
  delete outSegments;
  delete spectrum;
  delete triggers;
  if(FFL!=NULL) delete FFL;
  if(inject!=NULL){
    for(int c=0; c<(int)fChannels.size(); c++) delete inject[c];
    delete inject;
  }
  delete tile;
  delete GPlot;
  delete ChunkVect;
  delete WhiteChunkVect;
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
bool Omicron::InitSegments(Segments *aSeg){
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

  // input segment monitor
  // (try to optimize the update of inSegments)
  if(inSegments->GetLiveTime()&&aSeg->GetStart(0)>=inSegments->GetStart(inSegments->GetNsegments()-1))
    inSegments->Append(aSeg);
  else
    for(int s=0; s<aSeg->GetNsegments(); s++) inSegments->AddSegment(aSeg->GetStart(s),aSeg->GetEnd(s));
  
  // data structure
  if(!tile->SetSegments(aSeg)){
    cerr<<"Omicron::InitSegments: cannot initiate data segments."<<endl;
    return false;
  }

  // update channel list
  if(FFL!=NULL){
    if(!FFL->ExtractChannels(aSeg->GetStart(0))){
      cerr<<"Omicron::InitSegments: cannot update FFL info."<<endl;
      return false;
    }
  }
  
  if(fVerbosity>1){
    cout<<"\t- N segments       = "<<aSeg->GetNsegments()<<endl;
    cout<<"\t- Livetime         = "<<aSeg->GetLiveTime()<<endl;
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
    triggers[c]->SetOutputDirectory(outdir[c]);
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
  if(fVerbosity) cout<<"Omicron::NewChunk: call a new chunk..."<<endl;
  if(!tile->NewChunk()) return false;
  if(fVerbosity>1) cout<<"\t- chunk "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<" is loaded"<<endl;

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
    double *dvector_inj = FFL->GetData(dsize_inj, fInjChan[chanindex], tile->GetChunkTimeStart(), tile->GetChunkTimeEnd());

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
  if(aInVect[0]==aInVect[aInVectSize-1]){
    cerr<<"Omicron::Condition: input vector appears to be flat ("<<fChannels[chanindex]<<" "<<tile->GetChunkTimeStart()<<"-"<<tile->GetChunkTimeEnd()<<")"<<endl;
    return 3;
  }
    
  if(fVerbosity) cout<<"Omicron::Condition: condition data vector..."<<endl;


  // transform data vector
  if(fVerbosity>1) cout<<"\t- transform data vector..."<<endl;
  if(!triggers[chanindex]->Transform(aInVectSize, aInVect, tile->GetTimeRange()*triggers[chanindex]->GetWorkingFrequency(), ChunkVect)) return 4;

  // apply Tukey Window
  if(fVerbosity>1) cout<<"\t- apply Tukey window..."<<endl;
  int chunksize=tile->GetTimeRange()*triggers[chanindex]->GetWorkingFrequency();
  for(int i=0; i<chunksize; i++) ChunkVect[i] *= TukeyWindow[i];

  // Make spectrum
  if(fVerbosity>1) cout<<"\t- update spectrum..."<<endl;
  if(!spectrum->LoadData((tile->GetTimeRange()-tile->GetOverlapDuration())*triggers[chanindex]->GetWorkingFrequency(), ChunkVect, tile->GetOverlapDuration()/2)) return 5;

  // compute tiling power
  if(fVerbosity>1) cout<<"\t- compute tiling power..."<<endl;
  if(!tile->SetPower(spectrum)) return 6;

  // whiten data vector
  if(fVerbosity>1) cout<<"\t- whiten data vector..."<<endl;

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
  
  //locals
  double *dataRe;
  double *dataIm;

  // whitening
  if(fVerbosity>1) cout<<"\t- whiten chunk"<<endl;
  if(!Whiten(&dataRe, &dataIm)) return false;
    
  // project data
  if(fVerbosity>1) cout<<"\t- project chunk"<<endl;
  if(!tile->ProjectData(dataRe, dataIm, fTileDown)){
    delete dataRe; delete dataIm;
    return false;
  }

  // save whiten data for condition data products
  if(fOutProducts.find("white")!=string::npos){
    if(fVerbosity>1) cout<<"\t- save whitened data"<<endl;
    if(!offt->Backward(dataRe, dataIm)){// Back in time domain
      delete dataRe; delete dataIm;
      return false;
    }
    int chunksize=tile->GetTimeRange()*triggers[chanindex]->GetWorkingFrequency();
    for(int i=0; i<chunksize; i++) WhiteChunkVect[i]=2.0*offt->GetReOut(i)/(double)chunksize;
  }
  
  // not used anymore
  delete dataRe; delete dataIm;

  // save triggers
  if(fOutProducts.find("triggers")!=string::npos){
    if(fVerbosity>1) cout<<"\t- save triggers"<<endl;
    if(!tile->SaveTriggers(triggers[chanindex])) return false;
  }
  
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

  //*** WHITENED DATA
  if(fOutProducts.find("white")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write whitened data..."<<endl;
    SaveTS(true);// save whitened time series
    //if(!spectrum->LoadData((tile->GetCurrentChunkDuration()-tile->GetOverlapDuration())*triggers[chanindex]->GetWorkingFrequency(), CondChunkVect, tile->GetOverlapDuration()/2*triggers[chanindex]->GetWorkingFrequency())) return false;// compute spectrum excluding both chunk ends ov/2
    //SaveAPSD("ASD",true);// save whiten PSD
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

  // save info for html report
  if(fOutProducts.find("html")!=string::npos) chunkcenter.push_back(tile->GetChunkTimeCenter());

  //*** TRIGGERS
  if(fOutProducts.find("triggers")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write triggers "<<endl;
    if(!WriteTriggers()) return false;
  }

  // update monitoring segments
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
  if(!triggers[chanindex]->Write().compare("none")){
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
bool Omicron::Whiten(double **aDataRe, double **aDataIm){
////////////////////////////////////////////////////////////////////////////////////
  // fft-forward
  if(!offt->Forward(ChunkVect)){
    cerr<<"Omicron::Whiten: FFT-forward failed"<<endl;
    return false;
  }

  // NOTE: no FFT normalization when forwarding
  // the normalization (/N) is performed (by convention) after backwarding

  // get ffted data
  *aDataRe=offt->GetRe();
  *aDataIm=offt->GetIm();
  if(*aDataRe==NULL||*aDataIm==NULL){
    cerr<<"Omicron::Whiten: cannot retrieve FFT data"<<endl;
    return false;
  }
 
  int i=0;
 
  // zero-out below highpass frequency
  int n = (int)(tile->GetFrequencyMin()*tile->GetTimeRange());
  for(; i<n; i++){
    (*aDataRe)[i] = 0.0;
    (*aDataIm)[i] = 0.0;
  }
 
  // normalize data by the ASD
  n = tile->GetTimeRange()*triggers[chanindex]->GetWorkingFrequency()/2;
  double asdval;
  for(; i<n; i++){
    asdval=spectrum->GetPower((double)i/(double)tile->GetTimeRange());
    if(asdval<=0){
      cerr<<"Omicron::Whiten: could not retrieve power"<<endl;
      delete *aDataRe; delete *aDataIm;
      *aDataRe=NULL; *aDataIm=NULL;
      return false;
    }
    asdval=sqrt(asdval);
    (*aDataRe)[i] /= asdval;
    (*aDataIm)[i] /= asdval;

  }
 
  // zero-out negative frequencies
  n*=2;
  for(; i<n; i++){
    (*aDataRe)[i] = 0.0;
    (*aDataIm)[i] = 0.0;
  }
 
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveAPSD(const string type, const bool aWhite){
////////////////////////////////////////////////////////////////////////////////////

  // extract A/PSD 
  TGraph *GAPSD;
  if(!type.compare("ASD")) GAPSD = spectrum->GetASD(tile->GetFrequencyMin(),tile->GetFrequencyMax());
  else                     GAPSD = spectrum->GetPSD(tile->GetFrequencyMin(),tile->GetFrequencyMax());
  if(GAPSD==NULL) return;

  // extract sub-segment A/PSD
  TGraph **G;
  G = new TGraph* [spectrum->GetNSubSegmentsMax(0)+spectrum->GetNSubSegmentsMax(1)];
  if(!type.compare("ASD")){
    for(int i=0; i<spectrum->GetNSubSegmentsMax(0); i++) G[i] = spectrum->GetSubASD(0,i);
    for(int i=0; i<spectrum->GetNSubSegmentsMax(1); i++) G[spectrum->GetNSubSegmentsMax(0)+i] = spectrum->GetSubASD(1,i);
  }
  else{
    for(int i=0; i<spectrum->GetNSubSegmentsMax(0); i++) G[i] = spectrum->GetSubPSD(0,i);
    for(int i=0; i<spectrum->GetNSubSegmentsMax(1); i++) G[spectrum->GetNSubSegmentsMax(0)+i] = spectrum->GetSubPSD(1,i);
  }

  // cosmetics
  GPlot->SetLogx(1);
  GPlot->SetLogy(1);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);
  GPlot->Draw(GAPSD,"AP");
  for(int i=0; i<spectrum->GetNSubSegmentsMax(0)+spectrum->GetNSubSegmentsMax(1); i++){
    if(G[i]!=NULL) GPlot->Draw(G[i],"Lsame");
  }
  GAPSD->GetHistogram()->SetXTitle("Frequency [Hz]");
  if(type.compare("ASD")){
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
  if(aWhite) GAPSD->SetTitle(((string)GAPSD->GetTitle()+" (after whitening)").c_str());
  GAPSD->SetLineWidth(1);
  GAPSD->GetXaxis()->SetTitleOffset(1.1);
  GAPSD->GetXaxis()->SetLabelSize(0.045);
  GAPSD->GetYaxis()->SetLabelSize(0.045);
  GAPSD->GetXaxis()->SetTitleSize(0.045);
  GAPSD->GetYaxis()->SetTitleSize(0.045);
  //GAPSD->GetYaxis()->SetRangeUser(1e-25,1e-21);
  
  // set new name
  stringstream ss;
  ss<<type;
  if(aWhite) ss<<"white";
  ss<<"_"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter();
  GAPSD->SetName(ss.str().c_str());
  ss.str(""); ss.clear();

  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fpsd;
    ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_"<<type;
    if(aWhite) ss<<"white";
    ss<<".root";
    fpsd=new TFile((ss.str()).c_str(),"RECREATE");
    ss.str(""); ss.clear();
    fpsd->cd();
    GAPSD->Write();
    for(int i=0; i<spectrum->GetNSubSegmentsMax(0)+spectrum->GetNSubSegmentsMax(1); i++){
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
      ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter()<<"_"<<type;
      if(aWhite) ss<<"white";
      ss<<"."<<form[f];
      GPlot->Print(ss.str().c_str());
      ss.str(""); ss.clear();
    }
  }
  
  form.clear();
  delete GAPSD;
  for(int i=0; i<spectrum->GetNSubSegmentsMax(0)+spectrum->GetNSubSegmentsMax(1); i++)
    if(G[i]!=NULL) delete G[i];
  delete G;

  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveTS(const bool aWhite){
////////////////////////////////////////////////////////////////////////////////////

  // create graph 
  TGraph *GDATA = new TGraph(triggers[chanindex]->GetWorkingFrequency()*tile->GetTimeRange());
  if(GDATA==NULL) return;

  stringstream ss;
 
  // whitened data
  if(aWhite){
    ss<<"whitets_"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter();
    for(int i=0; i<triggers[chanindex]->GetWorkingFrequency()*tile->GetTimeRange(); i++) GDATA->SetPoint(i,(double)tile->GetChunkTimeStart()+(double)i/(double)(triggers[chanindex]->GetWorkingFrequency()),WhiteChunkVect[i]);
  }
  // conditioned data
  else{
    ss<<"ts_"<<fChannels[chanindex]<<"_"<<tile->GetChunkTimeCenter();
    for(int i=0; i<triggers[chanindex]->GetWorkingFrequency()*tile->GetTimeRange(); i++) GDATA->SetPoint(i,(double)tile->GetChunkTimeStart()+(double)i/(double)(triggers[chanindex]->GetWorkingFrequency()),ChunkVect[i]);
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

