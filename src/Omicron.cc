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

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------                OPTIONS                 --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(fVerbosity) cout<<"Omicron::Omicron: init options..."<<endl;

  // input
  fOptionFile=aOptionFile;

  // metadata field definition
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
  fOptionName.push_back("omicron_PARAMETER_SEGMENTDURATION"); fOptionType.push_back("i");
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

  // read option file
  ReadOptions();

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------                 OUTPUT                 --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // output directory
  maindir=fMaindir;
  for(int c=0; c<(int)fChannels.size(); c++) outdir.push_back(fMaindir);
  
  // plotting
  GPlot = new GwollumPlot ("Omicron",fOutStyle);
  fOutStyle = GPlot->GetCurrentStyle();

  // triggers
  if(fVerbosity) cout<<"Omicron::Omicron: init triggers..."<<endl;
  triggers = new MakeTriggers* [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){
    triggers[c] = new MakeTriggers(outdir[c],fChannels[c],fOutFormat,fVerbosity);
    status_OK*=triggers[c]->SetFrequencies(fSampleFrequency,fSampleFrequency,fFreqRange[0]);// default
    triggers[c]->SetClusterizeDt(fcldt);
    triggers[c]->SetProcessVersion(GetVersion());
    status_OK*=triggers[c]->InitUserMetaData(fOptionName,fOptionType);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],fFflFile);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[1],fChannels[c]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[2],triggers[c]->GetWorkingFrequency());
    if(fInjChan.size()) status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],fInjChan[c]);
    else status_OK*=triggers[c]->SetUserMetaData(fOptionName[3],"none");
    if(fInjChan.size()) status_OK*=triggers[c]->SetUserMetaData(fOptionName[4],fInjFact[c]);
    else status_OK*=triggers[c]->SetUserMetaData(fOptionName[4],0.0);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[5],fInjFilePat);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[6],fFreqRange[0]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[7],fFreqRange[1]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[8],fQRange[0]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[9],fQRange[1]);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[10],fChunkDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[11],fSegmentDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[12],fOverlapDuration);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[13],fMismatchMax);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[14],fSNRThreshold_trigger);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[15],fClusterAlgo);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[16],fcldt);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[17],fTileDown);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[18],fMaindir);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[19],fNTriggerMax);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[20],fVerbosity);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[21],fOutFormat);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[22],fOutProducts);
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[23],fOutStyle);
    triggers[c]->SetMprocessname("Omicron");
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------               INPUT DATA               --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // data FFL (if any)
  FFL=NULL; // no input data
  if(fFflFile.compare("none")){
    if(fVerbosity) cout<<"Omicron::Omicron: init data ffl..."<<endl;
    FFL = new ffl(fFflFile, fOutStyle, fVerbosity);
    status_OK*=FFL->DefineTmpDir(fMaindir);
    status_OK*=FFL->LoadFrameFile();
  }

  // software injection
  inject=NULL;
  if(fInjFilePat.compare("none")){
    if(fVerbosity) cout<<"Omicron::Omicron: init software injection..."<<endl;
    inject = new InjEct* [(int)fChannels.size()];
    for(int c=0; c<(int)fChannels.size(); c++){
      inject[c] = new InjEct(triggers[c],fInjFilePat,fVerbosity);
      status_OK*=inject[c]->GetStatus();
    }
  }
   
  // data sequence
  if(fVerbosity) cout<<"Omicron::Omicron: init data sequence..."<<endl;
  dataseq          = new Odata(fChunkDuration, fSegmentDuration, fOverlapDuration, fVerbosity);
  fChunkDuration   = dataseq->GetChunkDuration();
  fSegmentDuration = dataseq->GetSegmentDuration();
  fOverlapDuration = dataseq->GetOverlapDuration();

  // data Spectrum
  if(fVerbosity) cout<<"Omicron::Omicron: init data spectrum..."<<endl;
  int psdsize;
  if(fFreqRange[0]>=1.0)
    psdsize = triggers[0]->GetWorkingFrequency(); // 0.5Hz binning
  else
    psdsize=NextPowerOfTwo(10.0*floor((double)triggers[0]->GetWorkingFrequency()/fFreqRange[0]));// over-binning (>factor 20)
  spectrum = new Spectrum(triggers[0]->GetWorkingFrequency(),psdsize,0,fVerbosity);
  status_OK*=spectrum->GetStatus();

  // data container
  if(fVerbosity) cout<<"Omicron::Omicron: init data container..."<<endl;
  ChunkVect     = new double [fChunkDuration*triggers[0]->GetWorkingFrequency()];
  CondChunkVect = new double [fChunkDuration*triggers[0]->GetWorkingFrequency()];
  SegVect       = new double [fSegmentDuration*triggers[0]->GetWorkingFrequency()];
  TukeyWindow   = GetTukeyWindow(fSegmentDuration*triggers[0]->GetWorkingFrequency(),fOverlapDuration*triggers[0]->GetWorkingFrequency());
  offt          = new fft(fSegmentDuration*triggers[0]->GetWorkingFrequency(),"FFTW_ESTIMATE");
  dataRe        = new double* [dataseq->GetNSegments()];
  dataIm        = new double* [dataseq->GetNSegments()];
  for(int i=0; i<fChunkDuration*triggers[0]->GetWorkingFrequency(); i++)
    CondChunkVect[i]=0.0;

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------                INTERNALS               --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // process monitoring
  if(fVerbosity) cout<<"Omicron::Omicron: init process monitors..."<<endl;
  chanindex      = -1;
  timeoffset     = 0.0;
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
  
  // adjust plot time windows
  if(!fWindows.size()) fWindows.push_back(fSegmentDuration-fOverlapDuration);
  std::sort(fWindows.begin(), fWindows.end());


  // tiling
  if(fVerbosity) cout<<"Omicron::Omicron: init tiling..."<<endl;
  tile = new Otile(fSegmentDuration,fQRange[0],fQRange[1],fFreqRange[0],fFreqRange[1],triggers[0]->GetWorkingFrequency(),fMismatchMax,fOutStyle,fVerbosity);
  tile->SetSNRScale(fsnrscale);
  tile->SetSaveSelection(fSNRThreshold_map,fSNRThreshold_trigger,fNTriggerMax);

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
  delete dataRe;
  delete dataIm;
  delete spectrum;
  delete triggers;
  if(FFL!=NULL) delete FFL;
  if(inject!=NULL){
    for(int c=0; c<(int)fChannels.size(); c++) delete inject[c];
    delete inject;
  }
  delete dataseq;
  delete tile;
  delete GPlot;
  delete ChunkVect;
  delete CondChunkVect;
  delete SegVect;
  delete TukeyWindow;
  delete offt;

  mapcenter.clear();
  chunkstart.clear();
  chunkstop.clear();
  fOptionName.clear();
  fOptionType.clear();
  outdir.clear();
  fChannels.clear();
  fInjChan.clear();
  fInjFact.clear();
  fFreqRange.clear();
  fQRange.clear();
  fWindows.clear();
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

  // input segment monitor
  // (try to optimize the update of inSegments)
  if(inSegments->GetLiveTime()&&aSeg->GetStart(0)>=inSegments->GetStart(inSegments->GetNsegments()-1))
    inSegments->Append(aSeg);
  else
    for(int s=0; s<aSeg->GetNsegments(); s++) inSegments->AddSegment(aSeg->GetStart(s),aSeg->GetEnd(s));
  

  // data structure
  if(!dataseq->SetSegments(aSeg)){
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
  if(fVerbosity) cout<<"Omicron::NewChunk: load a new chunk..."<<endl;
  if(!dataseq->NewChunk()) return false;
  if(fVerbosity>1) cout<<"\t- chunk "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<" is loaded"<<endl;

  // save info for html report
  if(fOutProducts.find("html")!=string::npos){
    for(int s=0; s<dataseq->GetNSegments(); s++){
      mapcenter.push_back(dataseq->GetSegmentTimeStart(s)+dataseq->GetSegmentDuration()/2);
      if(!s){
	chunkstart.push_back(dataseq->GetChunkTimeStart()); 
	chunkstop.push_back(dataseq->GetChunkTimeEnd());
      }
    }
  }

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

  if(fVerbosity){
    if(fInjChan.size()) cout<<"Omicron::LoadData: load data vector and add injections..."<<endl;
    else cout<<"Omicron::LoadData: load data vector..."<<endl;
  }

  // get data vector
  *aDataVector = FFL->GetData(*aSize, fChannels[chanindex], dataseq->GetChunkTimeStart(), dataseq->GetChunkTimeEnd());

  // cannot retrieve data
  if(*aSize<=0){
    cerr<<"Omicron::LoadData: cannot retrieve data ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
    aDataVector=NULL; *aSize=0;
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
    cerr<<"Omicron::Condition: input vector is NULL ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
    return 1;
  }
  if(aInVectSize<=0){
    cerr<<"Omicron::Condition: input vector is empty ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
    return 2;
  }
  if(aInVect[0]==aInVect[aInVectSize-1]){
    cerr<<"Omicron::Condition: input vector appears to be flat ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
    return 3;
  }
    
  if(fVerbosity) cout<<"Omicron::Condition: condition data vector..."<<endl;

  // test native sampling (and update if necessary)
  int nativesampling = aInVectSize/(dataseq->GetCurrentChunkDuration());
  if(!triggers[chanindex]->SetNativeFrequency(nativesampling)){
    cerr<<"Omicron::Condition: incompatible native frequency ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
    return 4;
  }

  // transform data vector
  if(fVerbosity>1) cout<<"\t- transform data vector..."<<endl;
  if(!triggers[chanindex]->Transform(aInVectSize, aInVect, dataseq->GetCurrentChunkDuration()*triggers[chanindex]->GetWorkingFrequency(), ChunkVect)) return 5;

  // add software injections if any
  if(inject!=NULL){
    if(fVerbosity>1) cout<<"\t- perform injections..."<<endl;
    if(!inject[chanindex]->Inject(dataseq->GetCurrentChunkDuration()*triggers[chanindex]->GetWorkingFrequency(), ChunkVect, dataseq->GetChunkTimeStart())){
      cerr<<"Omicron::Condition: failed to inject ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
      return 6;

    }
  }
 
  // Make spectrum
  if(fVerbosity>1) cout<<"\t- make spectrum..."<<endl;
  if(!spectrum->LoadData(dataseq->GetCurrentChunkDuration()*triggers[chanindex]->GetWorkingFrequency(), ChunkVect)) return 7;

  // connect spectrum to tiling
  if(fVerbosity>1) cout<<"\t- connect spectrum to tiling..."<<endl;
  if(!tile->SetPower(spectrum)) return 8;

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
  int segsize=dataseq->GetSegmentDuration()*triggers[chanindex]->GetWorkingFrequency();
  int ovsize=dataseq->GetOverlapDuration()*triggers[chanindex]->GetWorkingFrequency();
  for(int s=0; s<dataseq->GetNSegments(); s++){
    if(fVerbosity>1) cout<<"\t- subsegment "<<dataseq->GetSegmentTimeStart(s)<<"-"<<dataseq->GetSegmentTimeEnd(s)<<endl;

    // fill segment vector (time-domain) and apply tukey window
    if(fVerbosity>2) cout<<"\t\t- fill subsegment"<<endl;
    for(int i=0; i<segsize; i++)
      SegVect[i] = ChunkVect[s*(segsize-ovsize)+i] * TukeyWindow[i];
    
    // get whiten data
    if(fVerbosity>2) cout<<"\t\t- whiten subsegment"<<endl;
    if(!Whiten(&(dataRe[s]), &(dataIm[s]))) return false;
    
    // project data
    if(fVerbosity>2) cout<<"\t\t- project subsegment"<<endl;
    if(!tile->ProjectData(dataRe[s], dataIm[s], fTileDown)){
      delete dataRe[s];
      delete dataIm[s];
      return false;
    }

    // save whiten data for condition data products
    if(fOutProducts.find("condition")!=string::npos){
      if(!offt->Backward(dataRe[s], dataIm[s])){// Back in time domain
	delete dataRe[s];
	delete dataIm[s];
	return false;
      }
      for(int i=ovsize/2; i<segsize-ovsize/2; i++) CondChunkVect[s*(segsize-ovsize)+i]=offt->GetReOut(i);// CHECKME: tukey window?
    }

    // not used anymore
    delete dataRe[s];
    delete dataIm[s];

    // write maps on disk
    double snr;
    if(fOutProducts.find("maps")!=string::npos){
      if(fVerbosity>2) cout<<"\t\t- write maps"<<endl;
      if(fOutProducts.find("html")==string::npos) // need to make thumbnails
	snr=tile->SaveMaps(outdir[chanindex],
			   fChannels[chanindex],
			   dataseq->GetSegmentTimeStart(s)+dataseq->GetSegmentDuration()/2,
			   fOutFormat,fWindows,false);
      else
	snr=tile->SaveMaps(outdir[chanindex],
			   fChannels[chanindex],
			   dataseq->GetSegmentTimeStart(s)+dataseq->GetSegmentDuration()/2,
			   fOutFormat,fWindows,true);
      if(snr>chan_mapsnrmax[chanindex]) chan_mapsnrmax[chanindex]=snr;
    }
        
    // save triggers
    if(fOutProducts.find("triggers")!=string::npos){
      if(fVerbosity>2) cout<<"\t\t- save triggers"<<endl;
      if(!tile->SaveTriggers(triggers[chanindex],
			     (double)(dataseq->GetCurrentOverlapDuration()-dataseq->GetOverlapDuration()/2),
			     (double)(dataseq->GetOverlapDuration()/2),
			     (double)(dataseq->GetSegmentTimeStart(s)+dataseq->GetSegmentDuration()/2))
	 ) return false;
    }
    
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
  
  //*** RAW TS
  if(fOutProducts.find("timeseries")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write time-series..."<<endl;
    SaveTS(false);
  }

  //*** CONDITIONNED DATA
  if(fOutProducts.find("condition")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write conditioned data..."<<endl;
    SaveTS(true);// save whitened time series
    spectrum->SetPadding(fOverlapDuration);
    if(!spectrum->LoadData(dataseq->GetCurrentChunkDuration()*triggers[chanindex]->GetWorkingFrequency(), CondChunkVect)) return false;// compute spectrum
    spectrum->SetPadding(0);
    SaveAPSD("PSD",true);// save whiten PSD
  }

  //*** TRIGGERS
  if(fOutProducts.find("triggers")!=string::npos){
    if(fVerbosity>1) cout<<"\t- write triggers "<<endl;

    // sort triggers
    triggers[chanindex]->SortTriggers();
    
    // clustering if any
    if(fClusterAlgo.compare("none")) triggers[chanindex]->Clusterize();
    
    // write triggers to disk
    if(!triggers[chanindex]->Write().compare("none")){
      cerr<<"Omicron::WriteOutput: triggers cannot be written to disk ("<<fChannels[chanindex]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")"<<endl;
      return false;
    }
  }

  chan_write_ctr[chanindex]++;
  outSegments[chanindex]->AddSegment(dataseq->GetChunkTimeStart()+fOverlapDuration/2,dataseq->GetChunkTimeEnd()-fOverlapDuration/2);// FIXME: known bug with the timing (last chunk)
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::PrintStatusInfo(void){
////////////////////////////////////////////////////////////////////////////////////

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
  int ssize=dataseq->GetSegmentDuration()*triggers[chanindex]->GetWorkingFrequency()/2;
  for(int i=icutoff; i<ssize; i++){
    asdval=spectrum->GetPower(i*(double)triggers[chanindex]->GetWorkingFrequency()/(double)(2*ssize));
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
void Omicron::SaveAPSD(const string type, const bool aCond){
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
  if(aCond) GAPSD->SetTitle(((string)GAPSD->GetTitle()+" (after conditioning)").c_str());
  GAPSD->SetLineWidth(2);
  GAPSD->GetXaxis()->SetLimits(fFreqRange[0],fFreqRange[1]);
  GAPSD->GetXaxis()->SetTitleOffset(1.1);
  GAPSD->GetXaxis()->SetLabelSize(0.045);
  GAPSD->GetYaxis()->SetLabelSize(0.045);
  GAPSD->GetXaxis()->SetTitleSize(0.045);
  GAPSD->GetYaxis()->SetTitleSize(0.045);
  
  // set new name
  stringstream ss;
  ss<<type;
  if(aCond) ss<<"cond";
  ss<<"_"<<fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd();
  GAPSD->SetName(ss.str().c_str());
  ss.str(""); ss.clear();

  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fpsd;
    ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"_"<<type;
    if(aCond) ss<<"cond";
    ss<<".root";
    fpsd=new TFile((ss.str()).c_str(),"RECREATE");
    ss.str(""); ss.clear();
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
      ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"_"<<type;
      if(aCond) ss<<"cond";
      ss<<"."<<form[f];
      GPlot->Print(ss.str().c_str());
      ss.str(""); ss.clear();
    }
  }
  
  form.clear();
  delete GAPSD;
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveTS(const bool aCond){
////////////////////////////////////////////////////////////////////////////////////

  // create graph 
  TGraph *GDATA = new TGraph(triggers[chanindex]->GetWorkingFrequency()*dataseq->GetCurrentChunkDuration());
  if(GDATA==NULL) return;

  stringstream ss;

  // condition data
  if(aCond){
    ss<<"tscond_"<<fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd();
    for(int i=0; i<triggers[chanindex]->GetWorkingFrequency()*dataseq->GetCurrentChunkDuration(); i++) GDATA->SetPoint(i,(double)dataseq->GetChunkTimeStart()+(double)i/(double)(triggers[chanindex]->GetWorkingFrequency())-timeoffset,CondChunkVect[i]);
  }
  // raw data
  else{
    ss<<"ts_"<<fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd();
    for(int i=0; i<triggers[chanindex]->GetWorkingFrequency()*dataseq->GetCurrentChunkDuration(); i++) GDATA->SetPoint(i,(double)dataseq->GetChunkTimeStart()+(double)i/(double)(triggers[chanindex]->GetWorkingFrequency())-timeoffset,ChunkVect[i]);
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
  if(aCond) GDATA->SetTitle((fChannels[chanindex]+": amplitude conditioned data time series").c_str());
  else GDATA->SetTitle((fChannels[chanindex]+": amplitude raw time series (high-passed)").c_str());
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
    if(aCond)
      ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"_tscond.root";
    else
      ss<<outdir[chanindex]<<"/"+fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"_ts.root";
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
    
    // give the time offset info in the title
    if(timeoffset){
      ss<<" centered at "<<fixed<<setprecision(3)<<timeoffset;
      GDATA->SetTitle(((string)GDATA->GetTitle()+ss.str()).c_str());
      ss.str(""); ss.clear();
    }

    // zoom
    double tcenter;
    if(timeoffset) tcenter=0.0;
    else tcenter=(double)(dataseq->GetChunkTimeStart())/2.0+(double)(dataseq->GetChunkTimeEnd())/2.0;
    for(int w=(int)fWindows.size()-1; w>=0; w--){
      GDATA->GetXaxis()->SetLimits(tcenter-(double)fWindows[w]/2.0,tcenter+(double)fWindows[w]/2.0);
      for(int f=0; f<(int)form.size(); f++){
	if(aCond)
	  ss<<outdir[chanindex]<<"/"<<fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"_tsconddt"<<fWindows[w]<<"."<<form[f];
	else
	  ss<<outdir[chanindex]<<"/"<<fChannels[chanindex]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"_tsdt"<<fWindows[w]<<"."<<form[f];
	GPlot->Print(ss.str().c_str());
	ss.str(""); ss.clear();
      }
    }
  }
  
  form.clear();
  delete GDATA;
  return;
}

