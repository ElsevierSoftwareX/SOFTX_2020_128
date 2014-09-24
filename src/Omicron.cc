//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

ClassImp(Omicron)

////////////////////////////////////////////////////////////////////////////////////
Omicron::Omicron(const string aOptionFile){ 
////////////////////////////////////////////////////////////////////////////////////
  PrintASCIIlogo();

  // init timer
  time ( &timer );
  timer_start=timer;
  ptm = gmtime ( &timer );  
  cout<<"#### Omicron::Omicron timer = "<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) ---> +"<<timer-timer_start<<"s ####"<<endl;

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
  if(fVerbosity) cout<<"Omicron::Omicron: init monitoring..."<<endl;
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
  FFL=NULL;
  if(fFflFile.compare("none")){
    if(fVerbosity) cout<<"Omicron::Omicron: init FFL..."<<endl;
    FFL = new ffl(fFflFile, fStyle, fVerbosity);
    status_OK*=FFL->DefineTmpDir(fMaindir);
    status_OK*=FFL->LoadFrameFile();
    
    // guess best sampling if not provided
    if(!fSampleFrequency){
      fSampleFrequency=2048;
      int sampling;
      for(int c=0; c<(int)fChannels.size(); c++){
	sampling=FFL->GetChannelSampling(fChannels[c]);
	if(sampling<=0) continue;
	if(sampling<fSampleFrequency) fSampleFrequency=sampling;
      }
    }
  }

  
  // init data
  if(fVerbosity) cout<<"Omicron::Omicron: init data container and procedures..."<<endl;
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
  spectrum = new Spectrum(fSampleFrequency,psdsize,fOverlapDuration,fVerbosity);
  status_OK*=spectrum->GetStatus();

  // output directory
  // will change for Process() and Scan()
  for(int c=0; c<(int)fChannels.size(); c++) fOutdir.push_back(fMaindir);
  fScandir="";// used for scans

  // init Triggers
  if(fVerbosity) cout<<"Omicron::Omicron: init Triggers..."<<endl;
  triggers    = new Triggers* [(int)fChannels.size()];
  string form = "";
  if(fOutFormat.find("xml")!=string::npos)  form+="xml";
  if(fOutFormat.find("txt")!=string::npos)  form+="txt";
  if(fOutFormat.find("root")!=string::npos) form+="root";
  if(!form.compare("")) form="root";
  for(int c=0; c<(int)fChannels.size(); c++){
    
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
    status_OK*=triggers[c]->SetUserMetaData(fOptionName[0],fMaindir);
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
  if(fVerbosity) cout<<"Omicron::Omicron: init Tile..."<<endl;
  tile = new Otile(fSegmentDuration,fOverlapDuration/2,fQRange[0],fQRange[1],fFreqRange[0],fFreqRange[1],fSampleFrequency,fMismatchMax,fSNRThreshold,fVerbosity);
  status_OK*=tile->GetStatus();

  // init Qmaps
  if(fVerbosity) cout<<"Omicron::Omicron: init Qmaps..."<<endl;
  Qmap = new TH2D* [tile->GetNQPlanes()];
  for(int q=0; q<tile->GetNQPlanes(); q++) Qmap[q]=NULL;
  Qmap_center=-1.0;
  MapChNumber=-1;

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
  
  if(!status_OK) cerr<<"Omicron::Omicron: Initialization failed"<<endl;
  cout<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
Omicron::~Omicron(void){
////////////////////////////////////////////////////////////////////////////////////
  time ( &timer );
  ptm = gmtime ( &timer );
  cout<<"#### Omicron::~Omicron timer = "<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) ---> +"<<timer-timer_start<<"s ####"<<endl;
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

  // timer
  if(fVerbosity){
    time ( &timer );
    ptm = gmtime ( &timer );
    cout<<"#### Omicron::Process timer = "<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) ---> +"<<timer-timer_start<<"s ####"<<endl;
  }

  // init segments
  if(!InitSegments(aSeg)){
    cerr<<"Omicron::Process: the requested segments cannot be initialized"<<endl;
    return false;
  }

  // locals
  int dsize;         // native data size
  double *dvector;   // data vector before resampling
  int res;

  // create trigger directories
  if(!MakeDirectories(0.0)){
    cerr<<"Omicron::Process: the directory structure cannot be created"<<endl;
    return false;
  }
  
  // loop over chunks
  while(NewChunk()){
    cout<<"Omicron::Process: ** chunk "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<endl;
    
    // loop over channels
    for(int c=0; c<(int)fChannels.size(); c++){
      cout<<"Omicron::Process: *  channel "<<fChannels[c]<<"..."<<endl;
      
      // one more chunk of data
      chunk_ctr[c]++;

      // get data vector
      dvector=FFL->GetData(dsize, fChannels[c], dataseq->GetChunkTimeStart(), dataseq->GetChunkTimeEnd());

      if(dsize<=0){
	cor_chunk_ctr[c]++;
	cerr<<"Omicron::Process: cannot retrieve data ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
	continue;
      }

      // condition data vector
      res=ConditionVector(c, dsize, dvector);
      if(res<0){
	cerr<<"Omicron::Process: conditioning failed ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
	delete dvector;// not needed anymore
	cor_data_ctr[c]++;
	return false;// fatal
      }
      if(res>0){
	cerr<<"Omicron::Process: conditioning failed ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
	delete dvector;// not needed anymore
	cor_data_ctr[c]++;
	continue;// skip channel
      }

      delete dvector;// not needed anymore

      // Save info if requested
      if(writetimeseries) SaveTS(c);
      if(writepsd)        SaveAPSD(c,"PSD");

      // get triggers above SNR threshold
      if(MakeTriggers(c)<0){
	cerr<<"Omicron::Process: cannot make triggers ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
	cor_data_ctr[c]++;
	continue;
      }
     
      // save triggers on disk
      if(!WriteTriggers(c)){
	cerr<<"Omicron::Process: This chunk of triggers is not saved ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
	max_chunk_ctr[c]++;
	continue; // move to next channel
      }
      	
    }

    // time for this chunk
    if(fVerbosity>1){ 
      time ( &timer );
      ptm = gmtime ( &timer );
      cout<<"#### Omicron::Process timer = "<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) ---> +"<<timer-timer_start<<"s ####"<<endl;
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
  if(fVerbosity){
    time ( &timer );
    ptm = gmtime ( &timer );
    cout<<"#### Omicron::Scan timer = "<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) ---> +"<<timer-timer_start<<"s ####"<<endl;
  }

  // locals
  int dsize;         // native data size
  double *dvector;   // data vector before resampling
  int res;
  
  // Segment to process
  if(fVerbosity) cout<<"Omicron::Scan: initiate data segments..."<<endl;
  Segments *ScanSeg = new Segments((int)aTimeCenter-fChunkDuration/2,(int)aTimeCenter+fChunkDuration/2);

  // data structure
  if(!InitSegments(ScanSeg)||!NewChunk()){
    cerr<<"Omicron::Scan: cannot initiate data segments."<<endl;
    delete ScanSeg;
    return false;
  }

  // create map directories
  if(!MakeDirectories(aTimeCenter)){
    cerr<<"Omicron::Scan: the directory structure cannot be created"<<endl;
    delete ScanSeg;
    return false;
  }

  // loop over channels
  for(int c=0; c<(int)fChannels.size(); c++){
    cout<<"Omicron::Process: **** channel "<<fChannels[c]<<"..."<<endl;
          
    // get data vector
    dvector=FFL->GetData(dsize, fChannels[c], dataseq->GetChunkTimeStart(), dataseq->GetChunkTimeEnd());

    if(dsize<=0){
      cor_chunk_ctr[c]++;
     cerr<<"Omicron::Scan: cannot retrieve data ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")."<<endl;
      continue;
    }

    // condition data vector
    res=ConditionVector(c, dsize, dvector);
    if(res<0){
      cerr<<"Omicron::Scan: conditioning failed ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
      delete dvector;// not needed anymore
      delete ScanSeg;
      cor_data_ctr[c]++;
      return false;// fatal
    }
    if(res>0){
      cerr<<"Omicron::Scan: conditioning failed ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
      delete dvector;// not needed anymore
      cor_data_ctr[c]++;
      continue;// skip channel
    }

    delete dvector;// not needed anymore

    // get maps
    if(!MakeMaps(c, aTimeCenter)){
      cerr<<"Omicron::Scan: cannot make maps ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
      cor_data_ctr[c]++;
      continue;
    }

    // apply SNR threshold (except the first channel)
    if(c&&Qmap_full[0]->GetMaximum()<fSNRThreshold){
      cout<<"Omicron::Scan: below SNR threshold "<<Qmap_full[0]->GetMaximum()<<"<"<<fSNRThreshold<<" ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
      max_chunk_ctr[c]++;// use this specific monitor here
      continue;
    }

    // save maps on disk
    if(!WriteMaps()){
      cerr<<"Omicron::Scan: cannot write maps for channel "<<fChannels[c]<<endl;
      cor_data_ctr[c]++;
      continue;
    }

    // save ASD on disk
    SaveAPSD(c,"ASD");
    
    // save data on disk
    SaveTS(c, aTimeCenter);

    // time for this channel
    if(fVerbosity>1){ 
      time ( &timer );
      ptm = gmtime ( &timer );
      cout<<"#### Omicron::Scan timer = "<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) ---> +"<<timer-timer_start<<"s ####"<<endl;
    }
  }
        
  delete ScanSeg;

  return true;
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

  // add requested segments (try to optimize the update of inSegments)
  if(inSegments->GetLiveTime()&&aSeg->GetStart(0)>=inSegments->GetStart(inSegments->GetNsegments()-1))
    inSegments->Append(aSeg);
  else{
    for(int s=0; s<aSeg->GetNsegments(); s++) inSegments->AddSegment(aSeg->GetStart(s),aSeg->GetEnd(s));
  }

  // data structure
  if(fVerbosity) cout<<"Omicron::InitSegments: initiate data segments..."<<endl;
  if(!dataseq->SetSegments(aSeg)){
    cerr<<"Omicron::InitSegments: cannot initiate data segments."<<endl;
    return false;
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
      if(fVerbosity>1) cout<<"                          - "<<fOutdir[c]<<endl;
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
      if(fVerbosity>1) cout<<"                          - "<<fOutdir[c]<<endl;
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
  if(fVerbosity) cout<<"Omicron::NewChunk: load a new chunk segment..."<<endl;
  if(!dataseq->NewChunk()){
    cerr<<"Omicron::NewChunk: no more chunk to load"<<endl;
    return false;
  }
  
  // size update
  // only the chunk size can change
  fChunkDuration=dataseq->GetChunkTimeEnd()-dataseq->GetChunkTimeStart();
  ChunkSize=fChunkDuration*fSampleFrequency;
  if(fVerbosity>1) cout<<"Omicron::NewChunk: chunk "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<" is loaded"<<endl;

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
  if((dataseq->GetChunkTimeEnd()-dataseq->GetChunkTimeStart())*fSampleFrequency!=ChunkSize){
    cerr<<"Omicron::ConditionVector: the chunk size is not compatible with the chunk duration"<<endl;
    return 2;
  }

  // Input data checks
  if(aInVect==NULL){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): input vector is NULL"<<endl;
    return 3;
  }
  if(aInVectSize<=0){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): input vector empty"<<endl;
    return 4;
  }
  if(aInVect[0]==aInVect[aInVectSize-1]){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): input vector is flat"<<endl;
    return 5;
  }
  
  // test native sampling (and update if necessary)
  int nativesampling = aInVectSize/(dataseq->GetChunkTimeEnd()-dataseq->GetChunkTimeStart());
  if(!sample[aChNumber]->SetFrequencies(nativesampling,fSampleFrequency,fFreqRange[0])){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): the Sample structure could not be parameterized"<<endl;
    return -2;
  }

  // transform data vector
  if(fVerbosity) cout<<"Omicron::ConditionVector: transform data vector..."<<endl;
  if(!sample[aChNumber]->Transform(aInVectSize, aInVect, ChunkSize, ChunkVect)){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): the transform failed"<<endl;
    return 6;
  }

  // Make spectrum
  if(fVerbosity) cout<<"Omicron::ConditionVector: make spectrum..."<<endl;
  if(!spectrum->LoadData(ChunkSize, ChunkVect)){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): the spectrum could not be estimated"<<endl;
    return 7;
  }

  // set new power for this chunk
  if(fVerbosity) cout<<"Omicron::ConditionVector: normalize tiling..."<<endl;
  if(!tile->SetPowerSpectrum(spectrum)){
    cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): the tiling could not be normalized"<<endl;
    return 8;
  }

  //loop over segments
  if(fVerbosity) cout<<"Omicron::ConditionVector: condition data segments..."<<endl;
  for(int s=0; s<dataseq->GetNSegments(); s++){
    	
    // fill segment vector (time-domain) and apply tukey window
    for(int i=0; i<SegmentSize; i++)
      SegVect[i] = ChunkVect[s*(SegmentSize-OverlapSize)+i] * TukeyWindow[i];
    
    // get conditioned data
    if(!Condition(&(dataRe[s]), &(dataIm[s]))){
      cerr<<"Omicron::ConditionVector ("<<fChannels[aChNumber]<<"): the conditioning failed"<<endl;
      return 9;
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
    delete dataRe[ss];
    delete dataIm[ss];
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

  // don't save if max flag
  if(triggers[aChNumber]->GetMaxFlag()){
    cerr<<"Omicron::WriteTriggers: channel "<<fChannels[aChNumber]<<" is maxed-out. This chunk of triggers is not saved"<<endl;
    triggers[aChNumber]->Reset();
    return false;
  }
      
  // save triggers for this chunk
  if(!triggers[aChNumber]->Write(fWriteMode, "default").compare("none")){
    cerr<<"Omicron::WriteTriggers: writing events failed for channel "<<fChannels[aChNumber]<<endl;
    return false;
  }

  // processed chunks
  else
    outSegments[aChNumber]->AddSegment(dataseq->GetChunkTimeStart()+fOverlapDuration/2,dataseq->GetChunkTimeEnd()-fOverlapDuration/2);
	
  return true;
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
  if(aTimeCenter-fWindowMax<dataseq->GetSegmentTimeStart(0)+fOverlapDuration/2){
    cerr<<"Omicron::MakeMaps: the requested central time "<<aTimeCenter<<" is not compatible with the data segmentation and the windowing -> you need to increase your segment"<<endl;
    return false;
  }
  if(aTimeCenter+fWindowMax>dataseq->GetSegmentTimeEnd(0)-fOverlapDuration/2){
    cerr<<"Omicron::MakeMaps: the requested central time "<<aTimeCenter<<" is not compatible with the data segmentation and the windowing -> you need to increase your segment"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::MakeMaps: make maps for channel "<<fChannels[aChNumber]<<" starting at "<<dataseq->GetChunkTimeStart()<<endl;

  Qmap_center=aTimeCenter;
  double toffset=Qmap_center-dataseq->GetSegmentTimeStart(0)-(double)fSegmentDuration/2.0;
  ostringstream tmpstream;

  // get Qmaps
  if(fVerbosity>1) cout<<"Omicron::MakeMaps: make Qmaps..."<<endl;
  for(int q=0; q<tile->GetNQPlanes(); q++){
    Qmap[q] = tile->GetMap(q, dataRe[0], dataIm[0], -toffset);
    if(Qmap[q]==NULL){
      cerr<<"Omicron::MakeMaps: maps for channel "<<fChannels[aChNumber]<<" are corrupted"<<endl;
      delete dataRe[0];
      delete dataIm[0];
      return false;
    }

    // map title
    tmpstream<<fChannels[aChNumber]<<": GPS="<<fixed<<setprecision(3)<<Qmap_center<<", Q="<<tile->GetQ(q);
    Qmap[q]->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
    
    // cosmetics
    //Qmap[q]->GetXaxis()->SetTitle("Time [s]");
    //Qmap[q]->GetYaxis()->SetTitle("Frequency [Hz]");
    //Qmap[q]->GetZaxis()->SetTitle("SNR [-]");
    Qmap[q]->GetXaxis()->SetTitleOffset(1.1);
    Qmap[q]->GetXaxis()->SetLabelSize(0.045);
    Qmap[q]->GetYaxis()->SetLabelSize(0.045);
    Qmap[q]->GetXaxis()->SetTitleSize(0.045);
    Qmap[q]->GetYaxis()->SetTitleSize(0.045);
    Qmap[q]->GetZaxis()->SetTitleSize(0.05);
  }
  delete dataRe[0];
  delete dataIm[0];

  MapChNumber=aChNumber;

  // create overall maps
  if(fVerbosity>1) cout<<"Omicron::MakeMaps: make overall maps..."<<endl;
  for(int w=0; w<(int)fWindows.size(); w++){
    tmpstream<<fChannels[aChNumber]<<": GPS="<<fixed<<setprecision(3)<<Qmap_center;
    Qmap_full[w]->Reset();
    Qmap_full[w]->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
  }
      
  //populate overall maps
  int bin_start, bin_stop, dummy, bin_start_t, bin_stop_t, bin_stop_f, bin_start_f;
  double content;
  for(int q=0; q<tile->GetNQPlanes(); q++){
    for(int bt=1; bt<=Qmap[q]->GetNbinsX(); bt++){
      for(int bf=1; bf<=Qmap[q]->GetNbinsY(); bf++){
      	content=Qmap[q]->GetBinContent(bt,bf);
	for(int w=0; w<(int)fWindows.size(); w++){
	  if(Qmap[q]->GetXaxis()->GetBinUpEdge(bt)<-(double)fWindows[w]/2.0) continue;
	  if(Qmap[q]->GetXaxis()->GetBinLowEdge(bt)>(double)fWindows[w]/2.0) continue;
	  if(Qmap[q]->GetYaxis()->GetBinUpEdge(bf)<Qmap_full[w]->GetYaxis()->GetBinLowEdge(1)) continue;
	  if(Qmap[q]->GetYaxis()->GetBinLowEdge(bf)>Qmap_full[w]->GetYaxis()->GetBinUpEdge(Qmap_full[w]->GetNbinsY())) continue;
	  bin_start = Qmap_full[w]->FindBin(Qmap[q]->GetXaxis()->GetBinLowEdge(bt),Qmap[q]->GetYaxis()->GetBinLowEdge(bf));
	  bin_stop = Qmap_full[w]->FindBin(Qmap[q]->GetXaxis()->GetBinUpEdge(bt),Qmap[q]->GetYaxis()->GetBinUpEdge(bf));
	  Qmap_full[w]->GetBinXYZ(bin_start,bin_start_t,bin_start_f,dummy);
	  Qmap_full[w]->GetBinXYZ(bin_stop,bin_stop_t,bin_stop_f,dummy);
	  for(int bbt=bin_start_t; bbt<=bin_stop_t; bbt++){// time-sweep the tile
	    for(int bbf=bin_start_f; bbf<=bin_stop_f; bbf++){// freq-sweep the tile
	      if(content>Qmap_full[w]->GetBinContent(bbt,bbf)) Qmap_full[w]->SetBinContent(bbt,bbf,content);
	    }
	  }
	}
      }
    }
  }
   
  // Q-map containing the loudest tile
  for(int w=0; w<(int)fWindows.size(); w++){
    loudest_qmap[w]=0;
    for(int q=1; q<tile->GetNQPlanes(); q++){
      Qmap[q]->GetXaxis()->SetRangeUser(-(double)fWindows[0]/2.0,(double)fWindows[0]/2.0);
      if(Qmap[q]->GetMaximum()>Qmap[loudest_qmap[w]]->GetMaximum()) loudest_qmap[w]=q;
      Qmap[q]->GetXaxis()->SetRange();
    }
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::WriteMaps(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::WriteMaps: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(MapChNumber<0){
    cerr<<"Omicron::WriteMaps: so far nothing has been mapped"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::WriteMaps: write maps for channel "<<fChannels[MapChNumber]<<endl;
  
  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fmap;
    if(first_save[MapChNumber]){
      fmap=new TFile((fOutdir[MapChNumber]+"/"+fChannels[MapChNumber]+"_data.root").c_str(),"RECREATE");
      first_save[MapChNumber]=false;
    }
    else fmap=new TFile((fOutdir[MapChNumber]+"/"+fChannels[MapChNumber]+"_data.root").c_str(),"UPDATE");
    fmap->cd();
    for(int q=0; q<tile->GetNQPlanes(); q++) Qmap[q]->Write();
    for(int w=0; w<(int)fWindows.size(); w++) Qmap_full[w]->Write();
    fmap->Close();
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
  if(!form.size()) return true;

  // summary tree
  TTree *MS = new TTree("mapsummary", "map summary");
  double map_gps, map_q, map_gpsmax, map_freqmax, map_snrmax;
  int map_win, map_nq, map_nwin, map_nformat, map_wstart, map_wstop;
  string map_format;
  MS->Branch("gps_center",&map_gps,"gps_center/D");
  MS->Branch("gps_whitening_start",&map_wstart,"gps_whitening_start/I");
  MS->Branch("gps_whitening_stop",&map_wstop,"gps_whitening_stop/I");
  MS->Branch("nQ",&map_nq,"nQ/I");
  MS->Branch("nwindow",&map_nwin,"nwindow/I");
  MS->Branch("Q",&map_q,"Q/D");
  MS->Branch("window",&map_win,"window/I");
  MS->Branch("gps_loudest",&map_gpsmax,"gps_loudest/D");
  MS->Branch("freq_loudest",&map_freqmax,"freq_loudest/D");
  MS->Branch("snr_loudest",&map_snrmax,"snr_loudest/D");
  MS->Branch("nformat",&map_nformat,"nformat/I");
  MS->Branch("format",&map_format);
  map_gps=Qmap_center;
  map_nwin=(int)fWindows.size();
  map_nq=tile->GetNQPlanes();
  map_nformat=(int)form.size();
  map_wstart=dataseq->GetChunkTimeStart();
  map_wstop=dataseq->GetChunkTimeEnd();

  // canvas style
  GPlot->SetLogx(0);
  GPlot->SetLogy(1);
  GPlot->SetLogz(1);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);

  // locals
  double center;
  ostringstream tmpstream;
  int xmax, ymax, zmax;

  // draw Qmaps
  for(int q=0; q<tile->GetNQPlanes(); q++){
    map_q=tile->GetQ(q);

    // plot
    Qmap[q]->GetZaxis()->SetRangeUser(1,50);
    GPlot->Draw(Qmap[q],"COLZ");
    
    // window resize
    center=(Qmap[q]->GetXaxis()->GetXmax()+Qmap[q]->GetXaxis()->GetXmin())/2.0;
    for(int w=0; w<(int)fWindows.size(); w++){
      map_win=fWindows[w];

      // zoom
      Qmap[q]->GetXaxis()->SetRangeUser(center-(double)fWindows[w]/2.0,center+(double)fWindows[w]/2.0);

      // get max bin
      Qmap[q]->GetMaximumBin(xmax, ymax, zmax);
         
      // loudest tile
      tmpstream<<"Loudest tile: GPS="<<fixed<<setprecision(3)<<Qmap_center+Qmap[q]->GetXaxis()->GetBinCenter(xmax)<<", f="<<Qmap[q]->GetYaxis()->GetBinCenter(ymax)<<" Hz, SNR="<<Qmap[q]->GetBinContent(xmax,ymax);
      GPlot->AddText(tmpstream.str(), 0.01,0.01,0.03);
      tmpstream.str(""); tmpstream.clear();
      map_gpsmax=Qmap_center+Qmap[q]->GetXaxis()->GetBinCenter(xmax);
      map_freqmax=Qmap[q]->GetYaxis()->GetBinCenter(ymax);
      map_snrmax=Qmap[q]->GetBinContent(xmax,ymax);

      // save qmaps
      for(int f=0; f<(int)form.size(); f++){
	tmpstream<<fOutdir[MapChNumber]<<"/"<<fChannels[MapChNumber]<<"_mapQ"<<q<<"_dt"<<fWindows[w]<<"."<<form[f];
	GPlot->Print(tmpstream.str());
	tmpstream.str(""); tmpstream.clear();
 	tmpstream<<fOutdir[MapChNumber]<<"/"<<fChannels[MapChNumber]<<"_mapQ"<<q<<"th_dt"<<fWindows[w]<<"."<<form[f];
	GPlot->SetGridx(0); GPlot->SetGridy(0);
	GPlot->Print(tmpstream.str(),0.5);
	GPlot->SetGridx(1); GPlot->SetGridy(1);
 	tmpstream.str(""); tmpstream.clear();
	map_format=form[f];
	
	// fill summary tree
	MS->Fill();
      }

    }
    Qmap[q]->GetXaxis()->UnZoom();
    Qmap[q]->GetZaxis()->UnZoom();
  }

  // draw overall maps
  for(int w=0; w<(int)fWindows.size(); w++){
    Qmap_full[w]->GetZaxis()->SetRangeUser(1,50);
    Qmap_full[w]->GetXaxis()->SetRangeUser(center-(double)fWindows[w]/2.0,center+(double)fWindows[w]/2.0);
    GPlot->Draw(Qmap_full[w],"COLZ");

    // get max bin
    Qmap_full[w]->GetMaximumBin(xmax, ymax, zmax);
         
    // loudest tile
    tmpstream<<"Loudest tile: GPS="<<fixed<<setprecision(3)<<Qmap_center+Qmap_full[w]->GetXaxis()->GetBinCenter(xmax)<<", f="<<Qmap_full[w]->GetYaxis()->GetBinCenter(ymax)<<" Hz, SNR="<<Qmap_full[w]->GetBinContent(xmax,ymax);
    GPlot->AddText(tmpstream.str(), 0.01,0.01,0.03);
    tmpstream.str(""); tmpstream.clear();

    // save maps
    for(int f=0; f<(int)form.size(); f++){
      tmpstream<<fOutdir[MapChNumber]<<"/"<<fChannels[MapChNumber]<<"_map_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();
      tmpstream<<fOutdir[MapChNumber]<<"/"<<fChannels[MapChNumber]<<"_mapth_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->SetGridx(0); GPlot->SetGridy(0);
      GPlot->Print(tmpstream.str(),0.5);
      GPlot->SetGridx(1); GPlot->SetGridy(1);
      tmpstream.str(""); tmpstream.clear();
    }
    Qmap_full[w]->GetXaxis()->UnZoom();
    Qmap_full[w]->GetZaxis()->UnZoom();
  }

  // save projections
  GPlot->SetLogx(0);
  GPlot->SetLogy(1);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);
  TH1D* proj;
  for(int w=0; w<(int)fWindows.size(); w++){
    Qmap_full[w]->GetMaximumBin(xmax, ymax, zmax);
    proj=Qmap_full[w]->ProjectionX("_pfx",ymax,ymax);
    tmpstream<<"SNR vs time at f = "<<fixed<<setprecision(3)<<Qmap_full[w]->GetYaxis()->GetBinCenter(ymax)<<"Hz";
    proj->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
    proj->GetXaxis()->SetTitle("Time [s]");
    proj->GetYaxis()->SetTitle("SNR");
    proj->GetXaxis()->SetTitleOffset(1.1);
    proj->GetXaxis()->SetLabelSize(0.045);
    proj->GetYaxis()->SetLabelSize(0.045);
    proj->GetXaxis()->SetTitleSize(0.045);
    proj->GetYaxis()->SetTitleSize(0.045);
    GPlot->SetLogx(0);
    GPlot->Draw(proj);
    for(int f=0; f<(int)form.size(); f++){
      tmpstream<<fOutdir[MapChNumber]<<"/"<<fChannels[MapChNumber]<<"_projt_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();
      tmpstream<<fOutdir[MapChNumber]<<"/"<<fChannels[MapChNumber]<<"_projtth_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->Print(tmpstream.str(),0.5);
      tmpstream.str(""); tmpstream.clear();
    }
    delete proj;

    proj=Qmap_full[w]->ProjectionY("_pfy",xmax,xmax);
    tmpstream<<"SNR vs frequency at GPS = "<<fixed<<setprecision(3)<<Qmap_full[w]->GetXaxis()->GetBinCenter(xmax)+Qmap_center;
    proj->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
    proj->GetXaxis()->SetTitle("Frequency [Hz]");
    proj->GetYaxis()->SetTitle("SNR");
    proj->GetXaxis()->SetTitleOffset(1.1);
    proj->GetXaxis()->SetLabelSize(0.045);
    proj->GetYaxis()->SetLabelSize(0.045);
    proj->GetXaxis()->SetTitleSize(0.045);
    proj->GetYaxis()->SetTitleSize(0.045);
    GPlot->SetLogx(1);
    GPlot->Draw(proj);
    for(int f=0; f<(int)form.size(); f++){
      tmpstream<<fOutdir[MapChNumber]<<"/"<<fChannels[MapChNumber]<<"_projf_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();
      tmpstream<<fOutdir[MapChNumber]<<"/"<<fChannels[MapChNumber]<<"_projfth_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->Print(tmpstream.str(),0.5);
      tmpstream.str(""); tmpstream.clear();
    }
    delete proj;
  }

  // save map summary
  TFile summary((fOutdir[MapChNumber]+"/"+fChannels[MapChNumber]+"_mapsummary.root").c_str(),"RECREATE");
  summary.cd();
  MS->Write();
  summary.Close();
  delete MS;

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
      fpsd=new TFile((fOutdir[c]+"/"+fChannels[c]+"_data.root").c_str(),"RECREATE");
      first_save[c]=false;
    }
    else fpsd=new TFile((fOutdir[c]+"/"+fChannels[c]+"_data.root").c_str(),"UPDATE");
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
      if(!type.compare("ASD")) ss<<fOutdir[c]<<"/"<<fChannels[c]<<"_asd_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"."<<form[f];
      else ss<<fOutdir[c]<<"/"<<fChannels[c]<<"_psd_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"."<<form[f];
      GPlot->Print(ss.str().c_str());
      ss.str(""); ss.clear();
      if(!type.compare("ASD")) ss<<fOutdir[c]<<"/"<<fChannels[c]<<"_asd_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"th."<<form[f];
      else ss<<fOutdir[c]<<"/"<<fChannels[c]<<"_psd_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd()<<"th."<<form[f];
      GPlot->Print(ss.str(),0.5);
      ss.str(""); ss.clear();
    }
  }
  
  form.clear();
  delete GAPSD;
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::SaveTS(const int c, double tcenter){
////////////////////////////////////////////////////////////////////////////////////

  TGraph *GDATA = new TGraph(ChunkSize);
  if(GDATA==NULL) return;

  stringstream ss;
  ss<<"ts_"<<fChannels[c]<<"_"<<dataseq->GetChunkTimeStart()<<"_"<<dataseq->GetChunkTimeEnd();
  GDATA->SetName(ss.str().c_str());
  ss.str(""); ss.clear();

  for(int i=0; i<ChunkSize; i++) GDATA->SetPoint(i,(double)dataseq->GetChunkTimeStart()+(double)i*(double)(dataseq->GetChunkTimeEnd()-dataseq->GetChunkTimeStart())/(double)ChunkSize-tcenter,ChunkVect[i]);
     
  // cosmetics
  GPlot->SetLogx(0);
  GPlot->SetLogy(0);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);
  GPlot->Draw(GDATA,"APL");
  //GDATA->GetXaxis()->SetTimeOffset(-toffset);
  GDATA->GetHistogram()->SetXTitle("Time [s]");
  GDATA->GetHistogram()->SetYTitle("Amplitude [?]");
  GDATA->SetTitle((fChannels[c]+": amplitude time series").c_str());
  GDATA->SetLineWidth(1);
  //GDATA->GetXaxis()->SetLimits(dataseq->GetChunkTimeStart()+toffset,dataseq->GetChunkTimeEnd()+toffset);
  GDATA->GetXaxis()->SetNoExponent();
  GDATA->GetXaxis()->SetTitleOffset(1.1);
  GDATA->GetXaxis()->SetLabelSize(0.045);
  GDATA->GetYaxis()->SetLabelSize(0.045);
  GDATA->GetXaxis()->SetTitleSize(0.045);
  GDATA->GetYaxis()->SetTitleSize(0.045);

  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fdata;
    if(first_save[c]){
      fdata=new TFile((fOutdir[c]+"/"+fChannels[c]+"_data.root").c_str(),"RECREATE");
      first_save[c]=false;
    }
    else fdata=new TFile((fOutdir[c]+"/"+fChannels[c]+"_data.root").c_str(),"UPDATE");
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
    
    if(!tcenter) tcenter=(double)(dataseq->GetChunkTimeStart()+dataseq->GetChunkTimeEnd())/2.0;
    else{
      ss<<fChannels[c]+": amplitude time series centred at "<<fixed<<setprecision(3)<<tcenter;
      tcenter=0;
      GDATA->SetTitle(ss.str().c_str());
      ss.str(""); ss.clear();
    }

    // zoom
    for(int w=(int)fWindows.size()-1; w>=0; w--){
      GDATA->GetXaxis()->SetLimits(tcenter-(double)fWindows[w]/2.0,tcenter+(double)fWindows[w]/2.0);
      for(int f=0; f<(int)form.size(); f++){
	ss<<fOutdir[c]<<"/"<<fChannels[c]<<"_ts_dt"<<fWindows[w]<<"."<<form[f];
	GPlot->Print(ss.str().c_str());
	ss.str(""); ss.clear();
	ss<<fOutdir[c]<<"/"<<fChannels[c]<<"_tsth_dt"<<fWindows[w]<<"."<<form[f];
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
void Omicron::PrintStatusInfo(void){
////////////////////////////////////////////////////////////////////////////////////

  
  cout<<"\n************* Omicron status info *************"<<endl;
  cout<<"requested start = "<<inSegments->GetStart(0)<<endl;
  cout<<"requested end   = "<<inSegments->GetEnd(inSegments->GetNsegments()-1)<<endl;
  cout<<"requested livetime   = "<<inSegments->GetLiveTime()<<"s"<<endl;
  for(int c=0; c<(int)fChannels.size(); c++){
    cout<<"\n*** "<<fChannels[c]<<endl;
    if(outSegments[c]->GetNsegments()){
      cout<<"start_out                  = "<<outSegments[c]->GetStart(0)<<endl;
      cout<<"end_out                    = "<<outSegments[c]->GetEnd(outSegments[c]->GetNsegments()-1)<<endl;
      cout<<"processed livetime         = "<<outSegments[c]->GetLiveTime()<<"s ("<<outSegments[c]->GetLiveTime()/inSegments->GetLiveTime()*100<<"%)"<<endl;
    }
    else{
      cout<<"start_out                    = -1"<<endl;
      cout<<"end_out                      = -1"<<endl;
      cout<<"processed livetime           = "<<"0s (0%)"<<endl;
    }
    cout<<"number of chunks             = "<<chunk_ctr[c]<<endl;
    cout<<"number of corrupted chunks   = "<<cor_chunk_ctr[c]<<endl;
    cout<<"number of unprocessed chunks = "<<cor_data_ctr[c]<<endl;
    cout<<"number of maxed-out chunks   = "<<max_chunk_ctr[c]<<endl;
  }
  cout<<"***********************************************\n"<<endl;

  return;
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
  cout<<"                                                               V1R3"<<endl;
  cout<<"############################################################################"<<endl;
  cout<<"############################################################################"<<endl;
  cout<<endl;
  cout<<endl;

  return;
}
