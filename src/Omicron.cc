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
  dataseq = new Odata(fChunkDuration, fSegmentDuration, fOverlapDuration, fVerbosity);
  status_OK*=dataseq->GetStatus();
  dataRe = new double* [dataseq->GetNSegments()];
  dataIm = new double* [dataseq->GetNSegments()];

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
  delete triggers;
  if(FFL!=NULL) delete FFL;
  delete dataseq;
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
  if(!dataseq->SetSegments(aSeg)){
    cerr<<"Omicron::Process: cannot initiate data segments."<<endl;
    return false;
  }

  // locals
  int dsize;         // native data size
  double *dvector;   // data vector before resampling
  int s, res;
  
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

      // process data vector
      res=ConditionVector(c, dsize, dvector);
      if(res<0){
	cerr<<"Omicron::Process: conditioning failed."<<endl;
	delete dvector;// not needed anymore
	return false;
      }
      if(res>0){
	cerr<<"Omicron::Process: conditioning failed."<<endl;
	delete dvector;// not needed anymore
	continue;
      }

      delete dvector;// not needed anymore

      // Save info if requested
      if(writetimeseries) SaveData(c,ChunkVect,dataseq->GetChunkTimeStart(),dataseq->GetChunkTimeEnd());
      if(writepsd)        SavePSD(c,dataseq->GetChunkTimeStart(),dataseq->GetChunkTimeEnd());

      // get triggers above SNR threshold
      for(s=0; s<dataseq->GetNSegments(); s++){
	if(!tile->GetTriggers(triggers[c],dataRe[s],dataIm[s],dataseq->GetSegmentTimeStart(s))){
	  cerr<<"Omicron::Process: could not get triggers for channel "<<fChannels[c]<<endl;
	  return false;
	}
	else{
	  triggers[c]->AddSegment(dataseq->GetSegmentTimeStart(s)+fOverlapDuration/2,dataseq->GetSegmentTimeEnd(s)-fOverlapDuration/2);
	  delete dataRe[s];
	  delete dataIm[s];
	  if(triggers[c]->GetMaxFlag()) break;
	}
      }
      for(int ss=s+1; ss<dataseq->GetNSegments(); ss++){
	delete dataRe[s];
	delete dataIm[s];
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
	outSegments[c]->AddSegment(dataseq->GetChunkTimeStart()+fOverlapDuration/2,dataseq->GetChunkTimeEnd()-fOverlapDuration/2);
	
    }
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
  ChunkSize=fChunkDuration*fSampleFrequency;
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
