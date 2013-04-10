//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

ClassImp(Omicron)

////////////////////////////////////////////////////////////////////////////////////
Omicron::Omicron(Segments *aSegments, const string aOptionFile){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // init
  status_OK=true;
  tiling_OK=false;
  
  // segments
  fSegments = aSegments;

  // read option file
  fOptionFile=aOptionFile;
  status_OK*=ReadOptions();

  // nulling structures
  tile=NULL;
  c_data[0] = NULL; c_data[1] = NULL;
  for(int c=0; c<(int)fChannels.size(); c++){
    odata[c]=NULL;
    triggers[c]=NULL;
  }

  // conditioned data vector
  if(status_OK){
    c_data[0] = new double [fSegmentDuration*fSampleFrequency/2]; // real part
    c_data[1] = new double [fSegmentDuration*fSampleFrequency/2]; // imaginary part
  }

  // init trigger objects
  if(status_OK){
    for(int c=0; c<(int)fChannels.size(); c++){
      triggers[c] = new Triggers(fOutdir[c],fChannels[c],fOutFormat,fVerbosity);// construct
      triggers[c]->SetNtriggerMax(fNtriggerMax);// maximum number of triggers per file
      if(!fClusterAlgo.empty()){
	status_OK*=triggers[c]->SetClustering(fClusterAlgo);// set clustering
	status_OK*=triggers[c]->SetClusterDeltaT(fcldt);// set dt
      }
    }
  }

  // init data objects
  if(status_OK){
    for(int c=0; c<(int)fChannels.size(); c++){
      odata[c] = new Odata(fFflFile,fChannels[c],fNativeFrequency[c],fSampleFrequency,fSegments,fChunkDuration,fSegmentDuration,fOverlapDuration,fFreqRange[0],fVerbosity);
      if(fInjChan.size()) status_OK*=odata[c]->SetInjectionChannel(fInjChan[c],fInjFact[c]);// add injection channel
      status_OK*=odata[c]->GetStatus();// update status
    }
  }

  if(!status_OK) cerr<<"Omicron::Omicron: Initialization failed"<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
Omicron::~Omicron(void){
////////////////////////////////////////////////////////////////////////////////////
  cout<<"delete Omicron"<<endl;

  if(fVerbosity) cout<<" -> delete tiling"<<endl;
  if(tile!=NULL) delete tile;
  if(fVerbosity>0) cout<<" -> delete data containers"<<endl;    
  if(c_data[0]!=NULL) delete c_data[0];
  if(c_data[1]!=NULL) delete c_data[1];
  for(int c=0; c<(int)fChannels.size(); c++){
    if(fVerbosity>0) cout<<" -> delete data for "<<fChannels[c]<<endl;
    if(odata[c]!=NULL) delete odata[c];
    if(fVerbosity>0) cout<<" -> delete triggers for "<<fChannels[c]<<endl;    
    if(triggers[c]!=NULL) delete triggers[c];
  }

}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::Process(){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::Process: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(!tiling_OK){
    cerr<<"Omicron::Process: the Otile object was not created or is corrupted"<<endl;
    return false;
  }

  // if requested, write whitened data
  TFile  *FWhite[NDATASTREAMS];
  TGraph *GWhite;
  if(writewhiteneddata){
    for(int c=0; c<(int)fChannels.size(); c++){
      FWhite[c] = new TFile((fOutdir[c]+"/white_"+fChannels[c]+".root").c_str(), "RECREATE");
    }
  }

  int psdsize;
  
  // looping status
  int keep_looping=0; 
  // 0 means keep calling new chunks of data
  // 1 means stop calling new chunks: end of data to process
  // 2 means the data chunk is corrupted -> move to the next one

  // loop over chunks
  while(keep_looping!=1){

    // loop over channels
    for(int c=0; c<(int)fChannels.size(); c++){
      if(fVerbosity) cout<<"processing channel "<<fChannels[c]<<"..."<<endl;
      
      // get new data chunk
      keep_looping=odata[c]->NewChunk();// get fresh data
      if(keep_looping==1) break;        // end of data to process
      if(keep_looping>1) continue;      // the data chunk is corrupted -> skip
  
      // write chunk info if requested
      if(writetimeseries) status_OK*=odata[c]->WriteTimeSeries(fOutdir[c]);
      if(writepsd) status_OK*=odata[c]->WritePSD(fOutdir[c]);

      //loop over segments
      for(int s=0; s<odata[c]->GetNSegments(); s++){

	// get conditioned data
	if(!odata[c]->GetConditionedData(s,c_data[0],c_data[1],&psd,psdsize)){
	  cerr<<"Omicron::Process: conditionned data are corrupted!"<<endl;
	  return false;
	}

	// set new power for this chunk
	if(!s) tile->SetPowerSpectrum(psd,psdsize);

	// save it if requested
	if(writewhiteneddata){
	  FWhite[c]->cd();
	  GWhite = new TGraph(fSegmentDuration*fSampleFrequency);
	  ostringstream graph_name;
	  graph_name<<"White_"<<odata[c]->GetSegmentTimeStart(s);
	  GWhite->SetName(graph_name.str().c_str());
	  for(int i=0; i<fSegmentDuration*fSampleFrequency/2; i++) 
	    GWhite->SetPoint(i,(double)i/(double)fSegmentDuration,fabs(c_data[0][i]));
	  GWhite->Write();
	  delete GWhite;
	}

	//get triggers
	cout<<" "<<fChannels[c]<<" Extracting triggers in "<<odata[c]->GetSegmentTimeStart(s)+fOverlapDuration/2<<"-"<<odata[c]->GetSegmentTimeStart(s)+fSegmentDuration-fOverlapDuration/2<<endl;
	if(!tile->GetTriggers(triggers[c],c_data[0],c_data[1],odata[c]->GetSegmentTimeStart(s))){
	  cerr<<"Omicron::Process: could not get trigger for channel "<<fChannels[c]
	      <<" in segment starting at "<<odata[c]->GetSegmentTimeStart(s)<<endl;
	  continue;
	}
	else
	  triggers[c]->AddSegment(odata[c]->GetSegmentTimeStart(s)+fOverlapDuration/2,odata[c]->GetSegmentTimeStart(s)+fSegmentDuration-fOverlapDuration/2);

	//stop getting triggers if max flag
	if(triggers[c]->GetMaxFlag()) break;
      }
      
      // don't save if max flag
      if(triggers[c]->GetMaxFlag()){
	cerr<<"Omicron::Process: channel "<<fChannels[c]<<" is maxed-out. This chunk is not saved"<<endl;
	triggers[c]->Reset();
	continue;
      }
      
      // save triggers
      if(!triggers[c]->Write("ALL", "default")){
	cerr<<"Omicron::Process: writing events failed for channel "<<fChannels[c]<<endl;
	return false;
      }
    }
  }

  // write whiten data
  if(writewhiteneddata)
    for(int c=0; c<(int)fChannels.size(); c++) FWhite[c]->Close();
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
int Omicron::ProcessOnline(const int aChNumber, FrVect *aVect){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::ProcessOnline: the Omicron object is corrupted"<<endl;
    return -1;
  }
  if(!tiling_OK){
    cerr<<"Omicron::ProcessOnline: the Otile object was not created or is corrupted"<<endl;
    return -2;
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
  if(writepsd) status_OK*=odata[aChNumber]->WritePSD(fOutdir[aChNumber]);

  // get conditioned data
  int psdsize;
  if(!odata[aChNumber]->GetConditionedData(0,c_data[0],c_data[1],&psd,psdsize)){
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
	
  // don't save if max flag
  if(triggers[aChNumber]->GetMaxFlag()){
    cerr<<"Omicron::ProcessOnline: channel "<<fChannels[aChNumber]<<" is maxed-out. This chunk is not saved"<<endl;
    triggers[aChNumber]->Reset();
    return 4;
  }

  // save triggers
  if(!triggers[aChNumber]->Write("ALL","default")){
    cerr<<"Omicron::ProcessOnline: writing events failed for channel "<<fChannels[aChNumber]<<endl;
    return 5;
  }
      
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::MakeTiling(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::MakeTiling: the Omicron object is corrupted"<<endl;
    return false;
  }
  tile = new Otile(fSegmentDuration,fOverlapDuration/2,fQRange[0],fQRange[1],fFreqRange[0],fFreqRange[1],fSampleFrequency,fMismatchMax,fSNRThreshold,fVerbosity);
  tiling_OK=tile->GetStatus();
  return tiling_OK;
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

  fOptions = new IO(fOptionFile.c_str());
  IO& io = *fOptions;

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
  //*****************************

  //***** List of channels *****
  fChannels.clear();
  io.GetOpt("DATA","CHANNELS", fChannels);
  if(fChannels.size()<=0){
    cerr<<"Omicron::ReadOptions: you must provide a list of channels: DATA/CHANNELS"<<endl;
    return false;
  }
  if(fChannels.size()>NDATASTREAMS){
    cerr<<"Omicron::ReadOptions: The number of channels cannot exceed "<<NDATASTREAMS<<endl;
    return false;
  }
  //*****************************

  //***** injection channels *****
  fInjChan.clear();
  io.GetOpt("INJECTION","CHANNELS", fInjChan);
  if(fInjChan.size()&&fInjChan.size()!=fChannels.size()){
    cerr<<"Omicron::ReadOptions: INJECTION/CHANNELS is inconsistent with the number of channels"<<endl;
    return false;
  }
  fInjFact.clear();
  io.GetOpt("INJECTION","FACTORS", fInjFact);
  if(fInjFact.size()!=fInjChan.size()){
    cerr<<"Omicron::ReadOptions: INJECTION/FACTORS is inconsistent with the number of channels"<<endl;
    return false;
  }
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
    LCF2FFL(fLcfFile,maindir+"/converted_lcf.ffl");
    fFflFile=maindir+"/converted_lcf.ffl";
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
  //*****************************
    
  //***** Native channel frequencies *****
  fNativeFrequency.clear();
  io.GetOpt("DATA","NATIVEFREQUENCY", fNativeFrequency);
  if(fNativeFrequency.size()!=fChannels.size()){
    cerr<<"Omicron::ReadOptions: Each channel must be given a native frequency: PARAMETER/NATIVEFREQUENCY"<<endl;
    return false;
  }
  //*****************************

  //***** Sampling frequency *****
  fSampleFrequency=-1;
  io.GetOpt("DATA","SAMPLEFREQUENCY", fSampleFrequency);
  if(fSampleFrequency<=0||fSampleFrequency>20000){
    cerr<<"Omicron::ReadOptions: Sampling frequency "<<fSampleFrequency<<" (PARAMETER/SAMPLEFREQUENCY) is not reasonable"<<endl;
    return false;
  }
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
  //*****************************

  //***** Q range *****
  fQRange.clear();
  io.GetOpt("PARAMETER","QRANGE", fQRange);
  if(fQRange.size()!=2||fQRange[0]>fQRange[1]){
    cerr<<"Omicron::ReadOptions: Q range (PARAMETER/QRANGE) is not correct"<<endl;
    return false;
  }
  //*****************************

  //***** timing *****
  fChunkDuration=288;
  fSegmentDuration=64;
  fOverlapDuration=8;
  io.GetOpt("PARAMETER","CHUNKDURATION", fChunkDuration);
  io.GetOpt("PARAMETER","BLOCKDURATION", fSegmentDuration);
  io.GetOpt("PARAMETER","OVERLAPDURATION", fOverlapDuration);
  if(fChunkDuration<4||fChunkDuration>1000||fChunkDuration%2){
    cerr<<"Omicron::ReadOptions: Chunk duration (PARAMETER/CHUNKDURATION) is not reasonable"<<endl;
    return false;
  }
  if(fFreqRange[0]<1.0/(double)fChunkDuration){
    cerr<<"Omicron::ReadOptions: Chunk duration (PARAMETER/CHUNKDURATION) is not reasonable given your frequency range"<<endl;
    return false;
  }
  if(fSegmentDuration<4||fSegmentDuration%2||fSegmentDuration>fChunkDuration){
    cerr<<"Omicron::ReadOptions: Block duration (PARAMETER/BLOCKDURATION) is not reasonable"<<endl;
    return false;
  }
  if (!(fSegmentDuration > 0 && !(fSegmentDuration & (fSegmentDuration - 1) ))){
    cerr<<"Omicron::ReadOptions: Block duration (PARAMETER/BLOCKDURATION) should be a power of 2"<<endl;
    return false;
  }
  if(fOverlapDuration<0||fOverlapDuration>fSegmentDuration/2||fOverlapDuration%2){
    cerr<<"Omicron::ReadOptions: Overlap duration (PARAMETER/OVERLAPDURATION) is not reasonable"<<endl;
    return false;
  }
  if(fSegmentDuration<2*fOverlapDuration){
    cerr<<"Omicron::ReadOptions: Block duration (PARAMETER/BLOCKDURATION) must be at least twice longer than Overlap duration"<<endl;
    return false;
  }
  if(fOverlapDuration < (int)ceil(16.0/fFreqRange[0])){
    cerr<<"Omicron::ReadOptions: Overlap duration must be larger than "<<(int)ceil(16.0/fFreqRange[0])<<"s to remove high-pass transients"<<endl;
    return false;
  }
  if((fChunkDuration-fOverlapDuration)%(fSegmentDuration-fOverlapDuration)){
    cerr<<"Omicron::ReadOptions: inconsistency in the segments structure, you should use:"<<endl;
    cerr<<"PARAMETER/CHUNKDURATION:   "<<(fChunkDuration-fOverlapDuration)/(fSegmentDuration-fOverlapDuration)*(fSegmentDuration-fOverlapDuration)+fOverlapDuration<<endl;
    cerr<<"PARAMETER/BLOCKDURATION:   "<<fSegmentDuration<<endl;
    cerr<<"PARAMETER/OVERLAPDURATION: "<<fOverlapDuration<<endl;
    return false;
  }
  //*****************************

  //***** maximum mismatch *****
  fMismatchMax=0.2;
  io.GetOpt("PARAMETER","MISMATCHMAX", fMismatchMax);
  if(fMismatchMax<=0||fMismatchMax>0.5){
    cerr<<"Omicron::ReadOptions: maximum mismatch (PARAMETER/MISMATCHMAX) is not reasonable"<<endl;
    return false;
  }
  //*****************************

  //***** SNR Threshold *****
  fSNRThreshold=8.0;
  io.GetOpt("TRIGGER","SNRTHRESHOLD", fSNRThreshold);
  //*****************************

  //***** set trigger limit *****
  fNtriggerMax=1000000;
  io.GetOpt("TRIGGER","NMAX", fNtriggerMax);
  //*****************************

  //***** set clustering *****
  fcldt=0.1;
  io.GetOpt("TRIGGER","CLUSTERING", fClusterAlgo);
  io.GetOpt("TRIGGER","CLUSTERDT", fcldt);
  //*****************************

  //***** set verbosity *****
  fVerbosity=0;
  io.GetOpt("OUTPUT","VERBOSITY", fVerbosity);
  //*****************************

  //***** set verbosity ***** 
  io.GetOpt("OUTPUT","FORMAT", fOutFormat);
  if(fOutFormat.compare("root")&&fOutFormat.compare("xml")&&fOutFormat.compare("txt")){
    cerr<<"Omicron::ReadOptions: possible output formats: root, xml or txt"<<endl;
    cerr<<"                      --> root format will be used"<<endl;
    fOutFormat="root";
  }
  //*****************************

  //***** set writing flags *****
  writepsd=0, writetimeseries=0, writewhiteneddata=0;
  io.GetOpt("OUTPUT","WRITEPSD", writepsd);
  io.GetOpt("OUTPUT","WRITETIMESERIES", writetimeseries);
  io.GetOpt("OUTPUT","WRITEWHITENEDDATA", writewhiteneddata);
  //*****************************

  // create output directories
  for(int c=0; c<(int)fChannels.size(); c++){
    fOutdir[c]=maindir+"/"+fChannels[c];
    system(("mkdir -p "+fOutdir[c]).c_str());
  }

  io.Dump(cout);

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
