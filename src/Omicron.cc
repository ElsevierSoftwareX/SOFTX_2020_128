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

  // data vector
  if(status_OK){
    c_data[0] = new double [fSegmentDuration*fSampleFrequency/2];
    c_data[1] = new double [fSegmentDuration*fSampleFrequency/2];
  }

  // trigger objects
  if(status_OK){
    for(int c=0; c<(int)fChannels.size(); c++){
      triggers[c] = new Triggers(fOutdir[c],fChannels[c],fOutFormat,fVerbosity);
      triggers[c]->SetNtriggerMax(fNtriggerMax);
      if(!fClusterAlgo.empty()){
	triggers[c]->SetClustering(fClusterAlgo);
	triggers[c]->SetClusterDeltaT(fcldt);
      }
    }
  }

  // Init data objects
  vector<string> s_chan;
  if(status_OK){
    for(int c=0; c<(int)fChannels.size(); c++){
      if(fChanSum.size()&&!fChanSum[0].compare(fChannels[c]))
	odata[c] = new Odata(fFflFile,fChanSum,fNativeFrequency[c],fSampleFrequency,fSegments,fChunkDuration,fSegmentDuration,fOverlapDuration,fFreqRange[0],fVerbosity);
      else{
	s_chan.push_back(fChannels[c]);
	odata[c] = new Odata(fFflFile,s_chan,fNativeFrequency[c],fSampleFrequency,fSegments,fChunkDuration,fSegmentDuration,fOverlapDuration,fFreqRange[0],fVerbosity);
      }
      s_chan.clear();
      status_OK*=odata[c]->GetStatus();
    }
  }

  if(!status_OK)
    cerr<<"Omicron::Omicron: Initialization failed"<<endl;
   
}

////////////////////////////////////////////////////////////////////////////////////
Omicron::~Omicron(void){
////////////////////////////////////////////////////////////////////////////////////
  cout<<"delete Omicron"<<endl;

  if(fVerbosity>0) cout<<" -> delete tiling"<<endl;
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
  
  // loop over chunks
  int keep_looping=0;
  while(keep_looping!=1){

    // loop over channels
    for(int c=0; c<(int)fChannels.size(); c++){
      if(fVerbosity>0) cout<<"processing channel "<<fChannels[c]<<"..."<<endl;
      
      // get new data chunk
      keep_looping=odata[c]->NewChunk();
      if(keep_looping==1) break;
      if(keep_looping>1) continue;
  
      // write chunk info if requested
      if(writetimeseries) status_OK*=odata[c]->WriteTimeSeries(fOutdir[c]);
      if(writepsd) status_OK*=odata[c]->WritePSD(fOutdir[c]);

      //loop over segments
      for(int s=0; s<odata[c]->GetNSegments(); s++){

	// get conditioned data
	if(!odata[c]->GetConditionedData(s,c_data[0],c_data[1])){
	  cerr<<"Omicron::Process: conditionned data are corrupted!"<<endl;
	  return false;
	}

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
  if(!odata[aChNumber]->GetConditionedData(0,c_data[0],c_data[1])){
    cerr<<"Omicron::ProcessOnline: conditionned data are corrupted!"<<endl;
    return 2;
  }
  
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
    cerr<<"Omicron::ReadOptions: option file cannot be found"<<endl;
    return false;
  }

  fOptions = new IO(fOptionFile.c_str());
  IO& io = *fOptions;

  // output directory
  string maindir;
  io.GetOpt("OUTPUT","DIRECTORY", maindir);
  if(maindir.empty()||!IsDirectory(maindir)){
    cerr<<"Omicron::ReadOptions: output directory cannot be found"<<endl;
    return false;
  }

  // List of channels
  fChannels.clear();
  io.GetOpt("DATA","CHANNELS", fChannels);
  if(fChannels.size()<=0){
    cerr<<"Omicron::ReadOptions: The list of channels (DATA/CHANNELS) has to be given"<<endl;
    return false;
  }
  if(fChannels.size()>NDATASTREAMS){
    cerr<<"Omicron::ReadOptions: The number of channels cannot exceed "<<NDATASTREAMS<<endl;
    return false;
  }

  // create output directories
  for(int c=0; c<(int)fChannels.size(); c++){
    fOutdir[c]=maindir+"/"+fChannels[c];
    system(("mkdir -p "+fOutdir[c]).c_str());
  }

  // channel sums
  fChanSum.clear();
  io.GetOpt("DATA","SUM", fChanSum);
  if(fChanSum.size()&&fChanSum.size()<3){
    cerr<<"Omicron::ReadOptions: there must be at least 2 channels to sum up"<<endl;
    return false;
  }
  if(fChanSum.size()){
    bool c_found=false;
    for(int c=0; c<(int)fChannels.size(); c++){
      if(!fChanSum[0].compare(fChannels[c])){
	c_found=true;
	break;
      }
    }
    if(!c_found){
      cerr<<"Omicron::ReadOptions: the SUM channel cannot be found"<<endl;
      return false;
    }
  }

  // ffl/cache file path
  online=false;
  io.GetOpt("DATA","FFL", fFflFile);
  io.GetOpt("DATA","CACHE", fCacheFile);
  if(!fFflFile.compare("ONLINE")){// online keyword
    online=true;
  }
  else if(!fFflFile.empty()){
    if(!IsTextFile(fFflFile)){
      cerr<<"Omicron::ReadOptions: FFL file (DATA/FFL) cannot be found"<<endl;
      return false;
    }
  }
  else{
    if(fCacheFile.empty()||!IsTextFile(fCacheFile)){
      cerr<<"Omicron::ReadOptions: FFL file (DATA/CACHE-FFL) cannot be found"<<endl;
      return false;
    }
    //convert cache to ffl
    if(!Cache2FFL(fCacheFile,maindir+"/convertedcache.ffl")){
      cerr<<"Omicron::ReadOptions: cannot convert cache file in ffl file"<<endl;
      return false;
    }
    fFflFile=maindir+"/convertedcache.ffl";
  }

  // Native channel frequencies
  fNativeFrequency.clear();
  io.GetOpt("DATA","NATIVEFREQUENCY", fNativeFrequency);
  if(fNativeFrequency.size()!=fChannels.size()){
    cerr<<"Omicron::ReadOptions: Each channel must be given a native frequency (PARAMETER/NATIVEFREQUENCY)"<<endl;
    return false;
  }

  // Sampling frequency (Def=4096)
  fSampleFrequency=-1;
  io.GetOpt("DATA","SAMPLEFREQUENCY", fSampleFrequency);
  if(fSampleFrequency<64||fSampleFrequency>20000){
    cerr<<"Omicron::ReadOptions: Sampling frequency (PARAMETER/SAMPLEFREQUENCY) is not reasonable"<<endl;
    return false;
  }

  // Frequency range
  fFreqRange.clear();
  io.GetOpt("PARAMETER","FREQUENCYRANGE", fFreqRange);
  if(fFreqRange.size()!=2||fFreqRange[0]>fFreqRange[1]){
    cerr<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) must be given"<<endl;
    return false;
  }
  if(fFreqRange[1]>fSampleFrequency/2){
    cout<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) goes beyond Nyquist frequency: "<<fFreqRange[1]<<">"<<fSampleFrequency/2<<endl;
    fFreqRange.pop_back(); fFreqRange.push_back(fSampleFrequency/2);
    return true;
  }

  // Q range
  fQRange.clear();
  io.GetOpt("PARAMETER","QRANGE", fQRange);
  if(fQRange.size()!=2||fQRange[0]>fQRange[1]){
    cerr<<"Omicron::ReadOptions: Q range (PARAMETER/QRANGE) must be given"<<endl;
    return false;
  }
 
  //***************** optional options *****************

  /* Segment structure
  
    PSD is computed over one chunk (median-mean)
    |________________chunk_________________|
    |                                   |____________chunk____________|
    |----seg----|                          |
             |----seg----|                 |
	              |----seg----|        |
		               |----seg----|
			                |--|
					overlap

    For each segment, triggers (t) are extracted as follows (overlap divided in 2):
    
    |---------seg---------|
       |ttttttttttttttt|
	               |
                    |---------seg---------|
		       |ttttttttttttttt|
  */
  fChunkDuration=288;
  fSegmentDuration=64;
  fOverlapDuration=8;
  io.GetOpt("PARAMETER","CHUNKDURATION", fChunkDuration);
  io.GetOpt("PARAMETER","BLOCKDURATION", fSegmentDuration);
  io.GetOpt("PARAMETER","OVERLAPDURATION", fOverlapDuration);
  if(fChunkDuration<4||fChunkDuration>500||fChunkDuration%2){
    cerr<<"Omicron::ReadOptions: Chunk duration (PARAMETER/CHUNKDURATION) is not reasonable"<<endl;
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
    cerr<<"Omicron::ReadOptions: Chunk duration (PARAMETER/OVERLAPDURATION) is not reasonable"<<endl;
    return false;
  }
  if((fChunkDuration-fOverlapDuration)%(fSegmentDuration-fOverlapDuration)){
    cerr<<"Omicron::ReadOptions: inconsistency in the segments structure, you should use:"<<endl;
    cerr<<"PARAMETER/CHUNKDURATION:   "<<(fChunkDuration-fOverlapDuration)/(fSegmentDuration-fOverlapDuration)*(fSegmentDuration-fOverlapDuration)+fOverlapDuration<<endl;
    cerr<<"PARAMETER/BLOCKDURATION:   "<<fSegmentDuration<<endl;
    cerr<<"PARAMETER/OVERLAPDURATION: "<<fOverlapDuration<<endl;
    return false;
  }
 
  // maximum mismatch
  fMismatchMax=0.2;
  io.GetOpt("PARAMETER","MISMATCHMAX", fMismatchMax);
  if(fMismatchMax<=0||fMismatchMax>0.5){
    cerr<<"Omicron::ReadOptions: maximum mismatch (PARAMETER/MISMATCHMAX) is not reasonable"<<endl;
    return false;
  }

  // SNR Threshold
  fSNRThreshold=5.0;
  io.GetOpt("TRIGGER","SNRTHRESHOLD", fSNRThreshold);
  
  // set trigger limit
  fNtriggerMax=1000000;
  io.GetOpt("TRIGGER","NMAX", fNtriggerMax);

  // set clustering
  fcldt=0.1;
  io.GetOpt("TRIGGER","CLUSTERING", fClusterAlgo);
  io.GetOpt("TRIGGER","CLUSTERDT", fcldt);

  // set verbosity
  fVerbosity=0;
  io.GetOpt("OUTPUT","VERBOSITY", fVerbosity);

  // set verbosity
  io.GetOpt("OUTPUT","FORMAT", fOutFormat);
  if(fOutFormat.compare("root")&&fOutFormat.compare("txt")){
    cerr<<"Omicron::ReadOptions: possible output formats: root or txt"<<endl;
    cerr<<"                      --> root format will be used"<<endl;
    fOutFormat="root";
  }


  // set writing flags
  writepsd=0, writetimeseries=0, writewhiteneddata=0;
  io.GetOpt("OUTPUT","WRITEPSD", writepsd);
  io.GetOpt("OUTPUT","WRITETIMESERIES", writetimeseries);
  io.GetOpt("OUTPUT","WRITEWHITENEDDATA", writewhiteneddata);


  io.Dump(cout);

  return true;
}
