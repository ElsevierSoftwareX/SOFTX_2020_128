//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

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
  if(!io->GetOpt("OUTPUT","DIRECTORY", fMaindir)){
    cout<<"Omicron::ReadOptions: No output directory (DATA/CHANNELS)                  --> set default: current"<<endl;
    fMaindir=".";
  }
  if(!IsDirectory(fMaindir)){
    cerr<<"Omicron::ReadOptions: output directory cannot be found: OUTPUT/DIRECTORY"<<endl;
    return false;
  }
  //*****************************

  //***** List of channels *****
  if(!io->GetAllOpt("DATA","CHANNELS", fChannels)){
    cout<<"Omicron::ReadOptions: No channel (DATA/CHANNELS)                           --> set default: V1:Pr_B1_ACp"<<endl;
    fChannels.push_back("V1:Pr_B1_ACp");
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
  if(io->GetOpt("DATA","FFL", fFflFile)) fFflFormat="ffl";
  else if(io->GetOpt("DATA","LCF", fFflFile)) fFflFormat="lcf";
  else{
    fFflFile="none";
    fFflFormat="none";
  }
  //*****************************

  //***** Sampling frequency *****
  if(!io->GetOpt("DATA","SAMPLEFREQUENCY", fSampleFrequency)){
    cout<<"Omicron::ReadOptions: No working sampling frequency (PARAMETER/SAMPLEFREQUENCY) --> set default: 2048 Hz"<<endl;
    fSampleFrequency=2048;// adjusted later
  }
  if(fSampleFrequency<16){
    cerr<<"Omicron::ReadOptions: the working sampling frequency (PARAMETER/SAMPLEFREQUENCY) should be at least 16Hz"<<endl;
    fSampleFrequency=2048;// adjusted later
    return false;
  }
  //*****************************

  //***** Frequency range *****
  fFreqRange.clear();
  if(!io->GetOpt("PARAMETER","FREQUENCYRANGE", fFreqRange)){
    cout<<"Omicron::ReadOptions: No search frequency range (PARAMETER/FREQUENCYRANGE) --> set default: 2-"<<fSampleFrequency/2<<"Hz"<<endl;
    fFreqRange.push_back(2);  fFreqRange.push_back(fSampleFrequency/2);
  }
  if(fFreqRange.size()!=2){
    cerr<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) is not correct"<<endl;
    return false;
  }
  //*****************************

  //***** Q range *****
  fQRange.clear();
  if(!io->GetOpt("PARAMETER","QRANGE", fQRange)){
    cout<<"Omicron::ReadOptions: No search Q range (PARAMETER/QRANGE)                 --> set default: 3.3167-100"<<endl;
    fQRange.push_back(3.3167);  fQRange.push_back(100);
  }
  if(fQRange.size()!=2){
    cerr<<"Omicron::ReadOptions: Q range (PARAMETER/QRANGE) is not correct"<<endl;
    return false;
  }
  //*****************************

  //***** timing *****
  if(!io->GetOpt("PARAMETER","CHUNKDURATION", fChunkDuration)){
    cout<<"Omicron::ReadOptions: No search chunk duration (PARAMETER/CHUNKDURATION)   --> set default: 512s"<<endl;
    fChunkDuration=512;
  }
  if(!io->GetOpt("PARAMETER","BLOCKDURATION", fSegmentDuration)){
    cout<<"Omicron::ReadOptions: No search block duration (PARAMETER/BLOCKDURATION)   --> set default: 64s"<<endl;
    fSegmentDuration=64;
  }
  if(!io->GetOpt("PARAMETER","OVERLAPDURATION", fOverlapDuration)){
    cout<<"Omicron::ReadOptions: No search block duration (PARAMETER/OVERLAPDURATION) --> set default: 8s"<<endl;
    fOverlapDuration=8;
  }
  //*****************************

  //***** maximum mismatch *****
  if(!io->GetOpt("PARAMETER","MISMATCHMAX", fMismatchMax)){
    cout<<"Omicron::ReadOptions: No mismatch (PARAMETER/MISMATCHMAX)                  --> set default: 0.2"<<endl;
    fMismatchMax=0.2;
  }
  //*****************************

  //***** scan windows *****
  if(!io->GetOpt("PARAMETER","WINDOWS", fWindows)){
    cout<<"Omicron::ReadOptions: No windows (PARAMETER/WINDOWS)                       --> set default: 2 8 32"<<endl;
    fWindows.push_back(2); fWindows.push_back(8); fWindows.push_back(32); 
  }
  //*****************************

  //***** SNR Threshold *****
  if(!io->GetOpt("TRIGGER","SNRTHRESHOLD", fSNRThreshold)){
    cout<<"Omicron::ReadOptions: No SNR threshold (PARAMETER/SNRTHRESHOLD)            --> set default: 8"<<endl;
    fSNRThreshold=8.0;
  }
  //*****************************
  
  //***** set trigger limit *****
  double ratemax;
  if(io->GetOpt("TRIGGER","RATEMAX", ratemax)) fNtriggerMax=(int)ceil(ratemax*fChunkDuration);
  else{
    if(!io->GetOpt("TRIGGER","NMAX", fNtriggerMax)){
      cout<<"Omicron::ReadOptions: No trigger limit (TRIGGER/NMAX)                      --> set default: 1,000,000"<<endl;
      fNtriggerMax=1000000;
    }
  }
  //*****************************

  //***** set clustering *****
  if(!io->GetOpt("TRIGGER","CLUSTERING", fClusterAlgo)){
    cout<<"Omicron::ReadOptions: No clustering (TRIGGER/CLUSTERING)                   --> set default: none"<<endl;
    fClusterAlgo="none";
  }
  if(!io->GetOpt("TRIGGER","CLUSTERDT", fcldt)) fcldt=0.1;
  //*****************************

  //***** set verbosity *****
  if(!io->GetOpt("OUTPUT","VERBOSITY", fVerbosity)) fVerbosity=1;
  //*****************************

  //***** set output format ***** 
  if(!io->GetOpt("OUTPUT","FORMAT", fOutFormat)){
    cout<<"Omicron::ReadOptions: No output format (OUTPUT/FORMAT)                     --> set default: root"<<endl;
    fOutFormat="root";
  }
  //*****************************

  //***** set plotting style ***** 
  if(!io->GetOpt("OUTPUT","PLOTSTYLE", fStyle)){
    cout<<"Omicron::ReadOptions: No plotting style (OUTPUT/PLOTSTYLE)                 --> set default: GWOLLUM"<<endl;
    fStyle="GWOLLUM";
  }
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

