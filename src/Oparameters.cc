//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

////////////////////////////////////////////////////////////////////////////////////
void Omicron::ReadOptions(void){
////////////////////////////////////////////////////////////////////////////////////

  // check that the option file exists
  if(!IsTextFile(fOptionFile)){
    cerr<<"Omicron::ReadOptions: option file "<<fOptionFile<<" cannot be found"<<endl;
    return;
  }
  

  // create parser
  IO *io = new IO(fOptionFile.c_str());


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------                  DATA                  --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //***** ffl file *****
  if(!io->GetOpt("DATA","FFL", fFflFile)&&!io->GetOpt("DATA","LCF", fFflFile)) fFflFile="none";
  //*****************************
  
  //***** List of channels *****
  if(!io->GetAllOpt("DATA","CHANNELS", fChannels)){
    cerr<<"Omicron::ReadOptions: a list of channels is required (DATA/CHANNELS)"<<endl;
    fChannels.push_back("M1:MISSING");
    status_OK=false;
  }
  //*****************************
  
  //***** Sampling frequency *****
  if(!io->GetOpt("DATA","SAMPLEFREQUENCY", fSampleFrequency)){
    cerr<<"Omicron::ReadOptions: a working sampling frequency is required (DATA/SAMPLEFREQUENCY)"<<endl;
    fSampleFrequency=16;
    status_OK=false;
  }
  if(fSampleFrequency<16){
    cerr<<"Omicron::ReadOptions: the working sampling frequency (DATA/SAMPLEFREQUENCY) should be at least 16Hz"<<endl;
    fSampleFrequency=16;
    status_OK=false;
  }
  //*****************************

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------                 OUTPUT                 --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //***** output directory *****
  if(!io->GetOpt("OUTPUT","DIRECTORY", fMaindir)){
    cerr<<"Omicron::ReadOptions: No output directory (OUTPUT/DIRECTORY)               --> set default: current"<<endl;
    fMaindir=".";
  }
  if(!IsDirectory(fMaindir)){
    cerr<<"Omicron::ReadOptions: output directory cannot be found: OUTPUT/DIRECTORY   --> set default: current"<<endl;
    fMaindir=".";
  }
  //*****************************
   
  //***** set verbosity *****
  if(!io->GetOpt("OUTPUT","VERBOSITY", fVerbosity)) fVerbosity=0;
  //*****************************

  //***** set output products *****
  if(!io->GetOpt("OUTPUT","PRODUCTS", fOutProducts)){
    cerr<<"Omicron::ReadOptions: No output products (OUTPUT/PRODUCTS)                 --> set default: triggers"<<endl;
    fOutProducts="triggers";
  }
  //*****************************

  //***** set output format ***** 
  if(!io->GetOpt("OUTPUT","FORMAT", fOutFormat)){
    cerr<<"Omicron::ReadOptions: No output format (OUTPUT/FORMAT)                     --> set default: root"<<endl;
    fOutFormat="root";
  }
  //*****************************

  //***** set output style *****
  if(!io->GetOpt("OUTPUT","STYLE", fOutStyle)){
    cerr<<"Omicron::ReadOptions: No output products (OUTPUT/STYLE)                    --> set default: GWOLLUM"<<endl;
    fOutStyle="GWOLLUM";
  }
  //*****************************



  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------                PARAMETER               --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //***** timing *****
  if(!io->GetOpt("PARAMETER","CHUNKDURATION", fChunkDuration)){
    cerr<<"Omicron::ReadOptions: No search chunk duration (PARAMETER/CHUNKDURATION)   --> set default: 512s"<<endl;
    fChunkDuration=0;// will be updated by Odata
  }
  if(!io->GetOpt("PARAMETER","SEGMENTDURATION", fSegmentDuration)){
    cerr<<"Omicron::ReadOptions: No search segment duration (PARAMETER/SEGMENTDURATION)   --> set default: 64s"<<endl;
    fSegmentDuration=0;// will be updated by Odata
  }
  if(!io->GetOpt("PARAMETER","OVERLAPDURATION", fOverlapDuration)){
    cerr<<"Omicron::ReadOptions: No search overlap duration (PARAMETER/OVERLAPDURATION) --> set default: 8s"<<endl;
    fOverlapDuration=0;// will be updated by Odata
  }
  //*****************************

  //***** Frequency range *****
  fFreqRange.clear();
  if(!io->GetOpt("PARAMETER","FREQUENCYRANGE", fFreqRange)){
    cerr<<"Omicron::ReadOptions: No search frequency range (PARAMETER/FREQUENCYRANGE) --> set default: 2-"<<fSampleFrequency/2<<" Hz"<<endl;
    fFreqRange.push_back(2); fFreqRange.push_back(fSampleFrequency/2);
  }
  if(fFreqRange[1]>(double)fSampleFrequency/2.0) fFreqRange[1]=(double)fSampleFrequency/2.0;
  if(fFreqRange.size()!=2||fFreqRange[0]>=fFreqRange[1]){
    cerr<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) is not correct --> set default: 2-"<<fSampleFrequency/2<<" Hz"<<endl;
    fFreqRange.clear();
    fFreqRange.push_back(2); fFreqRange.push_back(fSampleFrequency/2);
  }
  //*****************************

  //***** Q range *****
  fQRange.clear();
  if(!io->GetOpt("PARAMETER","QRANGE", fQRange)){
    cerr<<"Omicron::ReadOptions: No search Q range (PARAMETER/QRANGE)                 --> set default: 4-100"<<endl;
    fQRange.push_back(4);  fQRange.push_back(100);
  }
  if(fQRange.size()!=2||fQRange[0]>=fQRange[1]||fQRange[0]<=0.0){
    cerr<<"Omicron::ReadOptions: Q range (PARAMETER/QRANGE) is not correct            --> set default: 4-100"<<endl;
    fQRange.clear();
    fQRange.push_back(4);  fQRange.push_back(100);
  }
  //*****************************

  //***** maximum mismatch *****
  if(!io->GetOpt("PARAMETER","MISMATCHMAX", fMismatchMax)){
    cerr<<"Omicron::ReadOptions: No mismatch (PARAMETER/MISMATCHMAX)                  --> set default: 0.2"<<endl;
    fMismatchMax=0.2;
  }
  if(fMismatchMax<=0){
    cerr<<"Omicron::ReadOptions: mismatch (PARAMETER/MISMATCHMAX) is not correct      --> set default: 0.2"<<endl;
    fMismatchMax=0.2;
  }
  //*****************************
  
  //***** SNR Threshold *****
  if(!io->GetOpt("PARAMETER","SNRTHRESHOLD", fSNRThreshold)){
    cerr<<"Omicron::ReadOptions: No SNR threshold (PARAMETER/SNRTHRESHOLD)            --> set default: 8"<<endl;
    fSNRThreshold=8.0;
  }
  //*****************************
  
  //***** set clustering *****
  if(!io->GetOpt("PARAMETER","CLUSTERING", fClusterAlgo)){
    cerr<<"Omicron::ReadOptions: No clustering (PARAMETER/CLUSTERING)                 --> set default: none"<<endl;
    fClusterAlgo="none";
  }
  if(!io->GetOpt("PARAMETER","CLUSTERDT", fcldt)) fcldt=0.1;
  //*****************************
  
  //***** Down-tiling *****
  if(!io->GetOpt("PARAMETER","TILEDOWN", fTileDown)){
    cerr<<"Omicron::ReadOptions: No downtiling option (PARAMETER/TILEDOWN)            --> set default: NO"<<endl;
    fTileDown=0;
  }
  //*****************************
  
  //***** scan windows *****
  io->GetOpt("PARAMETER","WINDOWS", fWindows);
  //*****************************

  //***** vertical scale *****
  if(!io->GetOpt("PARAMETER","SNRSCALE", fsnrscale)){
    cerr<<"Omicron::ReadOptions: No snr scale option (PARAMETER/SNRSCALE)             --> set default: 50"<<endl;
    fsnrscale=50;
  }
  //*****************************


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------               INJECTIONS               --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //***** injection channels *****
  if(io->GetOpt("INJECTION","CHANNELS", fInjChan)){
    if(fInjChan.size()!=fChannels.size()){
      cerr<<"Omicron::ReadOptions: INJECTION/CHANNELS is inconsistent with the number of channels"<<endl;
      fInjChan.clear();
    }
    if(io->GetOpt("INJECTION","FACTORS", fInjFact)){
      if(fInjFact.size()!=fInjChan.size()){
	cerr<<"Omicron::ReadOptions: INJECTION/FACTORS is inconsistent with the number of channels"<<endl;
	fInjChan.clear();
	fInjFact.clear();
      }
    }
    else{
      for(int i=0; i<(int)fChannels.size(); i++) fInjFact.push_back(1.0);
    }
  }
  //*****************************

  //***** software injections *****
  if(!io->GetOpt("INJECTION","FILENAME", fInjFilePat)) fInjFilePat="none";
  //*****************************

  // dump options
  if(fVerbosity>1) io->Dump(cout);
  delete io;

  return;
}
