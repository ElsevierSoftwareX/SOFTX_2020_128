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
  //--------------                 OUTPUT                 --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //***** output directory *****
  if(!io->GetOpt("OUTPUT","DIRECTORY", fMaindir)){
    cerr<<"Omicron::ReadOptions: No output directory (OUTPUT/DIRECTORY)  --> set default: current"<<endl;
    fMaindir=".";
  }
  if(!IsDirectory(fMaindir)){
    cerr<<"Omicron::ReadOptions: output directory cannot be found (OUTPUT/DIRECTORY)  --> set default: current"<<endl;
    fMaindir=".";
  }
  //*****************************
   
  //***** verbosity *****
  if(!io->GetOpt("OUTPUT","VERBOSITY", fVerbosity)) fVerbosity=0;
  //*****************************

  //***** set output products *****
  if(!io->GetOpt("OUTPUT","PRODUCTS", fOutProducts)){
    cerr<<"Omicron::ReadOptions: No output products (OUTPUT/PRODUCTS)  --> set default: triggers"<<endl;
    fOutProducts="triggers";
  }
  //*****************************

  //***** set output format ***** 
  if(!io->GetOpt("OUTPUT","FORMAT", fOutFormat)){
    cerr<<"Omicron::ReadOptions: No output format (OUTPUT/FORMAT)  --> set default: root"<<endl;
    fOutFormat="root";
  }
  //*****************************

  //***** set output style *****
  string outstyle;
  if(!io->GetOpt("OUTPUT","STYLE", outstyle)){
    cerr<<"Omicron::ReadOptions: No output products (OUTPUT/STYLE)  --> set default: GWOLLUM"<<endl;
    outstyle="GWOLLUM";
  }
  GPlot = new GwollumPlot ("Omicron",outstyle);
  //*****************************

  //***** set output style *****
  if(!io->GetOpt("OUTPUT","NOLOGO", fNoLogo)) fNoLogo=false;
  //*****************************
  

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------                  DATA                  --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //***** ffl file *****
  string fflfile;
  if(io->GetOpt("DATA","FFL", fflfile)||io->GetOpt("DATA","LCF", fflfile))
    FFL = new ffl(fflfile, GPlot->GetCurrentStyle(), fVerbosity);
  else
    FFL=NULL;
  //*****************************
  
  //***** List of channels *****
  if(!io->GetAllOpt("DATA","CHANNELS", fChannels)){
    cerr<<"Omicron::ReadOptions: a list of channels is required (DATA/CHANNELS)"<<endl;
    fChannels.push_back("M1:MISSING");
    status_OK=false;
  }
  triggers = new MakeTriggers* [(int)fChannels.size()];// output triggers
  for(int c=0; c<(int)fChannels.size(); c++)
    triggers[c] = new MakeTriggers(fChannels[c],fVerbosity);
  //*****************************
  
  //***** Sampling frequency *****
  int sampling;
  if(!io->GetOpt("DATA","SAMPLEFREQUENCY", sampling)){
    cerr<<"Omicron::ReadOptions: a working sampling frequency is required (DATA/SAMPLEFREQUENCY)"<<endl;
    sampling=16;// not to crash
    status_OK=false;
  }
  if(sampling<16){
    cerr<<"Omicron::ReadOptions: the working sampling frequency (DATA/SAMPLEFREQUENCY) should be at least 16Hz"<<endl;
    sampling=16;
    status_OK=false;
  }
  for(int c=0; c<(int)fChannels.size(); c++)
    status_OK*=triggers[c]->SetFrequencies(sampling,sampling,0.0);
  //*****************************

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------                PARAMETER               --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //***** timing *****
  vector <int> timing;
  if(!io->GetOpt("PARAMETER", "TIMING",   timing)){
    timing.push_back(64); timing.push_back(4);
  }
  else if(timing.size()==1) timing.push_back(timing[0]/4);
  else;
  //*****************************
  
  //***** Frequency range *****
  vector <double> FRange;
  if(!io->GetOpt("PARAMETER","FREQUENCYRANGE", FRange)){
    cerr<<"Omicron::ReadOptions: No search frequency range (PARAMETER/FREQUENCYRANGE)  --> set default: 2-"<<triggers[0]->GetWorkingFrequency()/2<<" Hz"<<endl;
    FRange.push_back(2); FRange.push_back(triggers[0]->GetWorkingFrequency()/2);
  }
  if(FRange.size()!=2){// must be size 2
    cerr<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) is not correct  --> set default: 2-"<<triggers[0]->GetWorkingFrequency()/2<<" Hz"<<endl;
    FRange.clear();
    FRange.push_back(2); FRange.push_back(triggers[0]->GetWorkingFrequency()/2);
  }
  if(FRange[1]>(double)triggers[0]->GetWorkingFrequency()/2.0){// check for Nyquist
    FRange.pop_back();
    FRange.push_back(triggers[0]->GetWorkingFrequency()/2);
  }
  if(FRange[0]>=FRange[1]){// check the order
    cerr<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) is not correct  --> set default: 2-"<<triggers[0]->GetWorkingFrequency()/2<<" Hz"<<endl;
    FRange.clear();
    FRange.push_back(2); FRange.push_back(triggers[0]->GetWorkingFrequency()/2);
  }
  //*****************************

  //***** Q range *****
  vector <double> QRange;
  if(!io->GetOpt("PARAMETER","QRANGE", QRange)){
    cerr<<"Omicron::ReadOptions: No search Q range (PARAMETER/QRANGE)  --> set default: 4-100"<<endl;
    QRange.push_back(4); QRange.push_back(100);
  }
  if(QRange.size()!=2){// must be size 2
    cerr<<"Omicron::ReadOptions: Q range (PARAMETER/QRANGE) is not correct  --> set default: 4-100"<<endl;
    QRange.push_back(4); QRange.push_back(100);
  }
  if(QRange[0]>=QRange[1]){// check the order
    cerr<<"Omicron::ReadOptions: Q range (PARAMETER/QRANGE) is not correct  --> set default: 4-100"<<endl;
    QRange.push_back(4); QRange.push_back(100);
  }
  //*****************************

  //***** maximum mismatch *****
  double mmm;
  if(!io->GetOpt("PARAMETER","MISMATCHMAX", mmm)){
    cerr<<"Omicron::ReadOptions: No mismatch (PARAMETER/MISMATCHMAX)  --> set default: 0.25"<<endl;
    mmm=0.25;
  }
  tile = new Otile(timing[0],QRange[0],QRange[1],FRange[0],FRange[1],triggers[0]->GetWorkingFrequency(),mmm,GPlot->GetCurrentStyle(),fVerbosity);// tiling definition
  tile->SetOverlapDuration(timing[1]);
  QRange.clear(); FRange.clear();
  for(int c=0; c<(int)fChannels.size(); c++)
    status_OK*=triggers[c]->SetHighPassFrequency(tile->GetFrequencyMin());
  //*****************************
  
  //***** Tile selection *****
  vector <double> v;
  if(!io->GetOpt("PARAMETER","SNRTHRESHOLD", v)){
    cerr<<"Omicron::ReadOptions: No SNR threshold (PARAMETER/SNRTHRESHOLD)  --> set default: 7"<<endl;
    v.push_back(7.0); v.push_back(7.0);
  }
  if(v.size()==1) v.push_back(v[0]);
  tile->SetSaveSelection(v[1],v[0]);
  //*****************************
  
  //***** set spectrum *****
  double psdlength;
  if(!io->GetOpt("PARAMETER","PSDLENGTH", psdlength)) psdlength=tile->GetTimeRange()-tile->GetOverlapDuration();
  spectrum = new Spectrum* [(int)fChannels.size()];
  for(int c=0; c<(int)fChannels.size(); c++){
    if(tile->GetFrequencyMin()>1.0) // resolution = 0.5 Hz above 1 Hz
      spectrum[c] = new Spectrum(triggers[0]->GetWorkingFrequency(),psdlength,triggers[0]->GetWorkingFrequency(),fVerbosity);
    else // increase the resolution not to extrapolate the PSD.
      spectrum[c] = new Spectrum(2*(int)ceil((double)triggers[0]->GetWorkingFrequency()/tile->GetFrequencyMin()),psdlength,triggers[0]->GetWorkingFrequency(),fVerbosity);
  }
  //*****************************
  
  //***** set clustering *****
  double cldt=0.1;
  if(!io->GetOpt("PARAMETER","CLUSTERING", fClusterAlgo)) fClusterAlgo="none";
  if(!io->GetOpt("PARAMETER","CLUSTERDT", cldt)) cldt=0.1;
  for(int c=0; c<(int)fChannels.size(); c++) triggers[c]->SetClusterizeDt(cldt);
  //*****************************
    
  //***** plot windows *****
  if(!io->GetOpt("PARAMETER","WINDOWS", fWindows))
    fWindows.push_back(tile->GetTimeRange()-tile->GetOverlapDuration());
  //*****************************

  //***** vertical scale *****
  double snrscale;
  if(!io->GetOpt("PARAMETER","SNRSCALE", snrscale)){
    cerr<<"Omicron::ReadOptions: No snr scale option (PARAMETER/SNRSCALE)  --> set default: 50"<<endl;
    snrscale=50;
  }
  tile->SetSNRScale(fabs(snrscale));
  //*****************************


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------               INJECTIONS               --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //***** injection channels *****
  FFL_inject=NULL;
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
    if(io->GetOpt("INJECTION","FFL", fflfile)||io->GetOpt("INJECTION","LCF", fflfile)){
      FFL_inject = new ffl(fflfile, GPlot->GetCurrentStyle(), fVerbosity);
    }
    else
      FFL_inject=FFL;
  }
  //*****************************

  //***** software injections *****
  inject=NULL; string injfile;
  if(io->GetOpt("INJECTION","FILENAME", injfile)){
    inject = new InjEct* [(int)fChannels.size()];
    for(int c=0; c<(int)fChannels.size(); c++) inject[c] = new InjEct(triggers[c],injfile,fVerbosity);
  }
  //*****************************

  // dump options
  if(fVerbosity>1){
    cout<<"**********************************************"<<endl;
    cout<<"**********      PARAMETER FILE      **********"<<endl;
    cout<<"**********************************************"<<endl;
    io->Dump(cout);
    cout<<"**********************************************"<<endl;
  }
  delete io;
 
  return;
}
