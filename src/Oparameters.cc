/**
 * @file 
 * @brief Manage Omicron options.
 * @author Florent Robinet - <a href="mailto:florent.robinet@ijclab.in2p3.fr">florent.robinet@ijclab.in2p3.fr</a>
 */
#include "Oomicron.h"


/**
 * @brief Parse Omicron parameters.
 * @details Omicron parameters must be provided with an option file (see Omicron::Omicron()).
 * Each parameter is identified with a tag and a keyword:
 * @verbatim
TAG  KEYWORD  [PARAMETERS]
@endverbatim
 * The combination of tag/keyword/parameters is called an option.
 * @note For some options, multiple parameters can be used. They are separated by white spaces or tabs. 
 * @warning If an option is not provided, default values are used for the parameter.
 * 
 * Here we list all the options for Omicron.
 *
 * @section omicron_readoptions_output OUTPUT
 * 
 * @subsection omicron_readoptions_output_directory output directory 
 * @verbatim
OUTPUT  DIRECTORY  [PARAMETER]
@endverbatim
 * `[PARAMETER]` is the directory path (relative or absolute).
 * The current directory is used by default if this option is not provided or if the specified directory does not exist.
 *
 * @subsection omicron_readoptions_output_products output products
 * @verbatim
OUTPUT  PRODUCTS  [PARAMETERS]
@endverbatim
 * `[PARAMETERS]` is the list of products computed by Omicron.
 * Possible parameters are:
 * - `triggers`: Omicron tiles above a SNR threshold are saved in a file.
 * - `asd`: The amplitude spectral density function is saved in a file.
 * - `psd`: The power spectral density function is saved in a file.
 * - `html`: A html report is produced in the output directory specified with the @ref omicron_readoptions_output_directory "directory option".
 * - `timeseries`: The condition time series is saved in a file.
 * - `white`: The time series after whitening is saved in a file.
 * - `whitepsd`: The power spectral density after whitening is saved in a file.
 * - `mapsnr`: The snr spectrograms are saved in a file.
 * - `mapamplitude`: The amplitude spectrograms are saved in a file.
 * - `mapphase`: The phase spectrograms are saved in a file.
 *
 * By default, only `triggers` are produced.
 * The ouput file format is specified with the @ref omicron_readoptions_output_format "format option".
 *
 * @subsection omicron_readoptions_output_verbosity verbosity level
 * @verbatim
OUTPUT  VERBOSITY  [PARAMETER]
@endverbatim
 * `[PARAMETER]` is the verbosity level: 0, 1, 2, or 3.
 * By default = 0.
 *
 * @subsection omicron_readoptions_output_format file format
 * @verbatim
OUTPUT  FORMAT  [PARAMETERS]
@endverbatim
 * `[PARAMETERS]` is the list of file formats to save @ref omicron_readoptions_output_products "output products".
 * Supported formats: `root` (native), `hdf5` (for triggers only), and all usual graphical formats (`svg`, `gif`, `pdf`, `png`, `eps`...).
 * By default = `root`.
 *
 */
void Omicron::ReadOptions(const int aGpsRef, const bool aStrict){

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
    if(aStrict) status_OK=false;
    fMaindir=".";
  }
  if(!IsDirectory(fMaindir)){
    cerr<<"Omicron::ReadOptions: output directory "<<fMaindir<<" does not exist --> set default: current"<<endl;
    fMaindir=".";
  }
  //*****************************
   
  //***** verbosity *****
  if(!io->GetOpt("OUTPUT","VERBOSITY", fVerbosity)){
    cerr<<"Omicron::ReadOptions: No verbosity level (OUTPUT/VERBOSITY)  --> set default: 0"<<endl;
    fVerbosity=0;
  }
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
  //*****************************

  //***** set plot dimensions *****
  vector <int> dims;
  io->GetOpt("OUTPUT","PLOTDIMENSIONS", dims);
  //*****************************

  //***** set output style *****
  if(!io->GetOpt("OUTPUT","NOLOGO", fNoLogo)) fNoLogo=false;
  //*****************************
    

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------                  DATA                  --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //***** ffl file *****
  string fflfileopt;
  vector <string> fflfile;
  if(io->GetOpt("DATA","FFL", fflfileopt)||io->GetOpt("DATA","LCF", fflfileopt)){
    fflfile = SplitString(fflfileopt, ' ');
    FFL = new ffl(fflfile[0], outstyle, fVerbosity);
    FFL->SetName("mainffl");
    status_OK*=FFL->DefineTmpDir(fMaindir);

    if(fflfile.size()==1) status_OK*=FFL->LoadFrameFile(aGpsRef);
    else if(fflfile.size()==2) status_OK*=FFL->LoadFrameFile(aGpsRef, atoi(fflfile[1].c_str()));
    else status_OK*=FFL->LoadFrameFile(aGpsRef, atoi(fflfile[1].c_str()), atoi(fflfile[2].c_str()));
  }
  else
    FFL=NULL;
  //*****************************

  
  //***** Trigger buffer *****
  int bufsize;
  if(!io->GetOpt("PARAMETER","TRIGGERBUFFERSIZE", bufsize)){
    cerr<<"Omicron::ReadOptions: No trigger buffer (PARAMETER/TRIGGERBUFFERSIZE)"<<endl;
    bufsize=0;
  }
  
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

  //***** List of channels/streams  *****
  vector <string> channels; nchannels=0;
  vector <string> bl_channels;
  vector <string> channels_tmp;
  vector <string> channels_tmp2;
  int ndummy;
  bool bl;
  io->GetAllOpt("DATA","CHANNELS", channels_tmp);
  io->GetAllOpt("DATA","BLACKLISTEDCHANNELS", bl_channels);
  
  // select channels at GPS ref
  if(FFL!=NULL){
    for(int c1=0; c1<(int)channels_tmp.size(); c1++){// apply filters
      channels_tmp2=FFL->GetChannelList(channels_tmp[c1]);
      if(channels_tmp2.size()==0){
	cerr<<"Omicron::ReadOptions: No channels matching "<<channels_tmp[c1]<<" (@"<<aGpsRef<<") --> ignore"<<endl;
	if(aStrict) status_OK=false;
      }

      for(int c2=0; c2<(int)channels_tmp2.size(); c2++){

	// check channel validity
	if(FFL->GetChannelSampling(channels_tmp2[c2],ndummy)<16){// must be >=16 Hz
	  cerr<<"Omicron::ReadOptions: "<<channels_tmp2[c2]<<" sampling frequency is too low --> ignore"<<endl;
	  if(aStrict) status_OK=false;
	  continue;
	}
	if(FFL->GetChannelSampling(channels_tmp2[c2],ndummy)<sampling){// check working sampling frequency
	  cerr<<"Omicron::ReadOptions: "<<channels_tmp2[c2]<<" sampling frequency is below the working frequency --> ignore"<<endl;
	  if(aStrict) status_OK=false;
	  continue;
	}
	bl=false;
	for(int b=0; b<(int)bl_channels.size(); b++){
	  if(!fnmatch(bl_channels[b].c_str(), channels_tmp2[c2].c_str(),FNM_NOESCAPE)){ // black-listed
	    cerr<<"Omicron::ReadOptions: "<<channels_tmp2[c2]<<" is blacklisted --> ignore"<<endl;
	    bl=true;
	    break;
	  }
	}
	if(bl) continue;
	channels.push_back(channels_tmp2[c2]);
      }
    }
  }
  else{
    for(int c1=0; c1<(int)channels_tmp.size(); c1++){
      bl=false;
      for(int b=0; b<(int)bl_channels.size(); b++){
	if(!fnmatch(bl_channels[b].c_str(), channels_tmp[c1].c_str(),FNM_NOESCAPE)){ // black-listed
	  cerr<<"Omicron::ReadOptions: "<<channels_tmp[c1]<<" is blacklisted --> ignore"<<endl;
	  bl=true;
	  break;
	}
      }
      if(bl) continue;
      channels.push_back(channels_tmp[c1]);
    }
  }
  channels_tmp.clear();
  channels_tmp2.clear();

  if(channels.size()==0){
    cerr<<"Omicron::ReadOptions: the list of channels is empty (DATA/CHANNELS)"<<endl;
    channels.push_back("M1:MISSING");
    status_OK=false;
  }
  
  nchannels = (int)channels.size();
  triggers = new TriggerBuffer* [nchannels];
  for(int c=0; c<nchannels; c++){
    triggers[c] = new TriggerBuffer(bufsize,channels[c],fVerbosity);
    triggers[c]->SetDCRemoval(true);
    status_OK*=triggers[c]->SetFrequencies(sampling,sampling,0.0);
  }
  channels.clear();
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
  if(FRange.size()==1){// single value
    FRange.push_back(FRange[0]);
  }
  if(FRange.size()!=2){// must be size 2
    cerr<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) is not correct  --> set default: 2-"<<triggers[0]->GetWorkingFrequency()/2<<" Hz"<<endl;
    FRange.clear();
    if(aStrict) status_OK=false;
    FRange.push_back(2); FRange.push_back(triggers[0]->GetWorkingFrequency()/2);
  }
  if(FRange[1]>(double)triggers[0]->GetWorkingFrequency()/2.0){// check for Nyquist
    if(aStrict) status_OK=false;
    FRange.pop_back();
    FRange.push_back(triggers[0]->GetWorkingFrequency()/2);
  }
  if(FRange[0]>=FRange[1]){// check the order
    cerr<<"Omicron::ReadOptions: Frequency range (PARAMETER/FREQUENCYRANGE) is not correct  --> set default: 2-"<<triggers[0]->GetWorkingFrequency()/2<<" Hz"<<endl;
    if(aStrict) status_OK=false;
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
  if(QRange.size()==1){// single value
    QRange.push_back(QRange[0]);
  }
  if(QRange.size()!=2){// must be size 2
    cerr<<"Omicron::ReadOptions: Q range (PARAMETER/QRANGE) is not correct  --> set default: 4-100"<<endl;
    if(aStrict) status_OK=false;
    QRange.clear(); QRange.push_back(4); QRange.push_back(100);
  }
  if(QRange[0]>=QRange[1]){// check the order
    cerr<<"Omicron::ReadOptions: Q range (PARAMETER/QRANGE) is not correct  --> set default: 4-100"<<endl;
    if(aStrict) status_OK=false;
    QRange.clear(); QRange.push_back(4); QRange.push_back(100);
  }
  //*****************************

  //***** maximum mismatch *****
  double mmm;
  if(!io->GetOpt("PARAMETER","MISMATCHMAX", mmm)){
    cerr<<"Omicron::ReadOptions: No mismatch (PARAMETER/MISMATCHMAX)  --> set default: 0.25"<<endl;
    mmm=0.25;
  }
  tile = new Otile(timing[0],QRange[0],QRange[1],FRange[0],FRange[1],triggers[0]->GetWorkingFrequency(),mmm,outstyle,fVerbosity);// tiling definition
  if(dims.size()==2){
    tile->ResizePlot(dims[0],dims[1]);
  }
  tile->SetOverlapDuration(timing[1]);
  if(fOutProducts.find("mapsnr")!=string::npos) tile->SetMapFill("snr");
  else if(fOutProducts.find("mapamplitude")!=string::npos) tile->SetMapFill("amplitude");
  else if(fOutProducts.find("mapphase")!=string::npos) tile->SetMapFill("phase");
  else tile->SetMapFill("snr");
  QRange.clear(); FRange.clear();
  //*****************************
  
  //***** Tile selection *****
  vector <double> v;
  if(!io->GetOpt("PARAMETER","SNRTHRESHOLD", v)){
    cerr<<"Omicron::ReadOptions: No SNR threshold (PARAMETER/SNRTHRESHOLD)  --> set default: 7"<<endl;
    v.push_back(7.0); v.push_back(7.0);
  }
  if(v.size()==1) v.push_back(v[0]);
  tile->SetSNRThr(v[1],v[0]);
  //*****************************
  
  //***** set spectrum *****
  double psdlength;
  if(!io->GetOpt("PARAMETER","PSDLENGTH", psdlength)) psdlength=tile->GetTimeRange()-tile->GetOverlapDuration();
  spectrum1 = new Spectrum* [nchannels];
  spectrum2 = new Spectrum* [nchannels];
  if(tile->GetFrequencyMin()>1.0){ // resolution = 0.5 Hz above 1 Hz
    for(int c=0; c<nchannels; c++){
      spectrum1[c] = new Spectrum(triggers[0]->GetWorkingFrequency(),psdlength,triggers[0]->GetWorkingFrequency(),fVerbosity);
      spectrum2[c] = new Spectrum(triggers[0]->GetWorkingFrequency(),psdlength,triggers[0]->GetWorkingFrequency(),fVerbosity);
    }
    spectrumw = new Spectrum(triggers[0]->GetWorkingFrequency(),tile->GetTimeRange()-tile->GetOverlapDuration(),triggers[0]->GetWorkingFrequency(),0);
  }
  else{ // increase the resolution not to extrapolate the PSD.
    for(int c=0; c<nchannels; c++){
      spectrum1[c] = new Spectrum(2*NextPowerOfTwo((double)triggers[0]->GetWorkingFrequency()/tile->GetFrequencyMin()),psdlength,triggers[0]->GetWorkingFrequency(),fVerbosity);
      spectrum2[c] = new Spectrum(2*NextPowerOfTwo((double)triggers[0]->GetWorkingFrequency()/tile->GetFrequencyMin()),psdlength,triggers[0]->GetWorkingFrequency(),fVerbosity);
    }
    spectrumw = new Spectrum(2*NextPowerOfTwo((double)triggers[0]->GetWorkingFrequency()/tile->GetFrequencyMin()),tile->GetTimeRange()-tile->GetOverlapDuration(),triggers[0]->GetWorkingFrequency(),0);
  }
  //*****************************

  //***** set highpass filter *****
  double hplf;
  if(io->GetOpt("PARAMETER","HIGHPASS", hplf)){
    for(int c=0; c<nchannels; c++)
      status_OK*=triggers[c]->SetHighPassFrequency(hplf);
  }
  //*****************************
  
  //***** set clustering *****
  double cldt=0.1;
  if(!io->GetOpt("PARAMETER","CLUSTERING", fClusterAlgo)) fClusterAlgo="none";
  if(!io->GetOpt("PARAMETER","CLUSTERDT", cldt)) cldt=0.1;
  for(int c=0; c<(int)nchannels; c++) triggers[c]->SetClusterizeDt(cldt);
  //*****************************
    
  //***** plot windows *****
  if(!io->GetOpt("PARAMETER","WINDOWS", fWindows))
    fWindows.push_back(tile->GetTimeRange()-tile->GetOverlapDuration());
  //*****************************

  //***** map scale *****
  int mapvscale;
  if(!io->GetOpt("PARAMETER","MAPLOGSCALE", mapvscale)){
    cerr<<"Omicron::ReadOptions: No map vertical scale option (PARAMETER/MAPLOGSCALE)  --> set default: 1"<<endl;
    mapvscale=1;
  }
  tile->SetLogz(mapvscale);
  //*****************************

  //***** map vertical scale *****
  vector <double> mapvrange;
  if(!io->GetOpt("PARAMETER","MAPVRANGE", mapvrange)){
    if(!tile->GetMapFill().compare("snr")){
      cerr<<"Omicron::ReadOptions: No map vertical range option (PARAMETER/MAPVRANGE)  --> set default: 1-30"<<endl;
      mapvrange.push_back(1.0); mapvrange.push_back(30);
    }
    else if(!tile->GetMapFill().compare("amplitude")){
      cerr<<"Omicron::ReadOptions: No map vertical range option (PARAMETER/MAPVRANGE)  --> set default: automatic"<<endl;
      mapvrange.push_back(-1.0); mapvrange.push_back(-1.0);
    }
    else{
      cerr<<"Omicron::ReadOptions: No map vertical range option (PARAMETER/MAPVRANGE)  --> set default: -pi-pi"<<endl;
      mapvrange.push_back(-TMath::Pi()); mapvrange.push_back(-TMath::Pi());
    }
  }
  if(mapvrange.size()<2){
    if(mapvrange[0]>0) mapvrange.insert(mapvrange.begin(),mapvrange[0]/10);
    else mapvrange.push_back(mapvrange[0]);
  }
  tile->SetRangez(mapvrange[0],mapvrange[1]);
  //*****************************

  //***** fft plans *****
  if(!io->GetOpt("PARAMETER","FFTPLAN", fftplan)){
    cerr<<"Omicron::ReadOptions: No fftplan option (PARAMETER/FFTPLAN)  --> set default: FFTW_MEASURE"<<endl;
    fftplan="FFTW_MEASURE";
  }
  //*****************************

  //***** trigger max *****
  if(!io->GetOpt("PARAMETER","TRIGGERRATEMAX", fratemax)){
    cerr<<"Omicron::ReadOptions: No trigger rate limit option (PARAMETER/TRIGGERRATEMAX)  --> set default: 5000 Hz"<<endl;
    fratemax=5000.0;
  }
  //*****************************

  //***** set chirp mass *****
  vector <double> ch;
  if(!io->GetOpt("PARAMETER","CHIRP", ch)) tile->SetChirp(-1.0,-1.0);
  else{
    if(ch.size()==2) tile->SetChirp(ch[0],ch[1]);
    else if(ch.size()==1) tile->SetChirp(ch[0],-1.0);
    else tile->SetChirp(-1.0,-1.0);
  }

  //*****************************

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //--------------               INJECTIONS               --------------
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //***** injection channels *****
  fflfile.clear();
  FFL_inject=NULL;
  if(io->GetOpt("INJECTION","CHANNELS", fInjChan)){
    if((int)fInjChan.size()!=nchannels){
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
      for(int i=0; i<nchannels; i++) fInjFact.push_back(1.0);
    }
    if(io->GetOpt("INJECTION","FFL", fflfileopt)||io->GetOpt("INJECTION","LCF", fflfileopt)){
      fflfile = SplitString(fflfileopt, ' ');
      FFL_inject = new ffl(fflfile[0], outstyle, fVerbosity);
      FFL_inject->SetName("injffl");
      status_OK*=FFL_inject->DefineTmpDir(fMaindir);
      if(fflfile.size()==1) status_OK*=FFL_inject->LoadFrameFile(aGpsRef);
      else if(fflfile.size()==2) status_OK*=FFL_inject->LoadFrameFile(aGpsRef, atoi(fflfile[1].c_str()));
      else status_OK*=FFL_inject->LoadFrameFile(aGpsRef, atoi(fflfile[1].c_str()), atoi(fflfile[2].c_str()));
    }
    else
      FFL_inject=FFL;
  }
  //*****************************
  
  //***** software injections *****
  inject=NULL; string injfile;
  if(io->GetOpt("INJECTION","FILENAME", injfile)){
    inject = new InjEct* [nchannels];
    for(int c=0; c<nchannels; c++) inject[c] = new InjEct(triggers[c],injfile,fVerbosity);
  }
  //*****************************

  //***** sg injections *****
  oinj = new Oinject(tile->GetTimeRange());
  if(!io->GetOpt("INJECTION","SG", fsginj)) fsginj=0;

  vector <double> param;
  if(io->GetOpt("INJECTION","SGTIME", param)){
    if(param.size()==1) { param.push_back(param[0]); oinj->SetTimeRange(param[0],param[1]); }
    else if(param.size()==2) oinj->SetTimeRange(param[0],param[1]);
    else
      	cerr<<"Omicron::ReadOptions: INJECTION/SGTIME is incorrect"<<endl;
  }
  param.clear();
  if(io->GetOpt("INJECTION","SGFREQUENCY", param)){
    if(param.size()==1) { param.push_back(param[0]); oinj->SetFrequencyRange(param[0],param[1]); }
    else if(param.size()==2) oinj->SetFrequencyRange(param[0],param[1]);
    else
      	cerr<<"Omicron::ReadOptions: INJECTION/SGFREQUENCY is incorrect"<<endl;
  }
  param.clear();
  if(io->GetOpt("INJECTION","SGQ", param)){
    if(param.size()==1) { param.push_back(param[0]); oinj->SetQRange(param[0],param[1]); }
    else if(param.size()==2) oinj->SetQRange(param[0],param[1]);
    else
      	cerr<<"Omicron::ReadOptions: INJECTION/SGQ is incorrect"<<endl;
  }
  param.clear();
  if(io->GetOpt("INJECTION","SGAMPLITUDE", param)){
    if(param.size()==1) { param.push_back(param[0]); oinj->SetAmplitudeRange(param[0],param[1]); }
    else if(param.size()==2) oinj->SetAmplitudeRange(param[0],param[1]);
    else
      	cerr<<"Omicron::ReadOptions: INJECTION/SGAMPLITUDE is incorrect"<<endl;
  }
  param.clear();
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
