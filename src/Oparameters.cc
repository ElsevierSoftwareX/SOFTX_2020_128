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
 * Here is an example of an option file:
 * @include parameters.txt
 * 
 * Here we list all the options for Omicron.
 *
 * @section omicron_readoptions_output OUTPUT
 * 
 * @subsection omicron_readoptions_output_directory Output directory 
 * @verbatim
OUTPUT  DIRECTORY  [PARAMETER]
@endverbatim
 * `[PARAMETER]` is the directory path (relative or absolute).
 * The current directory is used by default if this option is not provided or if the specified directory does not exist.
 *
 * @subsection omicron_readoptions_output_products Output products
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
 * @subsection omicron_readoptions_output_format File format
 * @verbatim
OUTPUT  FORMAT  [PARAMETERS]
@endverbatim
 * `[PARAMETERS]` is the list of file formats to save @ref omicron_readoptions_output_products "output products".
 * Supported formats: `root` (native), `hdf5` (for triggers only), and all the usual graphical formats (`svg`, `gif`, `pdf`, `png`, `eps` and so on).
 * By default = `root`.
 *
 * @subsection omicron_readoptions_output_verbosity Verbosity level
 * @verbatim
OUTPUT  VERBOSITY  [PARAMETER]
@endverbatim
 * `[PARAMETER]` is the verbosity level: 0, 1, 2, or 3.
 * By default = 0.
 *
 * @subsection omicron_readoptions_output_style Output style
 * @verbatim
OUTPUT  STYLE  [PARAMETER]
@endverbatim
 * `[PARAMETER]` defines the css style for the web report (if requested with the @ref omicron_readoptions_output_products "output products" option). It also defines the graphical color palette for the plots. Several styles are supported:
 * - `GWOLLUM` (default): black background, blue style.
 * - `FIRE`: black background, color palette from red to yellow
 * - `STANDARD`: light background, rainbow color palette
 * - `PINK`: black background, pink style.
 *
 * @subsection omicron_readoptions_output_nologo No logo flag
 * @verbatim
OUTPUT  NOLOGO  [PARAMETER]
@endverbatim
 * If `[PARAMETER]` is set to a non-zero value, the omicron background logo does not appear on the web reports (if requested with the @ref omicron_readoptions_output_products "output products" option). By default, the logo is present (`[PARAMETER] = 0`).
 *
 * @subsection omicron_readoptions_output_plotdimensions Plot dimensions
 * @verbatim
OUTPUT  PLOTDIMENSIONS  [PARAMETERS]
@endverbatim
 * This option sets the graphical plot dimensions (if requested with the @ref omicron_readoptions_output_format "output format" option). Exactly two integer numbers should be provided with `[PARAMETERS]`: the width and the height measured in a number of pixels.
 *
 * @section omicron_readoptions_data DATA
 *
 * @subsection omicron_readoptions_data_ffl Frame file list
 * @verbatim
DATA  FFL  [PARAMETERS]
@endverbatim
 * or
 * @verbatim
DATA  LCF  [PARAMETERS]
@endverbatim
 * This option specifies the list of frame files to be processed by omicron. In the first option, it is provided as a frame file list (FFL). In the second option, it is provided as a LAL cache file (LCF). 
 * `[PARAMETERS]` must provide the path to the FFL or LCF file (absolute or relative). It can be complemented by 2 optional parameters which are used to load the FFL file multiple times in a row (for online applications):
 * - Number of retries to load the FFL,
 * - Time in [s] between two retries.
 *
 * By default, no FFL file is used.
 *
 * @subsection omicron_readoptions_data_samplefrequency Working sampling frequency
 * @verbatim
DATA  SAMPLEFREQUENCY  [PARAMETER]
@endverbatim
 * This option specifies the working sampling frequency of omicron in [Hz]. All channel time series are sampled to a common frequency before being processed. Only downsampling is supported. The working sampling frequency must be a power of 2 and must be at least 16 Hz.
 *
 * @subsection omicron_readoptions_data_channels List of channels
 * @verbatim
DATA  CHANNELS  [PARAMETERS]
@endverbatim
 * This option specifies the list of channels to process, e.g. `V1:Hrec_hoft_16384Hz V1:LSC_DARM`. Multiple channel names can be provided using wildcards. This option can also be used on multiple lines. This option is mandatory.
 *
 * @subsection omicron_readoptions_data_blacklistedchannels List of black-listed channels
 * @verbatim
DATA  BLACKLISTEDCHANNELS  [PARAMETERS]
@endverbatim
 * This option specifies the list of channels to exclude from the processing, e.g. `V1:Hrec_hoft_16384Hz V1:LSC_DARM`. Multiple channel names can be provided using wildcards. This option can also be used on multiple lines.
 *
 * @subsection omicron_readoptions_data_triggerbuffersize Trigger buffer size
 * @verbatim
DATA  TRIGGERBUFFERSIZE  [PARAMETER]
@endverbatim
 * This option is used for online applications. Triggers produced by omicron can be buffered. This parameter defines the size of the buffer.
 * By default, there is no buffer (size=0).
 *
 * @section omicron_readoptions_parameter PARAMETER
 *
 * @subsection omicron_readoptions_parameter_timing Analysis timing
 * @verbatim
PARAMETER  TIMING  [PARAMETERS]
@endverbatim
 * This option specifies the timing parameters to perform the analysis. It includes 2 parameters:
 * - The analysis window size in [s]. It must be a power of 2.
 * - Overlap duration [s]. Two consecutive analysis windows overlap by this amount.
 * By default, window size = 64 s and overlap = 4 s.
 *
 * @subsection omicron_readoptions_parameter_frequencyrange Frequency range
 * @verbatim
PARAMETER  FREQUENCYRANGE  [PARAMETERS]
@endverbatim
 * This option specifies the search frequency range with a lower bound and a higher bound, both in [Hz].
 * The higher bound must be compatible with the @ref omicron_readoptions_data_samplefrequency "sampling frequency" option.
 * By default, lower frequency = 2 Hz and higher frequency = sampling frequency /2.
 *
 * @subsection omicron_readoptions_parameter_qrange Q range
 * @verbatim
PARAMETER  QRANGE  [PARAMETERS]
@endverbatim
 * This option specifies the search Q range with a lower bound and a higher bound. The lower bound must be at least \f$\sqrt{11}\f$.
 * By default, lower Q = 4 and higher Q = 100.
 *
 * @subsection omicron_readoptions_parameter_mismatchmax Maximum mismatch between tiles
 * @verbatim
PARAMETER  MISMATCHMAX  [PARAMETER]
@endverbatim
 * This option specifies the maximum energy mismatch between time-frequency tiles. This value should be a number between 0 and 1.
 * By default, max mismatch = 0.25.
 *
 * @subsection omicron_readoptions_parameter_snrthreshold SNR thresholds
 * @verbatim
PARAMETER  SNRTHRESHOLD  [PARAMETERS]
@endverbatim
 * This option specifies the signal-to-noise ratio. There are 2 SNR threshold:
 * - The first SNR threshold applies to the individual time-frequency tiles to define a trigger.
 * - The second SNR threshold determines if spectogram plots must be generated. The plots are produced if at least one tile in the first plot time window is above this threshold. 
 * 
 * By default, both thresholds are set to 7.
 * @note only one `PARAMETER` can be provided. The same vlaue will then be given to both thresholds.
 *
 * @subsection omicron_readoptions_parameter_psdlength Power spectrum density length
 * @verbatim
PARAMETER  PSDLENGTH  [PARAMETER]
@endverbatim
 * This option specifies the duration [s] over which the power spectrum density is estimated.
 * By default, use the @ref omicron_readoptions_parameter_timing "analysis window duration" minus the overlap.
 *
 * @subsection omicron_readoptions_parameter_highpass High-pass filter
 * @verbatim
PARAMETER  HIGHPASS  [PARAMETER]
@endverbatim
 * This option specifies the frequency cutoff [Hz] at which to high-pass the data before processing.
 * By default, do not high-pass the data.
 *
 * @subsection omicron_readoptions_parameter_clustering Clustering algorithm
 * @verbatim
PARAMETER  CLUSTERING  [PARAMETER]
@endverbatim
 * This option selects the clustering method. Only one clustering method is currently supported: "TIME".
 * By default, no clustering.
 * @note This option is useless for triggers saved in a ROOT format. Clusters are nevered saved in ROOT files.
 *
 * @subsection omicron_readoptions_parameter_clusterdt Clustering time window
 * @verbatim
PARAMETER  CLUSTERDT  [PARAMETER]
@endverbatim
 * This option specifies the time window duration to cluster triggers in [s].
 * By default, = 0.1 s
 *
 * @subsection omicron_readoptions_parameter_windows Time windows for plots
 * @verbatim
PARAMETER  WINDOWS  [PARAMETERS]
@endverbatim
 * This option specifies the time window durations for plots. There is no limit on the number of parameters.
 * By default, on single value = @ref omicron_readoptions_parameter_timing "analysis window duration" minus the overlap.
 *
 * @subsection omicron_readoptions_parameter_maplogscale Log scale for SNR scales in plots
 * @verbatim
PARAMETER  MAPLOGSCALE  [PARAMETER]
@endverbatim
 * This option, when different from 0, activates the log scale when plotting the SNR.
 * By default, the log scale is used.
 *
 * @subsection omicron_readoptions_parameter_mapvrange Vertical range for spectrograms
 * @verbatim
PARAMETER  MAPVRANGE  [PARAMETERS]
@endverbatim
 * This option specifies the vertical range for spectrogram plots.
 *
 * @subsection omicron_readoptions_parameter_fftplan FFT plan
 * @verbatim
PARAMETER  FFTPLAN  [PARAMETER]
@endverbatim
 * This option specifies the plan to perform Fourier transforms with FFTW.
 * By default = "FFTW_MEASURE".
 *
 * @subsection omicron_readoptions_parameter_triggerratemax Maximum trigger rate
 * @verbatim
PARAMETER  TRIGGERRATEMAX [PARAMETER]
@endverbatim
 * This option specifies maximum trigger rate limit [Hz] when saving triggers to files. If the trigger rate is above this value (over the @ref omicron_readoptions_parameter_timing "analysis window"), the file is not saved.
 * By default = 5000 Hz.
 *
 * @subsection omicron_readoptions_parameter_chirp Newtonian chirp
 * @verbatim
PARAMETER  CHIRP [PARAMETERS]
@endverbatim
 * With this option, it is possible to draw a Newtonian chirp on top of omicron spectrograms. There can be one or two parameters:
 * - A chirp mass in solar masses: the chirp is drawn with an end time positioned at the center of the @ref omicron_readoptions_parameter_timing "analysis window".
 * - A chirp mass in solar masses + a GPS time [s]: the chirp is drawn with an end time positioned at the requested GPS time.
 *
 * By default: no chirp.
 *
 * @subsection omicron_readoptions_injection INJECTION
 *
 * @subsection omicron_readoptions_injection_channels Injection channel list 
 * @verbatim
INJECTION  CHANNELS [PARAMETERS]
@endverbatim
 * This option specifies the list of channels to be used as injection signals. There should be as many channels as listed in the @ref omicron_readoptions_data_channels "main channel list".
 *
 * @subsection omicron_readoptions_injection_factors Injection scaling factors
 * @verbatim
INJECTION  FACTORS [PARAMETERS]
@endverbatim
 * This option specifies the list of amplitude scaling factors used to inject signals. There should be as many factors as the number of @ref omicron_readoptions_injection_channels "injection channels".
 *
 * @subsection omicron_readoptions_injection_ffl Frame file list 
 * @verbatim
INJECTION  FFL [PARAMETERS]
@endverbatim
 * or
 * @verbatim
INJECTION  LCF [PARAMETERS]
@endverbatim 
 * If the @ref omicron_readoptions_injection_channels "injection channel(s)" are not included in the @ref omicron_readoptions_data_ffl "main ffl file", this option specifies an additional ffl file. `[PARAMETER]` is the relative or absolute path to the ffl/lcf file. This parameter can be complemented by 2 optional parameters which are used to load the FFL file multiple times in a row (for online applications):
 * - Number of retries to load the FFL,
 * - Time in [s] between two retries.
 *
 * @subsection omicron_readoptions_injection_filename Injection file
 * @verbatim
INJECTION  FILENAME [PARAMETER]
@endverbatim
 * Injections can also be performed through an injection file listing the source/waveform parameters. The injection file must be a ROOT file generated with the InjGen class. This option provides the file absolute/relative path to the injection file.
 *
 * @subsection omicron_readoptions_injection_sg Sine-Gauss injections
 * @verbatim
INJECTION  SG [PARAMETER]
@endverbatim
 * Sine-gauss injections can be performed by setting `[PARAMETER]` to a value different from 0. One waveform is injected in every @ref omicron_readoptions_parameter_timing "analysis window".
 *
 * @subsection omicron_readoptions_injection_sgtime Sine-Gauss injection time
 * @verbatim
INJECTION  SGTIME [PARAMETERS]
@endverbatim
* By default, sine-Gaussian waveforms are always injected at the center of the @ref omicron_readoptions_parameter_timing "analysis window". With this option, the time of the injection can be taken as a random value in a given time range. The time range is defined with 2 values: a negative and a positive range from the center (in seconds).
 *
 * @subsection omicron_readoptions_injection_sgfrequency Sine-Gauss injection frequency
 * @verbatim
INJECTION  SGFREQUENCY [PARAMETERS]
@endverbatim
* With this option, the frequency of the injection is taken as a random value in a given frequency range, following a logarithmic distribution. If only one value is provided, the injection frequency is fixed at that value.
 *
 * @subsection omicron_readoptions_injection_sgq Sine-Gauss injection Q
 * @verbatim
INJECTION  SGQ [PARAMETERS]
@endverbatim
* With this option, the quality factor of the injection is taken as a random value in a given range, following a logarithmic distribution. If only one value is provided, the injection Q is fixed at that value.
 *
 * @subsection omicron_readoptions_injection_sgamplitude Sine-Gauss amplitude
 * @verbatim
INJECTION  SGAMPLITUDE [PARAMETERS]
@endverbatim
* With this option, the amplitude of the injection is taken as a random value in a given range, following a logarithmic distribution. If only one value is provided, the injection amplitude is fixed at that value.
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
    cerr<<"Omicron::ReadOptions: No timing (PARAMETER/TIMING)  --> set default: 64s with 4s overlap"<<endl;
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
  if(!io->GetOpt("PARAMETER","PSDLENGTH", psdlength)){
    cerr<<"Omicron::ReadOptions: No PSD length (PARAMETER/PSDLENGTH)  --> set default: "<<tile->GetTimeRange()-tile->GetOverlapDuration()<<endl;
    psdlength=tile->GetTimeRange()-tile->GetOverlapDuration();
  }
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
