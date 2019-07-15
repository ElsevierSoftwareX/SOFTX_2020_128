//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "ReadTriggers.h"
#include "OmicronUtils.h"


using namespace std;

void PrintUsage(void){
  cerr<<endl;
  cerr<<"Usage:"<<endl;
  cerr<<endl;
  cerr<<"omicron-print.exe channel=[channel name] \\"<<endl;
  cerr<<"                  file=[trigger file pattern] \\"<<endl;
  cerr<<"                  gps-start=[GPS start] \\"<<endl;
  cerr<<"                  gps-end=[GPS end] \\"<<endl;
  cerr<<"                  snr-min=[minimum SNR] \\"<<endl;
  cerr<<"                  snr-max=[maximum SNR] \\"<<endl;
  cerr<<"                  freq-min=[minimum frequency] \\"<<endl;
  cerr<<"                  freq-max=[maximum frequency] \\"<<endl;
  cerr<<"                  print=[output type] \\"<<endl;
  cerr<<"                  cluster-dt=[cluster time window] \\"<<endl;
  cerr<<"                  print-time=[1/0] \\"<<endl;
  cerr<<"                  print-freq=[1/0] \\"<<endl;
  cerr<<"                  print-snr=[1/0] \\"<<endl;
  cerr<<"                  print-q=[1/0] \\"<<endl;
  cerr<<"                  print-amplitude=[1/0] \\"<<endl;
  cerr<<"                  print-phase=[1/0] \\"<<endl;
  cerr<<"                  print-time-start=[1/0] \\"<<endl;
  cerr<<"                  print-time-end=[1/0] \\"<<endl;
  cerr<<"                  print-freq-start=[1/0] \\"<<endl;
  cerr<<"                  print-freq-end=[1/0] \\"<<endl;
  cerr<<"                  print-duration=[1/0] \\"<<endl;
  cerr<<"                  print-bandwidth=[1/0]"<<endl;
 cerr<<endl;
  cerr<<"[channel name]            channel name used to retrieve centralized Omicron triggers"<<endl;
  cerr<<"[trigger file pattern]    file pattern to ROOT trigger files (GWOLLUM convention)"<<endl;
  cerr<<"[GPS start]               starting GPS time (integer only)"<<endl;
  cerr<<"[GPS end]                 stopping GPS time (integer only)"<<endl;
  cerr<<"[minimum SNR]             minimum SNR value"<<endl;
  cerr<<"[maximum SNR]             maximum SNR value"<<endl;
  cerr<<"[minimum frequency]       minimum frequency value [Hz]"<<endl;
  cerr<<"[maximum frequency]       maximum frequency value [Hz]"<<endl;
  cerr<<"[output type]             \"triggers\", \"clusters\" or \"segments\". By default, print=\"clusters\""<<endl;
  cerr<<"[cluster time window]     cluster time window [s]. By default, cluster-dt=0.1"<<endl;
  cerr<<endl;
  cerr<<"For \"print-\" options, 1=yes and 0=no."<<endl;
  cerr<<endl;
  return;
}

int main (int argc, char* argv[]){

  // number of arguments
  if(argc<2){
    PrintUsage();
    return 1;
  }

  // list of parameters + default
  string chname="";     // channel name
  string tfile_pat="";  // trigger file pattern
  int gps_start=-1;     // GPS start
  int gps_end=-1;       // GPS end
  double snrmin=-1.0;   // SNR min
  double snrmax=1.0e20; // SNR max
  double freqmin=-1.0;  // frequency min
  double freqmax=1.0e20;// frequency max
  string printtype="clusters";// print type
  double cluster_dt=0.1;// cluster time window
  bool ptime=true;
  bool pfreq=true;
  bool psnr=true;
  bool pq=false;
  bool pamp=false;
  bool pph=false;
  bool ptstart=false;
  bool ptend=false;
  bool pfstart=false;
  bool pfend=false;
  bool pduration=false;
  bool pbandwidth=false;

  // loop over arguments
  vector <string> sarg;
  for(int a=1; a<argc; a++){
    sarg=SplitString((string)argv[a],'=');
    if(sarg.size()!=2) continue;
    if(!sarg[0].compare("channel"))        chname=(string)sarg[1];
    if(!sarg[0].compare("file"))           tfile_pat=sarg[1];
    if(!sarg[0].compare("gps-start"))      gps_start=atoi(sarg[1].c_str());
    if(!sarg[0].compare("gps-end"))        gps_end=atoi(sarg[1].c_str());
    if(!sarg[0].compare("snr-min"))        snrmin=atof(sarg[1].c_str());
    if(!sarg[0].compare("snr-max"))        snrmax=atof(sarg[1].c_str());
    if(!sarg[0].compare("freq-min"))       freqmin=atof(sarg[1].c_str());
    if(!sarg[0].compare("freq-max"))       freqmax=atof(sarg[1].c_str());
    if(!sarg[0].compare("print"))          printtype=(string)sarg[1];
    if(!sarg[0].compare("cluster-dt"))     cluster_dt=atof(sarg[1].c_str());
    if(!sarg[0].compare("print-time"))     ptime=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("print-freq"))     pfreq=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("print-snr"))      psnr=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("print-q"))        pq=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("print-amplitude")) pamp=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("print-phase"))    pph=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("print-time-start")) ptstart=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("print-time-end")) ptend=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("print-freq-start")) pfstart=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("print-freq-end")) pfend=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("print-duration")) pduration=!!(atoi(sarg[1].c_str()));
    if(!sarg[0].compare("print-bandwidth")) pbandwidth=!!(atoi(sarg[1].c_str()));
  }

  // centralized trigger files
  if(!tfile_pat.compare("")){
    if(chname.compare("")&&gps_start>0&&gps_end>0) tfile_pat=GetOmicronFilePattern(chname,gps_start,gps_end);// get centralized trigger files
    else{
      cerr<<"Triggers must be specified, either with a channel name and a time range or with a file pattern"<<endl;
      cerr<<"Type omicron-print.exe for help"<<endl;
      return 1;
    }
  }
  stringstream tmpstream;

  // check file pattern
  if(!tfile_pat.compare("")){
    cerr<<"No trigger files"<<endl;
    return 2;
  }
  
  // triggers
  ReadTriggers *RT = new ReadTriggers(tfile_pat,"",0);
  if(!RT->GetStatus()||!RT->GetSegments()->GetLiveTime()) return 2;

  // use trigger time range
  if(gps_start<0) gps_start=RT->GetTimeMin();
  if(gps_end<0)   gps_end=RT->GetTimeMax();

  //*************** print segments
  if(!printtype.compare("segments")){

    Segments *S = new Segments(gps_start, gps_end);
    S->Intersect(RT->GetSegments());
    S->Dump();
    delete S;
  }

  //*************** print triggers
  else if(!printtype.compare("triggers")){
    
    RT->SetTriggerBranchStatus("*",false);
    RT->SetTriggerBranchStatus("tstart",true);
    RT->SetTriggerBranchStatus("tend",true);
    RT->SetTriggerBranchStatus("frequency",true);
    RT->SetTriggerBranchStatus("snr",true);
    
    // header + optimize speed
    cout<<"# raw triggers";
    cout<<endl;
    if(ptstart) cout<<"# starting time [GPS]"<<endl;
    if(ptime){
      RT->SetTriggerBranchStatus("time",true);
      cout<<"# peak time [GPS]"<<endl;
    }
    if(ptend) cout<<"# ending time [GPS]"<<endl;
    if(pduration) cout<<"# duration [s]"<<endl;
    if(pfstart){
      RT->SetTriggerBranchStatus("fstart",true);
      cout<<"# starting frequency [Hz]"<<endl;
    }
    if(pfreq) cout<<"# peak frequency [Hz]"<<endl;
    if(pfend){
      RT->SetTriggerBranchStatus("fend",true);
      cout<<"# ending frequency [Hz]"<<endl;
    }
    if(pbandwidth){
      RT->SetTriggerBranchStatus("fstart",true);
      RT->SetTriggerBranchStatus("fend",true);
      cout<<"# bandwidth [Hz]"<<endl;
    }
    if(pq){
      RT->SetTriggerBranchStatus("q",true);
      cout<<"# Q [-]"<<endl;
    }
    if(psnr) cout<<"# SNR [-]"<<endl;
    if(pamp){
      RT->SetTriggerBranchStatus("amplitude",true);
      cout<<"# amplitude [Hz^-1/2]"<<endl;
    }
    if(pph){
      RT->SetTriggerBranchStatus("phase",true);
      cout<<"# phase [rad]"<<endl;
    }

    for(int c=0; c<RT->GetNtriggers(); c++){
      // trigger selection
      if(RT->GetTriggerTimeEnd(c)<gps_start) continue;
      if(RT->GetTriggerTimeStart(c)>=gps_end) break;
      //if(RT->GetTriggerTime(c)<gps_start||RT->GetTriggerTime(c)>=gps_end) continue;
      if(RT->GetTriggerSNR(c)<snrmin||RT->GetTriggerSNR(c)>=snrmax) continue;
      if(RT->GetTriggerFrequency(c)<freqmin||RT->GetTriggerFrequency(c)>=freqmax) continue;
      
      // print
      if(ptstart)    cout<<fixed<<setprecision(4)<<RT->GetTriggerTimeStart(c)<<" ";
      if(ptime)      cout<<fixed<<setprecision(4)<<RT->GetTriggerTime(c)<<" ";
      if(ptend)      cout<<fixed<<setprecision(4)<<RT->GetTriggerTimeEnd(c)<<" ";
      if(pduration)  cout<<fixed<<setprecision(4)<<RT->GetTriggerDuration(c)<<" ";
      if(pfstart)    cout<<fixed<<setprecision(2)<<RT->GetTriggerFrequencyStart(c)<<" ";
      if(pfreq)      cout<<fixed<<setprecision(2)<<RT->GetTriggerFrequency(c)<<" ";
      if(pfend)      cout<<fixed<<setprecision(2)<<RT->GetTriggerFrequencyEnd(c)<<" ";
      if(pbandwidth) cout<<fixed<<setprecision(2)<<RT->GetTriggerBandWidth(c)<<" ";
      if(pq)         cout<<fixed<<setprecision(2)<<RT->GetTriggerQ(c)<<" ";
      if(psnr)       cout<<fixed<<setprecision(2)<<RT->GetTriggerSNR(c)<<" ";
      if(pamp)       cout<<scientific<<setprecision(4)<<RT->GetTriggerAmplitude(c)<<" ";
      if(pph)        cout<<fixed<<setprecision(4)<<RT->GetTriggerPhase(c)<<" ";
      cout<<endl;

    }
  }

  //*************** print clusters
  else{
    RT->SetClusterizeDt(cluster_dt);
    if(!RT->Clusterize(1)){ delete RT; return 3; }

    RT->SetClusterBranchStatus("*",false);
    RT->SetTriggerBranchStatus("*",false);
    RT->SetClusterBranchStatus("tstart",true);
    RT->SetClusterBranchStatus("tend",true);
    RT->SetClusterBranchStatus("frequency",true);
    RT->SetClusterBranchStatus("snr",true);

    // header + optimize speed
    cout<<"# time-clustered triggers, dt = "<<cluster_dt<<"s";
    cout<<endl;
    if(ptstart) cout<<"# starting time [GPS]"<<endl;
    if(ptime){
      RT->SetClusterBranchStatus("time",true);
      cout<<"# peak time [GPS]"<<endl;
    }
    if(ptend) cout<<"# ending time [GPS]"<<endl;
    if(pduration) cout<<"# duration [s]"<<endl;
    if(pfstart){
      RT->SetClusterBranchStatus("fstart",true);
      cout<<"# starting frequency [Hz]"<<endl;
    }
    if(pfreq) cout<<"# peak frequency [Hz]"<<endl;
    if(pfend){
      RT->SetClusterBranchStatus("fend",true);
      cout<<"# ending frequency [Hz]"<<endl;
    }
    if(pbandwidth){
      RT->SetClusterBranchStatus("fstart",true);
      RT->SetClusterBranchStatus("fend",true);
      cout<<"# bandwidth [Hz]"<<endl;
    }
    if(pq){
      RT->SetClusterBranchStatus("q",true);
      cout<<"# Q [-]"<<endl;
    }
    if(psnr) cout<<"# SNR [-]"<<endl;
    if(pamp){
      RT->SetClusterBranchStatus("amplitude",true);
      cout<<"# amplitude [Hz^-1/2]"<<endl;
    }
    if(pph){
      RT->SetClusterBranchStatus("phase",true);
      cout<<"# phase [rad]"<<endl;
    }

    for(int c=0; c<RT->GetNclusters(); c++){
      // cluster selection
      if(RT->GetClusterTimeEnd(c)<gps_start) continue;
      if(RT->GetClusterTimeStart(c)>=gps_end) break;
      //if(RT->GetClusterTime(c)<gps_start||RT->GetClusterTime(c)>=gps_end) continue;
      if(RT->GetClusterSNR(c)<snrmin||RT->GetClusterSNR(c)>=snrmax) continue;
      if(RT->GetClusterFrequency(c)<freqmin||RT->GetClusterFrequency(c)>=freqmax) continue;
      
      // print
      if(ptstart)    cout<<fixed<<setprecision(4)<<RT->GetClusterTimeStart(c)<<" ";
      if(ptime)      cout<<fixed<<setprecision(4)<<RT->GetClusterTime(c)<<" ";
      if(ptend)      cout<<fixed<<setprecision(4)<<RT->GetClusterTimeEnd(c)<<" ";
      if(pduration)  cout<<fixed<<setprecision(4)<<RT->GetClusterDuration(c)<<" ";
      if(pfstart)    cout<<fixed<<setprecision(2)<<RT->GetClusterFrequencyStart(c)<<" ";
      if(pfreq)      cout<<fixed<<setprecision(2)<<RT->GetClusterFrequency(c)<<" ";
      if(pfend)      cout<<fixed<<setprecision(2)<<RT->GetClusterFrequencyEnd(c)<<" ";
      if(pbandwidth) cout<<fixed<<setprecision(2)<<RT->GetClusterBandWidth(c)<<" ";
      if(pq)         cout<<fixed<<setprecision(2)<<RT->GetClusterQ(c)<<" ";
      if(psnr)       cout<<fixed<<setprecision(2)<<RT->GetClusterSNR(c)<<" ";
      if(pamp)       cout<<scientific<<setprecision(4)<<RT->GetClusterAmplitude(c)<<" ";
      if(pph)        cout<<fixed<<setprecision(4)<<RT->GetClusterPhase(c)<<" ";
      cout<<endl;
    }
  }


    
  delete RT;
  
  return 0;
}

