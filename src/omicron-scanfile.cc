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
  cerr<<"omicron-scanfile.exe channel=[channel name] \\"<<endl;
  cerr<<"                     file=[trigger file pattern] \\"<<endl;
  cerr<<"                     gps-start=[GPS start] \\"<<endl;
  cerr<<"                     gps-end=[GPS end]"<<endl;
  cerr<<endl;
  cerr<<"[channel name]            channel name used to retrieve centralized Omicron triggers"<<endl;
  cerr<<"[trigger file pattern]    file pattern to ROOT trigger files (GWOLLUM convention)"<<endl;
  cerr<<"[GPS start]               starting GPS time (integer only)"<<endl;
  cerr<<"[GPS end]                 stopping GPS time (integer only)"<<endl;
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

  // loop over arguments
  vector <string> sarg;
  for(int a=1; a<argc; a++){
    sarg=SplitString((string)argv[a],'=');
    if(sarg.size()!=2) continue;
    if(!sarg[0].compare("channel"))        chname=(string)sarg[1];
    if(!sarg[0].compare("file"))           tfile_pat=sarg[1];
    if(!sarg[0].compare("gps-start"))      gps_start=atoi(sarg[1].c_str());
    if(!sarg[0].compare("gps-end"))        gps_end=atoi(sarg[1].c_str());
  }

  // centralized trigger files
  if(!tfile_pat.compare("")){
    if(chname.compare("")&&gps_start>0&&gps_end>0) tfile_pat=GetOmicronFilePattern(chname,gps_start,gps_end);// get centralized trigger files
    else{
      cerr<<"Triggers must be specified, either with a channel name and a time range or with a file pattern"<<endl;
      cerr<<"Type omicron-scanfile.exe for help"<<endl;
      return 1;
    }
  }
  stringstream tmpstream;

  // check file pattern
  if(!tfile_pat.compare("")){
    cerr<<"No trigger files"<<endl;
    return 2;
  }

  // list of files
  string ls;
  ls = "ls -U "+tfile_pat;
  
  // file descriptor
  FILE *in;
  in = popen(ls.c_str(), "r");

  // loop over files
  char filename[100000];
  while(fscanf(in,"%s",filename) != EOF){
    
    TFile *tmp = new TFile(filename);
    if( !tmp || tmp->IsZombie() ){
      cout<<filename<<endl;
      continue;
    }
    if( tmp->ReadKeys() ) tmp->Close();
    else cout<<filename<<endl;
      
  }
  pclose(in);

  return 0;
}

