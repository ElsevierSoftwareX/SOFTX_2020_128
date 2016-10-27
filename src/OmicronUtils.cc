#include "OmicronUtils.h"

string GetOmicronFilePattern(const string aChannelName, const int aTimeStart, const int aTimeEnd){

  // trigger directory
  string trigdir = getenv("OMICRON_TRIGGERS");
  if(!IsDirectory(trigdir)) return "";
  
  // channel stream
  Streams *S = new Streams(aChannelName,0);
  
  // LIGO-Virgo directory
  stringstream lv_dir;
  lv_dir<< trigdir << "/" << S->GetNamePrefix() << "/" << S->GetNameSuffixUnderScore() << "_OMICRON";

  int g_start = aTimeStart/100000;
  int g_stop = aTimeEnd/100000;
  string filelist="";
  vector <string> vfilefrag;
  int fstart, fstop;
  
  // get list of files in the first directory
  stringstream sg;
  sg<<g_start;
  vector <string> vfiles=Glob((lv_dir.str()+"/"+sg.str()+"/"+S->GetNamePrefix()+"-*.root").c_str());
  sg.clear(); sg.str("");
  
  // loop over files and select the relevant ones
  for(int v=0; v<(int)vfiles.size(); v++){
    
    // get file fragments
    vfilefrag.clear();
    vfilefrag = SplitString(GetFileNameFromPath(vfiles[v]),'-');
    
    // check file naming convention (4 fragments)
    if(vfilefrag.size()!=4) continue;
    
    // file timing
    fstart = atoi(vfilefrag[2].c_str());
    fstop = fstart+atoi(vfilefrag[3].substr(0,vfilefrag[3].size()-5).c_str());
    
    // check file start time
    if(fstart>=aTimeEnd) break;
    
    // check file stop time
    if(fstop<aTimeStart) continue;
    
    // select file
    filelist+=vfiles[v]+" ";
  }
  vfiles.clear();
  vfilefrag.clear();

  
  // loop over intermediate GPS directory
  string idir;
  for(int g=g_start+1; g<g_stop; g+=1){
    sg<<g;
    idir=lv_dir.str()+"/"+sg.str();
    sg.clear(); sg.str("");
    if(IsDirectory(idir)) filelist+=(idir+"/"+S->GetNamePrefix()+"-*.root ");
  }

  // last directory
  if(g_stop!=g_start){

  // get list of files in the last directory
    sg<<g_stop;
    vfiles=Glob((lv_dir.str()+"/"+sg.str()+"/"+S->GetNamePrefix()+"-*.root").c_str());
    sg.clear(); sg.str("");

    // loop over files and select the relevant ones
    for(int v=0; v<(int)vfiles.size(); v++){
      
      // get file fragments
      vfilefrag.clear();
      vfilefrag = SplitString(GetFileNameFromPath(vfiles[v]),'-');
      
      // check file naming convention (4 fragments)
      if(vfilefrag.size()!=4) continue;
      
      // file timing
      fstart = atoi(vfilefrag[2].c_str());
      fstop = fstart+atoi(vfilefrag[3].substr(0,vfilefrag[3].size()-5).c_str());
      
      // check file start time
      if(fstart>=aTimeEnd) break;
      
      // check file stop time
      if(fstop<aTimeStart) continue;
      
      // select file
      filelist+=vfiles[v]+" ";
    }
    vfiles.clear();
    vfilefrag.clear();
  }

  return filelist;
}
