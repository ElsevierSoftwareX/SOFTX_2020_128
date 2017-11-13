#include "OmicronUtils.h"

string GetOmicronFilePattern(const string aChannelName, const int aTimeStart, const int aTimeEnd){

  // trigger directory
  string trigdir = getenv("OMICRON_TRIGGERS");

  // special case of hpss
  if(trigdir.find("hpss")!=string::npos) return GetOmicronFilePatternFromHpss(aChannelName, aTimeStart, aTimeEnd);
  
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

string GetOmicronFilePatternFromHpss(const string aChannelName, const int aTimeStart, const int aTimeEnd){
  string filelist="";

  // trigger directory
  string trigdir = getenv("OMICRON_TRIGGERS");

  // channel stream
  Streams *S = new Streams(aChannelName,0);

  // random id
  srand(time(NULL));
  int randid = rand();

  // tmp file
  stringstream ss;
  ss<<getenv("TMP")<<"/omf-"<<S->GetNameSuffixUnderScore()<<"-"<<aTimeStart<<"-"<<aTimeEnd<<"."<<randid<<".tmp";
  string tmpfile=ss.str();
  ss.clear(); ss.str("");

  // rfdir command
  ss<<"rfdir ${OMICRON_TRIGGERS}/"<<S->GetNamePrefix()<<"/"<<S->GetNameSuffixUnderScore()<<"_OMICRON | grep -v \"\\.\" | awk '{print $9}' | sort -u> "<<tmpfile;

  // dump list of channel sub-directories in a tmp file
  if(system(ss.str().c_str())) return filelist;
  ss.clear(); ss.str("");

  ReadAscii *R = new ReadAscii(tmpfile,"i");
  if(!R->GetNRow()) { delete R; return filelist; }

  // list of gps directories
  int stoproot = aTimeEnd / 100000;
  vector <int> gpsdir; int gps;
  for(int g=0; g<R->GetNRow(); g++){
    R->GetElement(gps,g,0);
    if(gps<=stoproot) gpsdir.push_back(gps);
  }

  delete R;
  
  // list relevant files in gps directories
  vector <string> vfilefrag; string filename;
  int fstart, fstop;
  for(int g=0; g<(int)gpsdir.size(); g++){
    ss<<"rfdir ${OMICRON_TRIGGERS}/"<<S->GetNamePrefix()<<"/"<<S->GetNameSuffixUnderScore()<<"_OMICRON/"<<gpsdir[g]<<" | grep OMICRON | awk '{print $9}' | sort -u> "<<tmpfile;
    if(system(ss.str().c_str())) return "";
    ss.clear(); ss.str("");

    R = new ReadAscii(tmpfile,"s");

    // select relevant files in that gps directories
    for(int gg=0; gg<R->GetNRow(); gg++){
      R->GetElement(filename,gg,0);
     
      // get file fragments
      vfilefrag.clear();
      vfilefrag = SplitString(filename,'-');
    
      // check file naming convention (4 fragments)
      if(vfilefrag.size()!=4) continue;
    
      // file timing
      fstart = atoi(vfilefrag[2].c_str());
      fstop = fstart+atoi(vfilefrag[3].substr(0,vfilefrag[3].size()-5).c_str());
    
      // check file start time
      if(fstart>=aTimeEnd) continue;
    
      // check file stop time
      if(fstop<aTimeStart) continue;

      // full file name
      ss<<"root://ccxroot:1999/${OMICRON_TRIGGERS}/"<<S->GetNamePrefix()<<"/"+S->GetNameSuffixUnderScore()<<"_OMICRON/"<<gpsdir[g]<<"/"<<filename;

      // select file
      filelist+=ss.str()+" ";
      ss.clear(); ss.str("");

    }

    system(("rm -f "+tmpfile).c_str());

    delete R;

  }

  return filelist;
}
