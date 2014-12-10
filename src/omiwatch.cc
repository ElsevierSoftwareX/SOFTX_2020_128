//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

using namespace std;

int main (int argc, char* argv[]){
  gErrorIgnoreLevel = 3000;

  // check the argument
  if(argc<3){
    cerr<<endl;
    cerr<<" usage:"<<endl; 
    cerr<<endl;
    cerr<<" omiscan [central time] [option file #1] [option file #2] ... [option file #N]"<<endl; 
    cerr<<endl;
    cerr<<" [central time]    =  GPS where the maps will be centred"<<endl; 
    cerr<<" [option file #i]  =  Omicron option file"<<endl; 
    cerr<<"                      For available options, please visit:"<<endl; 
    cerr<<"                      http://virgo.in2p3.fr/GWOLLUM/Friends/omicron.html"<<endl; 
    cerr<<endl;
    return 1;
  }
  
  if(!IsDirectory("$OMICRON_ONLINE_TRIGGERS")){
    cerr<<"omiwatch ERROR: there is no online trigger directory"<<endl;
    return 2;
  }

  vector <string> chandir;
  if(!ListDirectories("$OMICRON_ONLINE_TRIGGERS",chandir)){
    cerr<<"omiwatch ERROR: cannot scan the online trigger directory"<<endl;
    return 2;
  }

  ReadTriggerSegments *R;
  Segments *S;

  // loop over channels
  for(int c=0; c<(int)chandir.size(); c++){
    cout<<"\n**** testing "<<chandir[c]<<endl;

    R = new ReadTriggerSegments("${OMICRON_ONLINE_TRIGGERS}/"+chandir[c]+"/*.root");
    S = R->GetSegments();
 
    // test segment status
    if(!S->GetStatus()){
      cerr<<"omiwatch "<<chandir[c]<<": triggers segments are corrupted"<<endl;
      delete R;
      continue;
    }

    // test segment livetime
    if(!S->GetLiveTime()){
      cerr<<"omiwatch "<<chandir[c]<<": no livetime"<<endl;
      delete R;
      continue;
    }

    cout<<"OK"<<endl;
    delete R;
  }

  return 0;
}
  
