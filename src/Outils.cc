//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

////////////////////////////////////////////////////////////////////////////////////
void Omicron::PrintASCIIlogo(void){
////////////////////////////////////////////////////////////////////////////////////

  cout<<endl;
  cout<<endl;
  cout<<"############################################################################"<<endl;
  cout<<"############################################################################"<<endl;
  cout<<endl;
  cout<<"          ██████╗ ███╗   ███╗██╗ ██████╗██████╗  ██████╗ ███╗   ██╗"<<endl;
  cout<<"         ██╔═══██╗████╗ ████║██║██╔════╝██╔══██╗██╔═══██╗████╗  ██║"<<endl;
  cout<<"         ██║   ██║██╔████╔██║██║██║     ██████╔╝██║   ██║██╔██╗ ██║"<<endl;
  cout<<"         ██║   ██║██║╚██╔╝██║██║██║     ██╔══██╗██║   ██║██║╚██╗██║"<<endl;
  cout<<"         ╚██████╔╝██║ ╚═╝ ██║██║╚██████╗██║  ██║╚██████╔╝██║ ╚████║"<<endl;
  cout<<"          ╚═════╝ ╚═╝     ╚═╝╚═╝ ╚═════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝"<<endl;
  cout<<"                                                           "<<GetVersion()<<endl;
  cout<<"############################################################################"<<endl;
  cout<<"############################################################################"<<endl;
  cout<<endl;
  cout<<endl;

  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::PrintMessage(const string aMessage){
////////////////////////////////////////////////////////////////////////////////////
  time ( &timer );
  ptm = gmtime ( &timer );
  cout<<">>>>>>>>>> ("<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) [+"<<timer-timer_start<<"s]: "<<aMessage<<" <<<<<<<<<<"<<endl;
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::PrintStatusInfo(void){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK) return;
  if(!inSegments->GetNsegments()) return;
  
  cout<<"\n************* Omicron status info *************"<<endl;
  cout<<"requested start         = "<<(int)inSegments->GetStart(0)<<endl;
  cout<<"requested end           = "<<(int)inSegments->GetEnd(inSegments->GetNsegments()-1)<<endl;
  cout<<"requested livetime      = "<<(int)inSegments->GetLiveTime()<<"s"<<endl;
  cout<<"number of loaded chunks = "<<chunk_ctr<<endl;

  for(int c=0; c<nchannels; c++){
    cout<<"\n*** "<<triggers[c]->GetName()<<endl;
    cout<<"number of calls                = "<<chan_ctr[c]<<endl;
    cout<<"number of data calls           = "<<chan_data_ctr[c]<<endl;
    cout<<"number of conditioning calls   = "<<chan_cond_ctr[c]<<endl;
    cout<<"number of projection calls     = "<<chan_proj_ctr[c]<<endl;
    cout<<"number of write calls          = "<<chan_write_ctr[c]<<endl;
    if(outSegments[c]->GetNsegments()){
      cout<<"start_out           = "<<(int)outSegments[c]->GetStart(0)<<endl;
      cout<<"end_out             = "<<(int)outSegments[c]->GetEnd(outSegments[c]->GetNsegments()-1)<<endl;
      cout<<"trigger livetime    = "<<(int)outSegments[c]->GetLiveTime()<<"s ("<<outSegments[c]->GetLiveTime()/inSegments->GetLiveTime()*100<<"%)"<<endl;
    }
  }
  cout<<"***********************************************\n"<<endl;

  return;
}

////////////////////////////////////////////////////////////////////////////////////
double* Omicron::GetTukeyWindow(const int aSize, const int aFractionSize){
////////////////////////////////////////////////////////////////////////////////////

  int FracSize_2=aFractionSize/2;
  double *Window = new double [aSize];
  
  double factor = TMath::Pi()/(double)FracSize_2;

  for (int i=0; i<FracSize_2; i++)
    Window[i] = 0.5*(1+TMath::Cos(factor*(double)(i-FracSize_2)));
  for (int i=FracSize_2; i<aSize-FracSize_2; i++)
    Window[i] = 1.0;
  for (int i=aSize-FracSize_2; i<aSize; i++)
    Window[i] = 0.5*(1+TMath::Cos(factor*(double)(i-aSize+FracSize_2)));
 
  return Window;
}


////////////////////////////////////////////////////////////////////////////////////
bool Omicron::IsFlat(const int aInVectSize, double *aInVect){
////////////////////////////////////////////////////////////////////////////////////
  double val = aInVect[0];
  for(int i=0; i<aInVectSize; i++){
    if(aInVect[i]!=val) return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
Segments* Omicron::GetTriggerSegments(TH1D *aThr, const double aInfValue){
////////////////////////////////////////////////////////////////////////////////////
  Segments* empty = new Segments();
  if(!status_OK){
    cerr<<"Omicron::GetTriggerSegments: the Omicron object is corrupted"<<endl;
    return empty;
  }
  delete empty;

  return triggers[chanindex]->GetTriggerSegments(aThr,aInfValue);
}

////////////////////////////////////////////////////////////////////////////////////
Segments* Omicron::GetLastTriggerSegments(TH1D *aThr, const double aInfValue){
////////////////////////////////////////////////////////////////////////////////////
  Segments* empty = new Segments();
  if(!status_OK){
    cerr<<"Omicron::GetLasTriggerSegments: the Omicron object is corrupted"<<endl;
    return empty;
  }
  delete empty;

  return triggers[chanindex]->GetLastTriggerSegments(aThr,aInfValue);
}

////////////////////////////////////////////////////////////////////////////////////
string Omicron::GetColorCode(const double aSNRratio){
////////////////////////////////////////////////////////////////////////////////////
  if(aSNRratio<=0) return "";
  double inc = 0.1;

  if(aSNRratio>=17.0*inc) return colorcode[16];// saturation

  return colorcode[((int)floor(aSNRratio/inc))%17];
}

////////////////////////////////////////////////////////////////////////////////////
vector <string> Omicron::GetChannels(void){
////////////////////////////////////////////////////////////////////////////////////
  vector <string> chanlist;
  for(int c=0; c<nchannels; c++) chanlist.push_back(triggers[c]->GetName());
  return chanlist;
}

