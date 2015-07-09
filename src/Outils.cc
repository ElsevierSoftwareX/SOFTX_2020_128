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
string Omicron::GetColorCode(const double aSNRratio){
////////////////////////////////////////////////////////////////////////////////////
  if(aSNRratio<=0) return "";
  double inc = 0.1;

  if(aSNRratio>=17.0*inc) return colorcode[16];// saturation

  return colorcode[((int)floor(aSNRratio/inc))%17];
}
