//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Otile.h"

ClassImp(Otile)

////////////////////////////////////////////////////////////////////////////////////
Otile::Otile(const int aTimeRange,
	     const double aQMin, const double aQMax, 
	     const double aFrequencyMin, const double aFrequencyMax, 
	     const int aSampleFrequency, const double aMaximumMismatch, 
	     const int aVerbosity): GwollumPlot("otile","GWOLLUM"){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // save parameters
  fTimeRange=aTimeRange;
  fQMin=aQMin;
  fQMax=aQMax;
  fFrequencyMin=aFrequencyMin;
  fFrequencyMax=aFrequencyMax;
  fSampleFrequency=aSampleFrequency;
  fMaximumMismatch=aMaximumMismatch;
  fVerbosity=aVerbosity;
  
  // init
  fNumberOfTiles=0;
  
  // adjust Q min
  if(fQMin<sqrt(11.0)) fQMin=sqrt(11.0);
  if(fQMax<sqrt(11.0)) fQMax=sqrt(11.0);

  // check parameters
  status_OK=CheckParameters();
  fMismatchStep=2.0*sqrt(aMaximumMismatch/3.0);

  // compute Q values
  fQs = ComputeQs(fQMin,fQMax,fMaximumMismatch);
 
  // loop over planes
  qplanes = new Oqplane* [(int)fQs.size()];
  for(int q=0; q<(int)fQs.size(); q++){
    qplanes[q]=new Oqplane(fQs[q],fSampleFrequency,fTimeRange,fFrequencyMin,fFrequencyMax,fMismatchStep);
    fNumberOfTiles+=qplanes[q]->Ntiles;
  }
    
  if(aVerbosity){
    for(int q=0; q<(int)fQs.size(); q++){
      cout<<"Otile::Otile: plane #"<<q<<endl;
      if(aVerbosity>1) qplanes[q]->PrintParameters();
    }
  }
  
}

////////////////////////////////////////////////////////////////////////////////////
Otile::~Otile(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"Otile::~Otile"<<endl;
  fQs.clear();
  for(int p=0; p<(int)fQs.size(); p++) delete qplanes[p];
  delete qplanes;
}
 
/* 
////////////////////////////////////////////////////////////////////////////////////
bool Otile::GetTriggers(MakeTriggers *aTriggers, const int aTimeStart, const int aExtraTimePadMin){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Otile::GetTriggers: the Otile object is corrupted"<<endl;
    return false;
  }

  // get triggers for each Q plane
  for(int p=0; p<(int)fQs.size(); p++){
    if(!qplanes[p]->GetTriggers(aTriggers,aTimeStart,aExtraTimePadMin)){
      cerr<<"Otile::GetTriggers: the q-plane #"<<p<<" has failed to produce triggers"<<endl;
      return false;
    }
  }

  return true;
}
*/

/*
////////////////////////////////////////////////////////////////////////////////////
TH2D* Otile::GetMap(const int qindex, double *aDataRe, double *aDataIm, const double time_offset, const bool printamplitude){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Otile::GetMap: the Otile object is corrupted"<<endl;
    return NULL;
  }
  if(qindex<0||qindex>=(int)fQs.size()){
    cerr<<"Otile::GetMap: this q index is not allowed"<<endl;
    return NULL;
  }

  return qplanes[qindex]->GetMap(aDataRe,aDataIm,time_offset,printamplitude);
}
*/

////////////////////////////////////////////////////////////////////////////////////
bool Otile::SetPower(Spectrum *aSpec){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Otile::SetPower: the Otile object is corrupted"<<endl;
    return false;
  }
  if(!aSpec->GetStatus()){
    cerr<<"Otile::SetPower: the Spectrum object is corrupted"<<endl;
    return false;
  }

  // set power for each Q plane
  for(int p=0; p<(int)fQs.size(); p++){
    if(!qplanes[p]->SetPower(aSpec)){
      cerr<<"Otile::SetPower: cannot set power for plane #"<<p<<endl;
      return false;
    }
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::CheckParameters(void){
////////////////////////////////////////////////////////////////////////////////////  
  if(fTimeRange<1){
    cerr<<"Otile::CheckParameters: time range does not make any sense"<<endl;
    return false;
  }
  if(fQMax<fQMin){
    cerr<<"Otile::CheckParameters: the Q range does not make any sense"<<endl;
    return false;
  }
  if(fFrequencyMax<=fFrequencyMin){
    cerr<<"Otile::CheckParameters: the frequency range does not make any sense"<<endl;
    return false;
  }
  if(fMaximumMismatch>0.5){
    cerr<<"Otile::CheckParameters: maximum mismatch is not reasonable (has to be <=0.5)"<<endl;
    return false;
  }
  if(!IsPowerOfTwo(fSampleFrequency * fTimeRange)){
    cerr<<"Otile::CheckParameters: data length is not an integer power of two"<<endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
double Otile::GetQ(const int qindex){
////////////////////////////////////////////////////////////////////////////////////  
  if(!status_OK){
    cerr<<"Otile::GetQ: the Otile object is corrupted"<<endl;
    return -1.0;
  }
  if(qindex<0||qindex>=(int)fQs.size()){
    cerr<<"Otile::GetQ: this q index is not allowed"<<endl;
    return -1.0;
  }

  return fQs[qindex];
}

////////////////////////////////////////////////////////////////////////////////////
vector <double> Otile::ComputeQs(const double aQMin, const double aQMax, const double aMaximumMismatch){
////////////////////////////////////////////////////////////////////////////////////  
  
  // number of planes
  double QCumulativeMismatch = log(aQMax/aQMin)/sqrt(2.0);// cumulative mismatch across Q range
  double mismatchstep = 2.0*sqrt(aMaximumMismatch/3.0);
  int n = (int)ceil(QCumulativeMismatch/mismatchstep);
  if(n<=0) n=1;
  double Qmismatchstep = QCumulativeMismatch/(double)n;
  
  // compute Q values
  vector <double> qs;
  for(int i=0; i<n; i++) qs.push_back(aQMin * exp(sqrt(2.0) * (0.5+i) * Qmismatchstep));
  return qs;
}

////////////////////////////////////////////////////////////////////////////////////
void Otile::PrintInfo(void){
////////////////////////////////////////////////////////////////////////////////////  
  if(!status_OK){
    cerr<<"Otile::PrintInfo: the Otile object is corrupted"<<endl;
    return;
  }
  cout<<"\n--- Otile::PrintInfo ---"<<endl;
  cout<<" Time Range: "<<fTimeRange<<endl;
  cout<<" Q range: "<<fQMin<<" - "<<fQMax<<endl;
  cout<<" Frequency range: "<<fFrequencyMin<<" - "<<fFrequencyMax<<endl;
  cout<<" Sampling frequency: "<<fSampleFrequency<<endl;
  cout<<" Fractional energy loss due to mismatch: "<<fMaximumMismatch<<endl;
  cout<<" Mismatch step between tiles: "<<fMismatchStep<<endl;
  cout<<" Number of Q planes: "<<fQs.size()<<endl;
  cout<<" Total number of tiles: "<<fNumberOfTiles<<endl<<endl;
  
  /*
  cout<<" List of Qs: "<<endl;
  for(int p=0; p<(int)fQs.size(); p++) cout<<"   "<<p<<":"<<fQs[p]<<endl;
  cout<<" Effective frequency range: "<<endl;
  for(int p=0; p<(int)fQs.size(); p++) cout<<"   "<<p<<":"<<qplanes[p]->fFrequencyMin<<"-"<<qplanes[p]->fFrequencyMax<<endl;
  cout<<" Number of frequency rows: "<<endl;
  for(int p=0; p<(int)fQs.size(); p++) cout<<"   "<<p<<":"<<qplanes[p]->fNumberOfRows<<endl;
  cout<<" Number of tiles: "<<endl;
  for(int p=0; p<(int)fQs.size(); p++) cout<<"   "<<p<<":"<<qplanes[p]->fNumberOfTiles<<endl;
  */
  return;
  
}
