//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Otile.h"

ClassImp(Otile)

////////////////////////////////////////////////////////////////////////////////////
Otile::Otile(const int aTimeRange, const int aTimePad,
	     const double aQMin, const double aQMax, 
	     const double aFrequencyMin, const double aFrequencyMax, 
	     const int aSampleFrequency, const double aMaximumMismatch, 
	     const double aSNRThreshold, const int aVerbosity){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // save parameters
  fTimeRange=aTimeRange;
  fTimePad=aTimePad;
  fQMin=aQMin;
  fQMax=aQMax;
  fFrequencyMin=aFrequencyMin;
  fFrequencyMax=aFrequencyMax;
  fSampleFrequency=aSampleFrequency;
  fMaximumMismatch=aMaximumMismatch;
  fSNRThreshold=aSNRThreshold;
  fVerbosity=aVerbosity;
  
  // init
  status_OK=true;
  fNumberOfTiles=0;
  
  // adjust Q range
  if(fQMin<sqrt(11.0)) fQMin=sqrt(11.0);
  if(fQMin>=fQMax){
    cerr<<"Otile::Otile: the Q range does not make sense"<<endl;
    status_OK=false;
  }

  // number of planes
  double QCumulativeMismatch = log(fQMax/fQMin)/sqrt(2.0);// cumulative mismatch across Q range
  fMismatchStep = 2.0*sqrt(fMaximumMismatch/3.0);
  fNumberOfPlanes = (int)ceil(QCumulativeMismatch / fMismatchStep);
  if(fNumberOfPlanes<=0) fNumberOfPlanes=1;
  fQMismatchStep = QCumulativeMismatch / fNumberOfPlanes;
  
  // nulling
  for(int p=0; p<NQPLANEMAX; p++) qplanes[p]=NULL;

  // check parameters
  status_OK*=CheckParameters();

  // build Q-planes
  if(status_OK){
    fQs.clear();
    for(int i=0; i<fNumberOfPlanes; i++)
      fQs.push_back(fQMin * exp(sqrt(2.0) * (0.5+i) * fQMismatchStep));
    
    // create Q planes
    CreatePlanes();
  }
  
  if(!status_OK) cerr<<"Otile::Otile: initialization failed!"<<endl;
  if(fVerbosity>1) PrintInfo();
}

////////////////////////////////////////////////////////////////////////////////////
Otile::~Otile(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"Otile::~Otile"<<endl;
  fQs.clear();
  for(int p=0; p<fNumberOfPlanes; p++) if(qplanes[p]!=NULL) delete qplanes[p];
}
  
////////////////////////////////////////////////////////////////////////////////////
bool Otile::GetTriggers(Triggers *aTriggers, double *aDataRe, double *aDataIm, const int aTimeStart){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Otile::GetTriggers: the Otile object is corrupted"<<endl;
    return false;
  }

  // get triggers for each Q plane
  for(int p=0; p<fNumberOfPlanes; p++){
    if(!qplanes[p]->GetTriggers(aTriggers,aDataRe,aDataIm,aTimeStart)){
      cerr<<"Otile::GetTriggers: the q-plane "<<p<<" has failed to produce triggers"<<endl;
      return false;
    }
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
TH2D* Otile::GetMap(const int qindex, double *aDataRe, double *aDataIm, const double time_offset, const bool printamplitude){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Otile::GetMap: the Otile object is corrupted"<<endl;
    return NULL;
  }
  if(qindex<0||qindex>=fNumberOfPlanes){
    cerr<<"Otile::GetMap: this q index is not allowed"<<endl;
    return NULL;
  }

  return qplanes[qindex]->GetMap(aDataRe,aDataIm,time_offset,printamplitude);
}

////////////////////////////////////////////////////////////////////////////////////
bool Otile::SetPowerSpectrum(double *aPSD, const int aPSDsize){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Otile::SetPowerSpectrum: the Otile object is corrupted"<<endl;
    return false;
  }

  // set power for each Q plane
  for(int p=0; p<fNumberOfPlanes; p++){
    if(!qplanes[p]->SetPowerSpectrum(aPSD,aPSDsize)){
      cerr<<"Otile::SetPowerSpectrum: the q-plane "<<p<<" has failed"<<endl;
      return false;
    }
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
void Otile::CreatePlanes(void){
////////////////////////////////////////////////////////////////////////////////////

  // loop over planes
  for(int p=0; p<fNumberOfPlanes; p++){
    qplanes[p]=new Oqplane(fQs[p],fSampleFrequency,fTimeRange,fTimePad,fFrequencyMin,fFrequencyMax,fMismatchStep,fSNRThreshold);
    status_OK*=qplanes[p]->status_OK;
    fNumberOfTiles+=qplanes[p]->fNumberOfTiles;
  }
  return;
}


////////////////////////////////////////////////////////////////////////////////////
bool Otile::CheckParameters(void){
////////////////////////////////////////////////////////////////////////////////////  
  if(!status_OK) return false;

  if(fTimeRange<=0){
    cerr<<"Otile::CheckParameters: time range does not make any sense"<<endl;
    return false;
  }
  if(2*fTimePad>=fTimeRange){
    cerr<<"Otile::CheckParameters: time pads do not make any sense"<<endl;
    return false;
  }
  if(fQMax<=fQMin){
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
  if(fNumberOfPlanes > NQPLANEMAX){
    cerr<<"Otile::CheckParameters: the number of Q-planes is too large: "<<fNumberOfPlanes<<">"<<NQPLANEMAX<<endl;
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
  if(qindex<0||qindex>=fNumberOfPlanes){
    cerr<<"Otile::GetQ: this q index is not allowed"<<endl;
    return -1.0;
  }

  return fQs[qindex];
}

////////////////////////////////////////////////////////////////////////////////////
void Otile::PrintInfo(void){
////////////////////////////////////////////////////////////////////////////////////  

  cout<<"\n--- Otile::PrintInfo ---"<<endl;
  cout<<" Time Range: "<<fTimeRange<<endl;
  cout<<" Q range: "<<fQMin<<" - "<<fQMax<<endl;
  cout<<" Frequency range: "<<fFrequencyMin<<" - "<<fFrequencyMax<<endl;
  cout<<" Sampling frequency: "<<fSampleFrequency<<endl;
  cout<<" Fractional energy loss due to mismatch: "<<fMaximumMismatch<<endl;
  cout<<" Mismatch step between tiles: "<<fMismatchStep<<endl;
  cout<<" Number of Q planes: "<<fNumberOfPlanes<<endl;
  cout<<" Mismatch step between neighboring planes: "<<fQMismatchStep<<endl;
  cout<<" First Q: "<<fQs[0]<<endl;
  cout<<" Last Q: "<<fQs[fNumberOfPlanes-1]<<endl;
  cout<<" Total number of tiles: "<<fNumberOfTiles<<endl<<endl;
  
  cout<<" List of Qs: "<<endl;
  for(int p=0; p<fNumberOfPlanes; p++) cout<<"   "<<p<<":"<<fQs[p]<<endl;
  cout<<" Effective frequency range: "<<endl;
  for(int p=0; p<fNumberOfPlanes; p++) cout<<"   "<<p<<":"<<qplanes[p]->fFrequencyMin<<"-"<<qplanes[p]->fFrequencyMax<<endl;
  cout<<" Number of frequency rows: "<<endl;
  for(int p=0; p<fNumberOfPlanes; p++) cout<<"   "<<p<<":"<<qplanes[p]->fNumberOfRows<<endl;
  cout<<" Number of tiles: "<<endl;
  for(int p=0; p<fNumberOfPlanes; p++) cout<<"   "<<p<<":"<<qplanes[p]->fNumberOfTiles<<endl;

  return;
  
}
