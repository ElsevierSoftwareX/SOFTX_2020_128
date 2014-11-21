//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Oplot.h"

ClassImp(Oplot)

////////////////////////////////////////////////////////////////////////////////////
Oplot::Oplot(const string aPattern, const string aDirectory, const int aVerbose) : TriggerPlot (5, aPattern, aDirectory, aVerbose){ 
////////////////////////////////////////////////////////////////////////////////////

  // Default SNR thresholds
  SetSNRThresholds(TMath::Min(5.0,Msnrmin_stat),8.0,10.0,20.0);
  
  // Clusterize
  TriggerPlot::Clusterize("TIME",1);
  for(int s=0; s<4; s++) TriggerPlot::SetCollectionUseClusters(s,1);// unflagged
  TriggerPlot::SetCollectionUseClusters(4,2);// flagged

  // set plot style
  TriggerPlot::SetDateFormat(true);
  TriggerPlot::SetCollectionColor(0,9);
  TriggerPlot::SetCollectionColor(1,4);
  TriggerPlot::SetCollectionColor(2,kGreen+3);
  TriggerPlot::SetCollectionColor(3,2);
  TriggerPlot::SetCollectionColor(4,1);
  TriggerPlot::SetCollectionMarker(0,1);
  TriggerPlot::SetCollectionMarker(1,7);
  TriggerPlot::SetCollectionMarker(2,4,0.8);
  TriggerPlot::SetCollectionMarker(3,8,0.8);
  TriggerPlot::SetCollectionMarker(4,2,0.5);

  Eloud=NULL;

}

////////////////////////////////////////////////////////////////////////////////////
Oplot::~Oplot(void){
////////////////////////////////////////////////////////////////////////////////////
  if(ReadTriggerSegments::Verbose>1) cout<<"Oplot::~Oplot"<<endl;
  if(Eloud!=NULL) delete Eloud;
}

////////////////////////////////////////////////////////////////////////////////////
void Oplot::SetTimeRange(const int aTimeMin, const int aTimeMax){
////////////////////////////////////////////////////////////////////////////////////
  if(aTimeMin>=aTimeMax){
    cerr<<"Oplot::SetTimeRange: the time range does not make sense"<<endl;
    return;
  }

  // time binning
  int nbins;
  if(aTimeMax-aTimeMin<=3600)        nbins = (aTimeMax-aTimeMin)/60+1;
  else if(aTimeMax-aTimeMin<=100000) nbins = (aTimeMax-aTimeMin)/1000+1;
  else                               nbins = (aTimeMax-aTimeMin)/3600+1;

  for(int s=0; s<5; s++){

    // new time range
    TriggerPlot::GetCollectionSelection(s)->SetTimeRange(aTimeMin,aTimeMax);

    // adjust time binning
    TriggerPlot::GetCollectionSelection(s)->SetTimeResolution(nbins);
  }
    
  return;
}

////////////////////////////////////////////////////////////////////////////////////
void Oplot::SetSNRThresholds(const double aSNR0, const double aSNR1, const double aSNR2, const double aSNR3){
////////////////////////////////////////////////////////////////////////////////////
  if(aSNR0>aSNR1){
    cerr<<"Oplot::SetSNRThresholds: the snr thresholds should be increasingly sorted"<<endl;
    return;
  }
  if(aSNR1>aSNR2){
    cerr<<"Oplot::SetSNRThresholds: the snr thresholds should be increasingly sorted"<<endl;
    return;
  }
  if(aSNR2>aSNR3){
    cerr<<"Oplot::SetSNRThresholds: the snr thresholds should be increasingly sorted"<<endl;
    return;
  }
  
  // SNR thresholds
  snrthr[0]=aSNR0;
  snrthr[1]=aSNR1;
  snrthr[2]=aSNR2;
  snrthr[3]=aSNR3;
  snrthr[4]=aSNR0;

  if(ReadTriggerSegments::Verbose>1){
    cout<<"Oplot::SetSNRThresholds: SNR thresholds:"<<endl;
    for(int s=0; s<5; s++) cout<<"              "<<s+1<<"/ SNR > "<<snrthr[s]<<endl;
  }

  // SNR max
  double smax;
  if(TriggerPlot::GetCollectionSelection(0)->GetSNRMax()<50.0) smax=50.0;
  else if(TriggerPlot::GetCollectionSelection(0)->GetSNRMax()<100.0) smax=100.0;
  else if(TriggerPlot::GetCollectionSelection(0)->GetSNRMax()<1000.0) smax=1000.0;
  else if(TriggerPlot::GetCollectionSelection(0)->GetSNRMax()<10000.0) smax=10000.0;
  else smax=100000.0;

  // at least one decade
  if(smax/snrthr[3]<10) smax=10*snrthr[3];

  // apply SNR selection
  for(int s=0; s<5; s++) TriggerPlot::GetCollectionSelection(s)->SetSNRRange(snrthr[s],smax);
  
  // Set plot legends
  ostringstream tmpstream;
  for(int s=0; s<5; s++){
    tmpstream<<"SNR > "<<fixed<<setprecision(2)<<snrthr[s];
    if(s==4) tmpstream<<", flagged";
    TriggerPlot::SetCollectionLegend(s,tmpstream.str());
    tmpstream.str(""); tmpstream.clear();
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////
double Oplot::GetSNRThreshold(const int aCollIndex){
////////////////////////////////////////////////////////////////////////////////////
  if(aCollIndex<0||aCollIndex>=Ncoll){
    cerr<<"Oplot::GetSNRThreshold: Collection #"<<aCollIndex<<" does not exist"<<endl;
    return 0.0;
  }
  return snrthr[aCollIndex];
}

////////////////////////////////////////////////////////////////////////////////////
void Oplot::PrintLoudestEventMap(const string aFileName){
////////////////////////////////////////////////////////////////////////////////////
  if(!Gfreqtimeloud->GetN()){
    cerr<<"Oplot::PrintLoudestEventMap: no plot has been built"<<endl;
    return;
  }

  double x,y;
  Gfreqtimeloud->GetPoint(0,x,y);

  if(Eloud!=NULL) delete Eloud;
  Eloud = new EventMap(1,GetTriggerFiles(x-3,x+3),"",0);
  Eloud->BuildMap(0,x);
  Eloud->PrintMap(0);
  
  if(aFileName.compare("")) Eloud->Print(aFileName);
    
  return;
}

