//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Oplot.h"

ClassImp(Oplot)

////////////////////////////////////////////////////////////////////////////////////
Oplot::Oplot(const string aPattern, const string aDirectory, const int aVerbose) : TriggerPlot (4, aPattern, aDirectory, aVerbose){ 
////////////////////////////////////////////////////////////////////////////////////

  // SNR thresholds
  if(Msnrmin_stat<10000) snrthr[0]=Msnrmin_stat;
  else                   snrthr[0]=5.0;
  if(snrthr[0]<5)        snrthr[1]=5.0;
  else if(snrthr[0]<8)   snrthr[1]=8.0;
  else if(snrthr[0]<10)  snrthr[1]=10.0;
  else if(snrthr[0]<20)  snrthr[1]=20.0;
  else                   snrthr[1]=2.0*snrthr[0];
  if(snrthr[1]<8)        snrthr[2]=8.0;
  else if(snrthr[1]<10)  snrthr[2]=10.0;
  else if(snrthr[1]<20)  snrthr[2]=20.0;
  else                   snrthr[2]=2.0*snrthr[1];
  if(snrthr[2]<10)       snrthr[3]=10.0;
  else if(snrthr[2]<20)  snrthr[3]=20.0;
  else                   snrthr[3]=2.0*snrthr[2];
  if(ReadTriggerSegments::Verbose>1){
    cout<<"Oplot::Oplot: 4 SNR thresholds:"<<endl;
    for(int s=0; s<4; s++) cout<<"              "<<s+1<<"/ SNR > "<<snrthr[s]<<endl;
  }

  // SNR max
  double smax;
  if(TriggerPlot::GetCollectionSelection(0)->GetSNRMax()<50.0) smax=50.0;
  else if(TriggerPlot::GetCollectionSelection(0)->GetSNRMax()<100.0) smax=100.0;
  else if(TriggerPlot::GetCollectionSelection(0)->GetSNRMax()<1000.0) smax=1000.0;
  else if(TriggerPlot::GetCollectionSelection(0)->GetSNRMax()<10000.0) smax=10000.0;
  else smax=100000.0;
  if(smax<=TriggerPlot::GetCollectionSelection(0)->GetSNRMin()) smax=2*TriggerPlot::GetCollectionSelection(0)->GetSNRMin();

  // apply SNR selection
  for(int s=0; s<4; s++) TriggerPlot::GetCollectionSelection(s)->SetSNRRange(snrthr[s],smax);
  
  // Clusterize
  TriggerPlot::Clusterize("TIMEFREQ");
  for(int s=0; s<4; s++) TriggerPlot::SetCollectionUseClusters(s,true);

  // Set plot legends
  ostringstream tmpstream;
  for(int s=0; s<4; s++){
    tmpstream<<"SNR > "<<fixed<<setprecision(2)<<snrthr[s];
    TriggerPlot::SetCollectionLegend(s,tmpstream.str());
    tmpstream.str(""); tmpstream.clear();
  }

  // set plot style
  TriggerPlot::SetDateFormat(true);
  TriggerPlot::SetCollectionColor(0,9);
  TriggerPlot::SetCollectionColor(1,4);
  TriggerPlot::SetCollectionColor(2,1);
  TriggerPlot::SetCollectionColor(3,2);
  TriggerPlot::SetCollectionMarker(0,1);
  TriggerPlot::SetCollectionMarker(1,7);
  TriggerPlot::SetCollectionMarker(2,4,0.8);
  TriggerPlot::SetCollectionMarker(3,8,0.8);

}

////////////////////////////////////////////////////////////////////////////////////
Oplot::~Oplot(void){
////////////////////////////////////////////////////////////////////////////////////
  if(ReadTriggerSegments::Verbose>1) cout<<"Oplot::~Oplot"<<endl;
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

  for(int s=0; s<4; s++){

    // new time range
    TriggerPlot::GetCollectionSelection(s)->SetTimeRange(aTimeMin,aTimeMax);

    // adjust time binning
    TriggerPlot::GetCollectionSelection(s)->SetTimeResolution(nbins);
  }
    
  return;
}

