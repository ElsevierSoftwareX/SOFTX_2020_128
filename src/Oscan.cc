//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::Scan(const double aTimeCenter){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::Scan: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aTimeCenter<700000000){
    cerr<<"Omicron::Scan: the input time is not reasonable"<<endl;
    return false;
  }
  if(fChunkDuration!=fSegmentDuration){
    cerr<<"Omicron::Scan: this function can only be called if the chunk and segment durations are identical"<<endl;
    return false;
  }
  if(fVerbosity){
    time ( &timer );
    ptm = gmtime ( &timer );
    cout<<"#### Omicron::Scan timer = "<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) ---> +"<<timer-timer_start<<"s ####"<<endl;
  }

  // locals
  int dsize;          // native data size
  double *dvector;    // data vector before resampling
  int dsize_inj;      // native data size (inj)
  double *dvector_inj;// data vector before resampling (inj)
  int res;
  
  // Segment to process
  if(fVerbosity) cout<<"Omicron::Scan: initiate data segments..."<<endl;
  Segments *ScanSeg = new Segments((int)aTimeCenter-fChunkDuration/2,(int)aTimeCenter+fChunkDuration/2);

  // data structure
  if(!InitSegments(ScanSeg)||!NewChunk()){
    cerr<<"Omicron::Scan: cannot initiate data segments."<<endl;
    delete ScanSeg;
    return false;
  }

  // create map directories
  if(!MakeDirectories(aTimeCenter)){
    cerr<<"Omicron::Scan: the directory structure cannot be created"<<endl;
    delete ScanSeg;
    return false;
  }

  // loop over channels
  for(int c=0; c<(int)fChannels.size(); c++){
    cout<<"Omicron::Scan: **** channel "<<fChannels[c]<<"..."<<endl;
          
    // get data vector
    dvector=FFL->GetData(dsize, fChannels[c], dataseq->GetChunkTimeStart(), dataseq->GetChunkTimeEnd());

    if(dsize<=0){
      cor_chunk_ctr[c]++;
     cerr<<"Omicron::Scan: cannot retrieve data ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<"-"<<dataseq->GetChunkTimeEnd()<<")."<<endl;
      continue;
    }

    // get injection vector and inject it
    if(fInjChan.size()){
      dvector_inj=FFL->GetData(dsize_inj, fInjChan[c], dataseq->GetChunkTimeStart(), dataseq->GetChunkTimeEnd());
      if(dsize_inj<=0){
	cor_chunk_ctr[c]++;
	cerr<<"Omicron::Scan: cannot retrieve injection data ("<<fInjChan[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
	delete dvector;
	continue;
      }
      if(dsize_inj!=dsize){
	cor_chunk_ctr[c]++;
	cerr<<"Omicron::Scan: the sampling of the injection channel ("<<fInjChan[c]<<") is not the same as the sampling of the main channel --> skip chunk"<<endl;
	delete dvector;
	delete dvector_inj;
	continue;
      }
      for(int d=0; d<dsize; d++) dvector[d]+=(fInjFact[c]*dvector_inj[d]);
      delete dvector_inj;
    }

    // condition data vector
    res=ConditionVector(c, dsize, dvector);
    if(res<0){
      cerr<<"Omicron::Scan: conditioning failed ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
      delete dvector;
      delete ScanSeg;
      cor_data_ctr[c]++;
      return false;// fatal
    }
    if(res>0){
      cerr<<"Omicron::Scan: conditioning failed ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
      delete dvector;
      cor_data_ctr[c]++;
      continue;// skip channel
    }

    delete dvector;// not needed anymore

    // get maps
    if(!MakeMaps(c, aTimeCenter)){
      cerr<<"Omicron::Scan: cannot make maps ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
      cor_data_ctr[c]++;
      continue;
    }

    // apply SNR threshold (except the first channel)
    if(c&&Qmap_full[0]->GetMaximum()<fSNRThreshold){
      cout<<"Omicron::Scan: below SNR threshold "<<Qmap_full[0]->GetMaximum()<<"<"<<fSNRThreshold<<" ("<<fChannels[c]<<" "<<dataseq->GetChunkTimeStart()<<" "<<dataseq->GetChunkTimeEnd()<<")."<<endl;
      max_chunk_ctr[c]++;// use this specific monitor here
      continue;
    }

    // save maps on disk
    if(!WriteMaps(c)){
      cerr<<"Omicron::Scan: cannot write maps for channel "<<fChannels[c]<<endl;
      cor_data_ctr[c]++;
      continue;
    }

    // save ASD on disk
    SaveAPSD(c,"ASD");
    
    // save data on disk
    SaveTS(c, aTimeCenter);

    // time for this channel
    if(fVerbosity>1){ 
      time ( &timer );
      ptm = gmtime ( &timer );
      cout<<"#### Omicron::Scan timer = "<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) ---> +"<<timer-timer_start<<"s ####"<<endl;
    }
  }
        
  delete ScanSeg;

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::ScanTriggers(const double aTimeCenter){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::ScanTriggers: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aTimeCenter<700000000){
    cerr<<"Omicron::ScanTriggers: the input time is not reasonable"<<endl;
    return false;
  }
  if(!IsDirectory(fTrigDir)){
    cerr<<"Omicron::ScanTriggers: the trigger directory is required"<<endl;
    return false;
  }
  if(fVerbosity){
    time ( &timer );
    ptm = gmtime ( &timer );
    cout<<"#### Omicron::ScanTriggers timer = "<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) ---> +"<<timer-timer_start<<"s ####"<<endl;
  }

  // Segment to process
  if(fVerbosity) cout<<"Omicron::ScanTriggers: initiate data segments..."<<endl;
  int start= (int)aTimeCenter-fWindowMax/2-1;
  int stop= (int)aTimeCenter+fWindowMax/2+1;
  Segments *ScanSeg = new Segments(start,stop);

  // data structure
  if(!InitSegments(ScanSeg)){
    cerr<<"Omicron::ScanTriggers: cannot initiate data segments."<<endl;
    delete ScanSeg;
    return false;
  }

  // create map directories
  if(!MakeDirectories(aTimeCenter)){
    cerr<<"Omicron::ScanTriggers: the directory structure cannot be created"<<endl;
    delete ScanSeg;
    return false;
  }

  // locals
  EventMap *E;
  ReadTriggerMetaData *MD;
  string trigfiles;
  vector <double> qs;
  TH2D **qmap;
  double qmin, qmax, mmm;

  // loop over channels
  for(int c=0; c<(int)fChannels.size(); c++){
    cout<<"Omicron::ScanTriggers: **** channel "<<fChannels[c]<<"..."<<endl;

    // check trigger files
    MD = new ReadTriggerMetaData(fTrigDir+"/"+fChannels[c]+"/"+fChannels[c]+"_*.root","",0);
    if(!MD->GetNFiles()){
      cerr<<"Omicron::ScanTriggers: no trigger files --> skip channel"<<endl;
      delete MD;
      continue;
    }
  
    // select files of interest
    trigfiles = MD->GetTriggerFiles(start,stop);
    if(!trigfiles.compare("none")){
      cerr<<"Omicron::ScanTriggers: no trigger files for this time window --> skip channel"<<endl;
      delete MD;
      continue;
    }

    // get Qmin/max
    qmin=MD->GetMeta("omicron_PARAMETER_QMIN",(double)start);
    qmax=MD->GetMeta("omicron_PARAMETER_QMAX",(double)start);
    mmm=MD->GetMeta("omicron_PARAMETER_MISMATCHMAX",(double)start);
    if(qmin==-1.0e-20||qmax==-1.0e-20||mmm==-1.0e-20){
      cerr<<"Omicron::ScanTriggers: cannot retieve the Q parameters --> skip channel"<<endl;
      delete MD;
      continue;
    }
    delete MD;
    
    // get list of Qs
    qs=Otile::ComputeQs(qmin,qmax,mmm);
    qmap = new TH2D * [(int)qs.size()];
    
    // Map object
    E = new EventMap(trigfiles,"",2);
    E->SetMapTimeRange(0,stop-start);
    
    // get Q maps
    for(int q=0; q<(int)qs.size(); q++){
      E->SetMapQRange(0,0.99*qs[q],1.01*qs[q]);
      E->BuildMap(0,aTimeCenter);
      qmap[q] = E->GetMap(0);
    }
    
    // create overall maps
    MakeFullMaps(c,aTimeCenter,(int)qs.size(),qmap);

    // apply SNR threshold (except the first channel)

    // save maps on disk
    if(!WriteMaps(c,aTimeCenter,qs,qmap)){
      cerr<<"Omicron::Scan: cannot write maps for channel "<<fChannels[c]<<endl;
      cor_data_ctr[c]++;
      continue;
    }

    // time for this channel
    if(fVerbosity>1){ 
      time ( &timer );
      ptm = gmtime ( &timer );
      cout<<"#### Omicron::ScanTriggers timer = "<<setfill('0')<<setw(2)<<ptm->tm_hour<<":"<<setfill('0')<<setw(2)<<ptm->tm_min<<":"<<setfill('0')<<setw(2)<<ptm->tm_sec<<" (UTC) ---> +"<<timer-timer_start<<"s ####"<<endl;
    }

    // cleaning for next channel
    for(int q=0; q<(int)qs.size(); q++) delete qmap[q];
    delete qmap;
    delete E;
    qs.clear();
  }
        
  delete ScanSeg;

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::MakeMaps(const int aChNumber, const double aTimeCenter){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::MakeMaps: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::MakeMaps: channel number "<<aChNumber<<" does not exist"<<endl;
    return false;
  }
  if(aTimeCenter-fWindowMax<dataseq->GetSegmentTimeStart(0)+fOverlapDuration/2){
    cerr<<"Omicron::MakeMaps: the requested central time "<<aTimeCenter<<" is not compatible with the data segmentation and the windowing -> you need to increase your segment"<<endl;
    return false;
  }
  if(aTimeCenter+fWindowMax>dataseq->GetSegmentTimeEnd(0)-fOverlapDuration/2){
    cerr<<"Omicron::MakeMaps: the requested central time "<<aTimeCenter<<" is not compatible with the data segmentation and the windowing -> you need to increase your segment"<<endl;
    return false;
  }

  if(fVerbosity) cout<<"Omicron::MakeMaps: make maps for channel "<<fChannels[aChNumber]<<" starting at "<<dataseq->GetChunkTimeStart()<<endl;

  Qmap_center=aTimeCenter;
  double toffset=Qmap_center-dataseq->GetSegmentTimeStart(0)-(double)fSegmentDuration/2.0;
  ostringstream tmpstream;

  // get Qmaps
  if(fVerbosity>1) cout<<"Omicron::MakeMaps: make Qmaps..."<<endl;
  for(int q=0; q<tile->GetNQPlanes(); q++){
    Qmap[q] = tile->GetMap(q, dataRe[0], dataIm[0], -toffset);
    if(Qmap[q]==NULL){
      cerr<<"Omicron::MakeMaps: maps for channel "<<fChannels[aChNumber]<<" are corrupted"<<endl;
      delete dataRe[0];
      delete dataIm[0];
      return false;
    }

    // map title
    tmpstream<<fChannels[aChNumber]<<": GPS="<<fixed<<setprecision(3)<<Qmap_center<<", Q="<<tile->GetQ(q);
    Qmap[q]->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
    
    // cosmetics
    Qmap[q]->GetXaxis()->SetTitleOffset(1.1);
    Qmap[q]->GetXaxis()->SetLabelSize(0.045);
    Qmap[q]->GetYaxis()->SetLabelSize(0.045);
    Qmap[q]->GetXaxis()->SetTitleSize(0.045);
    Qmap[q]->GetYaxis()->SetTitleSize(0.045);
    Qmap[q]->GetZaxis()->SetTitleSize(0.05);
  }
  delete dataRe[0];
  delete dataIm[0];

  MapChNumber=aChNumber;

  // create overall maps
  MakeFullMaps(aChNumber,aTimeCenter,tile->GetNQPlanes(),Qmap);

  // Q-map containing the loudest tile
  for(int w=0; w<(int)fWindows.size(); w++){
    loudest_qmap[w]=0;
    for(int q=1; q<tile->GetNQPlanes(); q++){
      Qmap[q]->GetXaxis()->SetRangeUser(-(double)fWindows[0]/2.0,(double)fWindows[0]/2.0);
      if(Qmap[q]->GetMaximum()>Qmap[loudest_qmap[w]]->GetMaximum()) loudest_qmap[w]=q;
      Qmap[q]->GetXaxis()->SetRange();
    }
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////////
void Omicron::MakeFullMaps(const int aChNumber, const double aTimeCenter, const int aNq, TH2D **aQmap){
////////////////////////////////////////////////////////////////////////////////////

  // create overall maps
  if(fVerbosity>1) cout<<"Omicron::MakeFullMaps: make overall maps..."<<endl;
  ostringstream tmpstream;
  for(int w=0; w<(int)fWindows.size(); w++){
    tmpstream<<fChannels[aChNumber]<<": GPS="<<fixed<<setprecision(3)<<aTimeCenter;
    Qmap_full[w]->Reset();
    Qmap_full[w]->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
  }
      
  //populate overall maps
  int bin_start, bin_stop, dummy, bin_start_t, bin_stop_t, bin_stop_f, bin_start_f;
  double content;
  for(int q=0; q<aNq; q++){
    for(int bt=1; bt<=aQmap[q]->GetNbinsX(); bt++){
      for(int bf=1; bf<=aQmap[q]->GetNbinsY(); bf++){
      	content=aQmap[q]->GetBinContent(bt,bf);
	for(int w=0; w<(int)fWindows.size(); w++){
	  if(aQmap[q]->GetXaxis()->GetBinUpEdge(bt)<-(double)fWindows[w]/2.0) continue;
	  if(aQmap[q]->GetXaxis()->GetBinLowEdge(bt)>(double)fWindows[w]/2.0) continue;
	  if(aQmap[q]->GetYaxis()->GetBinUpEdge(bf)<Qmap_full[w]->GetYaxis()->GetBinLowEdge(1)) continue;
	  if(aQmap[q]->GetYaxis()->GetBinLowEdge(bf)>Qmap_full[w]->GetYaxis()->GetBinUpEdge(Qmap_full[w]->GetNbinsY())) continue;
	  bin_start = Qmap_full[w]->FindBin(aQmap[q]->GetXaxis()->GetBinLowEdge(bt),aQmap[q]->GetYaxis()->GetBinLowEdge(bf));
	  bin_stop = Qmap_full[w]->FindBin(aQmap[q]->GetXaxis()->GetBinUpEdge(bt),aQmap[q]->GetYaxis()->GetBinUpEdge(bf));
	  Qmap_full[w]->GetBinXYZ(bin_start,bin_start_t,bin_start_f,dummy);
	  Qmap_full[w]->GetBinXYZ(bin_stop,bin_stop_t,bin_stop_f,dummy);
	  for(int bbt=bin_start_t; bbt<=bin_stop_t; bbt++){// time-sweep the tile
	    for(int bbf=bin_start_f; bbf<=bin_stop_f; bbf++){// freq-sweep the tile
	      if(content>Qmap_full[w]->GetBinContent(bbt,bbf)) Qmap_full[w]->SetBinContent(bbt,bbf,content);
	    }
	  }
	}
      }
    }
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::WriteMaps(const int aChNumber, const double aTimeCenter, const vector <double>& aQ, TH2D **aQmap){
////////////////////////////////////////////////////////////////////////////////////
  if(!status_OK){
    cerr<<"Omicron::WriteMaps: the Omicron object is corrupted"<<endl;
    return false;
  }
  if(aChNumber<0||aChNumber>=(int)fChannels.size()){
    cerr<<"Omicron::WriteMaps: channel number "<<aChNumber<<" does not exist"<<endl;
    return false;
  }
  if(!aQ.size()&&MapChNumber!=aChNumber){
    cerr<<"Omicron::WriteMaps: the map in memory is not associated to channel "<<fChannels[aChNumber]<<endl;
    return false;
  }
  if(aQ.size()>0&&aQmap==NULL){
    cerr<<"Omicron::WriteMaps: the input maps are NULL for channel "<<fChannels[aChNumber]<<endl;
    return false;
  }
  if(aQ.size()>0){
    for(int q=0; q<(int)aQ.size(); q++){
      if(aQmap[q]==NULL){
	cerr<<"Omicron::WriteMaps: the input maps are NULL for channel "<<fChannels[aChNumber]<<endl;
	return false;
      }
    }
  }

  if(fVerbosity) cout<<"Omicron::WriteMaps: write maps for channel "<<fChannels[aChNumber]<<endl;

  // number of qmaps
  int nqs;
  if(aQ.size()>0) nqs=(int)aQ.size();
  else nqs=tile->GetNQPlanes();

  // pointers to Qmaps
  TH2D **map_tmp = new TH2D* [nqs];
  if(aQ.size()>0) for(int q=0; q<nqs; q++) map_tmp[q]=aQmap[q];
  else            for(int q=0; q<nqs; q++) map_tmp[q]=Qmap[q];
  
  // ROOT
  if(fOutFormat.find("root")!=string::npos){
    TFile *fmap;
    if(first_save[aChNumber]){
      fmap=new TFile((fOutdir[aChNumber]+"/"+fChannels[aChNumber]+"_data.root").c_str(),"RECREATE");
      first_save[aChNumber]=false;
    }
    else fmap=new TFile((fOutdir[aChNumber]+"/"+fChannels[aChNumber]+"_data.root").c_str(),"UPDATE");
    fmap->cd();
    for(int q=0; q<nqs; q++) map_tmp[q]->Write();
    for(int w=0; w<(int)fWindows.size(); w++) Qmap_full[w]->Write();
    fmap->Close();
  }

  // Graphix
  vector <string> form;
  if(fOutFormat.find("gif")!=string::npos) form.push_back("gif");
  if(fOutFormat.find("png")!=string::npos) form.push_back("png");
  if(fOutFormat.find("pdf")!=string::npos) form.push_back("pdf");
  if(fOutFormat.find("ps")!=string::npos)  form.push_back("ps");
  if(fOutFormat.find("xml")!=string::npos) form.push_back("xml");
  if(fOutFormat.find("eps")!=string::npos) form.push_back("eps"); 
  if(fOutFormat.find("jpg")!=string::npos) form.push_back("jpg"); 
  if(fOutFormat.find("svg")!=string::npos) form.push_back("svg"); 
  if(!form.size()) return true;

  // summary tree
  TTree *MS = new TTree("mapsummary", "map summary");
  double map_gps, map_q, map_gpsmax, map_freqmax, map_snrmax;
  int map_win, map_nq, map_nwin, map_nformat, map_wstart, map_wstop;
  string map_format;
  MS->Branch("gps_center",&map_gps,"gps_center/D");
  MS->Branch("gps_whitening_start",&map_wstart,"gps_whitening_start/I");
  MS->Branch("gps_whitening_stop",&map_wstop,"gps_whitening_stop/I");
  MS->Branch("nQ",&map_nq,"nQ/I");
  MS->Branch("nwindow",&map_nwin,"nwindow/I");
  MS->Branch("Q",&map_q,"Q/D");
  MS->Branch("window",&map_win,"window/I");
  MS->Branch("gps_loudest",&map_gpsmax,"gps_loudest/D");
  MS->Branch("freq_loudest",&map_freqmax,"freq_loudest/D");
  MS->Branch("snr_loudest",&map_snrmax,"snr_loudest/D");
  MS->Branch("nformat",&map_nformat,"nformat/I");
  MS->Branch("format",&map_format);
  if(aTimeCenter>0) map_gps=aTimeCenter;
  else map_gps=Qmap_center;
  map_nwin=(int)fWindows.size();
  map_nq=nqs;
  map_nformat=(int)form.size();
  map_wstart=dataseq->GetChunkTimeStart();
  map_wstop=dataseq->GetChunkTimeEnd();

  // canvas style
  GPlot->SetLogx(0);
  GPlot->SetLogy(1);
  GPlot->SetLogz(1);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);

  // locals
  double center;
  ostringstream tmpstream;
  int xmax, ymax, zmax;

  // draw Qmaps
  for(int q=0; q<nqs; q++){
    if(aQ.size()>0) map_q=aQ[q];
    else map_q=tile->GetQ(q);

    // plot
    map_tmp[q]->GetZaxis()->SetRangeUser(1,50);
    GPlot->Draw(map_tmp[q],"COLZ");
    
    // window resize
    center=(map_tmp[q]->GetXaxis()->GetXmax()+map_tmp[q]->GetXaxis()->GetXmin())/2.0;
    for(int w=0; w<(int)fWindows.size(); w++){
      map_win=fWindows[w];

      // zoom
      map_tmp[q]->GetXaxis()->SetRangeUser(center-(double)fWindows[w]/2.0,center+(double)fWindows[w]/2.0);

      // get max bin
      map_tmp[q]->GetMaximumBin(xmax, ymax, zmax);
         
      // loudest tile
      tmpstream<<"Loudest tile: GPS="<<fixed<<setprecision(3)<<Qmap_center+map_tmp[q]->GetXaxis()->GetBinCenter(xmax)<<", f="<<map_tmp[q]->GetYaxis()->GetBinCenter(ymax)<<" Hz, SNR="<<map_tmp[q]->GetBinContent(xmax,ymax);
      GPlot->AddText(tmpstream.str(), 0.01,0.01,0.03);
      tmpstream.str(""); tmpstream.clear();
      map_gpsmax=Qmap_center+map_tmp[q]->GetXaxis()->GetBinCenter(xmax);
      map_freqmax=map_tmp[q]->GetYaxis()->GetBinCenter(ymax);
      map_snrmax=map_tmp[q]->GetBinContent(xmax,ymax);

      // save qmaps
      for(int f=0; f<(int)form.size(); f++){
	tmpstream<<fOutdir[aChNumber]<<"/"<<fChannels[aChNumber]<<"_mapQ"<<q<<"_dt"<<fWindows[w]<<"."<<form[f];
	GPlot->Print(tmpstream.str());
	tmpstream.str(""); tmpstream.clear();
 	tmpstream<<fOutdir[aChNumber]<<"/"<<fChannels[aChNumber]<<"_mapQ"<<q<<"th_dt"<<fWindows[w]<<"."<<form[f];
	GPlot->SetGridx(0); GPlot->SetGridy(0);
	GPlot->Print(tmpstream.str(),0.5);
	GPlot->SetGridx(1); GPlot->SetGridy(1);
 	tmpstream.str(""); tmpstream.clear();
	map_format=form[f];
	
	// fill summary tree
	MS->Fill();
      }

    }
    map_tmp[q]->GetXaxis()->UnZoom();
    map_tmp[q]->GetZaxis()->UnZoom();
  }
  delete map_tmp;


  // draw overall maps
  for(int w=0; w<(int)fWindows.size(); w++){
    Qmap_full[w]->GetZaxis()->SetRangeUser(1,50);
    Qmap_full[w]->GetXaxis()->SetRangeUser(center-(double)fWindows[w]/2.0,center+(double)fWindows[w]/2.0);
    GPlot->Draw(Qmap_full[w],"COLZ");

    // get max bin
    Qmap_full[w]->GetMaximumBin(xmax, ymax, zmax);
         
    // loudest tile
    tmpstream<<"Loudest tile: GPS="<<fixed<<setprecision(3)<<Qmap_center+Qmap_full[w]->GetXaxis()->GetBinCenter(xmax)<<", f="<<Qmap_full[w]->GetYaxis()->GetBinCenter(ymax)<<" Hz, SNR="<<Qmap_full[w]->GetBinContent(xmax,ymax);
    GPlot->AddText(tmpstream.str(), 0.01,0.01,0.03);
    tmpstream.str(""); tmpstream.clear();

    // save maps
    for(int f=0; f<(int)form.size(); f++){
      tmpstream<<fOutdir[aChNumber]<<"/"<<fChannels[aChNumber]<<"_map_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();
      tmpstream<<fOutdir[aChNumber]<<"/"<<fChannels[aChNumber]<<"_mapth_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->SetGridx(0); GPlot->SetGridy(0);
      GPlot->Print(tmpstream.str(),0.5);
      GPlot->SetGridx(1); GPlot->SetGridy(1);
      tmpstream.str(""); tmpstream.clear();
    }
    Qmap_full[w]->GetXaxis()->UnZoom();
    Qmap_full[w]->GetZaxis()->UnZoom();
  }

  // save projections
  GPlot->SetLogx(0);
  GPlot->SetLogy(1);
  GPlot->SetGridx(1);
  GPlot->SetGridy(1);
  TH1D* proj;
  for(int w=0; w<(int)fWindows.size(); w++){
    Qmap_full[w]->GetMaximumBin(xmax, ymax, zmax);
    proj=Qmap_full[w]->ProjectionX("_pfx",ymax,ymax);
    tmpstream<<"SNR vs time at f = "<<fixed<<setprecision(3)<<Qmap_full[w]->GetYaxis()->GetBinCenter(ymax)<<"Hz";
    proj->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
    proj->GetXaxis()->SetTitle("Time [s]");
    proj->GetYaxis()->SetTitle("SNR");
    proj->GetXaxis()->SetTitleOffset(1.1);
    proj->GetXaxis()->SetLabelSize(0.045);
    proj->GetYaxis()->SetLabelSize(0.045);
    proj->GetXaxis()->SetTitleSize(0.045);
    proj->GetYaxis()->SetTitleSize(0.045);
    GPlot->SetLogx(0);
    GPlot->Draw(proj);
    for(int f=0; f<(int)form.size(); f++){
      tmpstream<<fOutdir[aChNumber]<<"/"<<fChannels[aChNumber]<<"_projt_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();
      tmpstream<<fOutdir[aChNumber]<<"/"<<fChannels[aChNumber]<<"_projtth_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->Print(tmpstream.str(),0.5);
      tmpstream.str(""); tmpstream.clear();
    }
    delete proj;

    proj=Qmap_full[w]->ProjectionY("_pfy",xmax,xmax);
    tmpstream<<"SNR vs frequency at GPS = "<<fixed<<setprecision(3)<<Qmap_full[w]->GetXaxis()->GetBinCenter(xmax)+Qmap_center;
    proj->SetTitle(tmpstream.str().c_str());
    tmpstream.str(""); tmpstream.clear();
    proj->GetXaxis()->SetTitle("Frequency [Hz]");
    proj->GetYaxis()->SetTitle("SNR");
    proj->GetXaxis()->SetTitleOffset(1.1);
    proj->GetXaxis()->SetLabelSize(0.045);
    proj->GetYaxis()->SetLabelSize(0.045);
    proj->GetXaxis()->SetTitleSize(0.045);
    proj->GetYaxis()->SetTitleSize(0.045);
    GPlot->SetLogx(1);
    GPlot->Draw(proj);
    for(int f=0; f<(int)form.size(); f++){
      tmpstream<<fOutdir[aChNumber]<<"/"<<fChannels[aChNumber]<<"_projf_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->Print(tmpstream.str());
      tmpstream.str(""); tmpstream.clear();
      tmpstream<<fOutdir[aChNumber]<<"/"<<fChannels[aChNumber]<<"_projfth_dt"<<fWindows[w]<<"."<<form[f];
      GPlot->Print(tmpstream.str(),0.5);
      tmpstream.str(""); tmpstream.clear();
    }
    delete proj;
  }

  // save map summary
  TFile summary((fOutdir[aChNumber]+"/"+fChannels[aChNumber]+"_mapsummary.root").c_str(),"RECREATE");
  summary.cd();
  MS->Write();
  summary.Close();
  delete MS;

  return true;
}
