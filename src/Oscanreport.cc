//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

////////////////////////////////////////////////////////////////////////////////////
bool Omicron::ScanReport(const string aScanDir){
////////////////////////////////////////////////////////////////////////////////////

  if(!status_OK){
    cerr<<"Omicron::ScanReport: the Omicron object is corrupted"<<endl;
    return false;
  }

  // report directory: default or user-defined
  string reportdir = aScanDir;
  if(!reportdir.compare("")) reportdir=fScandir;

  // check that the report directory exists
  if(!IsDirectory(reportdir)){
    cerr<<"Omicron::ScanReport: missing scan directory "<<reportdir<<endl;
    return false;
  }
  
  // scan channel directories
  vector <string> chanfull, chan;
  vector <bool> chanOK;
  if(!ListDirectories(reportdir,chanfull)){
    cerr<<"Omicron::ScanReport: no channel directories in "<<reportdir<<endl;
    return false;
  }

  // remove incomplete directories
  for(int c=0; c<(int)chanfull.size(); c++){
    if(!IsBinaryFile(reportdir+"/"+chanfull[c]+"/"+chanfull[c]+"_mapsummary.root")){
      chanOK.push_back(false);
    }
    else{
      chan.push_back(chanfull[c]);
      chanOK.push_back(true);
    }
  }
  if(!chan.size()){
    cerr<<"Omicron::ScanReport: no scan to display"<<endl;
    return false;
  }

  // get map configuration (first channel)
  ostringstream tmpstream;
  TTree *MS;
  double map_gpsref, map_q, map_gpsmax, map_freqmax, map_snrmax, snrtest;
  int map_win, map_winref, wintest, map_nq, map_nwin, map_nformat, mtest, map_wstart, map_wstop;
  string *map_format=0, winset_full, formatset_full;
  vector<double> qset;
  vector<int> winset;
  vector<string> formatset;
  TFile sumfile((reportdir+"/"+chan[0]+"/"+chan[0]+"_mapsummary.root").c_str());
  if(sumfile.IsZombie()){
    cerr<<"Omicron::ScanReport: no channel summary in "<<reportdir<<"/"<<chan[0]<<endl;
    return false;
  }
  MS=(TTree*)sumfile.Get("mapsummary");
  if(!MS->GetEntries()){
    cerr<<"Omicron::ScanReport: no scan material "<<reportdir<<"/"<<chan[0]<<endl;
    delete MS;
    sumfile.Close();   
    return false;
  }
  MS->SetBranchAddress("gps_center",&map_gpsref);
  MS->SetBranchAddress("gps_whitening_start",&map_wstart);
  MS->SetBranchAddress("gps_whitening_stop",&map_wstop);
  MS->SetBranchAddress("nQ",&map_nq);
  MS->SetBranchAddress("nwindow",&map_nwin);
  MS->SetBranchAddress("window",&map_winref);
  MS->SetBranchAddress("nformat",&map_nformat);
  MS->GetEntry(0);
  MS->SetBranchAddress("window",&map_win);
  tmpstream<<"[";
  for(int m=0; m<map_nwin*map_nformat; m+=map_nformat){
    MS->GetEntry(m);
    if(m) tmpstream<<",";
    winset.push_back(map_win);
    tmpstream<<"'"<<map_win<<"'";
  }
  tmpstream<<"]";
  winset_full=tmpstream.str();
  tmpstream.str(""); tmpstream.clear();
  MS->SetBranchAddress("format",&map_format);
  tmpstream<<"[";
  for(int m=0; m<map_nformat; m++){
    MS->GetEntry(m);
    if(m) tmpstream<<",";
    formatset.push_back(*map_format);
    tmpstream<<"'"<<*map_format<<"'";
  }
  tmpstream<<"]";
  formatset_full=tmpstream.str();
  tmpstream.str(""); tmpstream.clear();

  delete MS;
  sumfile.Close();

  // timing
  tm utc;
  GPSToUTC (&utc, (int)map_gpsref);

  // import material
  system(("cp -f ${GWOLLUM_DOC}/style.css "+reportdir).c_str());
  system(("cp -f ${GWOLLUM_DOC}/Pics/gwollum_logo_min_trans.gif "+reportdir).c_str());
  system(("cp -f ${OMICRON_HTML}/pics/omicronlogo_xxl.gif "+reportdir).c_str());
  system(("cp -f ${OMICRON_HTML}/template/comparison_mode.html "+reportdir).c_str());

  // index header
  ofstream report((reportdir+"/index.html").c_str());
  report<<"<html>"<<endl;
  report<<"<head>"<<endl;
  report<<"<title>OmiScan of "<<setprecision(3)<<fixed<<map_gpsref<<"</title>"<<endl;
  report<<"<link rel=\"stylesheet\" href=\"style.css\" type=\"text/css\">"<<endl;
  report<<"<link rel=\"icon\" type=\"image/x-icon\" href=\"icon.gif\" />"<<endl;
  report<<"<link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"icon.gif\" />"<<endl;
  report<<"<script type=\"text/javascript\">"<<endl;
  report<<"function showImage(channel, type, timerange, format) {"<<endl;
  report<<"  for (var dt in timerange) {"<<endl;
  report<<"    for (var form in format) {"<<endl;
  report<<"      var basename ="<<endl;
  report<<"        channel + \"_\" + type + \"_dt\" + timerange[dt];"<<endl;
  report<<"      var basenameth ="<<endl;
  report<<"        channel + \"_\" + type + \"th_dt\" + timerange[dt];"<<endl;
  report<<"      document.getElementById(\"a_\" + channel + \"_dt\" + timerange[dt] + \"_\" + format[form]).href ="<<endl;
  report<<"        \"./\" + channel + \"/\" + basename + \".\" + format[form];"<<endl;
  report<<"      document.getElementById(\"img_\" + channel + \"_dt\" + timerange[dt] + \"_\" + format[form]).src ="<<endl;
  report<<"        \"./\" + channel + \"/\" + basenameth + \".\" + format[form];"<<endl;
  report<<"    }"<<endl;
  report<<"  }"<<endl;
  report<<"}"<<endl;
  report<<"</script>"<<endl;
  report<<"<script type=\"text/javascript\">"<<endl;
  report<<"function showSingleImage(channel, type, timerange, format) {"<<endl;
  report<<"  for (var form in format) {"<<endl;
  report<<"    var basename ="<<endl;
  report<<"      channel + \"_\" + type;"<<endl;
  report<<"    var basenameth ="<<endl;
  report<<"      channel + \"_\" + type + \"th\";"<<endl;
  report<<"      document.getElementById(\"a_\" + channel + \"_dt\" + timerange + \"_\" + format[form]).href ="<<endl;
  report<<"        \"./\" + channel + \"/\" + basename + \".\" + format[form];"<<endl;
  report<<"      document.getElementById(\"img_\" + channel + \"_dt\" + timerange + \"_\" + format[form]).src ="<<endl;
  report<<"        \"./\" + channel + \"/\" + basenameth + \".\" + format[form];"<<endl;
  report<<"  }"<<endl;
  report<<"}"<<endl;
  report<<"</script>"<<endl;
  report<<"</head>"<<endl;
  report<<"<body class=\"omiscan\">"<<endl;
  
  // index summary
  report<<"<table><tr><td>"<<endl;
  report<<"<h1>OmiScan of "<<setprecision(3)<<fixed<<map_gpsref<<"</h1>"<<endl;
  report<<"<a href=\"./comparison_mode.html\">Switch to comparison mode</a>"<<endl;
  report<<"<pre>"<<endl;
  report<<"Central GPS time:    "<<setprecision(3)<<fixed<<map_gpsref<<endl;
  report<<"Central Date (UTC):  "<<asctime(&utc);
  report<<"Scan run Date (UTC): "<<asctime(ptm);
  report<<"</pre>"<<endl;
  report<<"</td></tr></table>"<<endl;

  // channel index
  report<<"<h2>Scanned channels:</h2>"<<endl;
  report<<"<table class=\"omiscan\">"<<endl;
  for(int c=0; c<(int)chanfull.size(); c++){
    if(c%4==0) report<<"  <tr>";
    if(chanOK[c]) report<<"<td><a href=\"#"<<chanfull[c]<<"\">"<<chanfull[c]<<"</a></td>";
    else report<<"<td>"<<chanfull[c]<<"</td>";
    if((c+1)%4==0) report<<"  </tr>"<<endl;
  }
  for(int c=0; c<((int)chanfull.size())%4; c++){
    report<<"<td></td>";
    if(c==((int)chanfull.size())%4-1) report<<"</tr>"<<endl;
  }
  report<<"</table>"<<endl;
  
  // loop over channel directories
  for(int c=0; c<(int)chan.size(); c++){
    
    // get summary tree
    TFile sumfile((reportdir+"/"+chan[c]+"/"+chan[c]+"_mapsummary.root").c_str());
    if(sumfile.IsZombie()){
      cerr<<"Omicron::ScanReport: no channel summary in "<<reportdir<<"/"<<chan[c]<<endl;
      continue;
    }
    MS=(TTree*)sumfile.Get("mapsummary");
    MS->SetBranchAddress("Q",&map_q);
    MS->SetBranchAddress("window",&map_win);
    MS->SetBranchAddress("gps_loudest",&map_gpsmax);
    MS->SetBranchAddress("freq_loudest",&map_freqmax);
    MS->SetBranchAddress("snr_loudest",&map_snrmax);

    // get loudest tile over the first window
    for(int m=0; m<MS->GetEntries(); m++){
      MS->GetEntry(m);
      if(!m){
	wintest=map_win;
	snrtest=map_snrmax;
	mtest=m;
	qset.push_back(map_q);	
	continue;
      }
      if(map_win!=wintest) continue;
      if(map_snrmax>snrtest){
	snrtest=map_snrmax;
	mtest=m;
      }
      if(map_q!=qset.back()) qset.push_back(map_q);
    }	
    MS->GetEntry(mtest);

    // print loudest tile
    report<<"<table class=\"omiscan\">"<<endl;
    report<<"  <tr><td colspan=\""<<map_nwin*map_nformat<<"\"><a name=\""<<chan[c]<<"\"></a><h2>"<<chan[c]<<"</h2>"<<endl;
    report<<"      Loudest tile in &plusmn;"<<wintest/2.0<<"s: GPS="<<setprecision(3)<<fixed<<map_gpsmax<<", f="<<map_freqmax<<"Hz, Q="<<map_q<<", <b>SNR="<<map_snrmax<<"</b></td></tr>"<<endl;
    
    // print time series menu    
    report<<"  <tr><td colspan=\""<<map_nwin*map_nformat<<"\">Time series: <a href=\"javascript:showImage('"<<chan[c]<<"', 'ts', "<<winset_full<<", "<<formatset_full<<");\">raw</a></td></tr>"<<endl;

    // print map menu    
    report<<"  <tr><td colspan=\""<<map_nwin*map_nformat<<"\">Maps: ";
    for(int q=0; q<(int)qset.size(); q++)
      report<<"<a href=\"javascript:showImage('"<<chan[c]<<"', 'mapQ"<<q<<"', "<<winset_full<<", "<<formatset_full<<");\">Q="<<setprecision(3)<<fixed<<qset[q]<<"</a> / ";
    report<<"<a href=\"javascript:showImage('"<<chan[c]<<"', 'map', "<<winset_full<<", "<<formatset_full<<");\">full</a></td></tr>"<<endl;
  
    // print other menu    
    report<<"  <tr><td colspan=\""<<map_nwin*map_nformat<<"\">Misc: ";
    report<<"<a href=\"javascript:showSingleImage('"<<chan[c]<<"', 'asd_"<<map_wstart<<"_"<<map_wstop<<"', '"<<winset[0]<<"', "<<formatset_full<<");\">ASD</a> / ";
    report<<"<a href=\"javascript:showImage('"<<chan[c]<<"', 'projt', "<<winset_full<<", "<<formatset_full<<");\">SNR(time)</a> / "<<endl;
    report<<"<a href=\"javascript:showImage('"<<chan[c]<<"', 'projf', "<<winset_full<<", "<<formatset_full<<");\">SNR(frequency)</a> / "<<endl;
    if(IsBinaryFile(chan[c]+"/"+chan[c]+"_data.root")) report<<"<a href=\"./"<<chan[c]<<"/"<<chan[c]<<"_data.root\">Scan data</a>"<<endl;
    report<<"</td></tr>"<<endl;

    // print images
    report<<"  <tr>"<<endl;
    for(int w=0; w<map_nwin; w++)
      for(int f=0; f<map_nformat; f++)
	report<<"    <td><a id=\"a_"<<chan[c]<<"_dt"<<winset[w]<<"_"<<formatset[f]<<"\" href=\"./"<<chan[c]<<"/"<<chan[c]<<"_map_dt"<<winset[w]<<"."<<formatset[f]<<"\"><img id=\"img_"<<chan[c]<<"_dt"<<winset[w]<<"_"<<formatset[f]<<"\" src=\"./"<<chan[c]<<"/"<<chan[c]<<"_mapth_dt"<<winset[w]<<"."<<formatset[f]<<"\" alt=\""<<chan[c]<<"_map\" /></a></td>"<<endl;
    report<<"  </tr>"<<endl;

    qset.clear();
    delete MS;
    sumfile.Close();
    report<<"</table>"<<endl;
    report<<"<hr />"<<endl;
  }

  // index footer
  report<<"<table valign=\"bottom\" align=\"center\">"<<endl;
  report<<"<tr>"<<endl;
  report<<"<td valign=\"middle\" align=\"center\">"<<endl;
  report<<"<p>"<<endl;
  report<<"Florent Robinet<br><a href=\"mailto:robinet@lal.in2p3.fr\">Contact: robinet@lal.in2p3.fr</a>"<<endl;
  report<<"</p>"<<endl;
  report<<"</td>"<<endl;
  report<<"</tr>"<<endl;
  report<<"</table>"<<endl;
  report<<"</body>"<<endl;
  report<<"</html>"<<endl;
  report.close();
  


  return true;
}

