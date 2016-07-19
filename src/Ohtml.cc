//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Omicron.h"

////////////////////////////////////////////////////////////////////////////////////
void Omicron::MakeHtml(void){
////////////////////////////////////////////////////////////////////////////////////

  // import material
  system(("cp -f ${OMICRON_HTML}/style/style."+GPlot->GetCurrentStyle()+".css "+maindir+"/style.css").c_str());
  system(("cp -f ${GWOLLUM_DOC}/Pics/gwollum_logo_min_trans.gif "+maindir+"/icon.gif").c_str());
  system(("cp -f ${OMICRON_HTML}/pics/led-*.gif "+maindir).c_str());
  if(!fNoLogo) system(("cp -f ${OMICRON_HTML}/pics/omicronlogo."+GPlot->GetCurrentStyle()+".gif "+maindir+"/logo.gif").c_str());
  else system(("rm -f "+maindir+"/logo.gif").c_str());
  system(("cp -f "+fOptionFile+" "+maindir+"/omicron.parameters.txt").c_str());

  // select format
  string form;
  if(fOutFormat.find("png")!=string::npos) form="png";
  else if(fOutFormat.find("gif")!=string::npos) form="gif";
  else if(fOutFormat.find("jpg")!=string::npos) form="jpg";
  else form="";

  // window set
  ostringstream tmpstream;
  tmpstream<<"[";
  for(int w=0; w<(int)fWindows.size(); w++){
    if(w) tmpstream<<",";
    tmpstream<<"'"<<fWindows[w]<<"'";
  }
  tmpstream<<"]";
  string windowset=tmpstream.str();
  tmpstream.clear(); tmpstream.str("");

  // total processing time
  time ( &timer );
  
  // index header & scripts
  ofstream report((maindir+"/index.html").c_str());
  report<<"<html>"<<endl;
  report<<"<head>"<<endl;
  report<<"<title>Omicron Report "<<(int)inSegments->GetStart(0)<<"-"<<(int)inSegments->GetEnd(inSegments->GetNsegments()-1)<<"</title>"<<endl;
  report<<"<link rel=\"stylesheet\" href=\"style.css\" type=\"text/css\" />"<<endl;
  report<<"<link rel=\"icon\" type=\"image/x-icon\" href=\"./icon.gif\" />"<<endl;
  report<<"<link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"./icon.gif\" />"<<endl;
  report<<"<script type=\"text/javascript\">"<<endl;
  report<<"function showImage(channel, channelconv, type, timerange, format) {"<<endl;
  report<<"  for (var dt in timerange) {"<<endl;
  report<<"    var basename ="<<endl;
  report<<"      channelconv + \"_\" + type + \"-\" + timerange[dt];"<<endl;
  report<<"    document.getElementById(\"a_\" + channelconv + \"-\" + timerange[dt]).href ="<<endl;
  report<<"      \"./\" + channel + \"/\" + basename + \".\" + format;"<<endl;
  report<<"    document.getElementById(\"img_\" + channelconv + \"-\" + timerange[dt]).src ="<<endl;
  report<<"      \"./\" + channel + \"/th\" + basename + \".\" + format;"<<endl;
  report<<"  }"<<endl;
  report<<"}"<<endl;
  report<<"</script>"<<endl;
  report<<"<script type=\"text/javascript\">"<<endl;
  report<<"function toggle(anId) {"<<endl;
  report<<"  node = document.getElementById(anId);"<<endl;
  report<<"  if (node.style.visibility==\"hidden\") {"<<endl;
  report<<"    node.style.visibility = \"visible\";"<<endl;
  report<<"    node.style.height = \"auto\";"<<endl;
  report<<"  }"<<endl;
  report<<"  else {"<<endl;
  report<<"    node.style.visibility = \"hidden\";"<<endl;
  report<<"    node.style.height = \"0\";"<<endl;
  report<<"  }"<<endl;
  report<<"}"<<endl;
  report<<"</script>"<<endl;
  report<<"</head>"<<endl;
  report<<"<body>"<<endl;
  report<<endl;

  // index title
  report<<"<h1>Omicron Report </h1>"<<endl;
  report<<"<hr />"<<endl;
  report<<endl;

  // index summary
  report<<"<h2>Summary</h2>"<<endl;
  tm utc; int gps;
  report<<"<table>"<<endl;
  gps=(int)inSegments->GetStart(0);
  GPSToUTC (&utc, gps);
  report<<"  <tr><td>Omicron version:</td><td>"<<GetVersion()<<"</td></tr>"<<endl;
  report<<"  <tr><td>Omicron run by:</td><td>"<<getenv("USER")<<"</td></tr>"<<endl;
  report<<"  <tr><td>Omicron processing time:</td><td>"<<(int)(timer-timer_start)/3600<<"h, "<<((int)(timer-timer_start)%3600)/60<<"min</td></tr>"<<endl;
  report<<"  <tr><td>Processing Date:</td><td>"<<asctime(ptm)<<" (UTC)</td></tr>"<<endl;
  report<<"  <tr><td>Requested start:</td><td>"<<gps<<" &rarr; "<<asctime(&utc)<<" (UTC)</td></tr>"<<endl;
  gps=(int)inSegments->GetEnd(inSegments->GetNsegments()-1);
  GPSToUTC (&utc, gps);
  report<<"  <tr><td>Requested stop:</td><td>"<<gps<<" &rarr; "<<asctime(&utc)<<" (UTC)</td></tr>"<<endl;
  report<<"  <tr><td>Requested livetime:</td><td>"<<(int)inSegments->GetLiveTime()<<" sec &rarr; "<<setprecision(3)<<fixed<<inSegments->GetLiveTime()/3600.0/24<<" days</td></tr>"<<endl;
  inSegments->Write(maindir+"/omicron.segments.txt");
  report<<"  <tr><td>Requested segments:</td><td><a href=\"./omicron.segments.txt\">omicron.segments.txt</a></td></tr>"<<endl;
  report<<"  <tr><td>Number of chunks:</td><td>"<<chunk_ctr<<"</td></tr>"<<endl;
  report<<"  <tr><td>Configuration:</td><td><a href=\"./omicron.parameters.txt\">omicron.parameters.txt</a></td></tr>"<<endl;
  report<<"</table>"<<endl;
  report<<"<hr />"<<endl;
  report<<endl;

  // search parameters
  report<<"<h2>Parameters</h2>"<<endl;
  report<<"<table>"<<endl;
  report<<"  <tr><td>Timing:</td><td> chunks of "<<tile->GetTimeRange()<<" sec, overlapping by "<<tile->GetOverlapDuration()<<" sec</td></tr>"<<endl;
  report<<"  <tr><td>Sampling frequency:</td><td>"<<triggers[0]->GetWorkingFrequency()<<" Hz</td></tr>"<<endl;
  report<<"  <tr><td>Frequency range:</td><td>"<<tile->GetFrequencyMin()<<" &rarr; "<<tile->GetFrequencyMax()<<" Hz</td></tr>"<<endl;
  report<<"  <tr><td>Q range:</td><td>"<<tile->GetQ(0)<<" &rarr; "<<tile->GetQ(tile->GetNQ()-1)<<"</td></tr>"<<endl;
  report<<"  <tr><td>Tiling maximal mismatch:</td><td>"<<tile->GetMismatchMax()*100.0<<" %</td></tr>"<<endl;
  if(fOutProducts.find("triggers")!=string::npos) report<<"  <tr><td>SNR threshold (triggers):</td><td>SNR &gt; "<<tile->GetSNRTriggerThr()<<"</td></tr>"<<endl;
  if(fOutProducts.find("map")!=string::npos) report<<"  <tr><td>SNR threshold (maps):</td><td>SNR &gt; "<<tile->GetSNRMapThr()<<"</td></tr>"<<endl;
  if(fClusterAlgo.compare("none")) report<<"  <tr><td>Trigger clustering:</td><td>"<<fClusterAlgo<<", dt = "<<triggers[0]->GetClusterizeDt()<<" sec</td></tr>"<<endl;
  else report<<"  <tr><td>Trigger clustering:</td><td>NONE</td></tr>"<<endl;
  report<<"</table>"<<endl;
  report<<"<hr />"<<endl;
  report<<endl;

  // Injection parameters
  report<<"<h2>Injections</h2>"<<endl;
  if(fsginj){
    report<<"<h3>Sinusoidal Gaussian waveforms:</h3>"<<endl;
    report<<"<table>"<<endl;
    report<<"  <tr><td>Frequency range:</td><td>"<<oinj->GetFrequencyMin()<<" &rarr; "<<oinj->GetFrequencyMax()<<" Hz</td></tr>"<<endl;
    report<<"  <tr><td>Q range:</td><td>"<<oinj->GetQMin()<<" &rarr; "<<oinj->GetQMax()<<"</td></tr>"<<endl;
    report<<"  <tr><td>Amplitude range:</td><td>"<<scientific<<oinj->GetAmplitudeMin()<<" &rarr; "<<oinj->GetAmplitudeMax()<<"</td></tr>"<<endl;
    report<<"  <tr><td>Time range /chunk center:</td><td>"<<oinj->GetTimeMin()<<" &rarr; "<<oinj->GetTimeMax()<<" s</td></tr>"<<endl;
    report<<"  <tr><td>Parameters:</td><td><a href=\"./sginjections.txt\">injection parameters</a></td></tr>"<<endl;
    report<<"</table>"<<endl;
  }
  if(FFL_inject!=NULL){
    report<<"<h3>Injection channels:</h3>"<<endl;
    report<<"<table>"<<endl;
    report<<"  <tr><td>FFL data:</td><td>"<<FFL_inject->GetInputFfl()<<"</td></tr>"<<endl;
    report<<"</table>"<<endl;
    report<<"<table>"<<endl;
    report<<"  <tr><th>Main channel</th><th>Injection channel</th><th>Injection factor</th></tr>"<<endl;
    for(int c=0; c<nchannels; c++) report<<"  <tr><td>"<<triggers[c]->GetName()<<"</td><td>"<<fInjChan[c]<<"</td><td>"<<fInjFact[c]<<"</td></tr>"<<endl;
    report<<"</table>"<<endl;
  }
  if(!fsginj&&FFL_inject==NULL)
    report<<"<p>No injection was performed</p>"<<endl;
  report<<"<hr />"<<endl;
  report<<endl;

  // channel index
  report<<"<h2>Channel index</h2>"<<endl;
  if(fOutProducts.find("map")!=string::npos){// color scale for glitchiness
    report<<"<table><tr>"<<endl;
    for(int c=0; c<17; c++) report<<"<td bgcolor=\""<<colorcode[c]<<"\"></td>"<<endl;
    report<<"<td>event strength</td>"<<endl;
    report<<"</tr></table>"<<endl;
  }
  report<<"<table class=\"omicronindex\">"<<endl;
  string colcode;
  for(int c=0; c<nchannels; c++){
    colcode="";
    if(tile->GetSNRMapThr()>0) colcode=GetColorCode((chan_mapsnrmax[c]-tile->GetSNRMapThr())/tile->GetSNRMapThr());
    if(!(c%6)) report<<"  <tr>"<<endl;
    if(colcode.compare("")) report<<"    <td style=\"border:2px solid "<<colcode<<"\"><a href=\"#"<<triggers[c]->GetNameConv()<<"\">"<<triggers[c]->GetName()<<"</a></td>"<<endl;
    else report<<"    <td><a href=\"#"<<triggers[c]->GetNameConv()<<"\">"<<triggers[c]->GetName()<<"</a></td>"<<endl;
    if(!((c+1)%6)) report<<"  </tr>"<<endl;
  }
  for(int c=0; c<6-(nchannels)%6; c++) report<<"    <td></td>"<<endl;
  report<<"  </tr>"<<endl;
  report<<"</table>"<<endl;
  report<<"<hr />"<<endl;
  report<<endl;

  //**** channel report *********
  string type_first="", led, led_h;
  for(int c=0; c<nchannels; c++){

    // select processing led
    if(chan_write_ctr[c]==chunk_ctr){ led="green"; led_h="Processing OK"; }
    else if(chan_data_ctr[c]!=chunk_ctr){ led="blue"; led_h="Data access failed"; }
    else{ led="red"; led_h="Processing failed"; }

    // write processed segments
    outSegments[c]->Write(outdir[c]+"/omicron.segments.txt");

    // processing report
    if(fOutProducts.find("map")!=string::npos&&chan_mapsnrmax[c]<tile->GetSNRMapThr()){
      report<<"<h2 class=\"off\"><img src=\"./led-"<<led<<".gif\" alt=\""<<led_h<<"\" title=\""<<led_h<<"\"/>&nbsp;"<<triggers[c]->GetName()<<" <a href=\"javascript:void(0)\" name=\""<<triggers[c]->GetNameConv()<<"\" onclick=\"toggle('id_"<<triggers[c]->GetNameConv()<<"')\">[click here to expand/hide]</a></h2>"<<endl;
      report<<"<div class=\"omicronchannel\" id=\"id_"<<triggers[c]->GetNameConv()<<"\" style=\"visibility:hidden;height:0;\">"<<endl;
    }
    else{
      report<<"<h2 class=\"on\"><img src=\"./led-"<<led<<".gif\" />&nbsp;"<<triggers[c]->GetName()<<" <a href=\"javascript:void(0)\" name=\""<<triggers[c]->GetNameConv()<<"\" onclick=\"toggle('id_"<<triggers[c]->GetNameConv()<<"')\">[click here to expand/hide]</a></h2>"<<endl;
      report<<"<div class=\"omicronchannel\" id=\"id_"<<triggers[c]->GetNameConv()<<"\" style=\"visibility:visible;height:auto;\">"<<endl;
    }
    report<<"Processing:"<<endl;
    report<<"  <table class=\"omicronsummary\">"<<endl;
    report<<"    <tr><td>Number of calls [load/data/condition/projection/write]:</td><td>"<<chan_ctr[c]<<"/"<<chan_data_ctr[c]<<"/"<<chan_cond_ctr[c]<<"/"<<chan_proj_ctr[c]<<"/"<<chan_write_ctr[c]<<"</td></tr>"<<endl;
    report<<"    <tr><td>Processed livetime:</td><td>"<<(int)outSegments[c]->GetLiveTime()<<" sec ("<<setprecision(3)<<fixed<<outSegments[c]->GetLiveTime()/inSegments->GetLiveTime()*100.0<<"%) &rarr; "<<setprecision(3)<<fixed<<outSegments[c]->GetLiveTime()/3600.0/24<<" days <a href=\"./"<<triggers[c]->GetName()<<"/omicron.segments.txt\">segments</a></td></tr>"<<endl;
    report<<"  </table>"<<endl;
    
    // output products
    report<<"Output:"<<endl;
    report<<"<table class=\"omicronsummary\">"<<endl;
      
    // loop over chunks
    for(int s=0; s<(int)chunkcenter.size(); s++){
      report<<"  <tr>"<<endl;
      report<<"    <td>"<<chunkcenter[s]<<":</td>"<<endl;

      // triggers
      if(fOutProducts.find("triggers")!=string::npos)
	report<<"    <td><a href=\"./"<<triggers[c]->GetName()<<"/"<<chunktfile[s]<<"\">Triggers</a></td>"<<endl;
      
      // maps
      if(fOutProducts.find("map")!=string::npos){
	tmpstream.clear(); tmpstream.str("");
	tmpstream<<outdir[c]<<"/"<<triggers[c]->GetNameConv()<<"_OMICRONMAP-"<<chunkcenter[s]<<"-"<<fWindows[0]<<"."<<form;
	if(IsBinaryFile(tmpstream.str())){
	  for(int q=0; q<=tile->GetNQ(); q++){
	    if(q) report<<"    <td><a href=\"javascript:showImage('"<<triggers[c]->GetName()<<"', '"<<triggers[c]->GetNameConv()<<"', 'OMICRONMAPQ"<<q-1<<"-"<<chunkcenter[s]<<"', "<<windowset<<", '"<<form<<"');\">mapQ="<<setprecision(1)<<fixed<<tile->GetQ(q-1)<<"</a></td>"<<endl;
	    else  report<<"    <td><a href=\"javascript:showImage('"<<triggers[c]->GetName()<<"', '"<<triggers[c]->GetNameConv()<<"', 'OMICRONMAP"<<"-"<<chunkcenter[s]<<"', "<<windowset<<", '"<<form<<"');\">Full map</a></td>"<<endl;
	  }
	}
	else{
	  for(int q=0; q<=tile->GetNQ(); q++){
	    if(q) report<<"    <td>mapQ="<<setprecision(1)<<fixed<<tile->GetQ(q-1)<<"</td>"<<endl;
 	    else  report<<"    <td>Full map</td>"<<endl;
	  }
	}
      }
      
      // ASD
      if(fOutProducts.find("asd")!=string::npos)
	report<<"    <td><a href=\"./"<<triggers[c]->GetName()<<"/"<<triggers[c]->GetNameConv()<<"_OMICRONASD-"<<chunkcenter[s]-tile->GetTimeRange()/2<<"-"<<tile->GetTimeRange()<<"."<<form<<"\" target=\"_blank\">ASD</a></td>"<<endl;
	
      // PSD
      if(fOutProducts.find("psd")!=string::npos)
	report<<"    <td><a href=\"./"<<triggers[c]->GetName()<<"/"<<triggers[c]->GetNameConv()<<"_OMICRONPSD-"<<chunkcenter[s]-tile->GetTimeRange()/2<<"-"<<tile->GetTimeRange()<<"."<<form<<"\" target=\"_blank\">PSD</a></td>"<<endl;
	
      // conditioned time-series
      if(fOutProducts.find("timeseries")!=string::npos)
	report<<"    <td><a href=\"javascript:showImage('"<<triggers[c]->GetName()<<"', '"<<triggers[c]->GetNameConv()<<"', 'OMICRONCONDTS"<<"-"<<chunkcenter[s]<<"', "<<windowset<<", '"<<form<<"');\">Conditioned data"<<"</a></td>"<<endl;
      
      // whitened time-series
      if(fOutProducts.find("white")!=string::npos)
	report<<"    <td><a href=\"javascript:showImage('"<<triggers[c]->GetName()<<"', '"<<triggers[c]->GetNameConv()<<"', 'OMICRONWHITETS"<<"-"<<chunkcenter[s]<<"', "<<windowset<<", '"<<form<<"');\">Whitened data"<<"</a></td>"<<endl;
      
    }
    
      /*
             
      // injection
      if(fsginj==1&&fOutProducts.find("injection")!=string::npos){
      report<<"    <td><a href=\"./"<<triggers[c]->GetName()<<"/"<<triggers[c]->GetName()<<"_"<<chunkstart[s]<<"_sginjection.txt\">Injection</a></td>"<<endl;
      if(!type_first.compare("")) type_first="injection";
      }
      report<<"  </tr>"<<endl;
      }
     
      */
    report<<"</table>"<<endl;

    // add window plots
    if(fOutProducts.find("map")!=string::npos){
      report<<"<table>"<<endl;
      report<<"  <tr>"<<endl;
      for(int w=0; w<(int)fWindows.size(); w++)
	report<<"    <td><a id=\"a_"<<triggers[c]->GetNameConv()<<"-"<<fWindows[w]<<"\" href=\"./"<<triggers[c]->GetName()<<"/"<<triggers[c]->GetNameConv()<<"_OMICRONMAP-"<<chunkcenter[0]<<"-"<<fWindows[w]<<"."<<form<<"\"><img id=\"img_"<<triggers[c]->GetNameConv()<<"-"<<fWindows[w]<<"\" src=\"./"<<triggers[c]->GetName()<<"/th"<<triggers[c]->GetNameConv()<<"_OMICRONMAP-"<<chunkcenter[0]<<"-"<<fWindows[w]<<"."<<form<<"\" alt=\""<<triggers[c]->GetName()<<" "<<fWindows[w]<<"s\" /></a></td>"<<endl;
      report<<"  </tr>"<<endl;
      report<<"</table>"<<endl;
    }
    else if(fOutProducts.find("white")!=string::npos){
      report<<"<table>"<<endl;
      report<<"  <tr>"<<endl;
      for(int w=0; w<(int)fWindows.size(); w++)
	report<<"    <td><a id=\"a_"<<triggers[c]->GetNameConv()<<"-"<<fWindows[w]<<"\" href=\"./"<<triggers[c]->GetName()<<"/"<<triggers[c]->GetNameConv()<<"_OMICRONWHITETS-"<<chunkcenter[0]<<"-"<<fWindows[w]<<"."<<form<<"\"><img id=\"img_"<<triggers[c]->GetNameConv()<<"-"<<fWindows[w]<<"\" src=\"./"<<triggers[c]->GetName()<<"/th"<<triggers[c]->GetNameConv()<<"_OMICRONWHITETS-"<<chunkcenter[0]<<"-"<<fWindows[w]<<"."<<form<<"\" alt=\""<<triggers[c]->GetName()<<" "<<fWindows[w]<<"s\" /></a></td>"<<endl;
      report<<"  </tr>"<<endl;
      report<<"</table>"<<endl;
    }
    else if(fOutProducts.find("timeseries")!=string::npos){
      report<<"<table>"<<endl;
      report<<"  <tr>"<<endl;
      for(int w=0; w<(int)fWindows.size(); w++)
	report<<"    <td><a id=\"a_"<<triggers[c]->GetNameConv()<<"-"<<fWindows[w]<<"\" href=\"./"<<triggers[c]->GetName()<<"/"<<triggers[c]->GetNameConv()<<"_OMICRONCONDTS-"<<chunkcenter[0]<<"-"<<fWindows[w]<<"."<<form<<"\"><img id=\"img_"<<triggers[c]->GetNameConv()<<"-"<<fWindows[w]<<"\" src=\"./"<<triggers[c]->GetName()<<"/th"<<triggers[c]->GetNameConv()<<"_OMICRONCONDTS-"<<chunkcenter[0]<<"-"<<fWindows[w]<<"."<<form<<"\" alt=\""<<triggers[c]->GetName()<<" "<<fWindows[w]<<"s\" /></a></td>"<<endl;
      report<<"  </tr>"<<endl;
      report<<"</table>"<<endl;
    }
    else;
    
    report<<"</div>"<<endl;
    report<<"<hr />"<<endl;
    report<<endl;
  }

  // index footer
  report<<"<table>"<<endl;
  report<<"  <tr><td>Author: Florent Robinet, <a href=\"mailto:robinet@lal.in2p3.fr\">robinet@lal.in2p3.fr</a></td></tr>"<<endl;
  report<<"</table>"<<endl;
  report<<"</body>"<<endl;
  report<<"</html>"<<endl;
  report.close();
  
  return;
}

