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
  if(!fNoLogo) system(("cp -f ${OMICRON_HTML}/pics/omicronlogo."+GPlot->GetCurrentStyle()+".gif "+maindir+"/logo.gif").c_str());
  else system(("rm -f "+maindir+"/logo.gif").c_str());
  system(("cp -f "+fOptionFile+" "+maindir+"/omicron.parameters.txt").c_str());

  // select web supported image format
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
  report<<"function showImage(channel, type, timerange, format) {"<<endl;
  report<<"  for (var dt in timerange) {"<<endl;
  report<<"    var basename ="<<endl;
  report<<"      channel + \"_\" + type + \"dt\" + timerange[dt];"<<endl;
  report<<"    var basenameth ="<<endl;
  report<<"      channel + \"_\" + type + \"dt\" + timerange[dt] + \"th\";"<<endl;
  report<<"    document.getElementById(\"a_\" + channel + \"_dt\" + timerange[dt]).href ="<<endl;
  report<<"      \"./\" + channel + \"/\" + basename + \".\" + format;"<<endl;
  report<<"    document.getElementById(\"img_\" + channel + \"_dt\" + timerange[dt]).src ="<<endl;
  report<<"      \"./\" + channel + \"/\" + basenameth + \".\" + format;"<<endl;
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

  // FIXME: add injections

  // search parameters
  report<<"<h2>Parameters</h2>"<<endl;
  report<<"<table>"<<endl;
  report<<"  <tr><td>Timing:</td><td> chunks of "<<tile->GetTimeRange()<<" sec, overlapping by "<<tile->GetOverlapDuration()<<" sec</td></tr>"<<endl;
  report<<"  <tr><td>Sampling frequency:</td><td>"<<triggers[0]->GetWorkingFrequency()<<" Hz</td></tr>"<<endl;
  report<<"  <tr><td>Frequency range:</td><td>"<<tile->GetFrequencyMin()<<" &rarr; "<<tile->GetFrequencyMax()<<" Hz</td></tr>"<<endl;
  report<<"  <tr><td>Q range:</td><td>"<<tile->GetQ(0)<<" &rarr; "<<tile->GetQ(tile->GetNQ()-1)<<"</td></tr>"<<endl;
  report<<"  <tr><td>Tiling maximal mismatch:</td><td>"<<tile->GetMismatchMax()*100<<" %</td></tr>"<<endl;
  if(fOutProducts.find("triggers")!=string::npos) report<<"  <tr><td>SNR threshold (triggers):</td><td>SNR &gt; "<<tile->GetSNRTriggerThr()<<"</td></tr>"<<endl;
  if(fOutProducts.find("maps")!=string::npos) report<<"  <tr><td>SNR threshold (maps):</td><td>SNR &gt; "<<tile->GetSNRMapThr()<<"</td></tr>"<<endl;
  if(fClusterAlgo.compare("none")) report<<"  <tr><td>Trigger clustering:</td><td>"<<fClusterAlgo<<", dt = "<<triggers[0]->GetClusterizeDt()<<" sec</td></tr>"<<endl;
  else report<<"  <tr><td>Trigger clustering:</td><td>NONE</td></tr>"<<endl;
  report<<"</table>"<<endl;
  report<<"<hr />"<<endl;
  report<<endl;

  // channel index
  report<<"<h2>Channel index</h2>"<<endl;
  if(fOutProducts.find("maps")!=string::npos){// color scale for glitchiness
    report<<"<table><tr>"<<endl;
    for(int c=0; c<17; c++) report<<"<td bgcolor=\""<<colorcode[c]<<"\"></td>"<<endl;
    report<<"<td>glitch strength</td>"<<endl;
    report<<"</tr></table>"<<endl;
  }
  report<<"<table class=\"omicronindex\">"<<endl;
  string colcode;
  for(int c=0; c<(int)fChannels.size(); c++){
    colcode="";
    if(tile->GetSNRMapThr()>0) colcode=GetColorCode((chan_mapsnrmax[c]-tile->GetSNRMapThr())/tile->GetSNRMapThr());
    if(!(c%9)) report<<"  <tr>"<<endl;
    if(colcode.compare("")) report<<"    <td style=\"border:2px solid "<<colcode<<"\"><a href=\"#"<<fChannels[c]<<"\">"<<fChannels[c]<<"</a></td>"<<endl;
    else report<<"    <td><a href=\"#"<<fChannels[c]<<"\">"<<fChannels[c]<<"</a></td>"<<endl;
    if(!((c+1)%9)) report<<"  </tr>"<<endl;
  }
  for(int c=0; c<9-((int)fChannels.size())%9; c++) report<<"    <td></td>"<<endl;
  report<<"  </tr>"<<endl;
  report<<"</table>"<<endl;
  report<<"<hr />"<<endl;
  report<<endl;

  //**** channel report *********
  string type_first="";
  for(int c=0; c<(int)fChannels.size(); c++){

    // processing report
    if(fOutProducts.find("maps")!=string::npos&&chan_mapsnrmax[c]<tile->GetSNRMapThr()){
      report<<"<h2 class=\"off\">"<<fChannels[c]<<" -- below threshold (SNR &lt; "<<tile->GetSNRMapThr()<<") <a href=\"javascript:void(0)\" name=\""<<fChannels[c]<<"\" onclick=\"toggle('id_"<<fChannels[c]<<"')\">[click here to expand/hide]</a></h2>"<<endl;
      report<<"<div class=\"omicronchannel\" id=\"id_"<<fChannels[c]<<"\" style=\"visibility:hidden;height:0;\">"<<endl;
    }
    else{
      report<<"<h2 class=\"on\">"<<fChannels[c]<<" <a href=\"javascript:void(0)\" name=\""<<fChannels[c]<<"\" onclick=\"toggle('id_"<<fChannels[c]<<"')\">[click here to expand/hide]</a></h2>"<<endl;
      report<<"<div class=\"omicronchannel\" id=\"id_"<<fChannels[c]<<"\" style=\"visibility:visible;height:auto;\">"<<endl;
    }
    report<<"Processing:"<<endl;
    report<<"  <table class=\"omicronsummary\">"<<endl;
    report<<"    <tr><td>Number of calls [load/data/condition/projection/write]:</td><td>"<<chan_ctr[c]<<"/"<<chan_data_ctr[c]<<"/"<<chan_cond_ctr[c]<<"/"<<chan_proj_ctr[c]<<"/"<<chan_write_ctr[c]<<"</td></tr>"<<endl;
    report<<"    <tr><td>Processed livetime:</td><td>"<<(int)outSegments[c]->GetLiveTime()<<" sec ("<<setprecision(3)<<fixed<<outSegments[c]->GetLiveTime()/inSegments->GetLiveTime()*100.0<<"%) &rarr; "<<setprecision(3)<<fixed<<outSegments[c]->GetLiveTime()/3600.0/24<<" days</td></tr>"<<endl;
    outSegments[c]->Write(outdir[c]+"/omicron.segments.txt");
    report<<"    <tr><td>Processed segments:</td><td><a href=\"./"<<fChannels[c]<<"/omicron.segments.txt\">omicron.segments.txt</a></td></tr>"<<endl;
    report<<"  </table>"<<endl;
    
    // Plot links
    if(form.compare("") &&
       (fOutProducts.find("timeseries")!=string::npos ||
	fOutProducts.find("maps")!=string::npos ||
	fOutProducts.find("white")!=string::npos ||
	fOutProducts.find("asd")!=string::npos ||
	fOutProducts.find("psd")!=string::npos)){

      report<<"Plots:"<<endl;
      report<<"<table class=\"omicronsummary\">"<<endl;
      
      // loop over chunks
      for(int s=0; s<(int)chunkcenter.size(); s++){
	report<<"  <tr>"<<endl;
	report<<"    <td>"<<chunkcenter[s]<<":</td>"<<endl;
	
	// maps
	if(fOutProducts.find("maps")!=string::npos){
	  for(int q=0; q<=tile->GetNQ(); q++){
	    tmpstream<<outdir[c]<<"/"<<fChannels[c]<<"_"<<chunkcenter[s]<<"_fullmapdt"<<fWindows[0]<<"."<<form;
	    if(!IsBinaryFile(tmpstream.str())){// check full map exists
	      report<<"    <td colspan=\""<<tile->GetNQ()+1<<"\">no maps</td>"<<endl; tmpstream.clear(); tmpstream.str("");
	      break;
	    }
	    tmpstream.clear(); tmpstream.str("");
	    if(q) report<<"    <td><a href=\"javascript:showImage('"<<fChannels[c]<<"', '"<<chunkcenter[s]<<"_mapQ"<<q-1<<"', "<<windowset<<", '"<<form<<"');\">mapQ="<<setprecision(1)<<fixed<<tile->GetQ(q-1)<<"</a></td>"<<endl;
	    else  report<<"    <td><a href=\"javascript:showImage('"<<fChannels[c]<<"', '"<<chunkcenter[s]<<"_fullmap', "<<windowset<<", '"<<form<<"');\">Full map</a></td>"<<endl;
	    if(!type_first.compare("")) type_first="fullmap";
	  }
	}

	// ASD
	if(fOutProducts.find("asd")!=string::npos){
	  report<<"    <td><a href=\"./"<<fChannels[c]<<"/"<<fChannels[c]<<"_"<<chunkcenter[s]<<"_ASD."<<form<<"\" target=\"_blank\">ASD</a></td>"<<endl;
	  if(!type_first.compare("")) type_first="";
	}

	// PSD
	if(fOutProducts.find("psd")!=string::npos){
	  report<<"    <td><a href=\"./"<<fChannels[c]<<"/"<<fChannels[c]<<"_"<<chunkcenter[s]<<"_PSD."<<form<<"\" target=\"_blank\">PSD</a></td>"<<endl;
	  if(!type_first.compare("")) type_first="";
	}

	// Spectral
	if(fOutProducts.find("spectral")!=string::npos){
	  report<<"    <td><a href=\"./"<<fChannels[c]<<"/"<<fChannels[c]<<"_"<<chunkcenter[s]<<"_spec."<<form<<"\" target=\"_blank\">Spectral</a></td>"<<endl;
	  if(!type_first.compare("")) type_first="";
	}

	// time-series
	if(fOutProducts.find("timeseries")!=string::npos){
	  report<<"    <td><a href=\"javascript:showImage('"<<fChannels[c]<<"', '"<<chunkcenter[s]<<"_ts', "<<windowset<<", '"<<form<<"');\">Conditionned time series</a></td>"<<endl;
	  if(!type_first.compare("")) type_first="ts";
	}
	
	// whitened data
	if(fOutProducts.find("white")!=string::npos){
	  report<<"    <td><a href=\"javascript:showImage('"<<fChannels[c]<<"', '"<<chunkcenter[s]<<"_whitets', "<<windowset<<", '"<<form<<"');\">Whitened time series</a></td>"<<endl;
	  if(!type_first.compare("")) type_first="whitets";
	}

	// triggers
	if(fOutProducts.find("triggers")!=string::npos){
	  report<<"    <td><a href=\"./"<<fChannels[c]<<"/"<<fChannels[c]<<"_"<<chunkcenter[s]<<".root\">Triggers</a></td>"<<endl;
	  if(!type_first.compare("")) type_first="triggers";
	}
	
	report<<"  </tr>"<<endl;
      }

      /*      
  
    // Map links
    if(form.compare("")&&fOutProducts.find("maps")!=string::npos){
      report<<"Maps:"<<endl;
      report<<"<table class=\"omicronsummary\">"<<endl;
      
      // loop over segments
      sfirst=-1;
      for(int s=0; s<(int)mapcenter.size(); s++){
	
	// check image presence
	tmpstream<<fChannels[c]<<"_"<<mapcenter[s]<<"_fullmapdt"<<fWindows[0];
	if(!IsBinaryFile(outdir[c]+"/"+tmpstream.str()+"."+form)){
	  tmpstream.clear(); tmpstream.str("");
	  continue;
	}
       	tmpstream.clear(); tmpstream.str("");
	if(sfirst<0) sfirst=s;

	// add links
	report<<"  <tr>"<<endl;
	for(int q=0; q<=tile->GetNQ(); q++){
	  if(q) report<<"    <td><a href=\"javascript:showImage('"<<fChannels[c]<<"', '"<<mapcenter[s]<<"_mapQ"<<q-1<<"', "<<windowset<<", '"<<form<<"');\">Q="<<setprecision(1)<<fixed<<tile->GetQ(q-1)<<"</a></td>"<<endl;
	  else  report<<"    <td>"<<mapcenter[s]<<":</td><td><a href=\"javascript:showImage('"<<fChannels[c]<<"', '"<<mapcenter[s]<<"_fullmap', "<<windowset<<", '"<<form<<"');\">Full map</a></td>"<<endl;
	}
	report<<"  </tr>"<<endl;
      }
      if(sfirst<0) report<<"<tr><td>missing</td></tr>"<<endl;
      */
      report<<"</table>"<<endl;

      // add window plots
      if(type_first.compare("")){
	report<<"<table>"<<endl;
	report<<"  <tr>"<<endl;
	for(int w=0; w<(int)fWindows.size(); w++)
	  report<<"    <td><a id=\"a_"<<fChannels[c]<<"_dt"<<fWindows[w]<<"\" href=\"./"<<fChannels[c]<<"/"<<fChannels[c]<<"_"<<chunkcenter[0]<<"_"<<type_first<<"dt"<<fWindows[w]<<"."<<form<<"\"><img id=\"img_"<<fChannels[c]<<"_dt"<<fWindows[w]<<"\" src=\"./"<<fChannels[c]<<"/"<<fChannels[c]<<"_"<<chunkcenter[0]<<"_"<<type_first<<"dt"<<fWindows[w]<<"th."<<form<<"\" alt=\""<<fChannels[c]<<" "<<type_first<<" dt="<<fWindows[w]<<"\" /></a></td>"<<endl;
	report<<"  </tr>"<<endl;
	report<<"</table>"<<endl;
      }
    }
     
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

