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
  
  vector <string> chan;
  if(!ListDirectories(reportdir,chan)){
    cerr<<"Omicron::ScanReport: no channel directories in "<<reportdir<<endl;
    return false;
  }

  // import material
  system(("cp -f ${GWOLLUM_DOC}/style.css "+reportdir).c_str());
  system(("cp -f ${GWOLLUM_DOC}/Pics/gwollum_logo_min_trans.gif "+reportdir).c_str());
  system(("cp -f ${OMICRON_HTML}/pics/omicronlogo_xxl.gif "+reportdir).c_str());
 
  // index header
  ofstream report((reportdir+"/index.html").c_str());
  report<<"<html>"<<endl;
  report<<"<head>"<<endl;
  report<<"<title>GWOLLUM: OmiScan</title>"<<endl;
  report<<"<link rel=\"stylesheet\" href=\"style.css\" type=\"text/css\">"<<endl;
  report<<"<link rel=\"icon\" type=\"image/x-icon\" href=\"icon.gif\" />"<<endl;
  report<<"<link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"icon.gif\" />"<<endl;
  report<<"</head>"<<endl;
  report<<"<body class=\"omiscan\">"<<endl;

  // index summary
  report<<"<table><tr><td>"<<endl;
  report<<"<h1>OmiScan of [GPS_CENTER]</h1>"<<endl;
  report<<"<a href=\"./comparison_mode.html\">Switch to comparison mode</a>"<<endl;
  report<<"<pre>"<<endl;
  report<<"Central GPS time:   "<<endl;
  report<<"Central Date (UTC): "<<endl;
  report<<"</pre>"<<endl;
  report<<"</td></tr></table>"<<endl;




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

