<html>
<head>
<title>OmiScan: your scan</title>
<link rel="stylesheet" href="style.css" type="text/css">
<link rel="icon" type="image/x-icon" href="icon.gif" />
<link rel="shortcut icon" type="image/x-icon" href="icon.gif" />
</head>

<body>
    
<!-- check arguments -->
   
<?php
$badNumber = "/([^0-9.])/";

// get values
$gpscenter = $_REQUEST['gpscenter'];
$mainchannel = $_REQUEST['mainchannel'];
$thr = $_REQUEST['thr'];
$type = $_REQUEST['type'];

// check gps value
if(empty($gpscenter)) die("Please specify a GPS value");
if(preg_match($badNumber, $gpscenter)) die("Your GPS value does not make any sense");
if($gpscenter<700000000) die("Your GPS value does not make any sense");

// check main channel
if(empty($mainchannel)) die("Please specify a main channel (h_4096Hz maybe?)");

die("Sorry, OmiScan can not be started from the web (yet!). Come back later");

// FIXME: I can't write scripts in any other directory
// output directory in /tmp
$outdir=$OMICRON."/scan";
$scriptdir="/opt/w3/DataAnalysis/VOD/output/script";

// command
$command=$GWOLLUM."/local/scripts/GetOmicronScan.sh -m ".$mainchannel." -d ".$outdir." -s ".$thr." ".$gpscenter." >> /dev/null";
$cleaningcommand1="find {$outdir}/ -maxdepth 1 -type d -mtime +3 -exec rm -r {} \;";

// random key for unicity
$randomkey=rand();

// command file	
$commandfile=$scriptdir."/omiscan.".$randomkey.".sh";

//echo exec("ls -l /data/procdata/bufferv16/Omicron/Omicron | grep scan");
//echo exec("pwd");
//echo exec("whoami");

// open command file and fill it
$fp=fopen($commandfile,'w') or die("can't open file {$commandfile}");
fwrite($fp,"#!/bin/bash\n");
fwrite($fp,$cleaningcommand1."\n");
fwrite($fp,$command."\n");
fwrite($fp,"exit 0\n");
fclose($fp);
chmod($commandfile,0777);

// redirect to existing page if it exists
if(file_exists("{$outdir}/{$gpscenter}/index.html")){
  header("Location: {$outdir}/{$gpscenter}/index.html");
  exit;
}
else die("Your Omiscan failed.");
?>  


  <!-- foot page -->
  <table valign="bottom" align="center">
    <tr>
      <td valign="bottom" align="right">
	<img src="./omicronlogo_s.gif" width="42">
      </td>
      <td valign="middle" align="left">
	<p>
	  Florent Robinet<br>
	  <a HREF="mailto:robinet@lal.in2p3.fr">Contact: robinet@lal.in2p3.fr</a>
	</p>
      </td>
    </tr>
  </table>
  
</body>
</html>
