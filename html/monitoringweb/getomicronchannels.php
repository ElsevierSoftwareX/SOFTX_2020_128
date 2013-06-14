<html>
<head>
<title>Omicron channels</title>
<link rel="stylesheet" href="style.css" type="text/css">
<link rel="icon" type="image/x-icon" href="icon.gif" />
<link rel="shortcut icon" type="image/x-icon" href="icon.gif" />
</head>

<body>

<h1>List of available channel reports</h1>

<?php
include "tconvert.php";
$badNumber = "/([^0-9.])/";

// set the default timezone to use: UTC
date_default_timezone_set('UTC');

// get values
$timemode = $_REQUEST['timemode'];
$gpsvalue = trim($_REQUEST['gpsvalue']);
$year = $_REQUEST['year'];
$month = $_REQUEST['month'];
$day = $_REQUEST['day'];
$hour = $_REQUEST['hour'];
$predefvalue = $_REQUEST['predefvalue'];
$fullday = $_REQUEST['fullday'];

$found=0;
   
if($timemode=="date"){// ----------- Date input
  if($hour=="all") $fullday="yes";
  else $fullday="no";
}
else{
  if($timemode=="gps"){// ----------- GPS time input
    if(empty($gpsvalue)) die("Please specify a GPS value");
    if(preg_match($badNumber, $gpsvalue)) die("Your GPS value does not make any sense");
    if($gpsvalue<700000000) die("Your GPS value does not make any sense");
    $unixtime=gps2unix($gpsvalue);
    $date_elmt = explode(" ",date("Y m d H",$unixtime));
  }
  else{// ----------- Predefined time input
    if($predefvalue=="latestday"){// latest day

      foreach (glob("./latestday/*.html") as $html_link) {
	
	// retrieve channel name
	$channel_path=explode("/",$html_link);
	$channel=$channel_path[count($channel_path)-1];
	$channel=substr($channel, 0, -5);
	echo "<a href=\"$html_link\">$channel</a><br>";
      }
      exit;
    }
    if($predefvalue=="latesthour"){// latest hour
      foreach (glob("./latesthour/*.html") as $html_link) {
	
	// retrieve channel name
	$channel_path=explode("/",$html_link);
	$channel=$channel_path[count($channel_path)-1];
	$channel=substr($channel, 0, -5);
	echo "<a href=\"$html_link\">$channel</a><br>";
      }
      exit;
    }
    elseif($predefvalue=="thishour"){// this hour
      $date_elmt = explode(" ",date("Y m d H"));
      $fullday="no";
    }
    else{// today
      $date_elmt = explode(" ",date("Y m d H"));
      $fullday="yes";
    }
  }
  
  // get date elements
  $year=$date_elmt[0];
  $month=$date_elmt[1];
  $day=$date_elmt[2];
  $hour=$date_elmt[3];
}
  
// list channels for this day
if($fullday=="yes"){
  
  foreach (glob("./{$year}/{$month}/{$day}/*.html") as $html_link) {
    
    // retrieve channel name
    $channel_path=explode("/",$html_link);
    $channel=$channel_path[count($channel_path)-1];
    $channel=substr($channel, 0, -5);
    $found=1;
    echo "<a href=\"$html_link\">$channel</a><br>";
  }
  
}
else{

  foreach (glob("./{$year}/{$month}/{$day}/{$hour}/*.html") as $html_link) {
    
    // retrieve channel name
    $channel_path=explode("/",$html_link);
    $channel=$channel_path[count($channel_path)-1];
    $channel=substr($channel, 0, -5);
    $found=1;
    echo "<a href=\"$html_link\">$channel</a><br>";
  }
}

if($found==0){
  echo"Sorry, there is no channel report available for your request: ";
  if($fullday=="yes") echo"{$year}-{$month}-{$day}<br/>";
  else                echo"{$year}-{$month}-{$day} {$hour}h<br/>";
  echo"<br/>";
 }

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
