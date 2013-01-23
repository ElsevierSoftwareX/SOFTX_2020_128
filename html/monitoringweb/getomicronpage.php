<html>
<head>
<title>Omicron: your request</title>
<link rel="stylesheet" href="style.css" type="text/css">
<link rel="icon" type="image/x-icon" href="icon.gif" />
<link rel="shortcut icon" type="image/x-icon" href="icon.gif" />
</head>

<body>
    
   <!-- check arguments from file -->
    
   <?php
   $badNumber = "/([^0-9.])/";

   // environment. FIXME: hard-coded
   $OMICRON = "/virgoDev/Omicron/v0r1"

   // get values
   $timemode = $_REQUEST['timemode'];
   $gpsvalue = trim($_REQUEST['gpsvalue']);
   $year = $_REQUEST['year'];
   $month = $_REQUEST['month'];
   $day = $_REQUEST['day'];
   $hour = $_REQUEST['hour'];
   $predefvalue = $_REQUEST['predefvalue'];
   $channel = trim($_REQUEST['channel']);
   $fullday = $_REQUEST['fullday'];
   
   // Only list available channels
   if($_REQUEST['getchannel'] == 'Get Available Channels'){
     $chandirs = array();
     
     // loop over channel directories
     if ($handle = opendir("{$OMICRONTRIGGERS}")) {
       while (false !== ($chandir = readdir($handle))) {
	 if(!strpos($chandir,'_')) continue;
	 $chandirs[] = $chandir;
       }
       closedir($handle);
     }
     sort($chandirs);
     echo "<pre>\n";
     print_r($chandirs);
     //foreach($chandir in $chandirs) echo "$chandir\n";
     echo "</pre>\n";
   }
   else{
     if($timemode=="date"){// ----------- Date input
       if($hour=="all") $fullday="yes";
       else $fullday="no";
     }
     else{
       if($timemode=="gps"){// ----------- GPS time input
	 if(empty($gpsvalue)) die("Please specify a GPS value");
	 if(preg_match($badNumber, $gpsvalue)) die("Your GPS value does not make any sense");
	 if($gpsvalue<700000000) die("Your GPS value does not make any sense");
	 $tstamp=exec("source ${OMICRON}/cmt/setup.csh; tconvert -f \"%Y %m %d %H\" {$gpsvalue}");
       }
       else{// ----------- Predefined time input
	 if($predefvalue=="latestday"){// latest day
	   if(file_exists("./latestday/{$channel}.html")){
	     header("Location: ./latestday/{$channel}.html");
	     exit;
	   }
	   else die("Channel {$channel} is not available for the last day");
	 }
	 if($predefvalue=="latesthour"){// latest hour
	   if(file_exists("./latesthour/{$channel}.html")){
	     header("Location: ./latesthour/{$channel}.html");
	     exit;
	   }
	   else die("Channel {$channel} is not available for the last hour");
	 }
	 elseif($predefvalue=="thishour"){// this hour
	   $gps_now=exec("source ${OMICRON}/cmt/setup.csh; tconvert now");
	   $tstamp=exec("source ${OMICRON}/cmt/setup.csh; tconvert -f \"%Y %m %d %H\" {$gps_now}");
	   $fullday="no";	 
	 }
	 else{// today
	   $gps_now=exec("source ${OMICRON}/cmt/setup.csh; tconvert now");
	   $tstamp=exec("source ${OMICRON}/cmt/setup.csh; tconvert -f \"%Y %m %d %H\" {$gps_now}");
	   $fullday="yes";
	 }
       }
       // convert gps into date elements
       $date_elmt = explode(" ", $tstamp);
       $year=$date_elmt[0];
       $month=$date_elmt[1];
       $day=$date_elmt[2];
       $hour=$date_elmt[3];
     }
     
     // redirect to existing page if it exists
     if($fullday=="yes"){
       if(file_exists("./{$year}/{$month}/{$day}/{$channel}.html")){
	 header("Location: ./{$year}/{$month}/{$day}/{$channel}.html");
	 exit;
       }
       else die("Your request ({$channel}: {$year}-{$month}-{$day}) does not lead to an existing web page.<br>Do you want to process your request anyway and create a web report on demand?<br>Well, be patient; this feature will be available later...");
     }
     else{
       if(file_exists("./{$year}/{$month}/{$day}/{$hour}/{$channel}.html")){
	 header("Location: ./{$year}/{$month}/{$day}/{$hour}/{$channel}.html");
	 exit;
       }
       else die("Your request ({$channel}: {$year}-{$month}-{$day} {$hour}h) does not lead to an existing web page.<br>Do you want to process your request anyway and create a web report on demand?<br>Well, be patient; this feature will be available later...");
     }

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
