bplist00�_WebMainResource�	
_WebResourceMIMEType_WebResourceTextEncodingName^WebResourceURL_WebResourceFrameName_WebResourceDataYtext/htmlUutf-8_file:///index.htmlPO6<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
  <meta http-equiv="Content-Style-Type" content="text/css">
  <title></title>
  <meta name="Generator" content="Cocoa HTML Writer">
  <meta name="CocoaVersion" content="2022.6">
  <style type="text/css">
    p.p1 {margin: 0.0px 0.0px 0.0px 0.0px; font: 12.0px Helvetica}
    p.p2 {margin: 0.0px 0.0px 0.0px 0.0px; font: 12.0px Helvetica; min-height: 14.0px}
  </style>
</head>
<body>
<p class="p1"># Galaxy simulation</p>
<p class="p2"><br></p>
<p class="p1">This repository contain some Python and shell scripts in order to do a simple simulation of a galaxy. We compute the trajectories of a thousand of stars. We use the PP method (for Particule-Particule) and we'll try after the implement a MP method (for Mesh-Particule). This repository contain :</p>
<p class="p2"><br></p>
<p class="p1">- **main.py** : This file is the main file with the main *for loop*. We compute the sum of all forces on each particules by the others. We install a cut-off in order to limit the computing time.</p>
<p class="p2"><br></p>
<p class="p1">- **anim.py** : This file allows to create and save the animation that we produce by the *FuncAnimation* Python module.</p>
<p class="p2"><br></p>
<p class="p1">- **run_simu.sh** : This shell script allows to run the *main.py* file and give some argument like the number of stars $$ N $$, the computing time $$ t $$ or the time integration $$ dt $$ for example.</p>
<p class="p2"><br></p>
<p class="p1">- **run_anim.sh** : This shell script allows to run the *anim.py* file and give argument to it.</p>
<p class="p2"><br></p>
<p class="p1">- **Results** : This file contains a MarkDown script which present the results.</p>
</body>
</html>
    ( > \ k � � � � � �                           �