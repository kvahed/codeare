<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0"
				xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="html" encoding="utf-8" indent="yes"/>
  
  <xsl:template match="br">
	<br/>
  </xsl:template>
  
  <xsl:template match="a">
	<a href="{@href}"><xsl:apply-templates/></a>
  </xsl:template>
  
  <xsl:template match="main"> 

    <xsl:text disable-output-escaping='yes'>&lt;?xml version="1.0" standalone="yes" ?></xsl:text>
    <xsl:text disable-output-escaping='yes'>&lt;!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN"
    "http://www.w3.org/TR/1998/REC-html40-19980424/loose.dtd"></xsl:text>
	<html>
	  <head>
		
		<meta name="google-site-verification" content="gOYMZkytTZ8mW9ekuk7pbY67hJf0UkdjSHYKE30z8Jo" />

		<meta http-equiv="content-type" content="text/html; charset=utf-8"/>
		<meta http-equiv="content-language" content="en" />
		
		<meta name="description" lang="en" content="" />
		<meta name="keywords" content="MRI, image reconstruction, image reconstruction library, codeare,
									   Compressed sensing, non-cartesian, SENSE, kt-Points, 
									   variable density spiral, k-space	trajectory, GRAPPA, parallel imaging, 
									   selective excitation, open-source, Non-Uniform FFT, NUFFT" /> 
		
		<link rel="stylesheet" href="screen.css" type="text/css" media="screen" />
				
		<script type="text/javascript" src="sh_main.js"></script>
		<script type="text/javascript" src="sh_cpp.js"></script>
		<script type="text/javascript" src="sh_xml.js"></script>
		<script type="text/javascript" src="sh_sh.js"></script>
		<script type="text/javascript" src="nav.js"></script>
		
		<link type="text/css" rel="stylesheet" href="sh_style.css"/>

		<title>codeare 1 - <xsl:value-of select="title" disable-output-escaping="yes"/></title>
		
	  </head>
	  <body onload="sh_highlightDocument(); navhl();">
		
		<div id="topborder"></div>
		
		<div id="sidebar">
		  <h1><a href="http://codeare.org/">codeare<em>1.2</em></a></h1>
		  
		  <ul id="nav">
			<li><a id="index" href="index.html">home</a></li>
			<li><a id="download" href="download.html">download</a></li>
			<li><a id="install" href="install.html">install</a></li>
			<li><a id="gettingstarted" href="gettingstarted.html">getting started</a></li>
			<li><a id="developerguide" href="developerguide.html">developer's guide</a></li>
			<li><a id="api" href="api/html/index.html" target="_api">api annotation</a></li>
			<li><a id="contact" href="contact.html">contact</a></li>
		  </ul>
		  <script type="text/javascript" language="JavaScript">
			navhl();
		  </script>
		  <div class="callout">
			<a href="download.html" >v1.2 released<br />
			<span class="date">on Mar 26 2014 </span>
			<br/>by kaveh vahedipour</a>
		  </div>
		  <br/>

		</div>
		
		<div id="content">

		  <xsl:value-of select="article" disable-output-escaping="yes"/>
		  
		</div>
		
		<div id="rightbar">
		  <div class="callout">
			<h1>visit also</h1> <br/>
			<a href="http://www.jemris.org" target="_jemris">jemris
			mri simulator</a> 
			<br/>
			open-source full-featured parallel tx/rx mri sequence and
			hardware simulator by <em>tony stoecker, kaveh vahedipour 
			and daniel pflugfelder</em>.
			<br/><br/> 
			<a href="http://www.drcmr.dk/bloch" target="_hansen">bloch
			  simulator for education in mri and nmr</a> 
			<br/>
			free educational mri and nmr sofware. multi- os/platform 
			without installation of software. remarkable and
			beautiful effort by <em>lars hanson</em> from drcmr, hvidovre, dk.
			<br/><br/> 
			<a href="http://gadgetron.sourceforge.net"
			   target="_mriunbound">gagdetron</a>
			<br/>
			a comparable project to codeare by <em>michael hansen</em> and
			<em>thomas sorensen</em>.
			<br/><br/> 
			<a href="http://www.ismrm.org/mri_unbound/"
			   target="_mriunbound"> mri unbound</a>
			<br/>
			a collaborative forum for mri data acquisition and image
			reconstruction maintained by the unbreakable <em>jim pipe</em>. 

		  </div>
		</div>
		
	  </body>
	</html>
  </xsl:template>
</xsl:stylesheet>
    
