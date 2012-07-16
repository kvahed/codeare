function navhl () {
	
	var loc = location.href;
	

	if        (loc.match(/download.html/g)) {
		var el = document.getElementById('download');
		el.setAttribute ('class', 'active');
	} else if (loc.match(/install.html/g)) {
		var el = document.getElementById('install');
		el.setAttribute ('class', 'active');
	} else if (loc.match(/gettingstarted.html/g)) {
		var el = document.getElementById('gettingstarted');
		el.setAttribute ('class', 'active');
	} else if (loc.match(/developerguide.html/g)) {
		var el = document.getElementById('developerguide');
		el.setAttribute ('class', 'active');
	} else if (loc.match(/api.html/g)) {
		var el = document.getElementById('api');
		el.setAttribute ('class', 'active');
	} else if (loc.match(/contact.html/g)) {
		var el = document.getElementById('contact');
		el.setAttribute ('class', 'active');
	} else  {
		var el = document.getElementById('index');
		el.setAttribute ('class', 'active');
	} 

}


