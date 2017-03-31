/*
 * Copyright (C) 2013 NASA PhoneSat
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA 
 */
"use strict";

//Return SVG color from Hue Saturation and Value
//http://en.wikipedia.org/wiki/HSL_and_HSV#Converting_to_RGB
function getSvgColor(H, S, V) {
	var C = V * S;
	var m = V - C;
	var HP = H / 60;
	var X = C * (1 - Math.abs(HP % 2 - 1));
	var R,G,B;
	if (HP < 1) {
		R = C;
		G = X;
		B = 0;
	} else if (HP < 2) {
		R = X;
		G = C;
		B = 0;
	} else if (HP < 3) {
		R = 0;
		G = C;
		B = X;
	} else if (HP < 4) {
		R = 0;
		G = X;
		B = C;
	} else if (HP < 5) {
		R = X;
		G = 0;
		B = C;
	} else {
		R = C;
		G = 0;
		B = X;
	}
	R += m;
	G += m;
	B += m;
	return 'rgb(' + Math.floor(R) + ',' + Math.floor(G) + ',' + Math.floor(B)
			+ ')';
}

// ////////////////////////////////////////////////////
// Latitude and longitude orbit plot on the earth map
// ////////////////////////////////////////////////////

sat_t.prototype.mapInit = function(map) {

	this.mapSvg = document.getElementById('map');
	if (!this.mapSvg) {
		return;
	}

	// Create satellite dot
	this.mapDot = L.circleMarker(
			[0,0],
			{
				radius: 6,
				stroke: false,
				fill: true,
				fillColor: this.color,
				fillOpacity: 1,

			}
	);

	// Create satellite footprint
	this.mapFootprint = L.polygon(
			[],
			{
				stroke: true,
				color: this.color,
				weight: 2,
				fill: true,
			}
	);

	// Create and draw satellite orbit
	{
		var current_orbit = [];
		var all_orbits = [];
		var t = new Date();
		/* TODO launch_time for each satellite
		if (t.getTime() < launch_time) {
			t.setTime(launch_time);
		}
		*/
		var g = "";
		var previous = 0;
		for ( var i = 0; i < COUNT; i++) {
			predict_calc(this, qth, Julian_Date(t));
			if (Math.abs(this.ssplon - previous) > 180) {
				// orbit crossing -PI, PI
				all_orbits.push(current_orbit);
				current_orbit = [];
			}
			current_orbit.push([this.ssplat, this.ssplon]);
			previous = this.ssplon;
			// Increase time for next point
			t.setTime(t.getTime() + STEP);
		}
		this.mapOrbit = L.multiPolyline(
				all_orbits,
				{
					color: this.color,
					weight: 4				
				}
		);
		
		// Create the layer
		var layer = L.layerGroup([this.mapDot, this.mapFootprint, this.mapOrbit]);
		layer.addTo(map);
		
		layers_control.addOverlay(layer, this.name);
	}
};

sat_t.prototype.mapRefresh = function() {
	if (!this.mapSvg) {
		return;
	}

	// Refresh satellite dot
	this.mapDot.setLatLng([this.ssplat, this.ssplon]);

	// Refresh satellite footprint
	{
		var azi;
		//var msx, msy, ssx, ssy;
		var ssplat, ssplon, beta, azimuth, num, dem;
		//var rangelon, rangelat, mlon;

		var geo = new geodetic_t();

		/* Range circle calculations.
		 * Borrowed from gsat 0.9.0 by Xavier Crehueras, EB3CZS
		 * who borrowed from John Magliacane, KD2BD.
		 * Optimized by Alexandru Csete and William J Beksi.
		 */
		ssplat = Radians(this.ssplat);
		ssplon = Radians(this.ssplon);
		beta = (0.5 * this.footprint) / xkmper;

		var gn = "", gp = "", g = "";
		var pos_overlap = false;
		var neg_overlap = false;
		var points = [];

		for (azi = 0; azi < 360; azi += 5) {
			azimuth = de2ra * azi;
			geo.lat = asin(sin(ssplat) * cos(beta) + cos(azimuth) * sin(beta)
					* cos(ssplat));
			num = cos(beta) - (sin(ssplat) * sin(geo.lat));
			dem = cos(ssplat) * cos(geo.lat);

			if (azi == 0 && (beta > pio2 - ssplat))
				geo.lon = ssplon + pi;

			else if (azi == 180 && (beta > pio2 + ssplat))
				geo.lon = ssplon + pi;

			else if (fabs(num / dem) > 1.0)
				geo.lon = ssplon;

			else {
				if ((180 - azi) >= 0)
					geo.lon = ssplon - arccos(num, dem);
				else
					geo.lon = ssplon + arccos(num, dem);
			}

			points.push([Degrees(geo.lat), Degrees(geo.lon)]);
		}

		this.mapFootprint.setLatLngs(points);
	}
};

// ///////////////////////////////////////////////////////
// Azimuth and elevation polar plot
// ///////////////////////////////////////////////////////

sat_t.prototype.polarInit = function() {
	this.orbitSvg = document.getElementById('polar');
	if (!this.orbitSvg) {
		return;
	}
	// Create polar dot
	this.polarDot = document.createElementNS("http://www.w3.org/2000/svg", "circle");
	this.polarDot.setAttributeNS(null, 'fill', this.color);
	this.polarDot.setAttributeNS(null, 'stroke', 'none');
	this.polarDot.setAttributeNS(null, 'r', '4');
};

sat_t.prototype.polarRefresh = function() {
	if (!this.orbitSvg) {
		return;
	}
	if (this.el > 0) {
		var coord = this.polarGetXY(this.az, this.el);
		this.polarDot.setAttributeNS(null, 'cx', coord.x);
		this.polarDot.setAttributeNS(null, 'cy', coord.y);
		if (!this.polarDot.parentNode) {
			this.orbitSvg.appendChild(this.polarDot);
		}
	} else {
		if (this.polarDot.parentNode) {
			this.orbitSvg.removeChild(this.polarDot);
		}
	}
};

sat_t.prototype.polarGetXY = function(az, el) {
	var ret = new Object();
	ret.x = (90 - el) * Math.sin(Radians(az));
	ret.y = (el - 90) * Math.cos(Radians(az));
	return ret;
};

pass_t.prototype.polarInit = function() {
	this.orbitSvg = document.getElementById('polar');
	if (!this.orbitSvg) {
		return;
	}
	this.polarOrbit = document.createElementNS("http://www.w3.org/2000/svg", "path");
	this.polarOrbit.setAttributeNS(null, 'fill', 'none');
	this.polarOrbit.setAttributeNS(null, 'stroke', this.sat.color);
	this.polarOrbit.setAttributeNS(null, 'stroke-opacity', '0.2');
	this.polarOrbit.setAttributeNS(null, 'stroke-width', '1');
	// Draw the orbit pass on the polar az/el plot
	var g = "";
	for ( var i = 0; i < this.details.length; i++) {
		var coord = this.sat.polarGetXY(this.details[i].az, this.details[i].el);
		if (i == 0) {
			// Start of line
			g += "M";
		} else {
			// Continue line
			g += " L";
		}
		g += coord.x + " " + coord.y;
	}
	this.polarOrbit.setAttribute("d", g);
	this.orbitSvg.appendChild(this.polarOrbit);
};

pass_t.prototype.polarRemove = function() {
	this.orbitSvg.removeChild(this.polarOrbit);
};

//removeChild

pass_t.prototype.polarEnable = function(enabled) {
	if (!this.orbitSvg) {
		return;
	}
	if (enabled) {
		this.polarOrbit.setAttributeNS(null, "stroke-opacity", '1');
		this.polarOrbit.setAttributeNS(null, 'stroke-width', '2');
	} else {
		this.polarOrbit.setAttributeNS(null, "stroke-opacity", '0.2');
		this.polarOrbit.setAttributeNS(null, 'stroke-width', '1');
	}
};

// ////////////////////////////////////////////////////////////////
// Periodic refresh loop
// ////////////////////////////////////////////////////////////////

function refresh() {
	var jt = Julian_Date(new Date());

	// Compute current satellites longitude and latitude
	for (var i = 0; i < sats.length; i++) {
		predict_calc(sats[i], qth, jt);
		sats[i].mapRefresh();
		sats[i].polarRefresh();
	}

	refresh_sat_info();
}

function refresh_sat_info() {
	// Refresh selected pass information
	if (displayed_pass) {
		if (sat_name) {
			sat_name.innerHTML = displayed_pass.sat.name;
		}
		if (sat_az) {
			sat_az.innerHTML = displayed_pass.sat.az.toFixed(1) + "&deg;";
		}
		if (sat_el) {
			sat_el.innerHTML = displayed_pass.sat.el.toFixed(1) + "&deg;";
		}
		if (sat_freq) {
			sat_freq.innerHTML = (displayed_pass.sat.freq * (1 - displayed_pass.sat.range_rate / 299792.4580))
					.toFixed(3)
					+ " MHz";
		}
		if (sat_alt) {
			sat_alt.innerHTML = displayed_pass.sat.alt.toFixed(1) + " km";
		}
	}
}

function ISODateString(d) {
	function pad(n) {
		return n < 10 ? '0' + n : n;
	}
	return d.getUTCFullYear() + pad(d.getUTCMonth() + 1)
			+ pad(d.getUTCDate()) + 'T' + pad(d.getUTCHours())
			+ pad(d.getUTCMinutes()) + pad(d.getUTCSeconds()) + 'Z';
}

var sat_name = document.getElementById("sat_name");
var sat_az = document.getElementById("sat_az");
var sat_el = document.getElementById("sat_el");
var sat_freq = document.getElementById("sat_freq");
var sat_alt = document.getElementById("sat_alt");
var polar_window = document.getElementById("polar_window");
var next_passes_window = document.getElementById("next_passes_window");
var timeline_mask = document.getElementById("timeline_mask");
var displayed_pass = null;
var selected_row = null;
var all_passes = [];

function hide_timeline() {
	timeline_mask.hidden = true;
}

function show_timeline() {
	timeline_mask.hidden = false;
}

function onGroundStationChanged() {
	var next_passes = document.getElementById("next_passes");
	if (next_passes) {
		for (var i in all_passes) {
			all_passes[i].polarRemove();
		}
		// Compute next passes
		all_passes = [];
		{
			var t = Julian_Date(new Date());
			for (var i = 0; i < sats.length; i++) {
				all_passes = all_passes.concat(get_passes (sats[i], qth, t, 15, 10));
			}
		}
		for (var i = 0; i < all_passes.length; i++) {
			all_passes[i].polarInit();
		}
		all_passes.sort(function(a,b){return a.aos - b.aos;});

		// Create the next passes table
		while(next_passes.rows.length > 1) {
			next_passes.deleteRow(1);
		}
		for (var i = 0; i < all_passes.length; i++) {
			var row = next_passes.insertRow(i+1);
			row.pass = all_passes[i];
			row.className = "clickable";
			row.onmouseover = function() {
				if (displayed_pass != this.pass) {
					if (displayed_pass) {
						displayed_pass.polarEnable(false);
					}
					displayed_pass = this.pass;
					displayed_pass.polarEnable(true);
					refresh_sat_info();
				}
			};
			row.onmouseout = function() {
				if (displayed_pass != selected_row.pass) {
					if (displayed_pass) {
						displayed_pass.polarEnable(false);
					}
					displayed_pass = selected_row.pass;
					displayed_pass.polarEnable(true);
					refresh_sat_info();
				}
			};
			row.onclick = function() {
				// next line is to change the displayed pass on touch-screen devices
				this.onmouseover();
				if (selected_row) {
					selected_row.className = "clickable";
				}
				this.className = "selected";
				selected_row = this;
			};
			row.insertCell(0).innerHTML = all_passes[i].sat.name;
			row.insertCell(1).innerHTML = Date_Time(all_passes[i].aos).toLocaleString();
			row.insertCell(2).innerHTML = new Number(all_passes[i].max_el).toFixed(1) + "&deg;";

			if (i == 0) {
				row.onclick();
			}
		}

		// Create the download link
		var dtstamp = ISODateString(new Date());
		var ical=
			"BEGIN:VCALENDAR\r\n"+
			"PRODID:-//phonesat.org///EN\r\n"+
			"VERSION:2.0\r\n"+
			"METHOD:REQUEST\r\n";
		for (var i = 0; i < all_passes.length; i++) {
			ical +=
				"BEGIN:VEVENT\r\n"+
				"DTSTAMP:" + dtstamp + "\r\n" +
				"UID:" + dtstamp + "." + i + "@phonesat.org\r\n" +
				"DTSTART:" + ISODateString(Date_Time(all_passes[i].aos)) + "\r\n" +
				"DTEND:" + ISODateString(Date_Time(all_passes[i].los)) + "\r\n" +
				"SUMMARY:PhoneSat radio signal tracking\r\n" +
				"ATTACH:http://www.phonesat.org/\r\n" +
				"LOCATION:" + qth.lat + "\\," + qth.lon + "\r\n" +
				"END:VEVENT\r\n";
		}
		ical += "END:VCALENDAR\r\n";
		document.getElementById("ical").setAttribute("href", "data:text/calendar;base64," + window.btoa(ical));
	}
}
