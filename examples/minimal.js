/*
 * Copyright (C) 2013 NASA PhoneSat
 *               2017 Fabian P. Schmidt <kerel-fs@gmx.de>
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

//function satellite_t(){}
// ////////////////////////////////////////////////////
// Latitude and longitude orbit plot on the earth map
// ////////////////////////////////////////////////////

Sat_t.prototype.mapInit = function(map) {

	this.mapSvg = document.getElementById('map');
	if (!this.mapSvg) {
		return;
	}

	predict_calc(this, qth, julian_date(new Date()));

	// Create satellite dot
	this.mapDot = L.circleMarker(
			[this.ssplat, this.ssplon],
			{
				radius: 6,
				stroke: false,
				fill: true,
				fillColor: 'rgb(0,0,0)',
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

		var g = "";
		var previous = 0;
		for ( var i = 0; i < COUNT; i++) {
			predict_calc(this, qth, julian_date(t));
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

Sat_t.prototype.mapRefresh = function() {
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

		var geo = new Geodetic_t();

		/* Range circle calculations.
		 * Borrowed from gsat 0.9.0 by Xavier Crehueras, EB3CZS
		 * who borrowed from John Magliacane, KD2BD.
		 * Optimized by Alexandru Csete and William J Beksi.
		 */
		ssplat = radians(this.ssplat);
		ssplon = radians(this.ssplon);
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

			points.push([degrees(geo.lat), degrees(geo.lon)]);
		}

		this.mapFootprint.setLatLngs(points);
	}
};

// ////////////////////////////////////////////////////////////////
// Periodic refresh loop
// ////////////////////////////////////////////////////////////////

function refresh() {
	var jt = julian_date(new Date());

	// Compute current satellites longitude and latitude
	for (var i = 0; i < sats.length; i++) {
		predict_calc(sats[i], qth, jt);
		sats[i].mapRefresh();
	}
}
