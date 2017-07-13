"use strict";

function get_orbit_type(sat) {
	var orbit = ORBIT_TYPE_UNKNOWN;

	if (geostationary(sat)) {
		orbit = ORBIT_TYPE_GEO;
	} else if (decayed(sat)) {
		orbit = ORBIT_TYPE_DECAYED;
	} else {
		orbit = ORBIT_TYPE_UNKNOWN;
	}

	return orbit;
}

function geostationary(sat) {
	if (fabs(sat.meanmo - 1.0027) < 0.0002)
		return true;
	else
		return false;
}

function decayed(sat) {
	if (sat.jul_epoch
			+ ((16.666666 - sat.meanmo) / (10.0 * fabs(sat.tle.xndt2o
					/ (twopi / xmnpda / xmnpda)))) < sat.jul_utc)
		return true;
	else
		return false;
}

/** 
 * Determine whether satellite ever reaches AOS.
 * 
 * @author John A. Magliacane, KD2BD
 * @author Alexandru Csete, OZ9AEC
 * 
 * @param {Sat_t} sat Pointer to satellite data.
 * @param {qth_t} qth
 * 
 * @return {number} true if the satellite will reach AOS, false otherwise.
 */
function has_aos(sat, qth)
{
     var lin = 0.0, sma = 0.0, apogee = 0.0;
     var retcode = false;

     /* FIXME */
     if (sat.meanmo == 0.0) {
          retcode = false;
     }
     else {

          /* xincl is already in RAD by select_ephemeris */
          lin = sat.tle.xincl;
          if (lin >= pio2)
               lin = pi - lin;

          sma = 331.25 * Math.exp(Math.log(1440.0/sat.meanmo) * (2.0/3.0));
          apogee = sma * (1.0 + sat.tle.eo) - xkmper;

          if ((acos(xkmper/(apogee+xkmper))+(lin)) > fabs(qth.lat*de2ra))
               retcode = true;
          else
               retcode = false;

     }

     return retcode;
}
