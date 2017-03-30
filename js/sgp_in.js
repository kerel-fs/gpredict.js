/* Converts the strings in a raw two-line element set  */
/* to their intended numerical values. No processing   */
/* of these values is done, e.g. from deg to rads etc. */
/* This is done in the select_ephemeris() function.    */

function Convert_Satellite_Data(tle_set, tle) {
	tle.sat_name = tle_set[0];
	tle.catnr = tle_set[1].substr(2, 5);
	tle.idesg = tle_set[1].substr(9, 8);
	tle.epoch = tle_set[1].substr(18, 14);
	tle.epoch_year = tle_set[1].substr(18, 2);
	tle.epoch_year += 2000;
	tle.epoch_day = tle_set[1].substr(20, 3);
	tle.epoch_fod = tle_set[1].substr(23, 9);
	tle.xndt2o = tle_set[1].substr(33, 10);
	tle.xndd6o = tle_set[1].substr(44, 1) + "." + tle_set[1].substr(45, 5)
			+ "E" + tle_set[1].substr(50, 2);
	tle.bstar = tle_set[1].substr(53, 1) + "." + tle_set[1].substr(54, 5) + "E"
			+ tle_set[1].substr(59, 2);
	tle.elset = tle_set[1].substr(64, 4);
	tle.xincl = tle_set[2].substr(8, 8);
	tle.xnodeo = tle_set[2].substr(17, 8);
	tle.eo = "." + tle_set[2].substr(26, 7);
	tle.omegao = tle_set[2].substr(34, 8);
	tle.xmo = tle_set[2].substr(43, 8);
	tle.xno = tle_set[2].substr(52, 10);
	tle.revnum = tle_set[2].substr(63, 5);
}

/*------------------------------------------------------------------*/

/* Selects the apropriate ephemeris type to be used */
/* for predictions according to the data in the TLE */
/* It also processes values in the tle set so that  */
/* they are apropriate for the sgp4/sdp4 routines   */
function select_ephemeris(sat) {
	var ao, xnodp, dd1, dd2, delo, temp, a1, del1, r1;

	sat.tle.xnodeo *= de2ra;
	sat.tle.omegao *= de2ra;
	sat.tle.xmo *= de2ra;
	sat.tle.xincl *= de2ra;
	temp = twopi / xmnpda / xmnpda;

	sat.meanmo = sat.tle.xno;
	sat.tle.xno = sat.tle.xno * temp * xmnpda;
	sat.tle.xndt2o *= temp;
	sat.tle.xndd6o = sat.tle.xndd6o * temp / xmnpda;
	sat.tle.bstar /= ae;

	dd1 = (xke / sat.tle.xno);
	dd2 = tothrd;
	a1 = pow(dd1, dd2);
	r1 = cos(sat.tle.xincl);
	dd1 = (1.0 - sat.tle.eo * sat.tle.eo);
	temp = ck2 * 1.5 * (r1 * r1 * 3.0 - 1.0) / pow(dd1, 1.5);
	del1 = temp / (a1 * a1);
	ao = a1
			* (1.0 - del1
					* (tothrd * 0.5 + del1 * (del1 * 1.654320987654321 + 1.0)));
	delo = temp / (ao * ao);
	xnodp = sat.tle.xno / (delo + 1.0);

	if (twopi / xnodp / xmnpda >= .15625) {
		sat.flags |= DEEP_SPACE_EPHEM_FLAG;
	} else {
		sat.flags &= ~DEEP_SPACE_EPHEM_FLAG;
	}
}
