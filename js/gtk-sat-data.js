
/** \brief Read TLE data for a given satellite into memory.
 *  \param catnum The catalog number of the satellite.
 *  \param sat Pointer to a valid Sat_t structure.
 *  \return 0 if successfull, 1 if an I/O error occurred,
 *          2 if the TLE data appears to be bad.
 *
 */
function gtk_sat_data_read_sat(rawtle, sat) {
	Convert_Satellite_Data(rawtle, sat.tle);

	sat.name = sat.tle.sat_name.trim();
	sat.nickname = sat.name;
	sat.flags = 0;

	select_ephemeris(sat);

	/* initialise variable fields */
	sat.jul_utc = 0.0;
	sat.tsince = 0.0;
	sat.az = 0.0;
	sat.el = 0.0;
	sat.range = 0.0;
	sat.range_rate = 0.0;
	sat.ra = 0.0;
	sat.dec = 0.0;
	sat.ssplat = 0.0;
	sat.ssplon = 0.0;
	sat.alt = 0.0;
	sat.velo = 0.0;
	sat.ma = 0.0;
	sat.footprint = 0.0;
	sat.phase = 0.0;
	sat.aos = 0.0;
	sat.los = 0.0;

	/* calculate satellite data at epoch */
	gtk_sat_data_init_sat(sat);
}

/** \brief Initialise satellite data.
 *  \param sat The satellite to initialise.
 *  \param qth Optional QTH info, use (0,0) if NULL.
 *
 * This function calculates the satellite data at t = 0, ie. epoch time
 * The function is called automatically by gtk_sat_data_read_sat.
 */
function gtk_sat_data_init_sat(sat) {
	var sat_geodetic = new geodetic_t();
	var jul_utc, age;

	jul_utc = Julian_Date_of_Epoch(sat.tle.epoch); // => tsince = 0.0
	sat.jul_epoch = jul_utc;

	/* execute computations */
	if (sat.flags & DEEP_SPACE_EPHEM_FLAG)
		SDP4(sat, 0.0);
	else
		SGP4(sat, 0.0);

	/* scale position and velocity to km and km/sec */
	Convert_Sat_State(sat.pos, sat.vel);

	/* get the velocity of the satellite */
	Magnitude(sat.vel);
	sat.velo = sat.vel.w;
	Calculate_LatLonAlt(jul_utc, sat.pos, sat_geodetic);

	while (sat_geodetic.lon < -pi)
		sat_geodetic.lon += twopi;

	while (sat_geodetic.lon > (pi))
		sat_geodetic.lon -= twopi;

	sat.ssplat = Degrees(sat_geodetic.lat);
	sat.ssplon = Degrees(sat_geodetic.lon);
	sat.alt = sat_geodetic.alt;
	sat.ma = Degrees(sat.phase);
	sat.ma *= 256.0 / 360.0;
	sat.footprint = 2.0 * xkmper * acos(xkmper / sat.pos.w);
	age = 0.0;
	sat.orbit = floor((sat.tle.xno * xmnpda / twopi + age * sat.tle.bstar * ae)
			* age + sat.tle.xmo / twopi)
			+ sat.tle.revnum - 1;

	/* orbit type */
	sat.otype = get_orbit_type(sat);
}
