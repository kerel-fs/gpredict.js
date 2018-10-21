/* exported Calculate_RADec_and_Obs */
/* Procedure Calculate_User_PosVel passes the user's geodetic position */
/* and the time of interest and returns the ECI position and velocity  */
/* of the observer. The velocity calculation assumes the geodetic      */
/* position is stationary relative to the earth's surface.             */
function Calculate_User_PosVel(_time, geodetic, obs_pos, obs_vel) {
	var c, sq, achcp;

	geodetic.theta = (ThetaG_JD(_time) + geodetic.lon) % twopi;
	c = 1 / Math.sqrt(1 + __f * (__f - 2) * Sqr(sin(geodetic.lat)));
	sq = Sqr(1 - __f) * c;
	achcp = (xkmper * c + geodetic.alt) * Math.cos(geodetic.lat);
	obs_pos.x = achcp * Math.cos(geodetic.theta);
	obs_pos.y = achcp * Math.sin(geodetic.theta);
	obs_pos.z = (xkmper * sq + geodetic.alt) * Math.sin(geodetic.lat);
	obs_vel.x = -mfactor * obs_pos.y;
	obs_vel.y = mfactor * obs_pos.x;
	obs_vel.z = 0;
	Magnitude(obs_pos);
	Magnitude(obs_vel);
}

/*------------------------------------------------------------------*/

/* Procedure Calculate_LatLonAlt will calculate the geodetic  */
/* position of an object given its ECI position pos and time. */
/* It is intended to be used to determine the ground track of */
/* a satellite.  The calculations  assume the earth to be an  */
/* oblate spheroid as defined in WGS '72.                     */
function Calculate_LatLonAlt(_time, pos, geodetic) {
	var r, e2, phi, c;

	geodetic.theta = AcTan(pos.y, pos.x);
	geodetic.lon = FMod2p(geodetic.theta - ThetaG_JD(_time));
	r = sqrt(Sqr(pos.x) + Sqr(pos.y));
	e2 = __f * (2 - __f);
	geodetic.lat = AcTan(pos.z, r);

	do {
		phi = geodetic.lat;
		c = 1 / sqrt(1 - e2 * Sqr(sin(phi)));
		geodetic.lat = AcTan(pos.z + xkmper * c * e2 * Math.sin(phi), r);
	} while (fabs(geodetic.lat - phi) >= 1E-10);

	geodetic.alt = r / Math.cos(geodetic.lat) - xkmper * c;

	if (geodetic.lat > pio2)
		geodetic.lat -= twopi;

}

/*------------------------------------------------------------------*/

/* The procedures Calculate_Obs and Calculate_RADec calculate         */
/* the *topocentric* coordinates of the object with ECI position,     */
/* {pos}, and velocity, {vel}, from location {geodetic} at {time}.    */
/* The {obs_set} returned for Calculate_Obs consists of azimuth,      */
/* elevation, range, and range rate (in that order) with units of     */
/* radians, radians, kilometers, and kilometers/second, respectively. */
/* The WGS '72 geoid is used and the effect of atmospheric refraction */
/* (under standard temperature and pressure) is incorporated into the */
/* elevation calculation; the effect of atmospheric refraction on     */
/* range and range rate has not yet been quantified.                  */

/* The {obs_set} for Calculate_RADec consists of right ascension and  */
/* declination (in that order) in radians.  Again, calculations are   */
/* based on *topocentric* position using the WGS '72 geoid and        */
/* incorporating atmospheric refraction.                              */

/**
 * 
 */
function Calculate_Obs(_time, pos, vel, geodetic, obs_set) {
	var sin_lat, cos_lat, sin_theta, cos_theta, el, azim, top_s, top_e, top_z;

	var obs_pos = new vector_t();
	var obs_vel = new vector_t();
	var range = new vector_t();
	var rgvel = new vector_t();

	Calculate_User_PosVel(_time, geodetic, obs_pos, obs_vel);

	range.x = pos.x - obs_pos.x;
	range.y = pos.y - obs_pos.y;
	range.z = pos.z - obs_pos.z;

	rgvel.x = vel.x - obs_vel.x;
	rgvel.y = vel.y - obs_vel.y;
	rgvel.z = vel.z - obs_vel.z;

	Magnitude(range);

	sin_lat = Math.sin(geodetic.lat);
	cos_lat = Math.cos(geodetic.lat);
	sin_theta = Math.sin(geodetic.theta);
	cos_theta = Math.cos(geodetic.theta);
	top_s = sin_lat * cos_theta * range.x + sin_lat * sin_theta * range.y
			- cos_lat * range.z;
	top_e = -sin_theta * range.x + cos_theta * range.y;
	top_z = cos_lat * cos_theta * range.x + cos_lat * sin_theta * range.y
			+ sin_lat * range.z;
	azim = atan(-top_e / top_s);
	if (top_s > 0)
		azim = azim + pi;
	if (azim < 0)
		azim = azim + twopi;
	el = ArcSin(top_z / range.w);
	obs_set.az = azim;
	obs_set.el = el;
	obs_set.range = range.w;

	obs_set.range_rate = Dot(range, rgvel) / range.w;

	if (obs_set.el >= 0)
		SetFlag(VISIBLE_FLAG);
	else {
		obs_set.el = el;
		ClearFlag(VISIBLE_FLAG);
	}
}

/*------------------------------------------------------------------*/

function Calculate_RADec_and_Obs(_time, pos, vel, geodetic, obs_set) {
	var phi, theta, sin_theta, cos_theta, sin_phi, cos_phi, az, el, Lxh, Lyh;
	var Lzh, Sx, Ex, Zx, Sy, Ey, Zy, Sz, Ez, Zz, Lx, Ly, Lz, cos_delta;
	var sin_alpha, cos_alpha;

	var obs = new obs_set_t();

	Calculate_Obs(_time, pos, vel, geodetic, obs);

	az = obs.az;
	el = obs.el;
	phi = geodetic.lat;
	theta = FMod2p(ThetaG_JD(_time) + geodetic.lon);
	sin_theta = Math.sin(theta);
	cos_theta = Math.cos(theta);
	sin_phi = Math.sin(phi);
	cos_phi = Math.cos(phi);
	Lxh = -cos(az) * Math.cos(el);
	Lyh = Math.sin(az) * Math.cos(el);
	Lzh = Math.sin(el);
	Sx = sin_phi * cos_theta;
	Ex = -sin_theta;
	Zx = cos_theta * cos_phi;
	Sy = sin_phi * sin_theta;
	Ey = cos_theta;
	Zy = sin_theta * cos_phi;
	Sz = -cos_phi;
	Ez = 0;
	Zz = sin_phi;
	Lx = Sx * Lxh + Ex * Lyh + Zx * Lzh;
	Ly = Sy * Lxh + Ey * Lyh + Zy * Lzh;
	Lz = Sz * Lxh + Ez * Lyh + Zz * Lzh;
	obs_set.dec = ArcSin(Lz);
	cos_delta = sqrt(1 - Sqr(Lz));
	sin_alpha = Ly / cos_delta;
	cos_alpha = Lx / cos_delta;
	obs_set.ra = AcTan(sin_alpha, cos_alpha);
	obs_set.ra = FMod2p(obs_set.ra);
}
