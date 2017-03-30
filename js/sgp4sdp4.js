
/** Type definitions **/

// orbit_type_t
var ORBIT_TYPE_UNKNOWN = 0;
var ORBIT_TYPE_LEO = 1;	/* !< Low Earth orbit, up to 1200 km. */
var ORBIT_TYPE_ICO = 2;	/* !< Intermediate Circular Orbit, up to 1400 km. */ 
var ORBIT_TYPE_GEO = 3;	/* !< Geostationary. */
var ORBIT_TYPE_GSO = 4;	/* !< Geosynchronuous. */
var ORBIT_TYPE_MOLNIYA = 5;
var ORBIT_TYPE_TUNDRA = 6;
var ORBIT_TYPE_POLAR = 7;
var ORBIT_TYPE_SUNSYNC = 8;
var ORBIT_TYPE_DECAYED = 9;

/** \brief Operational status of satellite. */
// op_stat_t
var OP_STAT_UNKNOWN = 0,
	OP_STAT_OPERATIONAL = 1,/* !< Operational [+] */
	OP_STAT_NONOP = 2,		/* !< Nonoperational [-] */
	OP_STAT_PARTIAL = 3,	/* !< Partially operational [P] */
	OP_STAT_STDBY = 4,		/* !< Backup/Standby [B] */
	OP_STAT_SPARE = 5,		/* !< Spare [S] */
	OP_STAT_EXTENDED = 6;	/* !< Extended Mission [X] */

/** \brief Two-line-element satellite orbital data.
 *  \ingroup sgpsdpif
 *  \bug doc incomplete.
 */
function tle_t() {
    this.epoch = 0.0;            /*!< Epoch Time in NORAD TLE format YYDDD.FFFFFFFF */
    this.epoch_year = 0; /*!< Epoch: year */
    this.epoch_day = 0;  /*!< Epoch: day of year */
    this.epoch_fod = 0.0;        /*!< Epoch: Fraction of day. */
    this.xndt2o = 0.0;           /*!< 1. time derivative of mean motion */
    this.xndd6o = 0.0;           /*!< 2. time derivative of mean motion */
    this.bstar = 0.0;            /*!< Bstar drag coefficient. */
    this.xincl = 0.0;            /*!< Inclination */
    this.xnodeo = 0.0;           /*!< R.A.A.N. */
    this.eo = 0.0;               /*!< Eccentricity */
    this.omegao = 0.0;           /*!< argument of perigee */
    this.xmo = 0.0;              /*!< mean anomaly */
    this.xno = 0.0;              /*!< mean motion */

    this.catnr = 0;            /*!< Catalogue Number.  */
    this.elset = 0;            /*!< Element Set number. */
    this.revnum = 0;           /*!< Revolution Number at epoch. */

    this.sat_name = "";     /*!< Satellite name string. */
    this.idesg = "";         /*!< International Designator. */
    this.status = OP_STAT_UNKNOWN;        /*!< Operational status. */

	/* values needed for squint calculations */
    this.xincl1 = 0.0;
    this.xnodeo1 = 0.0;
    this.omegao1 = 0.0;
}

/** \brief Geodetic position data structure.
 *  \ingroup sgpsdpif
 *
 * \bug It is uncertain whether the units are uniform across all functions.
 */
function geodetic_t() {
	this.lat = 0.0;    /*!< Lattitude [rad] */
	this.lon = 0.0;    /*!< Longitude [rad] */
	this.alt = 0.0;    /*!< Altitude [km]? */
	this.theta = 0.0;
}

/** \brief General three-dimensional vector structure.
 *  \ingroup sgpsdpif
 */
function vector_t() {
	this.x = 0.0;   /*!< X component */
	this.y = 0.0;   /*!< Y component */
	this.z = 0.0;   /*!< Z component */
	this.w = 0.0;   /*!< Magnitude */
}

/** \brief Bearing to satellite from observer
 *  \ingroup sgpsdpif
 */
function obs_set_t() {
	this.az = 0.0;            /*!< Azimuth [deg] */
	this.el = 0.0;            /*!< Elevation [deg] */
	this.range = 0.0;         /*!< Range [km] */
	this.range_rate = 0.0;    /*!< Velocity [km/sec] */
}

function obs_astro_t() {
	this.ra = 0.0;   /*!< Right Ascension [dec] */
	this.dec = 0.0;  /*!< Declination [dec] */
}

/* Common arguments between deep-space functions */
function deep_arg_t() {
	/* Used by dpinit part of Deep() */
	this.eosq = 0.0;
	this.sinio = 0.0;
	this.cosio = 0.0;
	this.betao = 0.0;
	this.aodp = 0.0;
	this.theta2 = 0.0;
	this.sing = 0.0;
	this.cosg = 0.0;
	this.betao2 = 0.0;
	this.xmdot = 0.0;
	this.omgdot = 0.0;
	this.xnodot = 0.0;
	this.xnodp = 0.0;

	/* Used by dpsec and dpper parts of Deep() */
	this.xll = 0.0;
	this.omgadf = 0.0;
	this.xnode = 0.0;
	this.em = 0.0;
	this.xinc = 0.0;
	this.xn = 0.0;
	this.t = 0.0;

	/* Used by thetg and Deep() */
	this.ds50 = 0.0;
}

/* static data for SGP4 and SDP4 */
function sgpsdp_static_t() {
	this.aodp = 0.0; this.aycof= 0.0; this.c1 = 0.0; this.c4 = 0.0; this.c5 = 0.0; this.cosio = 0.0; this.d2 = 0.0; this.d3 = 0.0; this.d4 = 0.0; this.delmo = 0.0; this.omgcof = 0.0;
	this.eta = 0.0; this.omgdot= 0.0; this.sinio= 0.0; this.xnodp= 0.0; this.sinmo= 0.0; this.t2cof= 0.0; this.t3cof= 0.0; this.t4cof= 0.0; this.t5cof = 0.0;
	this.x1mth2 = 0.0; this.x3thm1= 0.0; this.x7thm1= 0.0; this.xmcof= 0.0; this.xmdot= 0.0; this.xnodcf= 0.0; this.xnodot= 0.0; this.xlcof = 0.0;
}

/* static data for DEEP */
function deep_static_t() {
	this.thgr = 0.0; this.xnq = 0.0; this.xqncl = 0.0; this.omegaq = 0.0; this.zmol = 0.0; this.zmos = 0.0; this.savtsn = 0.0; this.ee2 = 0.0; this.e3 = 0.0; this.xi2 = 0.0;
	this.xl2 = 0.0; this.xl3 = 0.0; this.xl4 = 0.0; this.xgh2 = 0.0; this.xgh3 = 0.0; this.xgh4 = 0.0; this.xh2 = 0.0; this.xh3 = 0.0; this.sse = 0.0; this.ssi = 0.0; this.ssg = 0.0; this.xi3 = 0.0;
	this.se2 = 0.0; this.si2 = 0.0; this.sl2 = 0.0; this.sgh2 = 0.0; this.sh2 = 0.0; this.se3 = 0.0; this.si3 = 0.0; this.sl3 = 0.0; this.sgh3 = 0.0; this.sh3 = 0.0; this.sl4 = 0.0; this.sgh4 = 0.0;
	this.ssl = 0.0; this.ssh = 0.0; this.d3210 = 0.0; this.d3222 = 0.0; this.d4410 = 0.0; this.d4422 = 0.0; this.d5220 = 0.0; this.d5232 = 0.0; this.d5421 = 0.0;
	this.d5433 = 0.0; this.del1 = 0.0; this.del2 = 0.0; this.del3 = 0.0; this.fasx2 = 0.0; this.fasx4 = 0.0; this.fasx6 = 0.0; this.xlamo = 0.0; this.xfact = 0.0;
	this.xni = 0.0; this.atime = 0.0; this.stepp = 0.0; this.stepn = 0.0; this.step2 = 0.0; this.preep = 0.0; this.pl = 0.0; this.sghs = 0.0; this.xli = 0.0;
	this.d2201 = 0.0; this.d2211 = 0.0; this.sghl = 0.0; this.sh1 = 0.0; this.pinc = 0.0; this.pe = 0.0; this.shs = 0.0; this.zsingl = 0.0; this.zcosgl = 0.0;
	this.zsinhl = 0.0; this.zcoshl = 0.0; this.zsinil = 0.0; this.zcosil = 0.0;
}

/** \brief Satellite data structure
 *  \ingroup sgpsdpif
 *
 */
function sat_t() {
    this.name = "";
    this.nickname = "";
    this.website = "";
    this.tle = new tle_t();     /*!< Keplerian elements */
    this.flags = 0;   /*!< Flags for algo ctrl */
    this.sgps = new sgpsdp_static_t();
    this.dps = new deep_static_t();
    this.deep_arg = new deep_arg_t();
    this.pos = new vector_t();       /*!< Raw position and range */
    this.vel = new vector_t();       /*!< Raw velocity */

	/*** FIXME: REMOVE */
    this.bearing = new obs_set_t();   /*!< Az, El, range and vel */
    this.astro = new obs_astro_t();     /*!< Ra and Decl */
	/*** END */

	/* time keeping fields */
    this.jul_epoch = 0.0;
    this.jul_utc = 0.0;
    this.tsince = 0.0;
    this.aos = 0.0;         /*!< Next AOS. */
    this.los = 0.0;         /*!< Next LOS */

    this.az = 0.0;          /*!< Azimuth [deg] */
    this.el = 0.0;          /*!< Elevation [deg] */
    this.range = 0.0;       /*!< Range [km] */
    this.range_rate = 0.0;  /*!< Range Rate [km/sec] */
    this.ra = 0.0;          /*!< Right Ascension [deg] */
    this.dec = 0.0;         /*!< Declination [deg] */
    this.ssplat = 0.0;      /*!< SSP latitude [deg] */
    this.ssplon = 0.0;      /*!< SSP longitude [deg] */
    this.alt = 0.0;         /*!< altitude [km] */
    this.velo = 0.0;        /*!< velocity [km/s] */
    this.ma = 0.0;          /*!< mean anomaly */
    this.footprint = 0.0;   /*!< footprint */
    this.phase = 0.0;       /*!< orbit phase */
    this.meanmo = 0.0;      /*!< mean motion kept in rev/day */
    this.orbit = 0;       /*!< orbit number */
    this.otype = ORBIT_TYPE_UNKNOWN;       /*!< orbit type. */
}

/** Table of constant values **/
var de2ra    =1.74532925E-2;   /* Degrees to Radians */
var pi       =3.1415926535898; /* Pi */
var pio2     =1.5707963267949; /* Pi/2 */
var x3pio2   =4.71238898;      /* 3*Pi/2 */
var twopi    =6.2831853071796; /* 2*Pi  */
var e6a      =1.0E-6;
var tothrd   =6.6666667E-1;    /* 2/3 */
var xj2      =1.0826158E-3;    /* J2 Harmonic */
var xj3      =-2.53881E-6;      /* J3 Harmonic */   
var xj4      =-1.65597E-6;      /* J4 Harmonic */
var xke      =7.43669161E-2;
var xkmper   =6.378135E3;      /* Earth radius km */
var xmnpda   =1.44E3;          /* Minutes per day */
var ae       =1.0;
var ck2      =5.413079E-4;
var ck4      =6.209887E-7;
var __f      =3.352779E-3;
var ge       =3.986008E5; 
var __s__    =1.012229;
var qoms2t   =1.880279E-09;
var secday   =8.6400E4;        /* Seconds per day */
var omega_E  =1.0027379;
var omega_ER =6.3003879;
var zns      =1.19459E-5;
var c1ss     =2.9864797E-6;
var zes      =1.675E-2;
var znl      =1.5835218E-4;
var c1l      =4.7968065E-7;
var zel      =5.490E-2;
var zcosis   =9.1744867E-1;
var zsinis   =3.9785416E-1;
var zsings   =-9.8088458E-1;
var zcosgs   =1.945905E-1;
var zcoshs   =1;
var zsinhs   =0;
var q22      =1.7891679E-6;
var q31      =2.1460748E-6;
var q33      =2.2123015E-7;
var g22      =5.7686396;
var g32      =9.5240898E-1;
var g44      =1.8014998;
var g52      =1.0508330;
var g54      =4.4108898;
var root22   =1.7891679E-6;
var root32   =3.7393792E-7;
var root44   =7.3636953E-9;
var root52   =1.1428639E-7;
var root54   =2.1765803E-9;
var thdt     =4.3752691E-3;
var rho      =1.5696615E-1;
var mfactor  =7.292115E-5;
var __sr__   =6.96000E5;      /*Solar radius - kilometers (IAU 76)*/
var AU       =1.49597870E8;   /*Astronomical unit - kilometers (IAU 76)*/

/* Entry points of Deep() 
FIXME: Change to enu */
var dpinit   =1; /* Deep-space initialization code */
var dpsec    =2; /* Deep-space secular code        */
var dpper    =3; /* Deep-space periodic code       */

/* Carriage return and line feed */
var CR  =0x0A;
var LF  =0x0D;

/* Flow control flag definitions */
var ALL_FLAGS              =-1;
var SGP_INITIALIZED_FLAG   =0x000001;
var SGP4_INITIALIZED_FLAG  =0x000002;
var SDP4_INITIALIZED_FLAG  =0x000004;
var SGP8_INITIALIZED_FLAG  =0x000008;
var SDP8_INITIALIZED_FLAG  =0x000010;
var SIMPLE_FLAG            =0x000020;
var DEEP_SPACE_EPHEM_FLAG  =0x000040;
var LUNAR_TERMS_DONE_FLAG  =0x000080;
var NEW_EPHEMERIS_FLAG     =0x000100;
var DO_LOOP_FLAG           =0x000200;
var RESONANCE_FLAG         =0x000400;
var SYNCHRONOUS_FLAG       =0x000800;
var EPOCH_RESTART_FLAG     =0x001000;
var VISIBLE_FLAG           =0x002000;
var SAT_ECLIPSED_FLAG      =0x004000;


/* SGP4 */
/* This function is used to calculate the position and velocity */
/* of near-earth (period < 225 minutes) satellites. tsince is   */
/* time since epoch in minutes, tle is a pointer to a tle_t     */
/* structure with Keplerian orbital elements and pos and vel    */
/* are vector_t structures returning ECI satellite position and */
/* velocity. Use Convert_Sat_State() to convert to km and km/s.*/
function SGP4(sat, tsince) {
	var cosuk, sinuk, rfdotk, vx, vy, vz, ux, uy, uz, xmy, xmx, cosnok, sinnok;
	var cosik, sinik, rdotk, xinck, xnodek, uk, rk, cos2u, sin2u, u, sinu, cosu;
	var betal, rfdot, rdot, r, pl, elsq, esine, ecose, epw, cosepw, x1m5th;
	var xhdot1, tfour, sinepw, capu, ayn, xlt, aynl, xll, axn, xn, beta, xl, e;
	var a, tcube, delm, delomg, templ, tempe, tempa, xnode, tsq, xmp, omega;
	var xnoddf, omgadf, xmdf, a1, a3ovk2, ao, betao, betao2, c1sq, c2, c3, coef;
	var coef1, del1, delo, eeta, eosq, etasq, perige, pinvsq, psisq, qoms24, s4;
	var temp, temp1, temp2, temp3, temp4, temp5, temp6, theta2, theta4, tsi;

	var i;

	if (~sat.flags & SGP4_INITIALIZED_FLAG) {

		sat.flags |= SGP4_INITIALIZED_FLAG;

		a1 = pow(xke / sat.tle.xno, tothrd);
		sat.sgps.cosio = Math.cos(sat.tle.xincl);
		theta2 = sat.sgps.cosio * sat.sgps.cosio;
		sat.sgps.x3thm1 = 3 * theta2 - 1.0;
		eosq = sat.tle.eo * sat.tle.eo;
		betao2 = 1 - eosq;
		betao = Math.sqrt(betao2);
		del1 = 1.5 * ck2 * sat.sgps.x3thm1 / (a1 * a1 * betao * betao2);
		ao = a1
				* (1 - del1 * (0.5 * tothrd + del1 * (1 + 134.0 / 81.0 * del1)));
		delo = 1.5 * ck2 * sat.sgps.x3thm1 / (ao * ao * betao * betao2);
		sat.sgps.xnodp = sat.tle.xno / (1.0 + delo);
		sat.sgps.aodp = ao / (1.0 - delo);

		if ((sat.sgps.aodp * (1.0 - sat.tle.eo) / ae) < (220.0 / xkmper + ae))
			sat.flags |= SIMPLE_FLAG;
		else
			sat.flags &= ~SIMPLE_FLAG;

		s4 = __s__;
		qoms24 = qoms2t;
		perige = (sat.sgps.aodp * (1 - sat.tle.eo) - ae) * xkmper;
		if (perige < 156.0) {
			if (perige <= 98.0)
				s4 = 20.0;
			else
				s4 = perige - 78.0;
			qoms24 = pow((120.0 - s4) * ae / xkmper, 4);
			s4 = s4 / xkmper + ae;
		}
		pinvsq = 1.0 / (sat.sgps.aodp * sat.sgps.aodp * betao2 * betao2);
		tsi = 1.0 / (sat.sgps.aodp - s4);
		sat.sgps.eta = sat.sgps.aodp * sat.tle.eo * tsi;
		etasq = sat.sgps.eta * sat.sgps.eta;
		eeta = sat.tle.eo * sat.sgps.eta;
		psisq = fabs(1.0 - etasq);
		coef = qoms24 * pow(tsi, 4);
		coef1 = coef / pow(psisq, 3.5);
		c2 = coef1
				* sat.sgps.xnodp
				* (sat.sgps.aodp * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) + 0.75
						* ck2
						* tsi
						/ psisq
						* sat.sgps.x3thm1
						* (8.0 + 3.0 * etasq * (8 + etasq)));
		sat.sgps.c1 = c2 * sat.tle.bstar;
		sat.sgps.sinio = Math.sin(sat.tle.xincl);
		a3ovk2 = -xj3 / ck2 * pow(ae, 3);
		c3 = coef * tsi * a3ovk2 * sat.sgps.xnodp * ae * sat.sgps.sinio
				/ sat.tle.eo;
		sat.sgps.x1mth2 = 1.0 - theta2;
		sat.sgps.c4 = 2.0
				* sat.sgps.xnodp
				* coef1
				* sat.sgps.aodp
				* betao2
				* (sat.sgps.eta * (2.0 + 0.5 * etasq) + sat.tle.eo
						* (0.5 + 2.0 * etasq) - 2.0
						* ck2
						* tsi
						/ (sat.sgps.aodp * psisq)
						* (-3.0
								* sat.sgps.x3thm1
								* (1.0 - 2.0 * eeta + etasq
										* (1.5 - 0.5 * eeta)) + 0.75
								* sat.sgps.x1mth2
								* (2.0 * etasq - eeta * (1.0 + etasq))
								* Math.cos(2.0 * sat.tle.omegao)));
		sat.sgps.c5 = 2.0 * coef1 * sat.sgps.aodp * betao2
				* (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
		theta4 = theta2 * theta2;
		temp1 = 3.0 * ck2 * pinvsq * sat.sgps.xnodp;
		temp2 = temp1 * ck2 * pinvsq;
		temp3 = 1.25 * ck4 * pinvsq * pinvsq * sat.sgps.xnodp;
		sat.sgps.xmdot = sat.sgps.xnodp + 0.5 * temp1 * betao * sat.sgps.x3thm1
				+ 0.0625 * temp2 * betao
				* (13.0 - 78.0 * theta2 + 137.0 * theta4);
		x1m5th = 1.0 - 5.0 * theta2;
		sat.sgps.omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2
				* (7.0 - 114.0 * theta2 + 395.0 * theta4) + temp3
				* (3.0 - 36.0 * theta2 + 49.0 * theta4);
		xhdot1 = -temp1 * sat.sgps.cosio;
		sat.sgps.xnodot = xhdot1
				+ (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3
						* (3.0 - 7.0 * theta2)) * sat.sgps.cosio;
		sat.sgps.omgcof = sat.tle.bstar * c3 * Math.cos(sat.tle.omegao);
		sat.sgps.xmcof = -tothrd * coef * sat.tle.bstar * ae / eeta;
		sat.sgps.xnodcf = 3.5 * betao2 * xhdot1 * sat.sgps.c1;
		sat.sgps.t2cof = 1.5 * sat.sgps.c1;
		sat.sgps.xlcof = 0.125 * a3ovk2 * sat.sgps.sinio
				* (3.0 + 5.0 * sat.sgps.cosio) / (1.0 + sat.sgps.cosio);
		sat.sgps.aycof = 0.25 * a3ovk2 * sat.sgps.sinio;
		sat.sgps.delmo = pow(1.0 + sat.sgps.eta * Math.cos(sat.tle.xmo), 3);
		sat.sgps.sinmo = Math.sin(sat.tle.xmo);
		sat.sgps.x7thm1 = 7.0 * theta2 - 1.0;
		if (~sat.flags & SIMPLE_FLAG) {
			c1sq = sat.sgps.c1 * sat.sgps.c1;
			sat.sgps.d2 = 4.0 * sat.sgps.aodp * tsi * c1sq;
			temp = sat.sgps.d2 * tsi * sat.sgps.c1 / 3.0;
			sat.sgps.d3 = (17.0 * sat.sgps.aodp + s4) * temp;
			sat.sgps.d4 = 0.5 * temp * sat.sgps.aodp * tsi
					* (221.0 * sat.sgps.aodp + 31.0 * s4) * sat.sgps.c1;
			sat.sgps.t3cof = sat.sgps.d2 + 2.0 * c1sq;
			sat.sgps.t4cof = 0.25 * (3.0 * sat.sgps.d3 + sat.sgps.c1
					* (12.0 * sat.sgps.d2 + 10.0 * c1sq));
			sat.sgps.t5cof = 0.2 * (3.0 * sat.sgps.d4 + 12.0 * sat.sgps.c1
					* sat.sgps.d3 + 6.0 * sat.sgps.d2 * sat.sgps.d2 + 15.0
					* c1sq * (2.0 * sat.sgps.d2 + c1sq));
		}
	}

	/* Update for secular gravity and atmospheric drag. */
	xmdf = sat.tle.xmo + sat.sgps.xmdot * tsince;
	omgadf = sat.tle.omegao + sat.sgps.omgdot * tsince;
	xnoddf = sat.tle.xnodeo + sat.sgps.xnodot * tsince;
	omega = omgadf;
	xmp = xmdf;
	tsq = tsince * tsince;
	xnode = xnoddf + sat.sgps.xnodcf * tsq;
	tempa = 1.0 - sat.sgps.c1 * tsince;
	tempe = sat.tle.bstar * sat.sgps.c4 * tsince;
	templ = sat.sgps.t2cof * tsq;
	if (~sat.flags & SIMPLE_FLAG) {
		delomg = sat.sgps.omgcof * tsince;
		delm = sat.sgps.xmcof
				* (pow(1 + sat.sgps.eta * Math.cos(xmdf), 3) - sat.sgps.delmo);
		temp = delomg + delm;
		xmp = xmdf + temp;
		omega = omgadf - temp;
		tcube = tsq * tsince;
		tfour = tsince * tcube;
		tempa = tempa - sat.sgps.d2 * tsq - sat.sgps.d3 * tcube - sat.sgps.d4
				* tfour;
		tempe = tempe + sat.tle.bstar * sat.sgps.c5
				* (sin(xmp) - sat.sgps.sinmo);
		templ = templ + sat.sgps.t3cof * tcube + tfour
				* (sat.sgps.t4cof + tsince * sat.sgps.t5cof);
	}

	a = sat.sgps.aodp * pow(tempa, 2);
	e = sat.tle.eo - tempe;
	xl = xmp + omega + xnode + sat.sgps.xnodp * templ;
	beta = Math.sqrt(1.0 - e * e);
	xn = xke / pow(a, 1.5);

	axn = e * Math.cos(omega);
	temp = 1.0 / (a * beta * beta);
	xll = temp * sat.sgps.xlcof * axn;
	aynl = temp * sat.sgps.aycof;
	xlt = xl + xll;
	ayn = e * Math.sin(omega) + aynl;

	capu = FMod2p(xlt - xnode);
	temp2 = capu;

	i = 0;
	do {
		sinepw = Math.sin(temp2);
		cosepw = Math.cos(temp2);
		temp3 = axn * sinepw;
		temp4 = ayn * cosepw;
		temp5 = axn * cosepw;
		temp6 = ayn * sinepw;
		epw = (capu - temp4 + temp3 - temp2) / (1.0 - temp5 - temp6) + temp2;
		if (fabs(epw - temp2) <= e6a)
			break;
		temp2 = epw;
	} while (i++ < 10);

	ecose = temp5 + temp6;
	esine = temp3 - temp4;
	elsq = axn * axn + ayn * ayn;
	temp = 1.0 - elsq;
	pl = a * temp;
	r = a * (1.0 - ecose);
	temp1 = 1.0 / r;
	rdot = xke * Math.sqrt(a) * esine * temp1;
	rfdot = xke * Math.sqrt(pl) * temp1;
	temp2 = a * temp1;
	betal = Math.sqrt(temp);
	temp3 = 1.0 / (1.0 + betal);
	cosu = temp2 * (cosepw - axn + ayn * esine * temp3);
	sinu = temp2 * (sinepw - ayn - axn * esine * temp3);
	u = AcTan(sinu, cosu);
	sin2u = 2.0 * sinu * cosu;
	cos2u = 2.0 * cosu * cosu - 1.0;
	temp = 1.0 / pl;
	temp1 = ck2 * temp;
	temp2 = temp1 * temp;

	rk = r * (1.0 - 1.5 * temp2 * betal * sat.sgps.x3thm1) + 0.5 * temp1
			* sat.sgps.x1mth2 * cos2u;
	uk = u - 0.25 * temp2 * sat.sgps.x7thm1 * sin2u;
	xnodek = xnode + 1.5 * temp2 * sat.sgps.cosio * sin2u;
	xinck = sat.tle.xincl + 1.5 * temp2 * sat.sgps.cosio * sat.sgps.sinio
			* cos2u;
	rdotk = rdot - xn * temp1 * sat.sgps.x1mth2 * sin2u;
	rfdotk = rfdot + xn * temp1
			* (sat.sgps.x1mth2 * cos2u + 1.5 * sat.sgps.x3thm1);

	sinuk = Math.sin(uk);
	cosuk = Math.cos(uk);
	sinik = Math.sin(xinck);
	cosik = Math.cos(xinck);
	sinnok = Math.sin(xnodek);
	cosnok = Math.cos(xnodek);
	xmx = -sinnok * cosik;
	xmy = cosnok * cosik;
	ux = xmx * sinuk + cosnok * cosuk;
	uy = xmy * sinuk + sinnok * cosuk;
	uz = sinik * sinuk;
	vx = xmx * cosuk - cosnok * sinuk;
	vy = xmy * cosuk - sinnok * sinuk;
	vz = sinik * cosuk;

	sat.pos.x = rk * ux;
	sat.pos.y = rk * uy;
	sat.pos.z = rk * uz;
	sat.vel.x = rdotk * ux + rfdotk * vx;
	sat.vel.y = rdotk * uy + rfdotk * vy;
	sat.vel.z = rdotk * uz + rfdotk * vz;

	sat.phase = xlt - xnode - omgadf + twopi;
	if (sat.phase < 0)
		sat.phase += twopi;
	sat.phase = FMod2p(sat.phase);

	sat.tle.omegao1 = omega;
	sat.tle.xincl1 = xinck;
	sat.tle.xnodeo1 = xnodek;

}

function SDP4(sat, tsince) {
	var i;

	var a, axn, ayn, aynl, beta, betal, capu, cos2u, cosepw, cosik, cosnok;
	var cosu, cosuk, ecose, elsq, epw, esine, pl, theta4, rdot, rdotk, rfdot;
	var rfdotk, rk, sin2u, sinepw, sinik, sinnok, sinu, sinuk, tempe, templ;
	var tsq, u, uk, ux, uy, uz, vx, vy, vz, xinck, xl, xlt, xmam, xmdf, xmx;
	var xmy, xnoddf, xnodek, xll, a1, a3ovk2, ao, c2, coef, coef1, x1m5th;
	var xhdot1, del1, r, delo, eeta, eta, etasq, perige, psisq, tsi, qoms24, s4;
	var pinvsq, temp, tempa, temp1, temp2, temp3, temp4, temp5, temp6;

	if (~sat.flags & SDP4_INITIALIZED_FLAG) {

		sat.flags |= SDP4_INITIALIZED_FLAG;

		a1 = pow(xke / sat.tle.xno, tothrd);
		sat.deep_arg.cosio = Math.cos(sat.tle.xincl);
		sat.deep_arg.theta2 = sat.deep_arg.cosio * sat.deep_arg.cosio;
		sat.sgps.x3thm1 = 3.0 * sat.deep_arg.theta2 - 1.0;
		sat.deep_arg.eosq = sat.tle.eo * sat.tle.eo;
		sat.deep_arg.betao2 = 1.0 - sat.deep_arg.eosq;
		sat.deep_arg.betao = Math.sqrt(sat.deep_arg.betao2);
		del1 = 1.5 * ck2 * sat.sgps.x3thm1
				/ (a1 * a1 * sat.deep_arg.betao * sat.deep_arg.betao2);
		ao = a1
				* (1.0 - del1
						* (0.5 * tothrd + del1 * (1.0 + 134.0 / 81.0 * del1)));
		delo = 1.5 * ck2 * sat.sgps.x3thm1
				/ (ao * ao * sat.deep_arg.betao * sat.deep_arg.betao2);
		sat.deep_arg.xnodp = sat.tle.xno / (1.0 + delo);
		sat.deep_arg.aodp = ao / (1.0 - delo);

		s4 = __s__;
		qoms24 = qoms2t;
		perige = (sat.deep_arg.aodp * (1.0 - sat.tle.eo) - ae) * xkmper;
		if (perige < 156.0) {
			if (perige <= 98.0)
				s4 = 20.0;
			else
				s4 = perige - 78.0;
			qoms24 = pow((120.0 - s4) * ae / xkmper, 4);
			s4 = s4 / xkmper + ae;
		}
		pinvsq = 1.0 / (sat.deep_arg.aodp * sat.deep_arg.aodp
				* sat.deep_arg.betao2 * sat.deep_arg.betao2);
		sat.deep_arg.sing = Math.sin(sat.tle.omegao);
		sat.deep_arg.cosg = Math.cos(sat.tle.omegao);
		tsi = 1.0 / (sat.deep_arg.aodp - s4);
		eta = sat.deep_arg.aodp * sat.tle.eo * tsi;
		etasq = eta * eta;
		eeta = sat.tle.eo * eta;
		psisq = fabs(1.0 - etasq);
		coef = qoms24 * pow(tsi, 4);
		coef1 = coef / pow(psisq, 3.5);
		c2 = coef1
				* sat.deep_arg.xnodp
				* (sat.deep_arg.aodp
						* (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) + 0.75
						* ck2 * tsi / psisq * sat.sgps.x3thm1
						* (8.0 + 3.0 * etasq * (8.0 + etasq)));
		sat.sgps.c1 = sat.tle.bstar * c2;
		sat.deep_arg.sinio = Math.sin(sat.tle.xincl);
		a3ovk2 = -xj3 / ck2 * pow(ae, 3);
		sat.sgps.x1mth2 = 1.0 - sat.deep_arg.theta2;
		sat.sgps.c4 = 2.0
				* sat.deep_arg.xnodp
				* coef1
				* sat.deep_arg.aodp
				* sat.deep_arg.betao2
				* (eta * (2.0 + 0.5 * etasq) + sat.tle.eo * (0.5 + 2.0 * etasq) - 2.0
						* ck2
						* tsi
						/ (sat.deep_arg.aodp * psisq)
						* (-3.0
								* sat.sgps.x3thm1
								* (1.0 - 2.0 * eeta + etasq
										* (1.5 - 0.5 * eeta)) + 0.75
								* sat.sgps.x1mth2
								* (2.0 * etasq - eeta * (1.0 + etasq))
								* Math.cos(2.0 * sat.tle.omegao)));
		theta4 = sat.deep_arg.theta2 * sat.deep_arg.theta2;
		temp1 = 3.0 * ck2 * pinvsq * sat.deep_arg.xnodp;
		temp2 = temp1 * ck2 * pinvsq;
		temp3 = 1.25 * ck4 * pinvsq * pinvsq * sat.deep_arg.xnodp;
		sat.deep_arg.xmdot = sat.deep_arg.xnodp + 0.5 * temp1
				* sat.deep_arg.betao * sat.sgps.x3thm1 + 0.0625 * temp2
				* sat.deep_arg.betao
				* (13.0 - 78.0 * sat.deep_arg.theta2 + 137.0 * theta4);
		x1m5th = 1.0 - 5.0 * sat.deep_arg.theta2;
		sat.deep_arg.omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2
				* (7.0 - 114.0 * sat.deep_arg.theta2 + 395.0 * theta4) + temp3
				* (3.0 - 36.0 * sat.deep_arg.theta2 + 49.0 * theta4);
		xhdot1 = -temp1 * sat.deep_arg.cosio;
		sat.deep_arg.xnodot = xhdot1
				+ (0.5 * temp2 * (4.0 - 19.0 * sat.deep_arg.theta2) + 2.0
						* temp3 * (3.0 - 7.0 * sat.deep_arg.theta2))
				* sat.deep_arg.cosio;
		sat.sgps.xnodcf = 3.5 * sat.deep_arg.betao2 * xhdot1 * sat.sgps.c1;
		sat.sgps.t2cof = 1.5 * sat.sgps.c1;
		sat.sgps.xlcof = 0.125 * a3ovk2 * sat.deep_arg.sinio
				* (3.0 + 5.0 * sat.deep_arg.cosio) / (1.0 + sat.deep_arg.cosio);
		sat.sgps.aycof = 0.25 * a3ovk2 * sat.deep_arg.sinio;
		sat.sgps.x7thm1 = 7.0 * sat.deep_arg.theta2 - 1.0;

		Deep(dpinit, sat);
	}

	xmdf = sat.tle.xmo + sat.deep_arg.xmdot * tsince;
	sat.deep_arg.omgadf = sat.tle.omegao + sat.deep_arg.omgdot * tsince;
	xnoddf = sat.tle.xnodeo + sat.deep_arg.xnodot * tsince;
	tsq = tsince * tsince;
	sat.deep_arg.xnode = xnoddf + sat.sgps.xnodcf * tsq;
	tempa = 1.0 - sat.sgps.c1 * tsince;
	tempe = sat.tle.bstar * sat.sgps.c4 * tsince;
	templ = sat.sgps.t2cof * tsq;
	sat.deep_arg.xn = sat.deep_arg.xnodp;

	sat.deep_arg.xll = xmdf;
	sat.deep_arg.t = tsince;

	Deep(dpsec, sat);

	xmdf = sat.deep_arg.xll;
	a = pow(xke / sat.deep_arg.xn, tothrd) * tempa * tempa;
	sat.deep_arg.em = sat.deep_arg.em - tempe;
	xmam = xmdf + sat.deep_arg.xnodp * templ;

	sat.deep_arg.xll = xmam;

	Deep(dpper, sat);

	xmam = sat.deep_arg.xll;
	xl = xmam + sat.deep_arg.omgadf + sat.deep_arg.xnode;
	beta = Math.sqrt(1.0 - sat.deep_arg.em * sat.deep_arg.em);
	sat.deep_arg.xn = xke / pow(a, 1.5);

	axn = sat.deep_arg.em * Math.cos(sat.deep_arg.omgadf);
	temp = 1.0 / (a * beta * beta);
	xll = temp * sat.sgps.xlcof * axn;
	aynl = temp * sat.sgps.aycof;
	xlt = xl + xll;
	ayn = sat.deep_arg.em * Math.sin(sat.deep_arg.omgadf) + aynl;

	capu = FMod2p(xlt - sat.deep_arg.xnode);
	temp2 = capu;

	i = 0;
	do {
		sinepw = Math.sin(temp2);
		cosepw = Math.cos(temp2);
		temp3 = axn * sinepw;
		temp4 = ayn * cosepw;
		temp5 = axn * cosepw;
		temp6 = ayn * sinepw;
		epw = (capu - temp4 + temp3 - temp2) / (1.0 - temp5 - temp6) + temp2;
		if (fabs(epw - temp2) <= e6a)
			break;
		temp2 = epw;
	} while (i++ < 10);

	ecose = temp5 + temp6;
	esine = temp3 - temp4;
	elsq = axn * axn + ayn * ayn;
	temp = 1.0 - elsq;
	pl = a * temp;
	r = a * (1.0 - ecose);
	temp1 = 1.0 / r;
	rdot = xke * Math.sqrt(a) * esine * temp1;
	rfdot = xke * Math.sqrt(pl) * temp1;
	temp2 = a * temp1;
	betal = Math.sqrt(temp);
	temp3 = 1.0 / (1.0 + betal);
	cosu = temp2 * (cosepw - axn + ayn * esine * temp3);
	sinu = temp2 * (sinepw - ayn - axn * esine * temp3);
	u = AcTan(sinu, cosu);
	sin2u = 2.0 * sinu * cosu;
	cos2u = 2.0 * cosu * cosu - 1.0;
	temp = 1.0 / pl;
	temp1 = ck2 * temp;
	temp2 = temp1 * temp;

	rk = r * (1.0 - 1.5 * temp2 * betal * sat.sgps.x3thm1) + 0.5 * temp1
			* sat.sgps.x1mth2 * cos2u;
	uk = u - 0.25 * temp2 * sat.sgps.x7thm1 * sin2u;
	xnodek = sat.deep_arg.xnode + 1.5 * temp2 * sat.deep_arg.cosio * sin2u;
	xinck = sat.deep_arg.xinc + 1.5 * temp2 * sat.deep_arg.cosio
			* sat.deep_arg.sinio * cos2u;
	rdotk = rdot - sat.deep_arg.xn * temp1 * sat.sgps.x1mth2 * sin2u;
	rfdotk = rfdot + sat.deep_arg.xn * temp1
			* (sat.sgps.x1mth2 * cos2u + 1.5 * sat.sgps.x3thm1);

	sinuk = Math.sin(uk);
	cosuk = Math.cos(uk);
	sinik = Math.sin(xinck);
	cosik = Math.cos(xinck);
	sinnok = Math.sin(xnodek);
	cosnok = Math.cos(xnodek);
	xmx = -sinnok * cosik;
	xmy = cosnok * cosik;
	ux = xmx * sinuk + cosnok * cosuk;
	uy = xmy * sinuk + sinnok * cosuk;
	uz = sinik * sinuk;
	vx = xmx * cosuk - cosnok * sinuk;
	vy = xmy * cosuk - sinnok * sinuk;
	vz = sinik * cosuk;

	sat.pos.x = rk * ux;
	sat.pos.y = rk * uy;
	sat.pos.z = rk * uz;
	sat.vel.x = rdotk * ux + rfdotk * vx;
	sat.vel.y = rdotk * uy + rfdotk * vy;
	sat.vel.z = rdotk * uz + rfdotk * vz;

	sat.phase = xlt - sat.deep_arg.xnode - sat.deep_arg.omgadf + twopi;
	if (sat.phase < 0.0)
		sat.phase += twopi;
	sat.phase = FMod2p(sat.phase);

	sat.tle.omegao1 = sat.deep_arg.omgadf;
	sat.tle.xincl1 = sat.deep_arg.xinc;
	sat.tle.xnodeo1 = sat.deep_arg.xnode;
}

function Deep(ientry, sat) {

	var a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, ainv2, alfdp, aqnv, sgh;
	var sini2, sinis, sinok, sh, si, sil, day, betdp, dalf, bfact, c, cc;
	var cosis, cosok, cosq, ctem, f322, zx, zy, dbet, dls, eoc, eq, f2, f220;
	var f221, f3, f311, f321, xnoh, f330, f441, f442, f522, f523, f542, f543;
	var g200, g201, g211, pgh, ph, s1, s2, s3, s4, s5, s6, s7, se, sel, ses;
	var xls, g300, g310, g322, g410, g422, g520, g521, g532, g533, gam, sinq;
	var sinzf, sis, sl, sll, sls, stem, temp, temp1, x1, x2, x2li, x2omi, x3;
	var x4, x5, x6, x7, x8, xl, xldot, xmao, xnddt, xndot, xno2, xnodce, xnoi;
	var xomi, xpidot, z1, z11, z12, z13, z2, z21, z22, z23, z3, z31, z32, z33;
	var ze, zf, zm, zn, zsing, zsinh, zsini, zcosg, zcosh, zcosi;
	var delt = 0.0, ft = 0.0;

	switch (ientry) {
	case dpinit:
		sat.dps.thgr = ThetaG(sat.tle.epoch, sat.deep_arg);
		eq = sat.tle.eo;
		sat.dps.xnq = sat.deep_arg.xnodp;
		aqnv = 1.0 / sat.deep_arg.aodp;
		sat.dps.xqncl = sat.tle.xincl;
		xmao = sat.tle.xmo;
		xpidot = sat.deep_arg.omgdot + sat.deep_arg.xnodot;
		sinq = Math.sin(sat.tle.xnodeo);
		cosq = Math.cos(sat.tle.xnodeo);
		sat.dps.omegaq = sat.tle.omegao;
		sat.dps.preep = 0;

		day = sat.deep_arg.ds50 + 18261.5;
		if (day != sat.dps.preep) {
			sat.dps.preep = day;
			xnodce = 4.5236020 - 9.2422029E-4 * day;
			stem = Math.sin(xnodce);
			ctem = Math.cos(xnodce);
			sat.dps.zcosil = 0.91375164 - 0.03568096 * ctem;
			sat.dps.zsinil = Math.sqrt(1.0 - sat.dps.zcosil * sat.dps.zcosil);
			sat.dps.zsinhl = 0.089683511 * stem / sat.dps.zsinil;
			sat.dps.zcoshl = Math.sqrt(1.0 - sat.dps.zsinhl * sat.dps.zsinhl);
			c = 4.7199672 + 0.22997150 * day;
			gam = 5.8351514 + 0.0019443680 * day;
			sat.dps.zmol = FMod2p(c - gam);
			zx = 0.39785416 * stem / sat.dps.zsinil;
			zy = sat.dps.zcoshl * ctem + 0.91744867 * sat.dps.zsinhl * stem;
			zx = AcTan(zx, zy);
			zx = gam + zx - xnodce;
			sat.dps.zcosgl = Math.cos(zx);
			sat.dps.zsingl = Math.sin(zx);
			sat.dps.zmos = 6.2565837 + 0.017201977 * day;
			sat.dps.zmos = FMod2p(sat.dps.zmos);
		}

		sat.dps.savtsn = 1E20;
		zcosg = zcosgs;
		zsing = zsings;
		zcosi = zcosis;
		zsini = zsinis;
		zcosh = cosq;
		zsinh = sinq;
		cc = c1ss;
		zn = zns;
		ze = zes;
		zmo = sat.dps.zmos;
		xnoi = 1.0 / sat.dps.xnq;

		for (;;) {
			a1 = zcosg * zcosh + zsing * zcosi * zsinh;
			a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
			a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
			a8 = zsing * zsini;
			a9 = zsing * zsinh + zcosg * zcosi * zcosh;
			a10 = zcosg * zsini;
			a2 = sat.deep_arg.cosio * a7 + sat.deep_arg.sinio * a8;
			a4 = sat.deep_arg.cosio * a9 + sat.deep_arg.sinio * a10;
			a5 = -sat.deep_arg.sinio * a7 + sat.deep_arg.cosio * a8;
			a6 = -sat.deep_arg.sinio * a9 + sat.deep_arg.cosio * a10;
			x1 = a1 * sat.deep_arg.cosg + a2 * sat.deep_arg.sing;
			x2 = a3 * sat.deep_arg.cosg + a4 * sat.deep_arg.sing;
			x3 = -a1 * sat.deep_arg.sing + a2 * sat.deep_arg.cosg;
			x4 = -a3 * sat.deep_arg.sing + a4 * sat.deep_arg.cosg;
			x5 = a5 * sat.deep_arg.sing;
			x6 = a6 * sat.deep_arg.sing;
			x7 = a5 * sat.deep_arg.cosg;
			x8 = a6 * sat.deep_arg.cosg;
			z31 = 12 * x1 * x1 - 3 * x3 * x3;
			z32 = 24 * x1 * x2 - 6 * x3 * x4;
			z33 = 12 * x2 * x2 - 3 * x4 * x4;
			z1 = 3 * (a1 * a1 + a2 * a2) + z31 * sat.deep_arg.eosq;
			z2 = 6 * (a1 * a3 + a2 * a4) + z32 * sat.deep_arg.eosq;
			z3 = 3 * (a3 * a3 + a4 * a4) + z33 * sat.deep_arg.eosq;
			z11 = -6 * a1 * a5 + sat.deep_arg.eosq
					* (-24 * x1 * x7 - 6 * x3 * x5);
			z12 = -6 * (a1 * a6 + a3 * a5) + sat.deep_arg.eosq
					* (-24 * (x2 * x7 + x1 * x8) - 6 * (x3 * x6 + x4 * x5));
			z13 = -6 * a3 * a6 + sat.deep_arg.eosq
					* (-24 * x2 * x8 - 6 * x4 * x6);
			z21 = 6 * a2 * a5 + sat.deep_arg.eosq
					* (24 * x1 * x5 - 6 * x3 * x7);
			z22 = 6 * (a4 * a5 + a2 * a6) + sat.deep_arg.eosq
					* (24 * (x2 * x5 + x1 * x6) - 6 * (x4 * x7 + x3 * x8));
			z23 = 6 * a4 * a6 + sat.deep_arg.eosq
					* (24 * x2 * x6 - 6 * x4 * x8);
			z1 = z1 + z1 + sat.deep_arg.betao2 * z31;
			z2 = z2 + z2 + sat.deep_arg.betao2 * z32;
			z3 = z3 + z3 + sat.deep_arg.betao2 * z33;
			s3 = cc * xnoi;
			s2 = -0.5 * s3 / sat.deep_arg.betao;
			s4 = s3 * sat.deep_arg.betao;
			s1 = -15 * eq * s4;
			s5 = x1 * x3 + x2 * x4;
			s6 = x2 * x3 + x1 * x4;
			s7 = x2 * x4 - x1 * x3;
			se = s1 * zn * s5;
			si = s2 * zn * (z11 + z13);
			sl = -zn * s3 * (z1 + z3 - 14 - 6 * sat.deep_arg.eosq);
			sgh = s4 * zn * (z31 + z33 - 6);
			sh = -zn * s2 * (z21 + z23);
			if (sat.dps.xqncl < 5.2359877E-2)
				sh = 0;
			sat.dps.ee2 = 2 * s1 * s6;
			sat.dps.e3 = 2 * s1 * s7;
			sat.dps.xi2 = 2 * s2 * z12;
			sat.dps.xi3 = 2 * s2 * (z13 - z11);
			sat.dps.xl2 = -2 * s3 * z2;
			sat.dps.xl3 = -2 * s3 * (z3 - z1);
			sat.dps.xl4 = -2 * s3 * (-21 - 9 * sat.deep_arg.eosq) * ze;
			sat.dps.xgh2 = 2 * s4 * z32;
			sat.dps.xgh3 = 2 * s4 * (z33 - z31);
			sat.dps.xgh4 = -18 * s4 * ze;
			sat.dps.xh2 = -2 * s2 * z22;
			sat.dps.xh3 = -2 * s2 * (z23 - z21);

			if (sat.flags & LUNAR_TERMS_DONE_FLAG)
				break;

			sat.dps.sse = se;
			sat.dps.ssi = si;
			sat.dps.ssl = sl;
			sat.dps.ssh = sh / sat.deep_arg.sinio;
			sat.dps.ssg = sgh - sat.deep_arg.cosio * sat.dps.ssh;
			sat.dps.se2 = sat.dps.ee2;
			sat.dps.si2 = sat.dps.xi2;
			sat.dps.sl2 = sat.dps.xl2;
			sat.dps.sgh2 = sat.dps.xgh2;
			sat.dps.sh2 = sat.dps.xh2;
			sat.dps.se3 = sat.dps.e3;
			sat.dps.si3 = sat.dps.xi3;
			sat.dps.sl3 = sat.dps.xl3;
			sat.dps.sgh3 = sat.dps.xgh3;
			sat.dps.sh3 = sat.dps.xh3;
			sat.dps.sl4 = sat.dps.xl4;
			sat.dps.sgh4 = sat.dps.xgh4;
			zcosg = sat.dps.zcosgl;
			zsing = sat.dps.zsingl;
			zcosi = sat.dps.zcosil;
			zsini = sat.dps.zsinil;
			zcosh = sat.dps.zcoshl * cosq + sat.dps.zsinhl * sinq;
			zsinh = sinq * sat.dps.zcoshl - cosq * sat.dps.zsinhl;
			zn = znl;
			cc = c1l;
			ze = zel;
			zmo = sat.dps.zmol;
			sat.flags |= LUNAR_TERMS_DONE_FLAG;
		}

		sat.dps.sse = sat.dps.sse + se;
		sat.dps.ssi = sat.dps.ssi + si;
		sat.dps.ssl = sat.dps.ssl + sl;
		sat.dps.ssg = sat.dps.ssg + sgh - sat.deep_arg.cosio
				/ sat.deep_arg.sinio * sh;
		sat.dps.ssh = sat.dps.ssh + sh / sat.deep_arg.sinio;

		sat.flags &= ~RESONANCE_FLAG;
		sat.flags &= ~SYNCHRONOUS_FLAG;

		if (!((sat.dps.xnq < 0.0052359877) && (sat.dps.xnq > 0.0034906585))) {
			if ((sat.dps.xnq < 0.00826) || (sat.dps.xnq > 0.00924))
				return;
			if (eq < 0.5)
				return;
			sat.flags |= RESONANCE_FLAG;
			eoc = eq * sat.deep_arg.eosq;
			g201 = -0.306 - (eq - 0.64) * 0.440;
			if (eq <= 0.65) {
				g211 = 3.616 - 13.247 * eq + 16.290 * sat.deep_arg.eosq;
				g310 = -19.302 + 117.390 * eq - 228.419 * sat.deep_arg.eosq
						+ 156.591 * eoc;
				g322 = -18.9068 + 109.7927 * eq - 214.6334 * sat.deep_arg.eosq
						+ 146.5816 * eoc;
				g410 = -41.122 + 242.694 * eq - 471.094 * sat.deep_arg.eosq
						+ 313.953 * eoc;
				g422 = -146.407 + 841.880 * eq - 1629.014 * sat.deep_arg.eosq
						+ 1083.435 * eoc;
				g520 = -532.114 + 3017.977 * eq - 5740 * sat.deep_arg.eosq
						+ 3708.276 * eoc;
			} else {
				g211 = -72.099 + 331.819 * eq - 508.738 * sat.deep_arg.eosq
						+ 266.724 * eoc;
				g310 = -346.844 + 1582.851 * eq - 2415.925 * sat.deep_arg.eosq
						+ 1246.113 * eoc;
				g322 = -342.585 + 1554.908 * eq - 2366.899 * sat.deep_arg.eosq
						+ 1215.972 * eoc;
				g410 = -1052.797 + 4758.686 * eq - 7193.992 * sat.deep_arg.eosq
						+ 3651.957 * eoc;
				g422 = -3581.69 + 16178.11 * eq - 24462.77 * sat.deep_arg.eosq
						+ 12422.52 * eoc;
				if (eq <= 0.715)
					g520 = 1464.74 - 4664.75 * eq + 3763.64 * sat.deep_arg.eosq;
				else
					g520 = -5149.66 + 29936.92 * eq - 54087.36
							* sat.deep_arg.eosq + 31324.56 * eoc;
			}

			if (eq < 0.7) {
				g533 = -919.2277 + 4988.61 * eq - 9064.77 * sat.deep_arg.eosq
						+ 5542.21 * eoc;
				g521 = -822.71072 + 4568.6173 * eq - 8491.4146
						* sat.deep_arg.eosq + 5337.524 * eoc;
				g532 = -853.666 + 4690.25 * eq - 8624.77 * sat.deep_arg.eosq
						+ 5341.4 * eoc;
			} else {
				g533 = -37995.78 + 161616.52 * eq - 229838.2
						* sat.deep_arg.eosq + 109377.94 * eoc;
				g521 = -51752.104 + 218913.95 * eq - 309468.16
						* sat.deep_arg.eosq + 146349.42 * eoc;
				g532 = -40023.88 + 170470.89 * eq - 242699.48
						* sat.deep_arg.eosq + 115605.82 * eoc;
			}

			sini2 = sat.deep_arg.sinio * sat.deep_arg.sinio;
			f220 = 0.75 * (1 + 2 * sat.deep_arg.cosio + sat.deep_arg.theta2);
			f221 = 1.5 * sini2;
			f321 = 1.875 * sat.deep_arg.sinio
					* (1 - 2 * sat.deep_arg.cosio - 3 * sat.deep_arg.theta2);
			f322 = -1.875 * sat.deep_arg.sinio
					* (1 + 2 * sat.deep_arg.cosio - 3 * sat.deep_arg.theta2);
			f441 = 35 * sini2 * f220;
			f442 = 39.3750 * sini2 * sini2;
			f522 = 9.84375
					* sat.deep_arg.sinio
					* (sini2
							* (1 - 2 * sat.deep_arg.cosio - 5 * sat.deep_arg.theta2) + 0.33333333 * (-2
							+ 4 * sat.deep_arg.cosio + 6 * sat.deep_arg.theta2));
			f523 = sat.deep_arg.sinio
					* (4.92187512
							* sini2
							* (-2 - 4 * sat.deep_arg.cosio + 10 * sat.deep_arg.theta2) + 6.56250012 * (1 + 2 * sat.deep_arg.cosio - 3 * sat.deep_arg.theta2));
			f542 = 29.53125
					* sat.deep_arg.sinio
					* (2 - 8 * sat.deep_arg.cosio + sat.deep_arg.theta2
							* (-12 + 8 * sat.deep_arg.cosio + 10 * sat.deep_arg.theta2));
			f543 = 29.53125
					* sat.deep_arg.sinio
					* (-2 - 8 * sat.deep_arg.cosio + sat.deep_arg.theta2
							* (12 + 8 * sat.deep_arg.cosio - 10 * sat.deep_arg.theta2));
			xno2 = sat.dps.xnq * sat.dps.xnq;
			ainv2 = aqnv * aqnv;
			temp1 = 3 * xno2 * ainv2;
			temp = temp1 * root22;
			sat.dps.d2201 = temp * f220 * g201;
			sat.dps.d2211 = temp * f221 * g211;
			temp1 = temp1 * aqnv;
			temp = temp1 * root32;
			sat.dps.d3210 = temp * f321 * g310;
			sat.dps.d3222 = temp * f322 * g322;
			temp1 = temp1 * aqnv;
			temp = 2 * temp1 * root44;
			sat.dps.d4410 = temp * f441 * g410;
			sat.dps.d4422 = temp * f442 * g422;
			temp1 = temp1 * aqnv;
			temp = temp1 * root52;
			sat.dps.d5220 = temp * f522 * g520;
			sat.dps.d5232 = temp * f523 * g532;
			temp = 2 * temp1 * root54;
			sat.dps.d5421 = temp * f542 * g521;
			sat.dps.d5433 = temp * f543 * g533;
			sat.dps.xlamo = xmao + sat.tle.xnodeo + sat.tle.xnodeo
					- sat.dps.thgr - sat.dps.thgr;
			bfact = sat.deep_arg.xmdot + sat.deep_arg.xnodot
					+ sat.deep_arg.xnodot - thdt - thdt;
			bfact = bfact + sat.dps.ssl + sat.dps.ssh + sat.dps.ssh;
		} else {
			sat.flags |= RESONANCE_FLAG;
			sat.flags |= SYNCHRONOUS_FLAG;

			g200 = 1 + sat.deep_arg.eosq * (-2.5 + 0.8125 * sat.deep_arg.eosq);
			g310 = 1 + 2 * sat.deep_arg.eosq;
			g300 = 1 + sat.deep_arg.eosq * (-6 + 6.60937 * sat.deep_arg.eosq);
			f220 = 0.75 * (1 + sat.deep_arg.cosio) * (1 + sat.deep_arg.cosio);
			f311 = 0.9375 * sat.deep_arg.sinio * sat.deep_arg.sinio
					* (1 + 3 * sat.deep_arg.cosio) - 0.75
					* (1 + sat.deep_arg.cosio);
			f330 = 1 + sat.deep_arg.cosio;
			f330 = 1.875 * f330 * f330 * f330;
			sat.dps.del1 = 3 * sat.dps.xnq * sat.dps.xnq * aqnv * aqnv;
			sat.dps.del2 = 2 * sat.dps.del1 * f220 * g200 * q22;
			sat.dps.del3 = 3 * sat.dps.del1 * f330 * g300 * q33 * aqnv;
			sat.dps.del1 = sat.dps.del1 * f311 * g310 * q31 * aqnv;
			sat.dps.fasx2 = 0.13130908;
			sat.dps.fasx4 = 2.8843198;
			sat.dps.fasx6 = 0.37448087;
			sat.dps.xlamo = xmao + sat.tle.xnodeo + sat.tle.omegao
					- sat.dps.thgr;
			bfact = sat.deep_arg.xmdot + xpidot - thdt;
			bfact = bfact + sat.dps.ssl + sat.dps.ssg + sat.dps.ssh;
		}

		sat.dps.xfact = bfact - sat.dps.xnq;

		sat.dps.xli = sat.dps.xlamo;
		sat.dps.xni = sat.dps.xnq;
		sat.dps.atime = 0;
		sat.dps.stepp = 720;
		sat.dps.stepn = -720;
		sat.dps.step2 = 259200;
		return;

	case dpsec:
		sat.deep_arg.xll = sat.deep_arg.xll + sat.dps.ssl * sat.deep_arg.t;
		sat.deep_arg.omgadf = sat.deep_arg.omgadf + sat.dps.ssg
				* sat.deep_arg.t;
		sat.deep_arg.xnode = sat.deep_arg.xnode + sat.dps.ssh * sat.deep_arg.t;
		sat.deep_arg.em = sat.tle.eo + sat.dps.sse * sat.deep_arg.t;
		sat.deep_arg.xinc = sat.tle.xincl + sat.dps.ssi * sat.deep_arg.t;
		if (sat.deep_arg.xinc < 0) {
			sat.deep_arg.xinc = -sat.deep_arg.xinc;
			sat.deep_arg.xnode = sat.deep_arg.xnode + pi;
			sat.deep_arg.omgadf = sat.deep_arg.omgadf - pi;
		}
		if (~sat.flags & RESONANCE_FLAG)
			return;

		do {
			if ((sat.dps.atime == 0)
					|| ((sat.deep_arg.t >= 0) && (sat.dps.atime < 0))
					|| ((sat.deep_arg.t < 0) && (sat.dps.atime >= 0))) {
				if (sat.deep_arg.t >= 0)
					delt = sat.dps.stepp;
				else
					delt = sat.dps.stepn;

				sat.dps.atime = 0;
				sat.dps.xni = sat.dps.xnq;
				sat.dps.xli = sat.dps.xlamo;
			} else {
				if (fabs(sat.deep_arg.t) >= fabs(sat.dps.atime)) {
					if (sat.deep_arg.t > 0)
						delt = sat.dps.stepp;
					else
						delt = sat.dps.stepn;
				}
			}

			do {
				if (fabs(sat.deep_arg.t - sat.dps.atime) >= sat.dps.stepp) {
					sat.flags |= DO_LOOP_FLAG;
					sat.flags &= ~EPOCH_RESTART_FLAG;
				} else {
					ft = sat.deep_arg.t - sat.dps.atime;
					sat.flags &= ~DO_LOOP_FLAG;
				}

				if (fabs(sat.deep_arg.t) < fabs(sat.dps.atime)) {
					if (sat.deep_arg.t >= 0)
						delt = sat.dps.stepn;
					else
						delt = sat.dps.stepp;
					sat.flags |= (DO_LOOP_FLAG | EPOCH_RESTART_FLAG);
				}

				if (sat.flags & SYNCHRONOUS_FLAG) {
					xndot = sat.dps.del1
							* Math.sin(sat.dps.xli - sat.dps.fasx2)
							+ sat.dps.del2
							* Math.sin(2 * (sat.dps.xli - sat.dps.fasx4))
							+ sat.dps.del3
							* Math.sin(3 * (sat.dps.xli - sat.dps.fasx6));
					xnddt = sat.dps.del1
							* Math.cos(sat.dps.xli - sat.dps.fasx2) + 2
							* sat.dps.del2
							* Math.cos(2 * (sat.dps.xli - sat.dps.fasx4)) + 3
							* sat.dps.del3
							* Math.cos(3 * (sat.dps.xli - sat.dps.fasx6));
				} else {
					xomi = sat.dps.omegaq + sat.deep_arg.omgdot * sat.dps.atime;
					x2omi = xomi + xomi;
					x2li = sat.dps.xli + sat.dps.xli;
					xndot = sat.dps.d2201 * Math.sin(x2omi + sat.dps.xli - g22)
							+ sat.dps.d2211 * Math.sin(sat.dps.xli - g22)
							+ sat.dps.d3210
							* Math.sin(xomi + sat.dps.xli - g32)
							+ sat.dps.d3222
							* Math.sin(-xomi + sat.dps.xli - g32)
							+ sat.dps.d4410 * Math.sin(x2omi + x2li - g44)
							+ sat.dps.d4422 * Math.sin(x2li - g44)
							+ sat.dps.d5220
							* Math.sin(xomi + sat.dps.xli - g52)
							+ sat.dps.d5232
							* Math.sin(-xomi + sat.dps.xli - g52)
							+ sat.dps.d5421 * Math.sin(xomi + x2li - g54)
							+ sat.dps.d5433 * Math.sin(-xomi + x2li - g54);
					xnddt = sat.dps.d2201
							* Math.cos(x2omi + sat.dps.xli - g22)
							+ sat.dps.d2211
							* Math.cos(sat.dps.xli - g22)
							+ sat.dps.d3210
							* Math.cos(xomi + sat.dps.xli - g32)
							+ sat.dps.d3222
							* Math.cos(-xomi + sat.dps.xli - g32)
							+ sat.dps.d5220
							* Math.cos(xomi + sat.dps.xli - g52)
							+ sat.dps.d5232
							* Math.cos(-xomi + sat.dps.xli - g52)
							+ 2
							* (sat.dps.d4410 * Math.cos(x2omi + x2li - g44)
									+ sat.dps.d4422 * Math.cos(x2li - g44)
									+ sat.dps.d5421
									* Math.cos(xomi + x2li - g54) + sat.dps.d5433
									* Math.cos(-xomi + x2li - g54));
				}

				xldot = sat.dps.xni + sat.dps.xfact;
				xnddt = xnddt * xldot;

				if (sat.flags & DO_LOOP_FLAG) {
					sat.dps.xli = sat.dps.xli + xldot * delt + xndot
							* sat.dps.step2;
					sat.dps.xni = sat.dps.xni + xndot * delt + xnddt
							* sat.dps.step2;
					sat.dps.atime = sat.dps.atime + delt;
				}
			} while ((sat.flags & DO_LOOP_FLAG)
					&& (~sat.flags & EPOCH_RESTART_FLAG));
		} while ((sat.flags & DO_LOOP_FLAG) && (sat.flags & EPOCH_RESTART_FLAG));

		sat.deep_arg.xn = sat.dps.xni + xndot * ft + xnddt * ft * ft * 0.5;
		xl = sat.dps.xli + xldot * ft + xndot * ft * ft * 0.5;
		temp = -sat.deep_arg.xnode + sat.dps.thgr + sat.deep_arg.t * thdt;

		if (~sat.flags & SYNCHRONOUS_FLAG)
			sat.deep_arg.xll = xl + temp + temp;
		else
			sat.deep_arg.xll = xl - sat.deep_arg.omgadf + temp;

		return;

	case dpper:
		sinis = Math.sin(sat.deep_arg.xinc);
		cosis = Math.cos(sat.deep_arg.xinc);
		if (fabs(sat.dps.savtsn - sat.deep_arg.t) >= 30) {
			sat.dps.savtsn = sat.deep_arg.t;
			zm = sat.dps.zmos + zns * sat.deep_arg.t;
			zf = zm + 2 * zes * Math.sin(zm);
			sinzf = Math.sin(zf);
			f2 = 0.5 * sinzf * sinzf - 0.25;
			f3 = -0.5 * sinzf * Math.cos(zf);
			ses = sat.dps.se2 * f2 + sat.dps.se3 * f3;
			sis = sat.dps.si2 * f2 + sat.dps.si3 * f3;
			sls = sat.dps.sl2 * f2 + sat.dps.sl3 * f3 + sat.dps.sl4 * sinzf;
			sat.dps.sghs = sat.dps.sgh2 * f2 + sat.dps.sgh3 * f3 + sat.dps.sgh4
					* sinzf;
			sat.dps.shs = sat.dps.sh2 * f2 + sat.dps.sh3 * f3;
			zm = sat.dps.zmol + znl * sat.deep_arg.t;
			zf = zm + 2 * zel * Math.sin(zm);
			sinzf = Math.sin(zf);
			f2 = 0.5 * sinzf * sinzf - 0.25;
			f3 = -0.5 * sinzf * Math.cos(zf);
			sel = sat.dps.ee2 * f2 + sat.dps.e3 * f3;
			sil = sat.dps.xi2 * f2 + sat.dps.xi3 * f3;
			sll = sat.dps.xl2 * f2 + sat.dps.xl3 * f3 + sat.dps.xl4 * sinzf;
			sat.dps.sghl = sat.dps.xgh2 * f2 + sat.dps.xgh3 * f3 + sat.dps.xgh4
					* sinzf;
			sat.dps.sh1 = sat.dps.xh2 * f2 + sat.dps.xh3 * f3;
			sat.dps.pe = ses + sel;
			sat.dps.pinc = sis + sil;
			sat.dps.pl = sls + sll;
		}

		pgh = sat.dps.sghs + sat.dps.sghl;
		ph = sat.dps.shs + sat.dps.sh1;
		sat.deep_arg.xinc = sat.deep_arg.xinc + sat.dps.pinc;
		sat.deep_arg.em = sat.deep_arg.em + sat.dps.pe;

		if (sat.dps.xqncl >= 0.2) {
			ph = ph / sat.deep_arg.sinio;
			pgh = pgh - sat.deep_arg.cosio * ph;
			sat.deep_arg.omgadf = sat.deep_arg.omgadf + pgh;
			sat.deep_arg.xnode = sat.deep_arg.xnode + ph;
			sat.deep_arg.xll = sat.deep_arg.xll + sat.dps.pl;
		} else {
			sinok = Math.sin(sat.deep_arg.xnode);
			cosok = Math.cos(sat.deep_arg.xnode);
			alfdp = sinis * sinok;
			betdp = sinis * cosok;
			dalf = ph * cosok + sat.dps.pinc * cosis * sinok;
			dbet = -ph * sinok + sat.dps.pinc * cosis * cosok;
			alfdp = alfdp + dalf;
			betdp = betdp + dbet;
			sat.deep_arg.xnode = FMod2p(sat.deep_arg.xnode);
			xls = sat.deep_arg.xll + sat.deep_arg.omgadf + cosis
					* sat.deep_arg.xnode;
			dls = sat.dps.pl + pgh - sat.dps.pinc * sat.deep_arg.xnode * sinis;
			xls = xls + dls;
			xnoh = sat.deep_arg.xnode;
			sat.deep_arg.xnode = AcTan(alfdp, betdp);

			if (fabs(xnoh - sat.deep_arg.xnode) > pi) {
				if (sat.deep_arg.xnode < xnoh)
					sat.deep_arg.xnode += twopi;
				else
					sat.deep_arg.xnode -= twopi;
			}

			sat.deep_arg.xll = sat.deep_arg.xll + sat.dps.pl;
			sat.deep_arg.omgadf = xls - sat.deep_arg.xll
					- Math.cos(sat.deep_arg.xinc) * sat.deep_arg.xnode;
		}
		return;

	}

}

var Flags = 0;

function isFlagSet(flag) {
	return (Flags & flag);
}

function isFlagClear(flag) {
	return (~Flags & flag);
}

function SetFlag(flag) {
	Flags |= flag;
}

function ClearFlag(flag) {
	Flags &= ~flag;
}

// //////////////////////////////////////////////////////////////////////////////
