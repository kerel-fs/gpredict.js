/*
    Gpredict: Real-time satellite tracking and orbit prediction program

    Copyright (C)  2001-2009  Alexandru Csete, OZ9AEC.

    Authors: Alexandru Csete <oz9aec@gmail.com>

    Comments, questions and bugreports should be submitted via
    http://sourceforge.net/projects/gpredict/
    More details can be found at the project home page:

            http://gpredict.oz9aec.net/
 
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
  
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
  
    You should have received a copy of the GNU General Public License
    along with this program; if not, visit http://www.fsf.org/
*/

/**
 * Satellite pass info default constructor.
 *
 * @returns {pass_t}
 * @constructor
 */
function pass_t() {
	this.satname = "";	/*!< satellite name */
	this.sat = null;
    this.aos = 0.0;      /*!< AOS time in "jul_utc" */
    this.tca = 0.0;      /*!< TCA time in "jul_utc" */
    this.los = 0.0;      /*!< LOS time in "jul_utc" */
    this.max_el = 0.0;   /*!< Maximum elevation during pass */
    this.aos_az = 0.0;   /*!< Azimuth at AOS */
    this.los_az = 0.0;   /*!< Azimuth at LOS */
    this.orbit = 0;    /*!< Orbit number */
    this.maxel_az = 0.0; /*!< Azimuth at maximum elevation */
    this.details = new Array();/*!< List of pass_detail_t entries */
}

/**
 * Pass detail entry.
 *
 * @constructor
 *
 * In order to ensure maximum flexibility at a minimal effort, only the raw
 * position and velocity is calculated. Calculations of the "human readable"
 * parameters are the responsibility of the consumer. This way we can use the
 * same prediction engine for various consumers without having too much overhead
 * and complexity in the low level code.
 */
function pass_detail_t() {
	this.time = 0.0;   /*!< time in "jul_utc" */
    this.pos = new vector_t();    /*!< Raw unprocessed position at time */
    this.vel = new vector_t();    /*!< Raw unprocessed velocity at time */
    this.velo = 0.0;
    this.az = 0.0;
    this.el = 0.0;
    this.range = 0.0;
    this.range_rate = 0.0;
    this.lat = 0.0;
    this.lon = 0.0;
    this.alt = 0.0;
    this.ma = 0.0;
    this.phase = 0.0;
    this.footprint = 0.0;
    this.orbit = 0;
}

//////////////////////////////////////////////////////////////////////////////
// from predict-tools.c
//////////////////////////////////////////////////////////////////////////////

/**
 * SGP4SDP4 driver for doing AOS/LOS calculations.
 *
 * @param {Sat_t} sat Pointer to the satellite data.
 * @param {qth_t} qth Pointer to the QTH data.
 * @param {number} t The time for calculation (Julian Date)
 */
function predict_calc(sat, qth, t)
{
	var obs_set = new obs_set_t();
    var sat_geodetic = new geodetic_t();
    var obs_geodetic = new geodetic_t();
    var age = 0.0;

    obs_geodetic.lon = qth.lon * de2ra;
    obs_geodetic.lat = qth.lat * de2ra;
    obs_geodetic.alt = qth.alt / 1000.0;
    obs_geodetic.theta = 0;

    sat.jul_utc = t;
    sat.tsince = (sat.jul_utc - sat.jul_epoch) * xmnpda;

    /* call the norad routines according to the deep-space flag */
    if (sat.flags & DEEP_SPACE_EPHEM_FLAG)
        SDP4 (sat, sat.tsince);
    else
        SGP4 (sat, sat.tsince);

    Convert_Sat_State (sat.pos, sat.vel);

    /* get the velocity of the satellite */
    Magnitude (sat.vel);
    sat.velo = sat.vel.w;
    Calculate_Obs (sat.jul_utc, sat.pos, sat.vel, obs_geodetic, obs_set);
    Calculate_LatLonAlt (sat.jul_utc, sat.pos, sat_geodetic);

    while (sat_geodetic.lon < -pi)
        sat_geodetic.lon += twopi;

    while (sat_geodetic.lon > (pi))
        sat_geodetic.lon -= twopi;

    sat.az = Degrees (obs_set.az);
    sat.el = Degrees (obs_set.el);
    sat.range = obs_set.range;
    sat.range_rate = obs_set.range_rate;
    sat.ssplat = Degrees (sat_geodetic.lat);
    sat.ssplon = Degrees (sat_geodetic.lon);
    sat.alt = sat_geodetic.alt;
    sat.ma = Degrees (sat.phase);
    sat.ma *= 256.0/360.0;
    sat.phase = Degrees (sat.phase);

    /* same formulas, but the one from predict is nicer */
    //sat.footprint = 2.0 * xkmper * acos (xkmper/sat.pos.w);
    sat.footprint = 12756.33 * acos (xkmper / (xkmper+sat.alt));
    age = sat.jul_utc - sat.jul_epoch;
    sat.orbit = floor((sat.tle.xno * xmnpda/twopi +
                    age * sat.tle.bstar * ae) * age +
                    sat.tle.xmo/twopi) + sat.tle.revnum - 1;
}

/**
 * Find the AOS time of the next pass.
 * 
 * @author Alexandru Csete, OZ9AEC
 * @author John A. Magliacane, KD2BD
 * 
 * @param {Sat_t} sat Pointer to the satellite data.
 * @param {qth_t} qth Pointer to the QTH data.
 * @param {number} start The time where calculation should start.
 * @param {number} maxdt The upper time limit in days (0.0 = no limit)
 *
 * @returns {number} The time of the next AOS or 0.0 if the satellite has no AOS.
 *
 * This function finds the time of AOS for the first coming pass taking place
 * no earlier that start.
 * If the satellite is currently within range, the function first calls
 * find_los to get the next LOS time. Then the calculations are done using
 * the new start time.
 */
function find_aos (sat, qth, start, maxdt)
{
    var t = start;
    var aostime = 0.0;

    /* make sure current sat values are
        in sync with the time
    */
    predict_calc (sat, qth, start);

    /* check whether satellite has aos */
    if ((sat.otype == ORBIT_TYPE_GEO) || 
        (sat.otype == ORBIT_TYPE_DECAYED) ||
        !has_aos (sat, qth)) {
        return 0.0;
    }

    if (sat.el > 0.0)
        t = find_los (sat, qth, start, maxdt) + 0.014; // +20 min

    /* invalid time (potentially returned by find_los) */
    if (t < 0.1)
        return 0.0;

    /* update satellite data */
    predict_calc (sat, qth, t);

    /* use upper time limit */
    if (maxdt > 0.0) {

        /* coarse time steps */
        while ((sat.el < -1.0) && (t <= (start + maxdt))) {
            t -= 0.00035 * (sat.el * ((sat.alt / 8400.0) + 0.46) - 2.0);
            predict_calc (sat, qth, t);
        }

        /* fine steps */
        while ((aostime == 0.0) && (t <= (start + maxdt))) {
            if (fabs (sat.el) < 0.005) {
                aostime = t;
            }
            else {
                t -= sat.el * sqrt (sat.alt) / 530000.0;
                predict_calc (sat, qth, t);
            }
        }
    }
    /* don't use upper time limit */
    else {

        /* coarse time steps */
        while (sat.el < -1.0) {

            t -= 0.00035 * (sat.el * ((sat.alt / 8400.0) + 0.46) - 2.0);
            predict_calc (sat, qth, t);
        }

        /* fine steps */
        while (aostime == 0.0) {

            if (fabs (sat.el) < 0.005) {
                aostime = t;
            }
            else {
                t -= sat.el * sqrt (sat.alt) / 530000.0;
                predict_calc (sat, qth, t);
            }
        }
    }
    return aostime;
}


/**
 * Find the LOS time of the next pass.
 *
 * @author Alexandru Csete, OZ9AEC
 * @author John A. Magliacane, KD2BD
 * @param {Sat_t} sat Pointer to the satellite data.
 * @param {qth_t} qth Pointer to the QTH data.
 * @param {number} start The time where calculation should start.
 * @param {number} maxdt The upper time limit in days (0.0 = no limit)
 * @returns {number} The time of the next LOS or 0.0 if the satellite has no LOS.
 *
 * This function finds the time of LOS for the first coming pass taking place
 * no earlier that start.
 * If the satellite is currently out of range, the function first calls
 * find_aos to get the next AOS time. Then the calculations are done using
 * the new start time.
 * The function has a built-in watchdog to ensure that we don't end up in
 * lengthy loops.
 */
function find_los (sat, qth, start, maxdt)
{
    var t = start;
    var lostime = 0.0;

    predict_calc (sat, qth, start);

    /* check whether satellite has aos */
    if ((sat.otype == ORBIT_TYPE_GEO) || 
        (sat.otype == ORBIT_TYPE_DECAYED) ||
        !has_aos (sat, qth)) {
        return 0.0;
    }

    if (sat.el < 0.0)
        t = find_aos (sat, qth, start, maxdt) + 0.001; // +1.5 min

    /* invalid time (potentially returned by find_aos) */
    if (t < 0.01)
        return 0.0;

    /* update satellite data */
    predict_calc (sat, qth, t);

    /* use upper time limit */
    if (maxdt > 0.0) {

        /* coarse steps */
        while ((sat.el >= 1.0) && (t <= (start + maxdt))) {
            t += cos((sat.el - 1.0) * de2ra) * sqrt(sat.alt) / 25000.0;
            predict_calc (sat, qth, t);
        }

        /* fine steps */
        while ((lostime == 0.0) && (t <= (start + maxdt)))  {
            
            t += sat.el * sqrt(sat.alt)/502500.0;
            predict_calc (sat, qth, t);
            
            if (fabs(sat.el) < 0.005)
                lostime = t;
        }
    }

    /* don't use upper limit */
    else {
        /* coarse steps */
        while (sat.el >= 1.0) {
            t += cos((sat.el - 1.0) * de2ra) * sqrt(sat.alt) / 25000.0;
            predict_calc (sat, qth, t);
        }

        /* fine steps */
        while (lostime == 0.0) {
            t += sat.el * sqrt(sat.alt)/502500.0;
            predict_calc (sat, qth, t);
            
            if (fabs(sat.el) < 0.005)
                lostime = t;
        }
    }

    return lostime;
}


/**
 * Find AOS time of current pass.
 *
 * @param {Sat_t} sat The satellite to find AOS for.
 * @param {qth_t} qth The ground station.
 * @param {number} start Start time, prefereably now.
 * @returns {number} The time of the previous AOS or 0.0 if the satellite has no AOS.
 *
 * This function can be used to find the AOS time in the past of the
 * current pass.
 */
function find_prev_aos (sat, qth, start)
{
    var aostime = start;

    /* make sure current sat values are
        in sync with the time
    */
    predict_calc (sat, qth, start);

    /* check whether satellite has aos */
    if ((sat.otype == ORBIT_TYPE_GEO) || 
        (sat.otype == ORBIT_TYPE_DECAYED) ||
        !has_aos (sat, qth)) {

        return 0.0;
    }

    while (sat.el >= 0.0) {
        aostime -= 0.0005; // 0.75 min
        predict_calc (sat, qth, aostime);
    }

    return aostime;
}

/**
 * Predict upcoming passes starting now
 *
 * @param {Sat_t} sat Pointer to the satellite data.
 * @param {qth_t} qth Pointer to the observer data.
 * @param {number} maxdt The maximum number of days to look ahead.
 * @param {number} num The number of passes to predict.
 * @returns {Array} A singly linked list of pass_t structures or NULL if
 *          there was an error.
 *
 * This function simply wraps the get_passes function using the
 * current time as parameter.
 *
 * @ the data in sat will be corrupt (future) and must be refreshed
 *       by the caller, if the caller will need it later on (eg. if the caller
 *       is GtkSatList).
 */
function get_next_passes (sat, qth, maxdt, num)
{
    /* get the current time and call
        the get_pass function */
    var now = Julian_Date(new Date());
    return get_passes (sat, qth, now, maxdt, num);
}


/**
 * Predict first pass after a certain time.
 *
 * @param {Sat_t} sat Pointer to the satellite data.
 * @param {qth_t} qth Pointer to the location data.
 * @param {number} start Starting time.
 * @param {number} maxdt The maximum number of days to look ahead (0 for no limit).
 * @param {pass_t} pass
 *
 * @returns Pointer to a newly allocated pass_t structure or NULL if
 *          there was an error.
 *
 * This function will find the first upcoming pass with AOS no earlier than
 * t = start and no later than t = (start+maxdt).
 *
 * \note For no time limit use maxdt = 0.0
 *
 * \note the data in sat will be corrupt (future) and must be refreshed
 *       by the caller, if the caller will need it later on (eg. if the caller
 *       is GtkSatList).
 *
 * \note Prepending to a singly linked list is much faster than appending.
 *       Therefore, the elements are prepended whereafter the GSList is
 *       reversed
 */
function get_pass(sat, qth, start, maxdt, pass)
{
    var aos = 0.0;    /* time of AOS */
    var tca = 0.0;    /* time of TCA */
    var los = 0.0;    /* time of LOS */
    var dt = 0.0;     /* time diff */
    var step = 0.0;   /* time step */
    var t0 = start;
    var t;            /* current time counter */
    var tres = 0.0; /* required time resolution */
    var max_el = 0.0; /* maximum elevation */
    var detail;
    var done = false;
    var iter = 0;      /* number of iterations */
    /* FIXME: watchdog */

    /* get time resolution; sat-cfg stores it in seconds */
    tres = SAT_CFG_INT_PRED_RESOLUTION / 86400.0;

    /* loop until we find a pass with elevation > SAT_CFG_INT_PRED_MIN_EL
        or we run out of time
        FIXME: we should have a safety break
    */
    while (!done) {

        /* Find los of next pass or of current pass */
        los = find_los (sat, qth, t0, maxdt); // See if a pass is ongoing
        aos = find_aos (sat, qth, t0, maxdt);
        /* sat_log_log(SAT_LOG_LEVEL_MSG, "%s:%s:%d: found aos %f and los %f for t0=%f", */
        /*          __FILE__,  */
        /*          __FUNCTION__, */
        /*          __LINE__, */
        /*          aos, */
        /*          los,  */
        /*          t0); */
        if (aos > los) {
            // los is from an currently happening pass, find previous aos
            aos = find_prev_aos(sat, qth, t0);
        }

        /* aos = 0.0 means no aos */
        if (aos == 0.0) {
            done = true;
        }

        /* check whether we are within time limits;
            maxdt = 0 mean no time limit.
        */
        else if ((maxdt > 0.0) && (aos > (start + maxdt)) ) {
            done = true;
        }
        else {
            //los = find_los (sat, qth, aos + 0.001, maxdt); // +1.5 min later
            dt = los - aos;

            /* get time step, which will give us the max number of entries */
            step = dt / SAT_CFG_INT_PRED_NUM_ENTRIES;

            /* but if this is smaller than the required resolution
                we go with the resolution
            */
            if (step < tres)
                step = tres;

            pass = new pass_t();
            pass.sat = sat;
            pass.aos = aos;
            pass.los = los;
            pass.satname = sat.nickname;

            /* iterate over each time step */
            for (t = pass.aos; t <= pass.los; t += step) {

                /* calculate satellite data */
                predict_calc (sat, qth, t);

                /* in the first iter we want to store
                    pass.aos_az
                */
                if (t == pass.aos) {
                    pass.aos_az = sat.az;
                    pass.orbit = sat.orbit;
                }

                /* append details to sat.details */
                detail = new pass_detail_t();
                detail.time = t;
                detail.pos.x = sat.pos.x;
                detail.pos.y = sat.pos.y;
                detail.pos.z = sat.pos.z;
                detail.pos.w = sat.pos.w;
                detail.vel.x = sat.vel.x;
                detail.vel.y = sat.vel.y;
                detail.vel.z = sat.vel.z;
                detail.vel.w = sat.vel.w;
                detail.velo = sat.velo;
                detail.az = sat.az;
                detail.el = sat.el;
                detail.range = sat.range;
                detail.range_rate = sat.range_rate;
                detail.lat = sat.ssplat;
                detail.lon = sat.ssplon;
                detail.alt = sat.alt;
                detail.ma = sat.ma;
                detail.phase = sat.phase;
                detail.footprint = sat.footprint;
                detail.orbit = sat.orbit;

                pass.details.push(detail);

                /* store elevation if greater than the
                    previously stored one
                */
                if (sat.el > max_el) {
                    max_el = sat.el;
                    tca = t;
                    pass.maxel_az = sat.az;
                }

                /*     g_print ("TIME: %f\tAZ: %f\tEL: %f (MAX: %f)\n", */
                /*           t, sat.az, sat.el, max_el); */
            }

            /* calculate satellite data */
            predict_calc (sat, qth, pass.los);
            /* store los_az, max_el and tca */
            pass.los_az = sat.az;
            pass.max_el = max_el;
            pass.tca    = tca;

            /* check whether this pass is good */
            if (max_el >= SAT_CFG_INT_PRED_MIN_EL) {
                done = true;
            }
            else {
                done = false;
                t0 = los + 0.014; // +20 min
                pass = null;
            }

            iter++;
        }
    }

    return pass;
}




/**
 * Predict passes after a certain time.
 *
 * @param {Sat_t} sat
 * @param {qth_t} qth
 * @param {number} start
 * @param {number} maxdt
 * @param {number} num
 *
 * @returns {Array}
 *
 * This function calculates num upcoming passes with AOS no earlier
 * than t = start and not later that t = (start+maxdt). The function will
 *  repeatedly call get_pass until
 * the number of predicted passes is equal to num, the time has reached
 * limit or the get_pass function returns NULL.
 *
 * \note For no time limit use maxdt = 0.0
 *
 * \note the data in sat will be corrupt (future) and must be refreshed
 *       by the caller, if the caller will need it later on (eg. if the caller
 *       is GtkSatList).
 */
function get_passes (sat, qth, start, maxdt, num)
{
    var passes = new Array();
    var pass;
    var i;
    var t;

    /* if no number has been specified
        set it to something big */
    if (num == 0)
        num = 100;

    t = start;

    for (i = 0; i < num; i++) {
        pass = get_pass (sat, qth, t, maxdt);

        if (pass != null) {
            passes.push(pass);
            t = pass.los + 0.014; // +20 min

            /* if maxdt > 0.0 check whether we have reached t = start+maxdt
                if yes finish predictions
            */
            if ((maxdt > 0.0) && (t >= (start+maxdt))) {
                i = num;
            }
        }
        else {
            /* we can't get any more passes */
            i = num;
        }

    }

    return passes;
}
