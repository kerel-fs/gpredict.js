var gulp = require('gulp');
var rename = require('gulp-rename');
var uglify = require('gulp-uglify');
var concat = require('gulp-concat');

// Concatenate JS
gulp.task('concat:js', function(){
    return gulp.src([
	    'js/orbit_tools.js',
        'js/sgp_in.js',
        'js/sgp4sdp4.js',
	    'js/predict-tools.js',
	    'js/sgp_time.js',
	    'js/sgp_obs.js',
	    'js/math.js',
	    'js/gtk-sat-data.js',
	    'js/qth-data.js',
	    'js/sgp_math.js'])
    .pipe(concat('gpredict.js'))
    .pipe(gulp.dest('./dist/'));
});

// Minify JS
gulp.task('minify:js', function() {
  return gulp.src('dist/gpredict.js')
    .pipe(rename({ suffix: '.min' }))
    .pipe(uglify())
    .pipe(gulp.dest('./dist/'));
});

// Default
gulp.task('default', ['concat:js', 'minify:js']);
