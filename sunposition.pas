(* Lazarus + FPC 0.9.30 + 2.4.4. On Linux for ARM, PPC, SPARC, x86. Lazarus + F *)

(* Some of the maths might require 2.8.0 on SPARC for reliable operation.       *)

unit SunPosition;

(* This calculates the time of sunrise/sunset for a given location, ignoring    *)
(* horizon elevation etc. It is transcribed from the JavaScript code at http:// *)
(* www.srrb.noaa.gov/highlights/sunrise/sunrise.html which is a superset of the *)
(* frequently-cited http://www.srrb.noaa.gov/highlights/sunrise/program.txt.    *)
(*                                                                              *)
(* The object of the exercise is to be able to predict an impending sunrise or  *)
(* sunset, which can reasonably be expected to result in a significant change   *)
(* in temperature. Abrupt changes of temperature at other times are probably    *)
(* noteworthy, since they might indicate an imminent storm.                     *)

(* Transcribed to Pascal by Mark Morgan Lloyd.                                  *)
{$hints off}{$notes off}

{$mode objfpc}

interface

uses
  Classes, SysUtils; 

FUNCTION CalcJD(year, month, day: INTEGER): DOUBLE;

//***********************************************************************/
//* Name:    CalcJD							*/
//* Type:    Function							*/
//* Purpose: Julian day from calendar day				*/
//* Arguments:								*/
//*   year : 4 digit year						*/
//*   month: January = 1						*/
//*   day  : 1 - 31							*/
//* Return value:							*/
//*   The Julian day corresponding to the date				*/
//* Note:								*/
//*   Number is returned for start of day.  Fractional days should be	*/
//*   added later.							*/
//***********************************************************************/

FUNCTION CalcTimeJulianCent(jd: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    CalcTimeJulianCent						*/
//* Type:    Function							*/
//* Purpose: convert Julian Day to centuries since J2000.0.		*/
//* Arguments:								*/
//*   jd : the Julian Day to convert					*/
//* Return value:							*/
//*   the T value corresponding to the Julian Day			*/
//***********************************************************************/

FUNCTION CalcSolNoonUTC(t, longitude: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    CalcSolNoonUTC						*/
//* Type:    Function							*/
//* Purpose: calculate the Universal Coordinated Time (UTC) of solar	*/
//*		noon for the given day at the given location on earth	*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//*   longitude : longitude of observer in degrees			*/
//* Return value:							*/
//*   time in minutes from zero Z					*/
//***********************************************************************/

FUNCTION CalcSunriseUTC(JD, latitude, longitude: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    CalcSunriseUTC						*/
//* Type:    Function							*/
//* Purpose: calculate the Universal Coordinated Time (UTC) of sunrise	*/
//*			for the given day at the given location on earth*/
//* Arguments:								*/
//*   JD  : julian day							*/
//*   latitude : latitude of observer in degrees			*/
//*   longitude : longitude of observer in degrees			*/
//* Return value:							*/
//*   time in minutes from zero Z					*/
//***********************************************************************/

FUNCTION CalcSunsetUTC(JD, latitude, longitude: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    CalcSunsetUTC						*/
//* Type:    Function							*/
//* Purpose: calculate the Universal Coordinated Time (UTC) of sunset	*/
//*			for the given day at the given location on earth*/
//* Arguments:								*/
//*   JD  : julian day							*/
//*   latitude : latitude of observer in degrees			*/
//*   longitude : longitude of observer in degrees			*/
//* Return value:							*/
//*   time in minutes from zero Z					*/
//***********************************************************************/

FUNCTION CalcSunElevation(t, utcmins, latitude, longitude: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    CalcSunElevation						*/
//* Type:    Function							*/
//* Purpose: calculate the elevation (altitude) of the sun		*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//*   utcmins: UTC time in minutes                                      */
//*   latitude : latitude of observer in degrees			*/
//*   longitude : longitude of observer in degrees			*/
//* Return value:							*/
//*   sun's elevation/altitude in degrees				*/
//***********************************************************************/

FUNCTION CalcSunAzimuth(t, utcmins, latitude, longitude: DOUBLE): DOUBLE;

// ***** NOT TESTED *****

//***********************************************************************/
//* Name:    CalcSunAzimuth						*/
//* Type:    Function							*/
//* Purpose: calculate the azimuth of the sun		                */
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//*   utcmins: UTC time in minutes                                      */
//*   latitude : latitude of observer in degrees			*/
//*   longitude : longitude of observer in degrees			*/
//* Return value:							*/
//*   sun's azimuth in degrees				                */
//***********************************************************************/


implementation

USES Math;

//***********************************************************************/
//***********************************************************************/
//*									*/
//*This section contains subroutines used in calculating solar position */
//*									*/
//***********************************************************************/
//***********************************************************************/

FUNCTION radToDeg(angleRad: DOUBLE): DOUBLE;

// Convert radian angle to degrees

BEGIN
  RESULT:= 180.0 * angleRad / Pi
END { radToDeg } ;


function degToRad(angleDeg: DOUBLE): DOUBLE;

// Convert degree angle to radians

BEGIN
  RESULT:= Pi * angleDeg / 180.0
END { degToRad } ;


FUNCTION calcDayOfYear(mn, dy: INTEGER; lpyr: BOOLEAN): INTEGER;

//***********************************************************************/
//* Name:    calcDayOfYear						*/
//* Type:    Function							*/
//* Purpose: Finds numerical day-of-year from mn, day and lp year info  */
//* Arguments:								*/
//*   month: January = 1						*/
//*   day  : 1 - 31							*/
//*   lpyr : 1 if leap year, 0 if not					*/
//* Return value:							*/
//*   The numerical day of year						*/
//***********************************************************************/

VAR     k, doy: INTEGER;

BEGIN
  IF lpyr THEN
    k:= 1
  ELSE
    k:= 2;
  doy:= Floor((275 * mn) / 9) - k * Floor((mn + 9) / 12) + dy - 30;
  RESULT:= doy
END { calcDayOfYear } ;


FUNCTION calcDayOfWeek(juld: DOUBLE): STRING;

//***********************************************************************/
//* Name:    calcDayOfWeek						*/
//* Type:    Function							*/
//* Purpose: Derives weekday from Julian Day				*/
//* Arguments:								*/
//*   juld : Julian Day							*/
//* Return value:							*/
//*   String containing name of weekday					*/
//***********************************************************************/

VAR     A: INTEGER;
        DOW: STRING;

BEGIN
  A:= Trunc(Frac((juld + 1.5) / 7) * 7);
  CASE A OF
    0: DOW:= 'Sunday';
    1: DOW:= 'Monday';
    2: DOW:= 'Tuesday';
    3: DOW:= 'Wednesday';
    4: DOW:= 'Thursday';
    5: DOW:= 'Friday'
  ELSE
    DOW:= 'Saturday'
  END;
  RESULT:= DOW
END { calcDayOfWeek } ;


FUNCTION CalcJD(year, month, day: INTEGER): DOUBLE;

//***********************************************************************/
//* Name:    CalcJD							*/
//* Type:    Function							*/
//* Purpose: Julian day from calendar day				*/
//* Arguments:								*/
//*   year : 4 digit year						*/
//*   month: January = 1						*/
//*   day  : 1 - 31							*/
//* Return value:							*/
//*   The Julian day corresponding to the date				*/
//* Note:								*/
//*   Number is returned for start of day.  Fractional days should be	*/
//*   added later.							*/
//***********************************************************************/

VAR     A, B: INTEGER;
        JD: DOUBLE;

BEGIN
  IF month <= 2 THEN BEGIN
    year -= 1;
    month += 12
  END;
  A:= Floor(year / 100);
  B:= 2 - A + Floor(A / 4);
  JD:= Floor(365.25 * (year + 4716)) + Floor(30.6001 * (month + 1)) + day + B - 1524.5;
  RESULT:= JD
END { CalcJD } ;


TYPE    TMonthList= ARRAY[0..11] OF STRING[15];

CONST   monthName: TMonthList= ('January', 'February', 'March', 'April',
                                'May', 'June', 'July', 'August',
                                'September', 'October', 'November', 'December');
        monthAbbr: TMonthList= ('Jan', 'Feb', 'Mar', 'Apr',
                                'May', 'Jun', 'Jul', 'Aug',
                                'Sep', 'Oct', 'Nov', 'Dec');


FUNCTION calcDateFromJD(jd: DOUBLE): STRING;

//***********************************************************************/
//* Name:    calcDateFromJD						*/
//* Type:    Function							*/
//* Purpose: Calendar date from Julian Day				*/
//* Arguments:								*/
//*   jd   : Julian Day							*/
//* Return value:							*/
//*   String date in the form DD-MONTHNAME-YYYY				*/
//* Note:								*/
//***********************************************************************/

VAR     z, f, A, alpha, B, C, D, E, day, month, year: INTEGER;

BEGIN
  z:= Floor(jd + 0.5);
  f:= Round(jd + 0.5) - z;
  IF z < 2299161 THEN
    A:= z
  ELSE BEGIN
    alpha:= Floor((z - 1867216.25) / 36524.25);
    A:= z + 1 + alpha - Floor(alpha / 4.0)
  END;
  B:= A + 1524;
  C:= Floor((B - 122.1) / 365.25);
  D:= Floor(365.25 * C);
  E:= Floor((B - D) / 30.6001);
  day:= B - D - Floor(30.6001 * E) + f;
  IF E < 14 THEN
    month:= E - 1
  ELSE
    month:= E - 13;
  IF month > 2 THEN
    year:= C - 4716
  ELSE
    year:= C - 4715;
  RESULT:= IntToStr(day) + '-' + monthName[month - 1] + '-' + IntToStr(year)
END { calcDateFromJD } ;


FUNCTION calcDayFromJD(jd: DOUBLE): STRING;

//***********************************************************************/
//* Name:    calcDayFromJD						*/
//* Type:    Function							*/
//* Purpose: Calendar day (minus year) from Julian Day			*/
//* Arguments:								*/
//*   jd   : Julian Day							*/
//* Return value:							*/
//*   String date in the form DD-MONTH					*/
//***********************************************************************/

VAR     z, f, A, alpha, B, C, D, E, day, month, year: INTEGER;

BEGIN
  z:= Floor(jd + 0.5);
  f:= Round(jd + 0.5) - z;
  IF z < 2299161 THEN
    A:= z
  ELSE BEGIN
    alpha:= Floor((z - 1867216.25) / 36524.25);
    A:= z + 1 + alpha - Floor(alpha / 4.0)
  END;
  B:= A + 1524;
  C:= Floor((B - 122.1) / 365.25);
  D:= Floor(365.25 * C);
  E:= Floor((B - D) / 30.6001);
  day:= B - D - Floor(30.6001 * E) + f;
  IF E < 14 THEN
    month:= E - 1
  ELSE
    month:= E - 13;
  IF month > 2 THEN
    year:= C - 4716
  ELSE
    year:= C - 4715;
  IF day < 10 THEN
    RESULT:= '0' + IntToStr(day) + monthAbbr[month - 1]
  ELSE
    RESULT:= IntToStr(day) + monthAbbr[month - 1]
END { calcDayFromJD } ;


FUNCTION CalcTimeJulianCent(jd: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    CalcTimeJulianCent						*/
//* Type:    Function							*/
//* Purpose: convert Julian Day to centuries since J2000.0.		*/
//* Arguments:								*/
//*   jd : the Julian Day to convert					*/
//* Return value:							*/
//*   the T value corresponding to the Julian Day			*/
//***********************************************************************/

BEGIN
  RESULT:= (jd - 2451545.0) / 36525.0
END { CalcTimeJulianCent } ;


FUNCTION calcJDFromJulianCent(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcJDFromJulianCent					*/
//* Type:    Function							*/
//* Purpose: convert centuries since J2000.0 to Julian Day.		*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   the Julian Day corresponding to the t value			*/
//***********************************************************************/

BEGIN
  RESULT:= t * 36525.0 + 2451545.0
END { calcJDFromJulianCent } ;


FUNCTION calcGeomMeanLongSun(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calGeomMeanLongSun						*/
//* Type:    Function							*/
//* Purpose: calculate the Geometric Mean Longitude of the Sun		*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   the Geometric Mean Longitude of the Sun in degrees		*/
//***********************************************************************/

VAR     L0: DOUBLE;

BEGIN
  L0:= 280.46646 + t * (36000.76983 + 0.0003032 * t);
  WHILE L0 > 360.0 DO
    L0 -= 360.0;
  WHILE L0 < 0.0 DO
    L0 += 360.0;
  RESULT:= L0                   // in degrees
END { calcGeomMeanLongSun } ;


FUNCTION calcGeomMeanAnomalySun(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calGeomAnomalySun						*/
//* Type:    Function							*/
//* Purpose: calculate the Geometric Mean Anomaly of the Sun		*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   the Geometric Mean Anomaly of the Sun in degrees			*/
//***********************************************************************/

BEGIN
  RESULT:= 357.52911 + t * (35999.05029 - 0.0001537 * t) // in degrees
END { calcGeomMeanAnomalySun } ;


FUNCTION calcEccentricityEarthOrbit(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcEccentricityEarthOrbit					*/
//* Type:    Function							*/
//* Purpose: calculate the eccentricity of earth's orbit		*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   the unitless eccentricity						*/
//***********************************************************************/

BEGIN
  RESULT:= 0.016708634 - t * (0.000042037 + 0.0000001267 * t) // unitless
END { calcEccentricityEarthOrbit } ;


FUNCTION calcSunEqOfCenter(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcSunEqOfCenter						*/
//* Type:    Function							*/
//* Purpose: calculate the equation of center for the sun		*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   in degrees							*/
//***********************************************************************/

VAR     m, mrad, sinm, sin2m, sin3m: DOUBLE;

BEGIN
  m:= calcGeomMeanAnomalySun(t);
  mrad:= degToRad(m);
  sinm:= Sin(mrad);
  sin2m:= Sin(2 * mrad);
  sin3m:= Sin(3 * mrad);
  RESULT:= sinm * (1.914602 - t * (0.004817 + 0.000014 * t)) + sin2m *
        (0.019993 - 0.000101 * t) + sin3m * 0.000289 // in degrees
END { calcSunEqOfCenter } ;


FUNCTION calcSunTrueLong(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcSunTrueLong						*/
//* Type:    Function							*/
//* Purpose: calculate the true longitude of the sun			*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   sun's true longitude in degrees					*/
//***********************************************************************/

VAR     l0, c: DOUBLE;

BEGIN
  l0:= calcGeomMeanLongSun(t);
  c:= calcSunEqOfCenter(t);
  RESULT:= l0 + c               // in degrees
END { calcSunTrueLong } ;


FUNCTION calcSunTrueAnomaly(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcSunTrueAnomaly						*/
//* Type:    Function							*/
//* Purpose: calculate the true anamoly of the sun			*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   sun's true anamoly in degrees					*/
//***********************************************************************/

VAR     m, c: DOUBLE;

BEGIN
  m:= calcGeomMeanAnomalySun(t);
  c:= calcSunEqOfCenter(t);
  RESULT:= m + c                // in degrees
END { calcSunTrueAnomaly } ;


FUNCTION calcSunRadVector(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcSunRadVector						*/
//* Type:    Function							*/
//* Purpose: calculate the distance to the sun in AU			*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   sun radius vector in AUs						*/
//***********************************************************************/

VAR     v, e: DOUBLE;

BEGIN
  v:= calcSunTrueAnomaly(t);
  e:= calcEccentricityEarthOrbit(t);
  RESULT:= (1.000001018 * (1 - e * e)) / (1 + e * Cos(degToRad(v))) // in AUs
END { calcSunRadVector } ;


FUNCTION calcSunApparentLong(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcSunApparentLong					*/
//* Type:    Function							*/
//* Purpose: calculate the apparent longitude of the sun		*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   sun's apparent longitude in degrees				*/
//***********************************************************************/

VAR     o, omega: DOUBLE;

BEGIN
  o:= calcSunTrueLong(t);
  omega:= 125.04 - 1934.136 * t;
  RESULT:= o - 0.00569 - 0.00478 * Sin(degToRad(omega)) // in degrees
END { calcSunApparentLong } ;


FUNCTION calcMeanObliquityOfEcliptic(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcMeanObliquityOfEcliptic				*/
//* Type:    Function							*/
//* Purpose: calculate the mean obliquity of the ecliptic		*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   mean obliquity in degrees						*/
//***********************************************************************/

VAR     seconds: DOUBLE;

BEGIN
  seconds:= 21.448 - t * (46.8150 + t * (0.00059 - t * (0.001813)));
  RESULT:= 23.0 + (26.0 + (seconds / 60.0)) / 60.0 // in degrees
END { calcMeanObliquityOfEcliptic } ;


FUNCTION calcObliquityCorrection(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcObliquityCorrection					*/
//* Type:    Function							*/
//* Purpose: calculate the corrected obliquity of the ecliptic		*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   corrected obliquity in degrees					*/
//***********************************************************************/

VAR     e0, omega: DOUBLE;

BEGIN
  e0:= calcMeanObliquityOfEcliptic(t);
  omega:= 125.04 - 1934.136 * t;
  RESULT:= e0 + 0.00256 * Cos(degToRad(omega)) // in degrees
END { calcObliquityCorrection } ;


FUNCTION calcSunRtAscension(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcSunRtAscension						*/
//* Type:    Function							*/
//* Purpose: calculate the right ascension of the sun			*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   sun's right ascension in degrees					*/
//***********************************************************************/

VAR     e, lambda, tananum, tanadenom: DOUBLE;

BEGIN
  e:= calcObliquityCorrection(t);
  lambda:= calcSunApparentLong(t);
  tananum:= Cos(degToRad(e)) * Sin(degToRad(lambda));
  tanadenom:= Cos(degToRad(lambda));
  RESULT:= radToDeg(ArcTan2(tananum, tanadenom)) // in degrees
END { calcSunRtAscension } ;


FUNCTION calcSunDeclination(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcSunDeclination						*/
//* Type:    Function							*/
//* Purpose: calculate the declination of the sun			*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   sun's declination in degrees					*/
//***********************************************************************/

VAR     e, lambda, sint: DOUBLE;

BEGIN
  e:= calcObliquityCorrection(t);
  lambda:= calcSunApparentLong(t);
  sint:= Sin(degToRad(e)) * Sin(degToRad(lambda));
  RESULT:= radToDeg(ArcSin(sint))       // in degrees
END { calcSunDeclination } ;


FUNCTION calcEquationOfTime(t: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcEquationOfTime						*/
//* Type:    Function							*/
//* Purpose: calculate the difference between true solar time and mean	*/
//*		solar time						*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//* Return value:							*/
//*   equation of time in minutes of time				*/
//***********************************************************************/

VAR     epsilon, l0,  e, m, y, sin2l0, sinm, cos2l0, sin4l0, sin2m, Etime: DOUBLE;

BEGIN
  epsilon:= calcObliquityCorrection(t);
  l0:= calcGeomMeanLongSun(t);
  e:= calcEccentricityEarthOrbit(t);
  m:= calcGeomMeanAnomalySun(t);
  y:= Tan(degToRad(epsilon) / 2.0);
  y *= y;
  sin2l0:= Sin(2.0 * degToRad(l0));
  sinm:= Sin(degToRad(m));
  cos2l0:= Cos(2.0 * degToRad(l0));
  sin4l0:= Sin(4.0 * degToRad(l0));
  sin2m:= Sin(2.0 * degToRad(m));
  Etime:= y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0
                        - 0.5 * y * y * sin4l0 - 1.25 * e * e * sin2m;
  RESULT:= radToDeg(Etime) * 4.0        // in minutes of time
END { calcEquationOfTime } ;


FUNCTION calcHourAngleSunrise(lat, solarDec: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcHourAngleSunrise					*/
//* Type:    Function							*/
//* Purpose: calculate the hour angle of the sun at sunrise for the	*/
//*			latitude					*/
//* Arguments:								*/
//*   lat : latitude of observer in degrees				*/
//*	solarDec : declination angle of sun in degrees			*/
//* Return value:							*/
//*   hour angle of sunrise in radians					*/
//***********************************************************************/

VAR     latRad, sdRad, HAarg: DOUBLE;

BEGIN
  latRad:= degToRad(lat);
  sdRad:= degToRad(solarDec);
  HAarg:= Cos(degToRad(90.833)) / (Cos(latRad) * Cos(sdRad)) - Tan(latRad) * Tan(sdRad);
  RESULT:= ArcCos(HAarg)        // in radians
END { calcHourAngleSunrise } ;


FUNCTION calcHourAngleSunset(lat, solarDec: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    calcHourAngleSunset					*/
//* Type:    Function							*/
//* Purpose: calculate the hour angle of the sun at sunset for the	*/
//*			latitude					*/
//* Arguments:								*/
//*   lat : latitude of observer in degrees				*/
//*	solarDec : declination angle of sun in degrees			*/
//* Return value:							*/
//*   hour angle of sunset in radians					*/
//***********************************************************************/

VAR     latRad, sdRad, HAarg: DOUBLE;

BEGIN
  latRad:= degToRad(lat);
  sdRad:= degToRad(solarDec);
  HAarg:= Cos(degToRad(90.833)) / (Cos(latRad) * Cos(sdRad)) - Tan(latRad) * Tan(sdRad);
  RESULT:= -ArcCos(HAarg);              // in radians
END { calcHourAngleSunset } ;


FUNCTION CalcSolNoonUTC(t, longitude: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    CalcSolNoonUTC						*/
//* Type:    Function							*/
//* Purpose: calculate the Universal Coordinated Time (UTC) of solar	*/
//*		noon for the given day at the given location on earth	*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//*   longitude : longitude of observer in degrees			*/
//* Return value:							*/
//*   time in minutes from zero Z					*/
//***********************************************************************/

VAR     tnoon, eqTime, solNoonUTC, newt: DOUBLE;

BEGIN
// First pass uses approximate solar noon to calculate eqtime
  tnoon:= CalcTimeJulianCent(calcJDFromJulianCent(t) + longitude / 360.0);
  eqTime:= calcEquationOfTime(tnoon);
  solNoonUTC:= 720 + (longitude * 4) - eqTime; // min
  newt:= CalcTimeJulianCent(calcJDFromJulianCent(t) -0.5 + solNoonUTC / 1440.0);
  eqTime:= calcEquationOfTime(newt);
  solNoonUTC:= 720 + (longitude * 4) - eqTime; // min
  RESULT:= solNoonUTC
END { CalcSolNoonUTC } ;


FUNCTION CalcSunriseUTC(JD, latitude, longitude: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    CalcSunriseUTC						*/
//* Type:    Function							*/
//* Purpose: calculate the Universal Coordinated Time (UTC) of sunrise	*/
//*			for the given day at the given location on earth*/
//* Arguments:								*/
//*   JD  : julian day							*/
//*   latitude : latitude of observer in degrees			*/
//*   longitude : longitude of observer in degrees			*/
//* Return value:							*/
//*   time in minutes from zero Z					*/
//***********************************************************************/

VAR     t, noonmin, tnoon, eqTime, solarDec,  hourAngle, delta, timeDiff,
                                                timeUTC, newt: DOUBLE;
BEGIN
  t:= CalcTimeJulianCent(JD);

// *** Find the time of solar noon at the location, and use that declination

  noonmin:= CalcSolNoonUTC(t, longitude);
  tnoon:= CalcTimeJulianCent (JD + noonmin / 1440.0);

// *** First pass to approximate sunrise (using solar noon)

  eqTime:= calcEquationOfTime(tnoon);
  solarDec:= calcSunDeclination(tnoon);
  hourAngle:= calcHourAngleSunrise(latitude, solarDec);
  delta:= longitude - radToDeg(hourAngle);
  timeDiff:= 4 * delta;
  timeUTC:= 720 + timeDiff - eqTime;

// *** Second pass includes fractional jday in gamma calc

  newt:= CalcTimeJulianCent(calcJDFromJulianCent(t) + timeUTC / 1440.0);
  eqTime:= calcEquationOfTime(newt);
  solarDec:= calcSunDeclination(newt);
  hourAngle:= calcHourAngleSunrise(latitude, solarDec);
  delta:= longitude - radToDeg(hourAngle);
  timeDiff:= 4 * delta;
  timeUTC:= 720 + timeDiff - eqTime;
  RESULT:= timeUTC              // in minutes
END { CalcSunriseUTC } ;


FUNCTION CalcSunsetUTC(JD, latitude, longitude: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    CalcSunsetUTC						*/
//* Type:    Function							*/
//* Purpose: calculate the Universal Coordinated Time (UTC) of sunset	*/
//*			for the given day at the given location on earth*/
//* Arguments:								*/
//*   JD  : julian day							*/
//*   latitude : latitude of observer in degrees			*/
//*   longitude : longitude of observer in degrees			*/
//* Return value:							*/
//*   time in minutes from zero Z					*/
//***********************************************************************/

VAR     t, noonmin, tnoon, eqTime, solarDec, hourAngle, delta,
                                        timeDiff, timeUTC, newt: DOUBLE;
BEGIN
  t:= CalcTimeJulianCent(JD);

// *** Find the time of solar noon at the location, and use that declination

  noonmin:= CalcSolNoonUTC(t, longitude);
  tnoon:= CalcTimeJulianCent(JD + noonmin / 1440.0);

// First calculates sunrise and approx length of day

  eqTime:= calcEquationOfTime(tnoon);
  solarDec:= calcSunDeclination(tnoon);
  hourAngle:= calcHourAngleSunset(latitude, solarDec);
  delta:= longitude - radToDeg(hourAngle);
  timeDiff:= 4 * delta;
  timeUTC:= 720 + timeDiff - eqTime;

// first pass used to include fractional day in gamma calc

  newt:= CalcTimeJulianCent(calcJDFromJulianCent(t) + timeUTC / 1440.0);
  eqTime:= calcEquationOfTime(newt);
  solarDec:= calcSunDeclination(newt);
  hourAngle:= calcHourAngleSunset(latitude, solarDec);
  delta:= longitude - radToDeg(hourAngle);
  timeDiff:= 4 * delta;
  timeUTC:= 720 + timeDiff - eqTime;    // in minutes
  RESULT:= timeUTC
END { CalcSunsetUTC } ;


FUNCTION CalcSunElevation(t, utcmins, latitude, longitude: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    CalcSunElevation						*/
//* Type:    Function							*/
//* Purpose: calculate the elevation (altitude) of the sun		*/
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//*   utcmins: UTC time in minutes                                      */
//*   latitude : latitude of observer in degrees			*/
//*   longitude : longitude of observer in degrees			*/
//* Return value:							*/
//*   sun's elevation/altitude in degrees				*/
//***********************************************************************/

(* This is not an original NOAA function, but is extracted from the calcSun()   *)
(* function in http://www.srrb.noaa.gov/highlights/sunrise/azel.html so that    *)
(* the current twilight state can be computed. MarkMLl.                         *)

VAR     eqTime, solarDec, solarTimeFix, trueSolarTime, hourAngle, haRad,
                csz, zenith, exoatmElevation, refractionCorrection,
                te, solarZen: DOUBLE;

BEGIN
  eqTime:= calcEquationOfTime(t);
  solarDec:= calcSunDeclination(t);
  solarTimeFix:= eqTime - 4.0 * longitude (* + 60.0 * zone *) ;
  trueSolarTime:= utcmins + solarTimeFix;       // in minutes
  WHILE trueSolarTime > 1440 DO
    trueSolarTime -= 1440;
  hourAngle:= trueSolarTime / 4.0 - 180.0;
  // Thanks to Louis Schwarzmayr for finding our error,
  // and providing the following 2 lines to fix it:
  IF hourAngle < -180 THEN
    hourAngle += 360.0;
  haRad:= degToRad(hourAngle);
  csz:= Sin(degToRad(latitude)) *
  	Sin(degToRad(solarDec)) +
  	Cos(degToRad(latitude)) *
  	Cos(degToRad(solarDec)) * Cos(haRad);
  IF csz > 1.0 THEN
    csz:= 1.0
  ELSE
    IF csz < -1.0 THEN
      csz:= -1.0;
  zenith:= radToDeg(ArcCos(csz));

  exoatmElevation:= 90.0 - zenith;
  IF exoatmElevation > 85.0 THEN
    refractionCorrection:= 0.0
  ELSE BEGIN
    te:= Tan(degToRad(exoatmElevation));
    IF exoatmElevation > 5.0 THEN
      refractionCorrection:= 58.1 / te - 0.07 / (te*te*te) +
  			                0.000086 / (te*te*te*te*te)
    ELSE
      IF exoatmElevation > -0.575 THEN
  	refractionCorrection:= 1735.0 + exoatmElevation *
  			(-518.2 + exoatmElevation * (103.4 +
  			exoatmElevation * (-12.79 +
  			exoatmElevation * 0.711) ) )
      ELSE
  	refractionCorrection:= -20.774 / te;
    refractionCorrection:= refractionCorrection / 3600.0
  END;
  solarZen:= zenith - refractionCorrection;
  RESULT:= 90.0 - solarZen
END { CalcSunElevation } ;


FUNCTION CalcSunAzimuth(t, utcmins, latitude, longitude: DOUBLE): DOUBLE;

//***********************************************************************/
//* Name:    CalcSunAzimuth						*/
//* Type:    Function							*/
//* Purpose: calculate the azimuth of the sun		                */
//* Arguments:								*/
//*   t : number of Julian centuries since J2000.0			*/
//*   utcmins: UTC time in minutes                                      */
//*   latitude : latitude of observer in degrees			*/
//*   longitude : longitude of observer in degrees			*/
//* Return value:							*/
//*   sun's azimuth in degrees				                */
//***********************************************************************/

(* This is not an original NOAA function, but is extracted from the calcSun()   *)
(* function in http://www.srrb.noaa.gov/highlights/sunrise/azel.html for        *)
(* completeness. MarkMLl.                                                       *)

// ***** NOT TESTED *****

VAR     eqTime, solarDec, solarTimeFix, trueSolarTime, hourAngle, haRad,
                csz, zenith, azDenom, azRad, azimuth: DOUBLE;

BEGIN
  eqTime:= calcEquationOfTime(t);
  solarDec:= calcSunDeclination(t);
  solarTimeFix:= eqTime - 4.0 * longitude (* + 60.0 * zone *) ;
  trueSolarTime:= utcmins + solarTimeFix;       // in minutes
  WHILE trueSolarTime > 1440 DO
    trueSolarTime -= 1440;
  hourAngle:= trueSolarTime / 4.0 - 180.0;
  // Thanks to Louis Schwarzmayr for finding our error,
  // and providing the following 2 lines to fix it:
  IF hourAngle < -180 THEN
    hourAngle += 360.0;
  haRad:= degToRad(hourAngle);
  csz:= Sin(degToRad(latitude)) *
  	Sin(degToRad(solarDec)) +
  	Cos(degToRad(latitude)) *
  	Cos(degToRad(solarDec)) * Cos(haRad);
  IF csz > 1.0 THEN
    csz:= 1.0
  ELSE
    IF csz < -1.0 THEN
      csz:= -1.0;
  zenith:= radToDeg(ArcCos(csz));

  azDenom:= Cos(degToRad(latitude)) * Sin(degToRad(zenith));
  IF Abs(azDenom) > 0.001 THEN BEGIN
    azRad:= (( Sin(degToRad(latitude)) *
  		Cos(degToRad(zenith)) ) -
  		Sin(degToRad(solarDec))) / azDenom;
    IF Abs(azRad) > 1.0 THEN
      IF azRad < 0 THEN
  	azRad:= -1.0
      ELSE
  	azRad:= 1.0;
    Azimuth:= 180.0 - radToDeg(ArcCos(azRad));
    IF hourAngle > 0.0 THEN
      azimuth:= -azimuth
    ELSE
      IF latitude > 0.0 THEN
  	azimuth:= 180.0
      ELSE
  	azimuth:= 0.0
  END;
  IF azimuth < 0.0 THEN
    azimuth += 360.0;
  RESULT:= azimuth
END { CalcSunAzimuth } ;


end.

