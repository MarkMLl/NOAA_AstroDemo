(* Lazarus + FPC 0.9.30 + 2.4.4. On Linux for ARM, PPC, SPARC, x86. Lazarus + F *)

(* Some of the maths might require 2.8.0 on SPARC for reliable operation.       *)

unit MoonPhase;

(* This is a transcription of Perl's Astro::MoonPhase.pm, which I believe is    *)
(* derived from John Walker's moontool.c.                                       *)

// ABOUT THE ALGORITHMS:
//
// MoonPhase calculates information about the phase of the moon
// at a given time.
//
// The algorithms used in this program to calculate the positions of Sun and
// Moon as seen from the Earth are given in the book Practical Astronomy
// With  Your  Calculator  by  Peter  Duffett-Smith,   Second   Edition,
// Cambridge University Press, 1981.  Ignore the word "Calculator" in the
// title;  this  is  an  essential  reference  if  you're  interested  in
// developing  software  which  calculates  planetary  positions, orbits,
// eclipses, and  the  like.   If  you're  interested  in  pursuing  such
// programming, you should also obtain:
//
// Astronomical  Formulae for Calculators by Jean Meeus, Third Edition,
// Willmann-Bell, 1985.  A must-have.
//
// Planetary  Programs  and  Tables  from  -4000  to  +2800  by  Pierre
// Bretagnon  and Jean-Louis Simon, Willmann-Bell, 1986.  If you want the
// utmost  (outside  of  JPL)  accuracy  for  the  planets,  it's   here.
//
// Celestial BASIC by Eric Burgess, Revised Edition, Sybex, 1985.  Very
// cookbook oriented, and many of the algorithms are hard to dig  out  of
// the turgid BASIC code, but you'll probably want it anyway.
//
// Many of these references can be obtained from Willmann-Bell, P.O.  Box
// 35025,  Richmond,  VA 23235, USA.  Phone: (804) 320-7016.  In addition
// to their own publications, they stock most of the standard  references
// for mathematical and positional astronomy.
//
// LICENCE:
//
// This  program is in the public domain: "Do what thou wilt shall be the
// whole of the law".
//
// NOTE: the above is the license text as found in MoonPhase.pm. The original
// license from John Walker as found in moontool.c is as below:
//
// This  program is in the public domain: "Do what thou wilt shall be the
// whole of the law".  I'd appreciate  receiving  any  bug  fixes  and/or
// enhancements,  which  I'll  incorporate  in  future  versions  of  the
// program.  Please leave the original attribution information intact  so
// that credit and blame may be properly apportioned.
//
// AUTHORS:
//
// The moontool.c Release 2.0:
//
//     A Moon for the Sun
//     Designed and implemented by John Walker in December 1987,
//     revised and updated in February of 1988.
//
// Initial Perl transcription:
//
//     Raino Pikkarainen, 1998
//     raino.pikkarainen@saunalahti.fi
//
// The moontool.c Release 2.4:
//
//     Major enhancements by Ron Hitchens, 1989
//
// Revisions:
//
//     Brett Hamilton  http://simple.be/
//     Bug fix, 2003
//     Second transcription and bugfixes, 2004
//
//     Christopher J. Madsen  http://www.cjmweb.net/
//     Added phaselist function, March 2007

(* Transcribed to Pascal by Mark Morgan Lloyd.                                  *)
{$hints off}{$notes off}

{$mode objfpc}

interface

uses
  Classes, SysUtils; 

CONST      UnixNow= Int64(-1);
           UnixToday= Int64(-2);
           UnixAddDay= Int64(-3);
           SecondsPerDay= 3600 * 24;

FUNCTION CalcUS(year, month, day: INTEGER): Int64;
FUNCTION CalcUS: Int64;

(* Calculate the unix seconds (i.e. based on an epoch of 1970-01-01 00:00:00)   *)
(* for the indicated year, month (based 1) and day (based 1). If parameters are *)
(* omitted calculate for the start of the current day (i.e. 00:00:00.0).        *)

TYPE       TPhaseResult= RECORD
                           MoonPhase, MoonIllum, MoonAge, MoonDist,
                           MoonAng, SunDist, SunAng: DOUBLE;
                           seconds_since_1970: Int64
                         END;

FUNCTION Phase(seconds_since_1970: Int64= UnixNow): TPhaseResult;
FUNCTION Phase(seconds_since_1970: Int64= UnixNow): DOUBLE;

(* Return either an array containing a description of the Moon's current        *)
(* appearance, or a single number giving its phase angle. In either case if the *)
(* parameter is omitted it is assumed to be the current date and time.          *)
(*                                                                              *)
(* Description below uses Perl notation.                                        *)

// The argument is the time for which the phase is requested,
// expressed as a time returned by the unix time function. If $seconds_since_1970
// is omitted, it does phase(time).
//
// Return value in scalar context is $MoonPhase,
// the terminator phase angle as a percentage of a full circle (i.e., 0 to 1).
//
// Return values in array context:
//
// $MoonPhase:
//
// the terminator phase angle as a percentage of a full circle (i.e., 0 to 1)
//
// $MoonIllum:
//
// the illuminated fraction of the Moon's disc
//
// $MoonAge:
//
// the Moon's age in days and fraction
//
// $MoonDist:
//
// the distance of the Moon from the centre of the Earth
//
// $MoonAng:
//
// the angular diameter subtended by the Moon as seen by
// an observer at the centre of the Earth
//
// $SunDist:
//
// the distance from the Sun in km
//
// $SunAng:
//
// the angular size of Sun in degrees
//
// Example:
//
//    ( $MoonPhase,
//      $MoonIllum,
//      $MoonAge,
//      $MoonDist,
//      $MoonAng,
//      $SunDist,
//      $SunAng ) = phase();
//
//      print "MoonPhase  = $MoonPhase\n";
//      print "MoonIllum  = $MoonIllum\n";
//      print "MoonAge    = $MoonAge\n";
//      print "MoonDist   = $MoonDist\n";
//      print "MoonAng    = $MoonAng\n";
//      print "SunDist    = $SunDist\n";
//      print "SunAng     = $SunAng\n";
//
// could print something like this:
//
//      MoonPhase  = 0.598939375319023
//      MoonIllum  = 0.906458030827876
//      MoonAge    = 17.6870323368022
//      MoonDist   = 372479.357420033
//      MoonAng    = 0.534682403555093
//      SunDist    = 152078368.820205
//      SunAng     = 0.524434538105092
//
// Knowing the month from PhaseHunt() (below), we can determine that these
// values correspond to Sun Jul 12 11:37:01 1998 UTC.

TYPE    TPhaseHuntResult= ARRAY[0..4] OF Int64;

FUNCTION PhaseHunt(seconds_since_1970: Int64= UnixNow): TPhaseHuntResult;

(* Return an array containing the date and times of the five significant phases *)
(* surrounding the indicated date and time, or surrounding the current date and *)
(* time if the parameter is omitted.                                            *)
(*                                                                              *)
(* Description below uses Perl notation.                                        *)

// Finds time of phases of the moon which surround the given
// date.  Five phases are found, starting and ending with the
// new moons which bound the current lunation.
//
// The argument is the time, expressed as a time returned
// by the unix time function. If $seconds_since_1970
// is omitted, it does phasehunt(time).
//
// Example:
//
//     @phases = phasehunt();
//     print "New Moon      = ", scalar(localtime($phases[0])), "\n";
//     print "First quarter = ", scalar(localtime($phases[1])), "\n";
//     print "Full moon     = ", scalar(localtime($phases[2])), "\n";
//     print "Last quarter  = ", scalar(localtime($phases[3])), "\n";
//     print "New Moon      = ", scalar(localtime($phases[4])), "\n";
//
// could print something like this:
//
//     New Moon      = Wed Jun 24 06:51:47 1998
//     First quarter = Wed Jul  1 21:42:19 1998
//     Full moon     = Thu Jul  9 19:02:47 1998
//     Last quarter  = Thu Jul 16 18:15:18 1998
//     New Moon      = Thu Jul 23 16:45:01 1998
//
// These values appear to correspond to Sun Jul 12 11:37:01 1998 UTC, but the
// timezone used for display is not GMT.

TYPE   TPhaseListResult= RECORD
                           phase: Int64;
                           times: ARRAY OF Int64
                         END;

FUNCTION PhaseList(start: Int64= UnixToday; stop: Int64= UnixAddDay): TPhaseListResult;

(* Return the first significant phase of the Moon, plus its date and time       *)
(* together with the date and time of any others between the parameters. If the *)
(* first parameter is omitted it is assumed to be the start of the current day, *)
(* if the second parameter is omitted it is assumed to be 24 hours after the    *)
(* first parameter.                                                             *)
(*                                                                              *)
(* Description below uses Perl notation.                                        *)

// Finds times of all phases of the moon which occur on or after
// $start but before $stop.  Both the arguments and the return
// values are expressed as seconds since 1970 (like the unix time function
// returns).
//
// $phase is an integer indicating the phase of the moon at
// $times[0], as shown in this table:
//
//     0  New Moon
//     1  First quarter
//     2  Full Moon
//     3  Last quarter
//
// The remaining values in @times indicate subsequent phases of the
// moon (in ascending order by time).  If there are no phases of the moon
// between $start and $stop, phaselist returns the empty list.
//
// Example:
//
//     @name = ("New Moon", "First quarter", "Full moon", "Last quarter");
//     ($phase, @times) = phaselist($start, $stop);
//
//     while (@times) {
//       printf "%-14s= %s\n", $name[$phase], scalar localtime shift @times;
//       $phase = ($phase + 1) % 4;
//     }
//
// could produce the same output as the phasehunt example above (given
// the appropriate start & stop times).


implementation

USES DateUtils, Math;

FUNCTION CalcUS(year, month, day: INTEGER): Int64;

(* Calculate the unix seconds (i.e. based on an epoch of 1970-01-01 00:00:00)   *)
(* for the indicated year, month (based 1) and day (based 1). If parameters are *)
(* omitted calculate for the start of the current day (i.e. 00:00:00.0).        *)

BEGIN
  RESULT:= DateTimeToUnix(EncodeDateTime(year, month, day, 0, 0, 0, 0))
END { CalcUS } ;


FUNCTION CalcUS: Int64;

(* Calculate the unix seconds (i.e. based on an epoch of 1970-01-01 00:00:00)   *)
(* for the indicated year, month (based 1) and day (based 1). If parameters are *)
(* omitted calculate for the start of the current day (i.e. 00:00:00.0).        *)

BEGIN
  RESULT:= DateTimeToUnix(DateOf(Now))
END { CalcUS } ;


CONST   Epoch		= 2444238.5;		// 1980 January 0.0

// Constants defining the Sun's apparent orbit.

CONST   Elonge		= 278.833540;		// ecliptic longitude of the Sun at epoch 1980.0
        Elongp		= 282.596403;		// ecliptic longitude of the Sun at perigee
        Eccent		= 0.016718;		// eccentricity of Earth's orbit
        Sunsmax		= 1.495985e8;		// semi-major axis of Earth's orbit, km
        Sunangsiz	= 0.533128;		// sun's angular size, degrees, at semi-major axis distance

// Elements of the Moon's orbit, epoch 1980.0.

        Mmlong		= 64.975464;		// moon's mean longitude at the epoch
        Mmlongp		= 349.383063;		// mean longitude of the perigee at the epoch
        Mlnode		= 151.950429;		// mean longitude of the node at the epoch
        Minc		= 5.145396;		// inclination of the Moon's orbit
        Mecc		= 0.054900;		// eccentricity of the Moon's orbit
        Mangsiz		= 0.5181;		// moon's angular size at distance a from Earth
        Msmax		= 384401.0;		// semi-major axis of Moon's orbit in km
        Mparallax	= 0.9507;		// parallax at distance a from Earth
        Synmonth	= 29.53058868;		// synodic month (new Moon to new Moon)

// Handy mathematical functions.

// Assume that Pi is available from the FP runtimes, and that the standard Sign()
// and Floor() functions replace sgn() and floor(), see initialisation section
// for tests of compatible operation.


FUNCTION fixAngle(d: DOUBLE): DOUBLE;

{ return ($_[0] - 360.0 * (floor($_[0] / 360.0))); }	// fix angle

BEGIN
  RESULT:= d - 360.0 * (Floor(d / 360.0))
END { fixAngle } ;


FUNCTION toRad(d: DOUBLE): DOUBLE;

BEGIN
  RESULT:= d * (Pi / 180.0)                             // deg->rad
END { toRad } ;


FUNCTION toDeg(r: DOUBLE): DOUBLE;

BEGIN
  RESULT:= r * (180.0 / Pi)                             // rad->deg
END { toDeg } ;


FUNCTION dSin(d: DOUBLE): DOUBLE;

BEGIN
  RESULT:= Sin(toRad(d))                                // sin from deg
END { dSin } ;


FUNCTION dCos(d: DOUBLE): DOUBLE;

BEGIN
  RESULT:= Cos(toRad(d))                                // cos from deg
END { dCos } ;


FUNCTION jTime(t: Int64): DOUBLE;

// jtime - convert internal date and time to astronomical Julian
// time (i.e. Julian date plus day fraction)

BEGIN
  RESULT:= t;
  RESULT:= (RESULT / 86400.0) + 2440587.5	// (seconds /(seconds per day)) + julian date of epoch
END { jTime } ;


FUNCTION jDayToSecs(jday: DOUBLE): Int64;

// jdaytosecs - convert Julian date to a UNIX epoch

BEGIN
  RESULT:= Round((jday - 2440587.5) * 86400)   // (juliandate - jdate of unix epoch)*(seconds per julian day)
END { jDayToSecs } ;


PROCEDURE jYear(td: DOUBLE; VAR yy, mm: Int64; VAR dd: DOUBLE);

// jyear - convert Julian date to year, month, day, which are
// returned via integer pointers to integers

VAR     z, a, alpha, b, c, d, e: Int64;
        f: DOUBLE;

BEGIN
  td:= td + 0.5;                       // astronomical to civil
  z:= Floor(td);
  f:= td - z;
  IF z < 2299161 THEN
    a:= z
  ELSE BEGIN
    alpha:= Floor((z - 1867216.25) / 36524.25);
    a:= z + 1 + alpha - Floor(alpha / 4)
  END;
  b:= a + 1524;
  c:= Floor((b - 122.1) / 365.25);
  d:= Floor(365.25 * c);
  e:= Floor((b - d) / 30.6001);

(* Note here: both MoonPhase.pm and moontool.c lose fractional days at this     *)
(* point since dd is declared as an integer; Perl programmers might not in      *)
(* fact notice this due to loose type handling. I think this is a bug since it  *)
(* prevents using non-GMT time ranges. MarkMLl.                                 *)

  dd:= b - d - Floor(30.6001 * e) + f;
  IF e < 14 THEN
    mm:= e - 1
  ELSE
    mm:= e - 13;
  IF mm > 2 THEN
    yy:= c - 4716
  ELSE
    yy:= c - 4715
END { jYear } ;


PROCEDURE jYear(td: DOUBLE; VAR yy, mm, dd: Int64);

// jyear - convert Julian date to year, month, day, which are
// returned via integer pointers to integers

(* This is equivalent to the original jYear() which loses fractional days.      *)
(* MarkMLl.                                                                     *)

VAR     y, m: Int64;
        d: DOUBLE;

BEGIN
  jYear(td, y, m, d);
  yy:= y;
  mm:= m;
  dd:= Trunc(d)
END { jYear } ;


FUNCTION meanPhase(sdate, k: DOUBLE): DOUBLE;

//  meanphase  --  Calculates  time  of  the mean new Moon for a given
//                 base date.  This argument K to this function is the
//                 precomputed synodic month index, given by:
//
//                        K = (year - 1900) * 12.3685
//
//                 where year is expressed as a year and fractional year.

VAR     t, t2, t3: DOUBLE;

BEGIN
  t:= (sdate - 2415020.0) / 36525; // Time in Julian centuries from 1900 January 0.5
  t2:= t * t;                      // Square for frequent use
  t3:= t2 * t;                     // Cube for frequent use
  RESULT:= 2415020.75933 + Synmonth * k
         + 0.0001178 * t2
         - 0.000000155 * t3
         + 0.00033 * dsin(166.56 + 132.87 * t - 0.009173 * t2)
END { meanPhase } ;


FUNCTION truePhase(k, phase: DOUBLE): DOUBLE;

// truephase - given a K value used to determine the mean phase of the
// new moon, and a phase selector (0.0, 0.25, 0.5, 0.75),
// obtain the true, corrected phase time

VAR     t, t2, t3, pt, m, mprime, f: DOUBLE;
        apcor: BOOLEAN= FALSE;

BEGIN
  k += phase;				// add phase to new moon time
  t:= k / 1236.85;			// time in Julian centuries from 1900 January 0.5
  t2:= t * t;				// square for frequent use
  t3:= t2 * t;				// cube for frequent use

// mean time of phase

  pt:= 2415020.75933
   + Synmonth * k
   + 0.0001178 * t2
   - 0.000000155 * t3
   + 0.00033 * dsin(166.56 + 132.87 * t - 0.009173 * t2);

// Sun's mean anomaly

  m:= 359.2242
  + 29.10535608 * k
  - 0.0000333 * t2
  - 0.00000347 * t3;

// Moon's mean anomaly

  mprime:= 306.0253
  + 385.81691806 * k
  + 0.0107306 * t2
  + 0.00001236 * t3;

// Moon's argument of latitude

  f:= 21.2964
  + 390.67050646 * k
  - 0.0016528 * t2
  - 0.00000239 * t3;

// Corrections for New and Full Moon.

  IF (phase < 0.01) OR (abs(phase - 0.5) < 0.01) THEN BEGIN
     pt += (0.1734 - 0.000393 * t) * dsin(m)
        + 0.0021 * dsin(2 * m)
        - 0.4068 * dsin(mprime)
        + 0.0161 * dsin(2 * mprime)
        - 0.0004 * dsin(3 * mprime)
        + 0.0104 * dsin(2 * f)
        - 0.0051 * dsin(m + mprime)
        - 0.0074 * dsin(m - mprime)
        + 0.0004 * dsin(2 * f + m)
        - 0.0004 * dsin(2 * f - m)
        - 0.0006 * dsin(2 * f + mprime)
        + 0.0010 * dsin(2 * f - mprime)
        + 0.0005 * dsin(m + 2 * mprime);
    apcor:= TRUE
  END ELSE
    IF (abs(phase - 0.25) < 0.01) OR (abs(phase - 0.75) < 0.01) THEN BEGIN
      pt += (0.1721 - 0.0004 * t) * dsin(m)
        + 0.0021 * dsin(2 * m)
        - 0.6280 * dsin(mprime)
        + 0.0089 * dsin(2 * mprime)
        - 0.0004 * dsin(3 * mprime)
        + 0.0079 * dsin(2 * f)
        - 0.0119 * dsin(m + mprime)
        - 0.0047 * dsin(m - mprime)
        + 0.0003 * dsin(2 * f + m)
        - 0.0004 * dsin(2 * f - m)
        - 0.0006 * dsin(2 * f + mprime)
        + 0.0021 * dsin(2 * f - mprime)
        + 0.0003 * dsin(m + 2 * mprime)
        + 0.0004 * dsin(m - 2 * mprime)
        - 0.0003 * dsin(2 * m + mprime);
      IF (phase < 0.5) THEN            // First quarter correction.
  	pt += 0.0028 - 0.0004 * dcos(m) + 0.0003 * dcos(mprime)
      ELSE                              // Last quarter correction.
  	pt += -0.0028 + 0.0004 * dcos(m) - 0.0003 * dcos(mprime);
      apcor:= TRUE
    END;
  IF NOT apcor THEN BEGIN
    RAISE Exception.Create('truePhase() called with invalid phase selector (' + FloatToStr(phase) + ').')
  END;
  RESULT:= pt
END { truePhase } ;


FUNCTION PhaseHunt(seconds_since_1970: Int64= UnixNow): TPhaseHuntResult;

(* Return an array containing the date and times of the five significant phases *)
(* surrounding the indicated date and time, or surrounding the current date and *)
(* time if the parameter is omitted.                                            *)

// phasehunt - find time of phases of the moon which surround the current
// date.  Five phases are found, starting and ending with the
// new moons which bound the current lunation

VAR     sdate, adate, k1, k2, ks, nt1, nt2: DOUBLE;
        yy, mm, dd: Int64;

BEGIN
  CASE seconds_since_1970 OF
    UnixNow:    seconds_since_1970:= DateTimeToUnix(Now);
    UnixToday:  seconds_since_1970:= DateTimeToUnix(DateOf(Now));
    UnixAddDay: seconds_since_1970:= DateTimeToUnix(DateOf(Now)) + SecondsPerDay
  ELSE
  END;
  sdate:= jtime(seconds_since_1970);
  adate:= sdate - 45;
  jyear(adate, yy, mm, dd);
  k1:= floor((yy + ((mm - 1) * (1.0 / 12.0)) - 1900) * 12.3685);
  nt1:= meanphase(adate, k1);
  adate:= nt1;
  WHILE TRUE DO BEGIN
    adate += Synmonth;
    k2:= k1 + 1;
    nt2:= meanphase(adate, k2);
    IF (nt1 <= sdate) AND (nt2 > sdate) THEN
      BREAK;
    nt1:= nt2;
    k1:= k2
  END;
  RESULT[0]:= jdaytosecs(truephase(k1, 0.0));
  RESULT[1]:= jdaytosecs(truephase(k1, 0.25));
  RESULT[2]:= jdaytosecs(truephase(k1, 0.5));
  RESULT[3]:= jdaytosecs(truephase(k1, 0.75));
  RESULT[4]:= jdaytosecs(truephase(k2, 0.0))
END { PhaseHunt } ;


FUNCTION PhaseList(start: Int64= UnixToday; stop: Int64= UnixAddDay): TPhaseListResult;

(* Return the first significant phase of the Moon, plus its date and time       *)
(* together with the date and time of any others between the parameters. If the *)
(* first parameter is omitted it is assumed to be the start of the current day, *)
(* if the second parameter is omitted it is assumed to be 24 hours after the    *)
(* first parameter.                                                             *)

// phaselist - find time of phases of the moon between two dates
// times (in & out) are seconds_since_1970

VAR     sdate, edate, k, d: DOUBLE;
        yy, mm: Int64;
        phaseX4: INTEGER;

BEGIN
  CASE start OF
    UnixNow:    start:= DateTimeToUnix(Now);
    UnixToday:  start:= DateTimeToUnix(DateOf(Now));
    UnixAddDay: start:= DateTimeToUnix(DateOf(Now)) + SecondsPerDay
  ELSE
  END;
  CASE stop OF
    UnixNow:    stop:= DateTimeToUnix(Now);
    UnixToday:  stop:= DateTimeToUnix(DateOf(Now));
    UnixAddDay: stop:= start + SecondsPerDay
  ELSE
  END;
  RESULT.phase:= -1;
  SetLength(RESULT.times, 0);
  sdate:= jtime(start);
  edate:= jtime(stop);
  jyear(sdate, yy, mm, d);
  k:= floor((yy + ((mm - 1) * (1.0 / 12.0)) - 1900) * 12.3685) - 2;
  WHILE TRUE DO BEGIN
    k += 1.0;
    FOR phaseX4:= 0 TO 3 DO BEGIN
      d:= truephase(k, phaseX4 / 4.0);
      IF d >= edate THEN
        EXIT;
      IF d >= sdate THEN BEGIN
        IF RESULT.phase = -1 THEN
          RESULT.phase:= phaseX4;
        SetLength(RESULT.times, Length(RESULT.times) + 1);
        RESULT.times[Length(RESULT.times) - 1]:= jdaytosecs(d)
      END
    END
  END
END { PhaseList } ;


FUNCTION kepler(m, ecc: DOUBLE): DOUBLE;

// kepler - solve the equation of Kepler

CONST   EPSILON= 1e-6;

VAR     delta: DOUBLE;

BEGIN
  m:= torad(m);
  RESULT:= m;
  REPEAT
    delta:= RESULT - ecc * sin(RESULT) - m;
    RESULT -= delta / (1 - ecc * cos(RESULT))
  UNTIL abs(delta) <= EPSILON
END { kepler } ;


FUNCTION Phase(seconds_since_1970: Int64= UnixNow): TPhaseResult;

(* Return either an array containing a description of the Moon's current        *)
(* appearance, or a single number giving its phase angle. In either case if the *)
(* parameter is omitted it is assumed to be the current date and time.          *)

// phase - calculate phase of moon as a fraction:
//
// The argument is the time for which the phase is requested,
// expressed as a Julian date and fraction.  Returns the terminator
// phase angle as a percentage of a full circle (i.e., 0 to 1),
// and stores into pointer arguments the illuminated fraction of
// the Moon's disc, the Moon's age in days and fraction, the
// distance of the Moon from the centre of the Earth, and the
// angular diameter subtended by the Moon as seen by an observer
// at the centre of the Earth.

VAR     pdate: DOUBLE;
        Day, N, M, Ec, Lambdasun, ml, MM, MN, Ev, Ae, A3, MmP,
           mEc, A4, lP, V, lPP, NP, y, x, Lambdamoon, BetaM,
           MoonAge, MoonPhase,
           MoonDist, MoonDFrac, MoonAng, MoonPar,
           F, SunDist, SunAng: DOUBLE;

BEGIN
  CASE seconds_since_1970 OF
    UnixNow:    seconds_since_1970:= DateTimeToUnix(Now);
    UnixToday:  seconds_since_1970:= DateTimeToUnix(DateOf(Now));
    UnixAddDay: seconds_since_1970:= DateTimeToUnix(DateOf(Now)) + SecondsPerDay
  ELSE
  END;
  RESULT.seconds_since_1970:= seconds_since_1970;
  pdate:= jtime(seconds_since_1970);

// Calculation of the Sun's position.

  Day:= pdate - Epoch;                  // date within epoch
  N:= fixangle((360 / 365.2422) * Day); // mean anomaly of the Sun
  M:= fixangle(N + Elonge - Elongp);    // convert from perigee co-ordinates to epoch 1980.0
  Ec:= kepler(M, Eccent);               // solve equation of Kepler
  Ec:= sqrt((1 + Eccent) / (1 - Eccent)) * tan(Ec / 2);
  Ec:= 2 * todeg(arctan(Ec));             // true anomaly
  Lambdasun:= fixangle(Ec + Elongp);    // Sun's geocentric ecliptic longitude

// Orbital distance factor.

  F:= ((1 + Eccent * cos(torad(Ec))) / (1 - Eccent * Eccent));
  SunDist:= Sunsmax / F;                // distance to Sun in km
  SunAng:= F * Sunangsiz;               // Sun's angular size in degrees

// Calculation of the Moon's position.

// Moon's mean longitude.

  ml:= fixangle(13.1763966 * Day + Mmlong);

// Moon's mean anomaly.

  MM:= fixangle(ml - 0.1114041 * Day - Mmlongp);

// Moon's ascending node mean longitude.

  MN:= fixangle(Mlnode - 0.0529539 * Day);

// Evection.

  Ev:= 1.2739 * sin(torad(2 * (ml - Lambdasun) - MM));

// Annual equation.

  Ae:= 0.1858 * sin(torad(M));

// Correction term.

  A3:= 0.37 * sin(torad(M));

// Corrected anomaly.

  MmP:= MM + Ev - Ae - A3;

// Correction for the equation of the centre.

  mEc:= 6.2886 * sin(torad(MmP));

// Another correction term.

  A4:= 0.214 * sin(torad(2 * MmP));

// Corrected longitude.

  lP:= ml + Ev + mEc - Ae + A4;

// Variation.

  V:= 0.6583 * sin(torad(2 * (lP - Lambdasun)));

// True longitude.

  lPP:= lP + V;

// Corrected longitude of the node.

  NP:= MN - 0.16 * sin(torad(M));

// Y inclination coordinate.

  y:= sin(torad(lPP - NP)) * cos(torad(Minc));

// X inclination coordinate.

  x:= cos(torad(lPP - NP));

// Ecliptic longitude.

  Lambdamoon:= todeg(arctan2(y, x));
  Lambdamoon += NP;

// Ecliptic latitude.

  BetaM:= todeg(arcsin(sin(torad(lPP - NP)) * sin(torad(Minc))));

// Calculation of the phase of the Moon.

// Age of the Moon in degrees.

  MoonAge:= lPP - Lambdasun;

// Phase of the Moon.

  MoonPhase:= (1 - cos(torad(MoonAge))) / 2;

// Calculate distance of moon from the centre of the Earth.

  MoonDist:= (Msmax * (1 - Mecc * Mecc)) / (1 + Mecc * cos(torad(MmP + mEc)));

// Calculate Moon's angular diameter.

  MoonDFrac:= MoonDist / Msmax;
  MoonAng:= Mangsiz / MoonDFrac;

// Calculate Moon's parallax.

  MoonPar:= Mparallax / MoonDFrac;

  RESULT.MoonPhase:= fixangle(MoonAge) / 360.0;
  RESULT.MoonIllum:= MoonPhase;
  RESULT.MoonAge:= Synmonth * (fixangle(MoonAge) / 360.0);
  RESULT.MoonDist:= MoonDist;
  RESULT.MoonAng:= MoonAng;
  RESULT.SunDist:= SunDist;
  RESULT.SunAng:= SunAng
END { Phase } ;


FUNCTION Phase(seconds_since_1970: Int64= UnixNow): DOUBLE;

(* Return either an array containing a description of the Moon's current        *)
(* appearance, or a single number giving its phase angle. In either case if the *)
(* parameter is omitted it is assumed to be the current date and time.          *)

VAR     temp: TPhaseResult;

BEGIN
  temp:= Phase(seconds_since_1970);
  RESULT:= temp.MoonPhase
END { Phase } ;


CONST   test_datum = 900243421;                 (* Unix seconds                 *)
        test_jatum = 2451006.984039351809770;   (* Astronomers' Julian date     *)
        test_yy= 1998;                          (* Civil date (year, month,     *)
        test_mm= 7;                             (* day with fraction giving     *)
        test_dd= 12.484039351809770;            (* time in UTC).                *)


FUNCTION test_jTime: BOOLEAN;

(* Test conversion of a known date and time (used by the author of MoonPhase.pm *)
(* even though he hasn't documented it) from a unix-style seconds count to an   *)
(* astronomer's Julian date.                                                    *)

CONST   limit= 1.0e-10;

VAR     temp, diff: DOUBLE;

BEGIN
  temp:= jTime(test_datum);
  diff:= temp - test_jatum;
  RESULT:= Abs(diff) <= limit
END { test_jTime } ;


FUNCTION test_jDayToSecs: BOOLEAN;

(* Test conversion of a known date and time (used by the author of MoonPhase.pm *)
(* even though he hasn't documented it) from an astronomer's Julian date to a   *)
(* unix-style seconds count, expecting an error of no worse than one second.    *)

CONST   limit= 1.0 / 86400;

VAR     temp, diff: DOUBLE;

BEGIN
  temp:= jDayToSecs(test_jatum);
  diff:= temp - test_datum;
  RESULT:= Abs(diff) <= limit
END { test_jDayToSecs } ;


FUNCTION test_jYear: BOOLEAN;

(* Test conversion of a known date and time (used by the author of MoonPhase.pm *)
(* even though he hasn't documented it) from an astronomer's Julian date to a   *)
(* civil date comprising year, month and day-of-month. Note that MoonPhase.pm   *)
(* documents the day as being an integer, but this is wrong since it loses      *)
(* time-of-day information which callers assume to be present.                  *)

CONST   limit= 1.0 / 86400;

VAR     yy, mm: Int64;
        dd, diff: DOUBLE;

BEGIN
  jYear(test_jatum, yy, mm, dd);
  diff:= dd - test_dd;
  RESULT:= (yy = test_yy) AND (mm = test_mm) AND (Abs(diff) <= limit)
END { test_jYear } ;


FUNCTION test_Phase: BOOLEAN;

(* Get lunar age etc. corresponding to a known date and time (used by the       *)
(* author of MoonPhase.pm even though he hasn't documented it). Tolerate loss   *)
(* of a couple of significant digits compared with the original Perl.           *)

VAR     phaseResult: TPhaseResult;
        diff: DOUBLE;

BEGIN
  phaseResult:= Phase(test_datum);
  diff:= phaseResult.MoonPhase - 0.598939375319023;
  RESULT:= Abs(diff) <= 1.0e-13;
  diff:= phaseResult.MoonIllum - 0.906458030827876;
  RESULT:= RESULT AND (Abs(diff) <= 1.0e-13);
  diff:= phaseResult.MoonAge - 17.6870323368022;
  RESULT:= RESULT AND (Abs(diff) <= 1.0e-11);
  diff:= phaseResult.MoonDist - 372479.357420033;
  RESULT:= RESULT AND (Abs(diff) <= 1.0e-7);
  diff:= phaseResult.MoonAng - 0.534682403555093;
  RESULT:= RESULT AND (Abs(diff) <= 1.0e-13);
  diff:= phaseResult.SunDist - 152078368.820205;
  RESULT:= RESULT AND (Abs(diff) <= 1.0e-4);
  diff:= phaseResult.SunAng - 0.524434538105092;
  RESULT:= RESULT AND (Abs(diff) <= 1.0e-13)
END { test_Phase } ;


FUNCTION test_PhaseHunt: BOOLEAN;

(* Get the five phases surrounding a known date and time (used by the author of *)
(* MoonPhase.pm even though he hasn't documented it and hasn't specified what   *)
(* timezone his time output assumes).                                           *)

CONST   n0= 898660308;
        q1= 899318539;
        q2= 900000167;
        q3= 900602119;
        n1= 901201502;

VAR     phaseHuntResult: TPhaseHuntResult;
        diff: Int64;

BEGIN
  phaseHuntResult:= PhaseHunt(test_datum);
  diff:= phaseHuntResult[0] - n0;
  RESULT:= Abs(diff) <= 1;
  diff:= phaseHuntResult[1] - q1;
  RESULT:= RESULT AND (Abs(diff) <= 1);
  diff:= phaseHuntResult[2] - q2;
  RESULT:= RESULT AND (Abs(diff) <= 1);
  diff:= phaseHuntResult[3] - q3;
  RESULT:= RESULT AND (Abs(diff) <= 1);
  diff:= phaseHuntResult[4] - n1;
  RESULT:= RESULT AND (Abs(diff) <= 1)
END { test_PhaseHunt } ;


FUNCTION test_PhaseList(dst: BOOLEAN): BOOLEAN;

(* Referring to http://scienceworld.wolfram.com/astronomy/Lunation.html there   *)
(* should be a New Moon at 2002-06-10 23:47 UTC, with a precision of around     *)
(* three minutes. This is a particularly useful test case since if a UK         *)
(* daylight saving time correction is applied it should shift to 2002-06-11     *)
(* 00:47 BST, i.e. there is a consequent date change.                           *)

VAR     dstOffset, start, stop, diff: Int64;
        temp: TPhaseListResult;

BEGIN
  IF dst THEN
    dstOffset:= 3600
  ELSE
    dstOffset:= 0;
  start:= CalcUS(2002, 06, 10) - dstOffset;
  stop:= CalcUS(2002, 06, 11) - dstOffset;
  temp:= PhaseList(start, stop);
  RESULT:= FALSE;
  IF dst THEN                           (* New Moon shifted to tomorrow         *)
    IF (temp.phase = -1) AND (Length(temp.times) = 0) THEN
      RESULT:= TRUE
    ELSE BEGIN END
  ELSE BEGIN
    IF temp.phase <> 0 THEN
      EXIT;
    diff:= temp.times[0] - 1023752886;
    IF Abs(diff) <= 1 THEN
      RESULT:= TRUE
  END
END { test_PhaseList } ;


INITIALIZATION

(* Test that standard Pascal functions conform to the original explicit         *)
(* functions.                                                                   *)

  Assert(Sign(-2) = -1);
  Assert(Sign(0) = 0);
  Assert(Sign(2) = 1);
  Assert(Floor(-1.5) = -2);
  Assert(Floor(0.0) = 0);
  Assert(Floor(1.5) = 1);

(* More advanced test cases, typically against unix's date command etc.         *)

  Assert(CalcUS(1970, 1, 1) = 0);
  Assert(test_jTime);
  Assert(test_jDayToSecs);
  Assert(test_jYear);

(* Tests similar to those at the end of the Perl module (i.e. in the POD        *)
(* content) using a date and time deduced from the second Perl example (i.e. in *)
(* late June or July 1998), since the author has not specified this.            *)

  Assert(test_Phase);
  Assert(test_PhaseHunt);

(* Test a specific data point from an external source, without and with a 1     *)
(* hour DST correction.                                                         *)

  Assert(test_PhaseList(FALSE));
  Assert(test_PhaseList(TRUE))
end.
