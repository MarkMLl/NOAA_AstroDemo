/*

    A Moon for the Sun  (and now the Penguin)

    Release 2.0

    Designed and implemented by John Walker in December 1987,
    revised and updated in February of 1988.

    Make with:

    cc -O moontool.c -o moontool -lm -lsuntool -lsunwindow -lpixrect

    Adding  appropriate  floating  point  options  to your hardware.  This
    program is a SunView tool which displays, as the  icon  for  a  closed
    window,  the  current phase of the Moon.  A subtitle in the icon gives
    the age of the Moon in days  and  hours.   If  called  with  the  "-t"
    switch,  it  rapidly  increments  forward  through time to display the
    cycle of phases.

    If you open the window, additional information is displayed  regarding
    the  Moon.   The  information  is  generally  accurate  to  within ten
    minutes.

    The algorithms used in this program to calculate the positions Sun and
    Moon as seen from the Earth are given in the book "Practical Astronomy
    With  Your  Calculator"  by  Peter  Duffett-Smith,   Second   Edition,
    Cambridge University Press, 1981.  Ignore the word "Calculator" in the
    title;  this  is  an  essential  reference  if  you're  interested  in
    developing  software  which  calculates  planetary  positions, orbits,
    eclipses, and  the  like.   If  you're  interested  in  pursuing  such
    programming, you should also obtain:

    "Astronomical  Formulae for Calculators" by Jean Meeus, Third Edition,
    Willmann-Bell, 1985.  A must-have.

    "Planetary  Programs  and  Tables  from  -4000  to  +2800"  by  Pierre
    Bretagnon  and Jean-Louis Simon, Willmann-Bell, 1986.  If you want the
    utmost  (outside  of  JPL)  accuracy  for  the  planets,  it's   here.

    "Celestial BASIC" by Eric Burgess, Revised Edition, Sybex, 1985.  Very
    cookbook oriented, and many of the algorithms are hard to dig  out  of
    the turgid BASIC code, but you'll probably want it anyway.

    Many of these references can be obtained from Willmann-Bell, P.O.  Box
    35025,  Richmond,  VA 23235, USA.  Phone: (804) 320-7016.  In addition
    to their own publications, they stock most of the standard  references
    for mathematical and positional astronomy.

    This program was written by:

       John Walker
       Autodesk, Inc.
       2320 Marinship Way
       Sausalito, CA  94965
       (415) 332-2344 Ext. 829

       Usenet: {sun!well}!acad!kelvin

    This  program is in the public domain: "Do what thou wilt shall be the
    whole of the law".  I'd appreciate  receiving  any  bug  fixes  and/or
    enhancements,  which  I'll  incorporate  in  future  versions  of  the
    program.  Please leave the original attribution information intact  so
    that credit and blame may be properly apportioned.

*/


#include "moon.h"
#include <math.h>


/*  Astronomical constants  */

#define epoch       2444238.5      /* 1980 January 0.0 */

/*  Constants defining the Sun's apparent orbit  */

#define elonge      278.833540     /* Ecliptic longitude of the Sun
                                      at epoch 1980.0 */
#define elongp      282.596403     /* Ecliptic longitude of the Sun at
                                      perigee */
#define eccent      0.016718       /* Eccentricity of Earth's orbit */
#define sunsmax     1.495985e8     /* Semi-major axis of Earth's orbit, km */
#define sunangsiz   0.533128       /* Sun's angular size, degrees, at
                                      semi-major axis distance */

/*  Elements of the Moon's orbit, epoch 1980.0  */

#define mmlong      64.975464      /* Moon's mean lonigitude at the epoch */
#define mmlongp     349.383063     /* Mean longitude of the perigee at the
                                      epoch */
#define mlnode      151.950429     /* Mean longitude of the node at the
                                      epoch */
#define minc        5.145396       /* Inclination of the Moon's orbit */
#define mecc        0.054900       /* Eccentricity of the Moon's orbit */
#define mangsiz     0.5181         /* Moon's angular size at distance a
                                      from Earth */
#define msmax       384401.0       /* Semi-major axis of Moon's orbit in km */
#define mparallax   0.9507         /* Parallax at distance a from Earth */
#define synmonth    29.53058868    /* Synodic month (new Moon to new Moon) */
#define lunatbase   2423436.0      /* Base date for E. W. Brown's numbered
                                      series of lunations (1923 January 16) */

/*  Properties of the Earth  */

#define earthrad    6378.16        /* Radius of Earth in kilometres */

#define PI 3.14159265358979323846  /* Assume not near black hole nor in
                                      Tennessee */

/* Prototypes */
void drawmoon (double ph);
void update_p1 (void);
long jdate (struct tm *t);
double jtime (struct tm *t);
void jyear (double td, long *yy, int *mm, int *dd);
void jhms (double j, int *h, int *m, int *s);
double meanphase (double sdate, double k);
double truephase (double k, double phase);
void phasehunt (double sdate, double phases[5]);
double kepler (double m, double ecc);
double phase (double pdate, double *pphase, double *mage, double *dist,
              double *angdia, double *sudist, double *suangdia);
void format_time (char *s, size_t sz, const struct tm *t, int short_form);
void jd2tm (double jd, struct tm *t);
const char *units (int un_type);


/*  Handy mathematical functions  */

#define sgn(x) (((x) < 0) ? -1 : ((x) > 0 ? 1 : 0))       /* Extract sign */
#define abs(x) ((x) < 0 ? (-(x)) : (x))                   /* Absolute val */
#define fixangle(a) ((a) - 360.0 * (floor((a) / 360.0)))  /* Fix angle    */
#define torad(d) ((d) * (PI / 180.0))                     /* Deg->Rad     */
#define todeg(d) ((d) * (180.0 / PI))                     /* Rad->Deg     */
#define dsin(x) (sin(torad((x))))                         /* Sin from deg */
#define dcos(x) (cos(torad((x))))                         /* Cos from deg */
#define plur(x) (x == 1 ? 0 : 1)

static int testmode = 0;      /* Rapid warp through time for debugging */
int running;                  /* True if the clock should be advance. */
static unsigned char moonmask[(64 * 64) / 8];
static unsigned char moonilast[(64 * 64) / 8];
static const unsigned char *moon_orig;
extern int localetimes_on;


/* First page of information. */
time_t p1_t, p1_t_disp, p1_t_diff;     /* Current time, time being displayed,
                                             and difference between them. */
struct tm p1_gm, p1_loc;               /* Time in GMT and local time zone. */
double p1_jd,                          /* Julian date. */
       p1_aom,                         /* Age of moon, in days. */
       p1_cphase,                      /* Percentage of moon visible. */
       p1_cdist,
       p1_cangdia,
       p1_csund,
       p1_csuang;
int p1_aom_d, p1_aom_h, p1_aom_m;      /* Age of moon (days, hours, mins) */
double p1_phasejd[5];                  /* Julian dates of phases. */
struct tm p1_phasetm[5];               /* Broken down times of phases. */
int p1_lunation;                       /* Lunation number of last new moon. */


/*  Main program  */
int
main (int argc, char *argv[])
{
#ifdef HAVE_SETLOCALE
   setlocale (LC_ALL, "");
#endif
   bindtextdomain (PACKAGE, LOCALEDIR);
   textdomain (PACKAGE);

   make_gtk_win (&argc, &argv, testmode);
   load_config();
   decode_args (argc, argv);

   /* Set up the moon mask. */
   moon_orig = get_orig_mask();
   memcpy (moonmask, moon_orig, sizeof (moonmask));

   /* Put values in initially before the clock tick is started. */
   p1_t_diff = 0;
   ringgg();

   start_timer (testmode ? 125 : 1000);

   do_gtk_stuff();
   return 0;
}


/*  DRAWMOON  --  Construct icon for moon, given phase of moon.  */

void
drawmoon (double ph)
{
   int i, j, lx, rx;
   int lb[8];            /* 8 bits in each one. */
   double cp, xscale;

   memset (moonmask, 0, sizeof (moonmask));

   xscale = cos(2 * PI * ph);
   for (i = 0; i < 24; i++)
   {
      for (j = 0; j < 8; j++)
         lb[j] = 0;
      cp = 24.0 * cos(asin(i / 24.0));
      if (ph < 0.5)
      {
         rx = 32 + cp;
         lx = 32 + xscale * cp;
      }
      else
      {
         lx = 33 - cp;
         rx = 33 - xscale * cp;
      }
      for (j = lx; j <= rx; j++)
         lb[j >> 3] |= (0x80 >> (7 - (j & 0x7)));
      for (j = 0; j < 8; j++)
         moonmask[(32 + i) * 8 + j] = moonmask[(32 - i) * 8 + j] = lb[j];
   }
}


/*  RINGGG  --  Refresh info displayed in window. */

int
ringgg (void)
{
   int wclosed;
   static double nptime = 0.0; /* Next new moon time */
   static int updyet = 0;      /* Update interval when window closed */
   static int firstime = 1;    /* Calculate text page first time */
   char amsg[32], tbuf[80];
   double p;
   int i;

#define CUPDINT  120               /* Update the icon every CUPDINT seconds
                                      when the window is iconic. */

   /* If the window is closed, only update the icon every
      two minutes. */
   wclosed = is_window_iconified();
   if (wclosed && (--updyet > 0) && !testmode)
      return 1;

   updyet = CUPDINT;
   update_p1();
   time(&p1_t);
   p1_t_disp = p1_t + p1_t_diff;
   p1_gm = *gmtime(&p1_t_disp);
   p1_jd = jtime(&p1_gm);

   p = phase(p1_jd, &p1_cphase, &p1_aom, &p1_cdist, &p1_cangdia,
             &p1_csund, &p1_csuang);
   drawmoon(p);
   p1_aom_d = (int) p1_aom;
   p1_aom_h = (int) (24 * (p1_aom - floor(p1_aom)));
   p1_aom_m = (int) (1440 * (p1_aom - floor(p1_aom))) % 60;
   if (p1_aom_d == 0)
      sprintf(amsg, _("%dh %dm"), p1_aom_h, p1_aom_m);
   else
      sprintf(amsg, _("%dd %dh"), p1_aom_d, p1_aom_h);
   update_icon_text (amsg);

   /* Only update icon if it changed (this eliminates gratuitous
      flashing of the icon on-screen). */

   if (memcmp(moonilast, moonmask, sizeof (moonmask)) != 0) {
      memcpy(moonilast, moonmask, sizeof (moonmask));
      update_moon (moonmask);
   }

   /* If we're iconic, there's nothing more to do. */

   if (wclosed && !firstime)
      return 1;

   /* Update textual information for open window. */

   firstime = 0;
   sprintf(tbuf, "%.5f", p1_jd);
   prt(0, tbuf);
   if (testmode)
      jd2tm (p1_jd, &p1_gm);
   format_time (tbuf, sizeof (tbuf), &p1_gm, 0);
   prt(1, tbuf);
   p1_loc = *localtime(&p1_t_disp);
   format_time (tbuf, sizeof (tbuf), &p1_loc, 0);
   prt(2, tbuf);
   update_discord (p1_t_disp);
   sprintf(tbuf, _("%.2f%%   (0%% = New, 100%% = Full)"),
           p1_cphase * 100.0);
   prt(5, tbuf);

   /* Information about the Moon */

   sprintf(tbuf, "%d %s, %d %s, %d %s.",
           p1_aom_d, units (U_DAY + plur (p1_aom_d)),
           p1_aom_h, units (U_HOUR + plur (p1_aom_h)),
           p1_aom_m, units (U_MINUTE + plur (p1_aom_m)));
   prt(4, tbuf);
   sprintf(tbuf, "%ld %s, %.1f %s.",
           (long) p1_cdist, units (U_KM + plur (p1_cdist)),
           p1_cdist / earthrad,
           units (U_EARTH_RAD + plur (p1_cdist / earthrad)));
   prt(6, tbuf);
   sprintf(tbuf, "%.4f %s.", p1_cangdia, units (U_DEGREE + plur (p1_cangdia)));
   prt(7, tbuf);

   /* Edit information about the Sun */

   sprintf(tbuf, "%.0f %s, %.3f %s.",
           p1_csund, units (U_KM + plur (p1_csund)),
           p1_csund / sunsmax,
           units (U_ASTR_UNIT + plur (p1_csund / sunsmax)));
   prt(9, tbuf);
   sprintf(tbuf, "%.4f %s.", p1_csuang,
           units (U_DEGREE + plur (p1_csuang)));
   prt(10, tbuf);

   /* Calculate times of phases of this lunation.  This is sufficiently
      time-consuming that we only do it once a month. */

   if (p1_jd > nptime)
   {
      phasehunt(p1_jd + 0.5, p1_phasejd);
      p1_lunation = floor(((p1_phasejd[0] + 7) - lunatbase) / synmonth) + 1;
      sprintf(tbuf, _("Lunation %d"), p1_lunation);
      prt(-1, tbuf);
      sprintf(tbuf, _("Lunation %d"), p1_lunation + 1);
      prt(-2, tbuf);
      for (i = 0; i < 5; i++)
      {
         jd2tm (p1_phasejd[i], &p1_phasetm[i]);
         format_time (tbuf, sizeof (tbuf), &p1_phasetm[i], 1);
         prt(12+i, tbuf);
      }
   }

   return 1;
}


/* Update the first page of information in the p1_* variables. */
void
update_p1 (void)
{
}


/*  JDATE  --  Convert internal GMT date and time to Julian day
               and fraction.  */

long
jdate (struct tm *t)
{
   long c, m, y;

   y = t->tm_year + 1900;
   m = t->tm_mon + 1;
   if (m > 2)
      m = m - 3;
   else
   {
      m = m + 9;
      y--;
   }
   c = y / 100L;              /* Compute century */
   y -= 100L * c;
   return t->tm_mday + (c * 146097L) / 4 + (y * 1461L) / 4 +
          (m * 153L + 2) / 5 + 1721119L;
}


/* JTIME --    Convert internal GMT date and time to astronomical Julian
               time (i.e. Julian date plus day fraction, expressed as
               a double).  */
/* Algorithm as given in Meeus, Astronomical Algorithms, Chapter 7, page 61 */
double
jtime (struct tm *t)
{
   long year;
   int mon, mday, hour, min, sec;
   int a, b, m;
   long y;

   year = t->tm_year + 1900;
   mon = t->tm_mon;
   mday = t->tm_mday;
   hour = t->tm_hour;
   min = t->tm_min;
   sec = t->tm_sec;

#ifdef PARANOID
    assert(mon  >= 0 && mon  < 12);
    assert(mday >  0 && mday < 32);
    assert(hour >= 0 && hour < 24);
    assert(min  >= 0 && min  < 60);
    assert(sec  >= 0 && sec  < 60);
#endif

    m = mon + 1;
    y = year;

   if (m <= 2) {
      y--;
      m += 12;
   }

   /* Determine whether date is in Julian or Gregorian calendar based on
      canonical date of calendar reform. */

   if ((year < 1582) || ((year == 1582) && ((mon < 9) ||
       (mon == 9 && mday < 5)))) {
      b = 0;
   } else {
      a = ((int) (y / 100));
      b = 2 - a + (a / 4);
   }

   return (((long) (365.25 * (y + 4716))) + ((int) (30.6001 * (m + 1))) +
           mday + b - 1524.5) +
           ((sec + 60L * (min + 60L * hour)) / 86400.0);
}

/*  JYEAR  --  Convert Julian date to year, month, day, which are
               returned via integer pointers to integers.  */

void
jyear (double td, long *yy, int *mm, int *dd)
{
   double z, f, a, alpha, b, c, d, e;

   td += 0.5;
   z = floor(td);
   f = td - z;

   if (z < 2299161.0) {
      a = z;
   } else {
      alpha = floor((z - 1867216.25) / 36524.25);
      a = z + 1 + alpha - floor(alpha / 4);
   }

   b = a + 1524;
   c = floor((b - 122.1) / 365.25);
   d = floor(365.25 * c);
   e = floor((b - d) / 30.6001);

   *dd = (int) (b - d - floor(30.6001 * e) + f);
   *mm = (int) ((e < 14) ? (e - 1) : (e - 13));
   *yy = (long) ((*mm > 2) ? (c - 4716) : (c - 4715));
}

/*  JHMS  --  Convert Julian time to hour, minutes, and seconds.  */

void
jhms (double j, int *h, int *m, int *s)
{
   long ij;

   j += 0.5;                  /* Astronomical to civil */
   ij = (long) (((j - floor(j)) * 86400.0) + 0.5);  // Round to nearest second
   *h = (int) (ij / 3600L);
   *m = (int) ((ij / 60L) % 60L);
   *s = (int) (ij % 60L);
}


/*  MEANPHASE  --  Calculates  time  of  the mean new Moon for a given
                   base date.  This argument K to this function is the
                   precomputed synodic month index, given by:

                          K = (year - 1900) * 12.3685

                   where year is expressed as a year and fractional year.  */

double
meanphase (double sdate, double k)
{
    double t, t2, t3, nt1;

    /* Time in Julian centuries from 1900 January 0.5 */
    t = (sdate - 2415020.0) / 36525;
    t2 = t * t;                       /* Square for frequent use */
    t3 = t2 * t;                      /* Cube for frequent use */

    nt1 = 2415020.75933 + synmonth * k
            + 0.0001178 * t2
            - 0.000000155 * t3
            + 0.00033 * dsin(166.56 + 132.87 * t - 0.009173 * t2);

    return nt1;
}


/*  TRUEPHASE  --  Given a K value used to determine the
                   mean phase of the new moon, and a phase
                   selector (0.0, 0.25, 0.5, 0.75), obtain
                   the true, corrected phase time.  */

double
truephase (double k, double phase)
{
        double t, t2, t3, pt, m, mprime, f;
        int apcor = 0;

        k += phase;                /* Add phase to new moon time */
        t = k / 1236.85;           /* Time in Julian centuries from
                                      1900 January 0.5 */
        t2 = t * t;                /* Square for frequent use */
        t3 = t2 * t;               /* Cube for frequent use */
        pt = 2415020.75933         /* Mean time of phase */
             + synmonth * k
             + 0.0001178 * t2
             - 0.000000155 * t3
             + 0.00033 * dsin(166.56 + 132.87 * t - 0.009173 * t2);

        m = 359.2242               /* Sun's mean anomaly */
            + 29.10535608 * k
            - 0.0000333 * t2
            - 0.00000347 * t3;
        mprime = 306.0253          /* Moon's mean anomaly */
            + 385.81691806 * k
            + 0.0107306 * t2
            + 0.00001236 * t3;
        f = 21.2964                /* Moon's argument of latitude */
            + 390.67050646 * k
            - 0.0016528 * t2
            - 0.00000239 * t3;
        if ((phase < 0.01) || (abs(phase - 0.5) < 0.01)) {

           /* Corrections for New and Full Moon */

           pt +=     (0.1734 - 0.000393 * t) * dsin(m)
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
           apcor = 1;
        } else if ((abs(phase - 0.25) < 0.01 || (abs(phase - 0.75) < 0.01))) {
           pt +=     (0.1721 - 0.0004 * t) * dsin(m)
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
           if (phase < 0.5)
              /* First quarter correction */
              pt += 0.0028 - 0.0004 * dcos(m) + 0.0003 * dcos(mprime);
           else
              /* Last quarter correction */
              pt += -0.0028 + 0.0004 * dcos(m) - 0.0003 * dcos(mprime);
           apcor = 1;
        }
        if (!apcor) {
           fprintf(stderr, "TRUEPHASE called with invalid phase selector.\n");
           abort();
        }
        return pt;
}

/*  PHASEHUNT  --  Find time of phases of the moon which surround
                   the current date.  Five phases are found, starting
                   and ending with the new moons which bound the
                   current lunation.  */

void
phasehunt (double sdate, double phases[5])
{
    double adate, k1, k2, nt1, nt2;
    long yy;
    int mm, dd;

    adate = sdate - 45;

    jyear(adate, &yy, &mm, &dd);
    k1 = floor((yy + ((mm - 1) * (1.0 / 12.0)) - 1900) * 12.3685);

    adate = nt1 = meanphase(adate, k1);
    while (1) {
        adate += synmonth;
        k2 = k1 + 1;
        nt2 = meanphase(adate, k2);
        if (nt1 <= sdate && nt2 > sdate)
            break;
        nt1 = nt2;
        k1 = k2;
    }
    phases[0] = truephase(k1, 0.0);
    phases[1] = truephase(k1, 0.25);
    phases[2] = truephase(k1, 0.5);
    phases[3] = truephase(k1, 0.75);
    phases[4] = truephase(k2, 0.0);
}


/*  KEPLER  --  Solve the equation of Kepler.  */
double
kepler (double m, double ecc)
{
        double e, delta;
#define EPSILON 1E-6

        e = m = torad(m);
        do {
           delta = e - ecc * sin(e) - m;
           e -= delta / (1 - ecc * cos(e));
        } while (abs(delta) > EPSILON);
        return e;
}


/*  PHASE  --  Calculate phase of moon as a fraction:

        The argument is the time for which the phase is requested,
        expressed as a Julian date and fraction.  Returns the terminator
        phase angle as a percentage of a full circle (i.e., 0 to 1),
        and stores into pointer arguments the illuminated fraction of
        the Moon's disc, the Moon's age in days and fraction, the
        distance of the Moon from the centre of the Earth, and the
        angular diameter subtended by the Moon as seen by an observer
        at the centre of the Earth.

*/

double phase(
   double pdate,                      /* Date for which to calculate phase */
   double *pphase,                    /* Illuminated fraction */
   double *mage,                      /* Age of moon in days */
   double *dist,                      /* Distance in kilometres */
   double *angdia,                    /* Angular diameter in degrees */
   double *sudist,                    /* Distance to Sun */
   double *suangdia)                  /* Sun's angular diameter */
{

        double Day, N, M, Ec, Lambdasun, ml, MM, MN, Ev, Ae, A3, MmP,
               mEc, A4, lP, V, lPP, NP, y, x, Lambdamoon, BetaM,
               MoonAge, MoonPhase,
               MoonDist, MoonDFrac, MoonAng, MoonPar,
               F, SunDist, SunAng;

        /* Calculation of the Sun's position */

        Day = pdate - epoch;                    /* Date within epoch */
        N = fixangle((360 / 365.2422) * Day);   /* Mean anomaly of the Sun */
        M = fixangle(N + elonge - elongp);      /* Convert from perigee
                                                co-ordinates to epoch 1980.0 */
        Ec = kepler(M, eccent);                 /* Solve equation of Kepler */
        Ec = sqrt((1 + eccent) / (1 - eccent)) * tan(Ec / 2);
        Ec = 2 * todeg(atan(Ec));               /* True anomaly */
        Lambdasun = fixangle(Ec + elongp);      /* Sun's geocentric ecliptic
                                                   longitude */
        /* Orbital distance factor */
        F = ((1 + eccent * cos(torad(Ec))) / (1 - eccent * eccent));
        SunDist = sunsmax / F;            /* Distance to Sun in km */
        SunAng = F * sunangsiz;           /* Sun's angular size in degrees */


        /* Calculation of the Moon's position */

        /* Moon's mean longitude */
        ml = fixangle(13.1763966 * Day + mmlong);

        /* Moon's mean anomaly */
        MM = fixangle(ml - 0.1114041 * Day - mmlongp);

        /* Moon's ascending node mean longitude */
        MN = fixangle(mlnode - 0.0529539 * Day);

        /* Evection */
        Ev = 1.2739 * sin(torad(2 * (ml - Lambdasun) - MM));

        /* Annual equation */
        Ae = 0.1858 * sin(torad(M));

        /* Correction term */
        A3 = 0.37 * sin(torad(M));

        /* Corrected anomaly */
        MmP = MM + Ev - Ae - A3;

        /* Correction for the equation of the centre */
        mEc = 6.2886 * sin(torad(MmP));

        /* Another correction term */
        A4 = 0.214 * sin(torad(2 * MmP));

        /* Corrected longitude */
        lP = ml + Ev + mEc - Ae + A4;

        /* Variation */
        V = 0.6583 * sin(torad(2 * (lP - Lambdasun)));

        /* True longitude */
        lPP = lP + V;

        /* Corrected longitude of the node */
        NP = MN - 0.16 * sin(torad(M));

        /* Y inclination coordinate */
        y = sin(torad(lPP - NP)) * cos(torad(minc));

        /* X inclination coordinate */
        x = cos(torad(lPP - NP));

        /* Ecliptic longitude */
        Lambdamoon = todeg(atan2(y, x));
        Lambdamoon += NP;

        /* Ecliptic latitude */
        BetaM = todeg(asin(sin(torad(lPP - NP)) * sin(torad(minc))));

        /* Calculation of the phase of the Moon */

        /* Age of the Moon in degrees */
        MoonAge = lPP - Lambdasun;

        /* Phase of the Moon */
        MoonPhase = (1 - cos(torad(MoonAge))) / 2;

        /* Calculate distance of moon from the centre of the Earth */
        MoonDist = (msmax * (1 - mecc * mecc)) /
                   (1 + mecc * cos(torad(MmP + mEc)));

        /* Calculate Moon's angular diameter */
        MoonDFrac = MoonDist / msmax;
        MoonAng = mangsiz / MoonDFrac;

        /* Calculate Moon's parallax */
        MoonPar = mparallax / MoonDFrac;

        *pphase = MoonPhase;
        *mage = synmonth * (fixangle(MoonAge) / 360.0);
        *dist = MoonDist;
        *angdia = MoonAng;
        *sudist = SunDist;
        *suangdia = SunAng;
        return (fixangle(MoonAge) / 360.0);
}


/* Format a time either in the traditional moontool format or in the way the
      current locale favours, depending on the menu option. */
void
format_time (char *s, size_t sz, const struct tm *t, int short_form)
{
   char buf[40];

   if (localetimes_on)
   {
      const char *frmt = "%c";
      strftime(s, sz, frmt, t);
   }
   else
   {
      strftime (buf, sizeof (buf), "%B", t);
      if (short_form)
         sprintf(s, "%02d:%02d UTC %d %s %d",
                 t->tm_hour, t->tm_min,
                 t->tm_mday, buf, t->tm_year + 1900);
      else
         sprintf(s, "%02d:%02d:%02d %d %s %d",
                 t->tm_hour, t->tm_min, t->tm_sec,
                 t->tm_mday, buf, t->tm_year + 1900);
   }
}


/* Fill a `tm' struct with time info for a given Julian date. */
void
jd2tm (double jd, struct tm *t)
{
   long y;
   int mon, d, h, min, s;

   jyear(jd, &y, &mon, &d);
   jhms(jd, &h, &min, &s);

   t->tm_sec = s;
   t->tm_min = min;
   t->tm_hour = h;
   t->tm_mday = d;
   t->tm_mon = mon - 1;
   t->tm_year = y - 1900;
   t->tm_wday = t->tm_yday = t->tm_isdst = 0;
   mktime (t);
}


/* Return a locale-specific string representation of the name of certain units
      of measurement. */
const char *
units (int un_type)
{
   const char *s;

   switch (un_type)
   {
      case U_DAY:  s = _("day");  break;
      case U_DAYS:  s = _("days");  break;
      case U_HOUR:  s = _("hour");  break;
      case U_HOURS:  s = _("hours");  break;
      case U_MINUTE:  s = _("minute");  break;
      case U_MINUTES:  s = _("minutes");  break;
      case U_DEGREE:  s = _("degree");  break;
      case U_DEGREES:  s = _("degrees");  break;
      case U_KM:  s = _("kilometre");  break;
      case U_KMS:  s = _("kilometres");  break;
      case U_ASTR_UNIT:  s = _("astronomical unit");  break;
      case U_ASTR_UNITS:  s = _("astronomical units");  break;
      case U_EARTH_RAD:  s = _("Earth radius");  break;
      case U_EARTH_RADS:  s = _("Earth radii");  break;
      default:
         assert (0);
   }

   return s;
}
