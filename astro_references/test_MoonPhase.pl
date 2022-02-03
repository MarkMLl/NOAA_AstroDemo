#!/usr/bin/perl -w

use MoonPhase;
use Time::Local;

# The worked examples in MoonPhase.pm don't state what date and time they're
# for, searching indicates that it's Sun Jul 12 11:37:01 1998 UTC. The
# "correct" values in the tests below are from an initial Perl run of the
# unmodified module on a Linux x86 system.

$datum = 900243421;
$jatum = 2451006.984039351809770;

printf "\nTest jtime(), jdaytosecs() and jyear(). Correct values are in parentheses():\n\n";

printf "jtime(%d) -> %17.15f (%17.15f, %+8.3e)\n",
		$datum, Astro::MoonPhase::jtime($datum), $jatum, Astro::MoonPhase::jtime($datum) - $jatum; 

printf "jdaytosecs(%17.15f) -> %17.15f (%17.15f, %+8.3e)\n",
		$jatum, Astro::MoonPhase::jdaytosecs($jatum), $datum, Astro::MoonPhase::jdaytosecs($jatum) - $datum;

{
  my ($yy, $mm, $dd);
  Astro::MoonPhase::jyear($jatum, \$yy, \$mm, \$dd);
  printf "jyear(%17.15f) -> %d %d %17.15f (1998 7 12.484039351809770)\n", $jatum, $yy, $mm, $dd;
}

printf "\nTest Astro::MoonPhase:phase(%d) for %s, significant digits should be the same on x86:\n\n", time, scalar(localtime(900243421));

   ( $MoonPhase,
     $MoonIllum,
     $MoonAge,
     $MoonDist,
     $MoonAng,
     $SunDist,
     $SunAng ) = Astro::MoonPhase::phase(900243421);

     $xMoonPhase = 0.598939375319023;
     $xMoonIllum = 0.906458030827876;
     $xMoonAge = 17.6870323368022;
     $xMoonDist = 372479.357420033;
     $xMoonAng = 0.534682403555093;
     $xSunDist = 152078368.820205;
     $xSunAng = 0.524434538105092;

     printf "MoonPhase  = %17.15f (%17.15f, %+8.3e)\n", $MoonPhase, $xMoonPhase, $MoonPhase - $xMoonPhase;
     printf "MoonIllum  = %17.15f (%17.15f, %+8.3e)\n", $MoonIllum, $xMoonIllum, $MoonIllum - $xMoonIllum;
     printf "MoonAge    = %17.15f (%17.15f, %+8.3e)\n", $MoonAge, $xMoonAge, $MoonAge - $xMoonAge;
     printf "MoonDist   = %17.15f (%17.15f, %+8.3e)\n", $MoonDist, $xMoonDist, $MoonDist - $xMoonDist;
     printf "MoonAng    = %17.15f (%17.15f, %+8.3e)\n", $MoonAng, $xMoonAng, $MoonAng - $xMoonAng;
     printf "SunDist    = %17.15f (%17.15f, %+8.3e)\n", $SunDist, $xSunDist, $SunDist - $xSunDist;
     printf "SunAng     = %17.15f (%17.15f, %+8.3e)\n", $SunAng, $xSunAng, $SunAng - $xSunAng;

# Referring to http://scienceworld.wolfram.com/astronomy/Lunation.html there
# should be a New Moon at 2002-06-10 23:47 UTC, with a precision of around
# three minutes. This is a particularly useful test case since if a UK
# daylight saving time correction is applied it should shift to 2002-06-11
# 00:47 BST, i.e. there is a consequent date change.

# Note that the second and third tests will only work properly if jyear()
# preserves fractional days, MoonPhase.pm states that days are returned as
# integers.
 
$start = timegm(0, 0, 0, 10, 5, 2002);
$stop = timegm(0, 0, 0, 11, 5, 2002);
printf "\nTest Astro::MoonPhase::phaseList(%d, %d), expecting New Moon at 2002-06-10 23:47 UTC:\n\n", $start, $stop;

    @name = ("New Moon", "First quarter", "Full moon", "Last quarter");
    ($phase, @times) = Astro::MoonPhase::phaselist($start, $stop);

    while (@times) {
      printf "%-14s= %s\n", $name[$phase], scalar localtime shift @times;
      $phase = ($phase + 1) % 4;
    }

$start = timegm(0, 0, 0, 10, 5, 2002) - 3600;
$stop = timegm(0, 0, 0, 11, 5, 2002) - 3600;
printf "\nTest Astro::MoonPhase::phaseList(%d, %d) for the same day with a simulated DST offset,\nthere should NOT be output in this case:\n\n", $start, $stop;

    @name = ("New Moon", "First quarter", "Full moon", "Last quarter");
    ($phase, @times) = Astro::MoonPhase::phaselist($start, $stop);

    while (@times) {
      printf "%-14s= %s\n", $name[$phase], scalar localtime 3600 + shift @times;
      $phase = ($phase + 1) % 4;
    }

$start = timegm(0, 0, 0, 11, 5, 2002) - 3600;
$stop = timegm(0, 0, 0, 12, 5, 2002) - 3600;
printf "\nTest Astro::MoonPhase::phaseList(%d, %d) for the next day with a simulated DST offset:\n\n", $start, $stop;

    @name = ("New Moon", "First quarter", "Full moon", "Last quarter");
    ($phase, @times) = Astro::MoonPhase::phaselist($start, $stop);

    while (@times) {
      printf "%-14s= %s\n", $name[$phase], scalar localtime 3600 + shift @times;
      $phase = ($phase + 1) % 4;
    }

printf "\nTest Astro::MoonPhase:phase(%d) for current date and time:\n\n", time;

   ( $MoonPhase,
     $MoonIllum,
     $MoonAge,
     $MoonDist,
     $MoonAng,
     $SunDist,
     $SunAng ) = Astro::MoonPhase::phase();

     print "MoonPhase  = $MoonPhase\n";
     print "MoonIllum  = $MoonIllum\n";
     print "MoonAge    = $MoonAge\n";
     print "MoonDist   = $MoonDist\n";
     print "MoonAng    = $MoonAng\n";
     print "SunDist    = $SunDist\n";
     print "SunAng     = $SunAng\n";

printf "\nTest Astro::MoonPhase::phasehunt(%d) for current date and time (UTC):\n\n", time;

    @phases = Astro::MoonPhase::phasehunt();
    print "New Moon      = ", scalar(gmtime($phases[0])), "\n";
    print "First quarter = ", scalar(gmtime($phases[1])), "\n";
    print "Full moon     = ", scalar(gmtime($phases[2])), "\n";
    print "Last quarter  = ", scalar(gmtime($phases[3])), "\n";
    print "New Moon      = ", scalar(gmtime($phases[4])), "\n";

$start = $phases[0];
$stop = $phases[4] + 1;
printf "\nTest Astro::MoonPhase::phaseList(%d, %d) for current lunation (UTC):\n\n", $start, $stop;

    @name = ("New Moon", "First quarter", "Full moon", "Last quarter");
    ($phase, @times) = Astro::MoonPhase::phaselist($start, $stop);

    while (@times) {
      my ($scratch, $sod, $secs);
      $scratch = int(shift @times);
      $sod = int($scratch / 86400) * 86400;
      $secs = $scratch - $sod;
      printf "%-14s= %s (%d, %d+%d)\n", $name[$phase], scalar gmtime $scratch, $scratch, $sod, $secs;
      $phase = ($phase + 1) % 4;
    }

print "\n";
0;

