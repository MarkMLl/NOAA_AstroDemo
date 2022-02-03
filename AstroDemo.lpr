program AstroDemo;

(* This demonstrates a Pascal transcript of the code underlying the sunrise     *)
(* etc. calculations at https://gml.noaa.gov/grad/solcalc/sunrise.html Note     *)
(* that that page displays equation of time etc. as minutes.decimal while this  *)
(* program uses minutes:seconds.                                MarkMLl         *)

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Classes, SysUtils, DateUtils, SunPosition, MoonPhase, IniFilesAbout;

const
  locationLongitude= 0.0;               (* Royston Priory being used as demo    *)
  locationLatitude= 52.0;               (* location.                            *)

var
  utc: TDateTime;
  unixSecs: longword;
  year, month, day, utcMins: word;
  jd, jc, sunriseUTC, noonUTC, sunsetUTC, elevation: double;
  recedingPhases, approachingPhases: TPhaseListResult;


function NowUTC(): TDateTime;

begin
{$if declared(UTC_Now) }
  result := UTC_Now()
{$else                 }
  result := Now()
{$endif                }
end { NowUtc } ;


function minsToTime(m: double): string;

var
  hours, mins, secs, msecs: integer;

begin
  result := '';
  hours := Trunc(m / 60);
  m := m - (hours * 60); // Residual integer minutes
  mins := Trunc(m);
  m := (m - mins) * 60; // Residual integer seconds
  secs := Trunc(m);
  mSecs := Trunc((m - secs) * 1000);
{$if declared(IsoFormatDateTime) }
  result := IsoFormatDateTime(EncodeTimeInterval(hours, mins, secs, mSecs), IsoTimeOnly, 0)
{$else                           }
  result := FormatDateTime('hh:nn:ss.zzz', EncodeTimeInterval(hours, mins, secs, mSecs))
{$endif                          }
end { minsToTime } ;


function dateTimeToDateAndTime(dt: TDateTime): string;

begin
{$if declared(IsoFormatDateTime) }
  result := IsoFormatDateTime(dt, IsoDateTime, 0)
{$else                           }
  result := DateToISO8601(dt, true)
{$endif                          }
end { dateTimeToDateAndTime } ;


function secsToDateAndTime(us: int64): string;

begin
{$if declared(IsoFormatDateTime) }
  result := IsoFormatDateTime(UnixToDateTime(us), IsoDateTime, 0)
{$else                           }
  result := DateToISO8601(UnixToDateTime(us), true)
{$endif                          }
end { SecsToDateAndTime } ;


begin
  utc := NowUTC();
  unixSecs := DateTimeToUnix(utc, true);
  DecodeDate(utc, year, month, day); (* Assume one-based                   *)
  utcMins := MinuteOfTheDay(utc);
  jd := CalcJD(year, month, day);
  jc := calcTimeJulianCent(jd);

(* Julian day checked against http://www.fourmilab.ch/documents/calendar. The   *)
(* result of these functions are in minutes relative to the start of the day.   *)

  sunriseUTC := calcSunriseUTC(jd, locationLatitude, locationLongitude);
  noonUTC := calcSolNoonUTC(jc, locationLongitude);
  sunsetUTC := calcSunsetUTC(jd, locationLatitude, locationLongitude);

(* Sunrise and Sunset (above) appear to have a constant "fudge factor" to       *)
(* to compensate for radius and refraction. Elevation (below) appears to apply  *)
(* a polynomial when the Sun is within 5 degrees of the horizon.                *)

  elevation := CalcSunElevation(jc, utcMins, locationLatitude, locationLongitude);

(* Report on the most recent New Moon and on the next phase.                    *)

  recedingPhases := PhaseList(unixSecs - Round(29.75 * 86400), unixSecs);
  approachingPhases := PhaseList(unixSecs, unixSecs + Round(29.75 * 86400));

(* Output the accumulated information.                                          *)

  WriteLn('Current minute:       ', minsToTime(utcMins));
  WriteLn('UTC Sunrise:          ', minsToTime(sunriseUTC));
  WriteLn('UTC Noon:             ', minsToTime(noonUTC));
  WriteLn('UTC Sunset:           ', minsToTime(sunsetUTC));
  Write  ('Solar elevation:      ', elevation:8:5, '°');
  if elevation >= 0.0 then
    WriteLn
  else
    case -Trunc(elevation / 6.0) of
      0: WriteLn(' (Civil twilight)');
      1: WriteLn(' (Nautical twilight)');
      2: WriteLn(' (Astronomical twilight)')
    otherwise
      WriteLn
    end;
  WriteLn('Start of Lunar Month: ', SecsToDateAndTime(recedingPhases.times[Length(recedingPhases.times) - approachingPhases.phase]                        ));
  WriteLn('Next Lunar phase:  Q', approachingPhases.phase, ' ',  SecsToDateAndTime(approachingPhases.times[0]))
end.

