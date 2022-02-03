(* Lazarus+FPC 0.9.24+2.2.0 to 0.9.30+2.6.0. On Linux for ARM, PPC, SPARC, x86. *)

unit IniFilesExtended;

(* Add stuff-that-I've-found-useful to TIniFile.                                *)

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, IniFiles;

TYPE    TIniState= (IniNotLoaded, IniLoadedLocked, IniLoaded, IniRunning, IniSaved, IniSavedLocked);

        TIniFileEx= CLASS(TIniFile)
        PRIVATE
          FIniState: TIniState;
        PUBLIC
	  FUNCTION ReadBoolean(CONST section, ident: STRING; default: BOOLEAN): BOOLEAN;
          PROCEDURE WriteBoolean(CONST section, ident: STRING; value: BOOLEAN);
	  FUNCTION ReadBooleanQualified(CONST section, ident, qualifier: STRING;
                                                default: BOOLEAN): BOOLEAN;
          FUNCTION ReadDegrees(CONST section, ident: STRING; default: DOUBLE): DOUBLE;
          FUNCTION ReadDegrees(CONST section, ident: STRING; default: STRING): DOUBLE;
          PROCEDURE WriteDegrees(CONST section, ident: STRING; value: DOUBLE);

          FUNCTION ReadStringAlt(CONST section, sectionAlt, ident, default: STRING): STRING;
          FUNCTION ReadIntegerAlt(CONST section, sectionAlt, ident: STRING; default: INTEGER): INTEGER;
          FUNCTION ReadFloatAlt(CONST section, sectionAlt, ident: STRING; default: DOUBLE): DOUBLE;
          FUNCTION ReadBooleanAlt(CONST section, sectionAlt, ident: STRING; default: BOOLEAN): BOOLEAN;
          FUNCTION ReadDegreesAlt(CONST section, sectionAlt, ident: STRING; default: DOUBLE): DOUBLE;
          FUNCTION ReadDegreesAlt(CONST section, sectionAlt, ident: STRING; default: STRING): DOUBLE;

          PROCEDURE WriteString3(CONST section, ident: STRING; index: INTEGER; CONST field, value: STRING);
          PROCEDURE WriteInteger3(CONST section, ident: STRING; index: INTEGER; CONST field: STRING; value: INTEGER);
          PROCEDURE WriteFloat3(CONST section, ident: STRING; index: INTEGER; CONST field: STRING; value: DOUBLE);
          PROCEDURE WriteBoolean3(CONST section, ident: STRING; index: INTEGER; CONST field: STRING; value: BOOLEAN);
          PROCEDURE WriteDegrees3(CONST section, ident: STRING; index: INTEGER; CONST field: STRING; value: DOUBLE);

          FUNCTION ReadStringAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field, default: STRING): STRING;
          FUNCTION ReadIntegerAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: INTEGER): INTEGER;
          FUNCTION ReadFloatAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: DOUBLE): DOUBLE;
          FUNCTION ReadBooleanAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: BOOLEAN): BOOLEAN;
          FUNCTION ReadDegreesAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: DOUBLE): DOUBLE;
          FUNCTION ReadDegreesAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: STRING): DOUBLE;
(*$IFDEF LCL *)
          FUNCTION ReadColour(CONST section, ident: STRING; default: LONGINT): LONGINT;
          PROCEDURE WriteColour(CONST section, ident: STRING; colour: LONGINT);
          FUNCTION ReadColourAlt(CONST section, sectionAlt, ident: STRING; default: LONGINT): LONGINT;
          PROCEDURE WriteColour3(CONST section, ident: STRING; index: INTEGER; CONST field: STRING; value: LONGINT);
          FUNCTION ReadColourAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: LONGINT): LONGINT;
(*$ENDIF *)
          PROPERTY IniState: TIniState READ FIniState WRITE FIniState;
        END;


implementation

USES StrUtils (*$IFDEF LCL *) , Graphics (*$ENDIF *) ;


FUNCTION TIniFileEx.ReadBoolean(CONST section, ident: STRING; default: BOOLEAN): BOOLEAN;

(* Make sure that TRUE, FALSE, YES and NO are recognised as Boolean values as   *)
(* well as numeric 1 and 0.					                *)

(* Liberated from Dialarm.                                                      *)

VAR	temp: STRING;

BEGIN
  ReadBoolean:= default;
  IF default THEN
    temp:= 'TRUE'
  ELSE
    temp:= 'FALSE';
  temp:= Trim(ReadString(section, ident, temp));
  IF UpperCase(temp) = 'TRUE' THEN
    ReadBoolean:= TRUE;
  IF UpperCase(temp) = 'FALSE' THEN
    ReadBoolean:= FALSE;
  IF UpperCase(temp) = 'YES' THEN
    ReadBoolean:= TRUE;
  IF UpperCase(temp) = 'NO' THEN
    ReadBoolean:= FALSE;
  IF temp = '1' THEN
    ReadBoolean:= TRUE;
  IF temp = '0' THEN
    ReadBoolean:= FALSE
END; (* TIniFileEx.ReadBool*)


PROCEDURE TIniFileEx.WriteBoolean(CONST section, ident: STRING; value: BOOLEAN);

(* Write a Boolean, respecting any existing style.                              *)

TYPE    Tstyle= (n, tf, yn, ctf, cyn);

VAR	temp: STRING;
        style: Tstyle;

BEGIN
  style:= n;
  temp:= Trim(ReadString(section, ident, ''));
  IF temp <> '' THEN
    CASE temp[1] OF
      't', 'f': style:= tf;
      'y', 'n': style:= yn;
      'T', 'F': style:= ctf;
      'Y', 'N': style:= cyn
    ELSE
    END;
  CASE style OF
    tf:  IF value THEN
           WriteString(section, ident, 'true')
         ELSE
           WriteString(section, ident, 'false');
    yn:  IF value THEN
           WriteString(section, ident, 'yes')
         ELSE
           WriteString(section, ident, 'no');
    ctf: IF value THEN
           WriteString(section, ident, 'True')
         ELSE
           WriteString(section, ident, 'False');
    cyn: IF value THEN
           WriteString(section, ident, 'Yes')
         ELSE
           WriteString(section, ident, 'No')
  ELSE
    WriteBool(section, ident, value)
  END
END { TIniFileEx.WriteBoolean } ;


FUNCTION TIniFileEx.ReadBooleanQualified(CONST section, ident, qualifier: STRING;
                                                default: BOOLEAN): BOOLEAN;

(* This is similar to ReadBoolean() except that it expects to find an entry in  *)
(* the .ini file like:						                *)
(*									        *)
(*	[Location]							        *)
(*	allowScriptedShellOps = fork, cp, mv				        *)
(*									        *)
(* and returns a Boolean depending upon the presence or otherwise of the        *)
(* qualifier parameter (case insensitive) in the string. Leading and trailing   *)
(* parentheses are ignored, embedded commas are treated as whitespace, "any",   *)
(* "all" and "none" are handled as expected.		                        *)

(* Liberated from Dialarm.                                                      *)

VAR	temp: STRING;
	index: INTEGER;

BEGIN
  ReadBooleanQualified:= default;
  IF default THEN
    temp:= 'any'
  ELSE
    temp:= 'none';
  temp:= Trim(ReadString(section, ident, temp));
  IF temp[1] = '(' THEN
    temp[1]:= ' ';
  IF temp[Length(temp)] = ')' THEN
    temp[Length(temp)]:= ' ';
  temp:= Trim(ReadString(section, ident, temp));
  FOR index:= 1 TO Length(temp) DO
    IF temp[index] <= ' ' THEN
      temp[index]:= ',';		(* Use , as separator		        *)
  temp:= ',' + UpperCase(temp) + ',';	(* Telomeres at ends		        *)
  index:= Pos(',ALL,', temp);
  IF index <> 0 THEN BEGIN		(* Convert "all" to "any"	        *)
    temp[index + 2]:= 'N';
    temp[index + 3]:= 'Y'
  END;
  IF Pos(',NONE,', temp) > 0 THEN
    ReadBooleanQualified:= FALSE
  ELSE
    IF Pos(',ANY,', temp) > 0 THEN
      ReadBooleanQualified:= TRUE
    ELSE
      ReadBooleanQualified:= Pos(UpperCase(Qualifier), temp) > 0
END; (* TIniFileEx.ReadBoolQualified *)


FUNCTION TIniFileEx.ReadDegrees(CONST section, ident: STRING; default: DOUBLE): DOUBLE;

(* Read a floating point number formatted as degrees, minutes and seconds. The  *)
(* degrees field is separated from the remainder by anything which isn't a      *)
(* digit to accomodate a degree symbol in any form, the minutes and seconds     *)
(* fields may be identified by ' and " respectively possibly with whitespace    *)
(* but these may be omitted.                                                    *)

BEGIN
  RESULT:= ReadDegrees(section, ident, FloatToStr(default))
END { TIniFileEx.ReadDegrees }  ;


FUNCTION parseDegrees(value: STRING; OUT ok: BOOLEAN): DOUBLE;

(* Parse two integers and a float as degrees, minutes and seconds. On error     *)
(* return with ok set FALSE and an indeterminate result.                        *)

VAR     i: INTEGER;
        degrees, minutes, seconds: DOUBLE;


  (* Parse the longest possible number from the start of the string, deleting   *)
  (* what we've parsed or ORing false into the ok parameter if there's text but *)
  (* it's not a number. No text at all is valid, and interpreted as 0.0.        *)
  //
  function parseFP(var value: string; var ok: boolean): double;

  type  charSet= set of CHAR;

  const digits: charSet= ['0'..'9', '.', '-', 'E', 'e'];
        bad= -99999999;

  var   scratch: string= '';
        ok2: boolean;

  begin
    result := 0.0;
    if value = '' then
      exit;
    while (value <> '') and (value[1] in digits) do begin
      scratch += value[1];
      Delete(value, 1, 1)
    end;
    value := TrimLeft(value);
    result := StrToFloatDef(scratch, bad);
    if result = bad then
      ok := ok and false
  end { parseFP } ;


BEGIN
  ok:= TRUE;
  result := 0.0;

(* Assuming that we don't know the codepage etc., assume that anything outside  *)
(* the standard ASCII range represents a degree symbol and obliterate it.       *)

  for i := 1 to Length(value) do
    if (value[i] < ' ') or (value[i] >= #$7f) then
      value[i] := ' ';
  value := Trim(DelSpace1(value));

(* Assume that we might have embedded ASCII ' and " characters identifying      *)
(* minutes and seconds. Parse degrees, but be prepared to decide that we've     *)
(* really got minutes or seconds.                                               *)

  seconds := 0.0;
  minutes := 0.0;
  degrees := parseFP(value, ok);
  if value <> '' then
    case value[1] of
      '''': begin
              minutes := degrees;
              degrees := 0.0;
              Delete(value, 1, 1)
            end;
      '"':  begin
              seconds := degrees;
              minutes := 0.0;
              degrees := 0.0;
              Delete(value, 1, 1)
            end
    otherwise
    end;
  value := TrimLeft(value);

(* Parse minutes, but be prepared to admit that we've really got seconds. We    *)
(* have to be careful here to not overwrite an existing minutes or seconds      *)
(* value if parseFP() sees an empty parameter.                                  *)

  if value <> '' then begin
    minutes := parseFP(value, ok);
    if value <> '' then
      case value[1] of
        '''': Delete(value, 1, 1);
        '"':  begin
                seconds := minutes;
                minutes := 0.0;
                Delete(value, 1, 1)
              end
      otherwise
      end;
    value := TrimLeft(value)
  end;

(* Parse seconds, avoiding overwrite. Expect nothing more than " to be left.    *)

  if value <> '' then
    seconds := parseFP(value, ok);
  ok := ok and ((value = '') or (value = '"'));
  result := degrees + (minutes / 60.0) + (seconds / 3600.0);
  while result > 360.0 do
    result -= 360.0;
  while result < -360.0 do
    result += 360.0
END { parseDegrees } ;


TYPE    DegreeConversionException= CLASS(Exception);

FUNCTION parseDegrees(const value: STRING): DOUBLE;

(* Parse two integers and a float as degrees, minutes and seconds. On error     *)
(* raise an exception.                                                          *)

VAR     ok: BOOLEAN;

BEGIN
  RESULT:= parseDegrees(value, ok);
  IF NOT ok THEN
    RAISE DegreeConversionException.Create('Bad degrees-minutes-seconds conversion')
END { parseDegrees } ;


procedure testParseDegrees;

const   a= 45.0 + 1/60 + 1/3600;        (* Should be evaluated by compiler      *)
        d= 0.1 / 3600;                  (* 0.1 second of arc                    *)

var     b, c: double;

begin
  b := 45.0 + 1/60 + 1/3600;            (* Should be exact match with a, but    *)
  Assert(Abs(a - b) < (d / 100));       (* allow slop for cross compilers.      *)
  c:= parseDegrees('45 1'' 1"');
  Assert(Abs(b - c) < d);
  c:= parseDegrees('45 1 1');
  Assert(Abs(b - c) < d);
  c := parseDegrees('45 1.0167');
  Assert(Abs(b - c) < d);
  c := parseDegrees('45 0 61');
  Assert(Abs(b - c) < d);
  c := parseDegrees('45 61"');
  Assert(Abs(b - c) < d);
  c := parseDegrees('45.016944');
  Assert(Abs(b - c) < d);
  b := parseDegrees('1');
  c := parseDegrees('361');
  Assert(b = c);
  b := parseDegrees('-1');
  c := parseDegrees('-361');
  Assert(b = c);
  b := parseDegrees('1');
  c := parseDegrees('-1');
  Assert(b <> c)
end { testParseDegrees } ;


FUNCTION TIniFileEx.ReadDegrees(CONST section, ident: STRING; default: STRING): DOUBLE;

(* Read a floating point number formatted as degrees, minutes and seconds. The  *)
(* degrees field is separated from the remainder by anything which isn't a      *)
(* digit to accomodate a degree symbol in any form, the minutes and seconds     *)
(* fields may be identified by ' and " respectively possibly with whitespace    *)
(* but these may be omitted.                                                    *)

CONST   bad= 'XXXXXXXX';

VAR     value: STRING;
        ok: BOOLEAN;

BEGIN
  value:= ReadString(section, ident, bad);
  IF value <> bad THEN BEGIN
    RESULT:= parseDegrees(value, ok);
    IF NOT ok THEN
      RESULT:= parseDegrees(default)
  END ELSE
    RESULT:= parseDegrees(default)
END { TIniFileEx.ReadDegrees }  ;


PROCEDURE TIniFileEx.WriteDegrees(CONST section, ident: STRING; value: DOUBLE);

(* Write a floating point number formatted as degrees, minutes and seconds with *)
(* decimals if necessary, e.g. 123 45' 6.78".                                   *)

VAR     degrees, minutes: INTEGER;

BEGIN
  degrees:= Trunc(value);
  value:= (value - degrees) * 60.0;     (* Now minutes                          *)
  minutes:= Trunc(value);
  value:= (value - minutes) * 60.0;     (* Now seconds                          *)
  WriteString(section, ident, IntToStr(degrees) + ' ' + IntToStr(minutes) +
                                                ''' ' + FloatToStr(value) + '"')
END { TIniFileEx.WriteDegrees } ;


FUNCTION TIniFileEx.ReadStringAlt(CONST section, sectionAlt, ident, default: STRING): STRING;

(* Similar to the standard function but tries to read the alternative before    *)
(* falling back to a hardcoded value.                                           *)

BEGIN
  RESULT:= ReadString(section, ident, ReadString(sectionAlt, ident, default))
END { TIniFileEx.ReadStringAlt } ;


FUNCTION TIniFileEx.ReadIntegerAlt(CONST section, sectionAlt, ident: STRING; default: INTEGER): INTEGER;

(* Similar to the standard function but tries to read the alternative before    *)
(* falling back to a hardcoded value.                                           *)

BEGIN
  RESULT:= ReadInteger(section, ident, ReadInteger(sectionAlt, ident, default))
END { TIniFileEx.ReadIntegerAlt } ;


FUNCTION TIniFileEx.ReadFloatAlt(CONST section, sectionAlt, ident: STRING; default: DOUBLE): DOUBLE;

(* Similar to the standard function but tries to read the alternative before    *)
(* falling back to a hardcoded value.                                           *)

BEGIN
  RESULT:= ReadFloat(section, ident, ReadFloat(sectionAlt, ident, default))
END { TIniFileEx.ReadFloatAlt } ;


FUNCTION TIniFileEx.ReadBooleanAlt(CONST section, sectionAlt, ident: STRING; default: BOOLEAN): BOOLEAN;

(* Similar to the standard function but tries to read the alternative before    *)
(* falling back to a hardcoded value.                                           *)

BEGIN
  RESULT:= ReadBoolean(section, ident, ReadBoolean(sectionAlt, ident, default))
END { TIniFileEx.ReadBooleanAlt } ;


FUNCTION TIniFileEx.ReadDegreesAlt(CONST section, sectionAlt, ident: STRING; default: DOUBLE): DOUBLE;

(* Similar to the standard function but tries to read the alternative before    *)
(* falling back to a hardcoded value.                                           *)

BEGIN
  RESULT:= ReadDegrees(section, ident, ReadDegrees(sectionAlt, ident, default))
END { TIniFileEx.ReadDegreesAlt } ;


FUNCTION TIniFileEx.ReadDegreesAlt(CONST section, sectionAlt, ident: STRING; default: STRING): DOUBLE;

(* Similar to the standard function but tries to read the alternative before    *)
(* falling back to a hardcoded value.                                           *)

BEGIN
  RESULT:= ReadDegrees(section, ident, ReadDegrees(sectionAlt, ident, default))
END { TIniFileEx.ReadDegreesAlt } ;


FUNCTION makeIdent3(CONST ident: STRING; index: INTEGER; CONST field: STRING): STRING;

(* Assemble the triplet of ident, index and field into an identifier,  assuming *)
(* that 'index' >= 0 and 'field' <> ''.                                         *)

BEGIN
  RESULT:= ident;
  IF index >= 0 THEN
    RESULT:= RESULT + '[' + IntToStr(index) + ']';
  IF field <> '' THEN
    RESULT:= RESULT + '.' + field
END { makeIdent3 } ;


PROCEDURE TIniFileEx.WriteString3(CONST section, ident: STRING; index: INTEGER; CONST field, value: STRING);

(* Similar to the standard WriteX methods except that the identifier is         *)
(* formatted as ident[index].field assuming that 'index' >= 0 and 'field' <> '' *)

BEGIN
  WriteString(section, makeIdent3(ident, index, field), value)
END { TIniFileEx.WriteString3 } ;


PROCEDURE TIniFileEx.WriteInteger3(CONST section, ident: STRING; index: INTEGER; CONST field: STRING; value: INTEGER);

(* Similar to the standard WriteX methods except that the identifier is         *)
(* formatted as ident[index].field assuming that 'index' >= 0 and 'field' <> '' *)

BEGIN
  WriteInteger(section, makeIdent3(ident, index, field), value)
END { TIniFileEx.WriteInteger3 } ;


PROCEDURE TIniFileEx.WriteFloat3(CONST section, ident: STRING; index: INTEGER; CONST field: STRING; value: DOUBLE);

(* Similar to the standard WriteX methods except that the identifier is         *)
(* formatted as ident[index].field assuming that 'index' >= 0 and 'field' <> '' *)

BEGIN
  WriteFloat(section, makeIdent3(ident, index, field), value)
END { TIniFileEx.WriteFloat3 } ;


PROCEDURE TIniFileEx.WriteBoolean3(CONST section, ident: STRING; index: INTEGER; CONST field: STRING; value: BOOLEAN);

(* Similar to the standard WriteX methods except that the identifier is         *)
(* formatted as ident[index].field assuming that 'index' >= 0 and 'field' <> '' *)

BEGIN
  WriteBoolean(section, makeIdent3(ident, index, field), value)
END { TIniFileEx.WriteBoolean3 } ;


PROCEDURE TIniFileEx.WriteDegrees3(CONST section, ident: STRING; index: INTEGER; CONST field: STRING; value: DOUBLE);

(* Similar to the standard WriteX methods except that the identifier is         *)
(* formatted as ident[index].field assuming that 'index' >= 0 and 'field' <> '' *)

BEGIN
  WriteDegrees(section, makeIdent3(ident, index, field), value)
END { TIniFileEx.WriteDegrees3 } ;


FUNCTION TIniFileEx.ReadStringAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field, default: STRING): STRING;

(* Similar to the standard ReadXAlt methods except that the identifier is       *)
(* formatted as ident[index].field assuming that 'index' >= 0 and 'field' <> '' *)

BEGIN
  RESULT:= ReadStringAlt(section, sectionAlt, makeIdent3(ident, index, field), default)
END { TIniFileEx.ReadStringAlt3 } ;


FUNCTION TIniFileEx.ReadIntegerAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: INTEGER): INTEGER;

(* Similar to the standard ReadXAlt methods except that the identifier is       *)
(* formatted as ident[index].field assuming that 'index' >= 0 and 'field' <> '' *)

BEGIN
  RESULT:= ReadIntegerAlt(section, sectionAlt, makeIdent3(ident, index, field), default)
END { TIniFileEx.ReadIntegerAlt3 } ;


FUNCTION TIniFileEx.ReadFloatAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: DOUBLE): DOUBLE;

(* Similar to the standard ReadXAlt methods except that the identifier is       *)
(* formatted as ident[index].field assuming that 'index' >= 0 and 'field' <> '' *)

BEGIN
  RESULT:= ReadFloatAlt(section, sectionAlt, makeIdent3(ident, index, field), default)
END { TIniFileEx.ReadFloatAlt3 } ;


FUNCTION TIniFileEx.ReadBooleanAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: BOOLEAN): BOOLEAN;

(* Similar to the standard ReadXAlt methods except that the identifier is       *)
(* formatted as ident[index].field, assuming that 'index' >= 0 and 'field' <> ''*)

BEGIN
  RESULT:= ReadBooleanAlt(section, sectionAlt, makeIdent3(ident, index, field), default)
END { TIniFileEx.ReadBooleanAlt3 } ;


FUNCTION TIniFileEx.ReadDegreesAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: DOUBLE): DOUBLE;

(* Similar to the standard ReadXAlt methods except that the identifier is       *)
(* formatted as ident[index].field, assuming that 'index' >= 0 and 'field' <> ''*)

BEGIN
  RESULT:= ReadDegreesAlt(section, sectionAlt, makeIdent3(ident, index, field), default)
END { TIniFileEx.ReadDegreesAlt3 } ;


FUNCTION TIniFileEx.ReadDegreesAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: STRING): DOUBLE;

(* Similar to the standard ReadXAlt methods except that the identifier is       *)
(* formatted as ident[index].field, assuming that 'index' >= 0 and 'field' <> ''*)

BEGIN
  RESULT:= ReadDegreesAlt(section, sectionAlt, makeIdent3(ident, index, field), default)
END { TIniFileEx.ReadDegreesAlt3 } ;


(*$IFDEF LCL *)

FUNCTION TIniFileEx.ReadColour(CONST section, ident: STRING; default: LONGINT): LONGINT;

(* Read a colour, applying a 'cl' prefix if necessary.                          *)

VAR     str: STRING;

BEGIN
  str:= ReadString(section, ident, '');
  IF Pos('cl', str) <> 1 THEN
    str:= 'cl' + str;
  IF NOT IdentToColor(str, RESULT) THEN
    RESULT:= default
END { TIniFileEx.ReadColour } ;


PROCEDURE TIniFileEx.WriteColour(CONST section, ident: STRING; colour: LONGINT);

(* Write a colour, removing the 'cl' prefix and assuming that the parameter is  *)
(* valid (i.e. not an arbitrary number which cannot be a real colour).          *)

VAR     str: STRING;

BEGIN
  IF NOT ColorToIdent(colour, str) THEN
    EXIT;
  Delete(str, 1, 2);
  WriteString(section, ident, str)
END { TIniFileEx.WriteColour } ;

FUNCTION colourToIdentDef(colour: LONGINT; CONST default: STRING): STRING;

(* Convert a numeric colour to a string, defaulting as given.                   *)

BEGIN
  IF NOT ColorToIdent(colour, RESULT) THEN
    RESULT:= default;
  IF Pos('cl', RESULT) = 1 THEN
    Delete(RESULT, 1, 2)
END { colourToIdentDef } ;


FUNCTION identToColourDef(str: STRING; default: LONGINT): LONGINT;

(* Convert a named colour to a numeric, defaulting as given.                    *)

BEGIN
  IF Pos('cl', str) <> 1 THEN
    str:= 'cl' + str;
  IF NOT IdentToColor(str, RESULT) THEN
    RESULT:= default
END { identToColourDef } ;


FUNCTION TIniFileEx.ReadColourAlt(CONST section, sectionAlt, ident: STRING; default: LONGINT): LONGINT;

(* Similar to the standard function but tries to read the alternative before    *)
(* falling back to a hardcoded value.                                           *)

BEGIN
  RESULT:= ReadColour(section, ident, ReadColour(sectionAlt, ident, default))
END { TIniFileEx.ReadColourAlt } ;


PROCEDURE TIniFileEx.WriteColour3(CONST section, ident: STRING; index: INTEGER; CONST field: STRING; value: LONGINT);

(* Similar to the standard WriteX methods except that the identifier is         *)
(* formatted as ident[index].field assuming that 'index' >= 0 and 'field' <> '' *)

BEGIN
  WriteColour(section, makeIdent3(ident, index, field), value)
END { TIniFileEx.WriteColour3 } ;


FUNCTION TIniFileEx.ReadColourAlt3(CONST section, sectionAlt, ident: STRING; index: INTEGER; CONST field: STRING; default: LONGINT): LONGINT;

(* Similar to the standard ReadXAlt methods except that the identifier is       *)
(* formatted as ident[index].field, assuming that 'index' >= 0 and 'field' <> ''*)

BEGIN
  RESULT:= ReadColourAlt(section, sectionAlt, makeIdent3(ident, index, field), default)
END { TIniFileEx.ReadColourAlt3 } ;

(*$ENDIF *)


begin
  testParseDegrees
end.

