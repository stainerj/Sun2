package com.example.sun;

import java.text.DecimalFormat;
import java.util.Scanner;

public class SunCalculator4 {

    public static Scanner keyboard = new Scanner(System.in);//input stream

    private static DecimalFormat df2 = new DecimalFormat("#0.00");
    private static DecimalFormat dfw = new DecimalFormat("##");
    private static DecimalFormat dfm = new DecimalFormat("00");

    public static double year = 0;
    public static double month = 0;
    public static double day = 0;
    public static double timeDay = 0;
    public static double julianDayRaw = 0;//Julian Day number not accounting for time of day. (Ends .5 - exact time at midnight)
    public static double julianDay = 0;//Julian Day number exact
    public static double meanAnomDeg = 0;//the mean anomaly in degrees
    public static double trueAnomaly = 0;//true anomaly
    public static double eclipticLongSun = 0;//ecliptic longitude of the sun
    public static double rASun = 0;//right ascension of sun
    public static double decSun = 0;//declination of sun
    public static double latitude = 0;//earth latitude, north positive   54.6 (nl 52)
    public static double longitude = 0;//earth longitude, east positive   -5.9 (nl 5)
    public static double jTransit = 0; //j (time of transit, measured in terms of Julian date)
    public static double jRise = 0;//sunrise time (Julian date)
    public static double jSet = 0;//sunset time (Julian date)

    public static final double J2K = 2451545;//J2000
    public static final double ELOP = 102.9373;//ecliptic longitude of perihelion
    public static final double OOTE = 23.4393;//obliquity of the ecliptic
    public static final double JZERO = 0.0009;
    public static final double JONE = 0.0053;
    public static final double JTWO = -0.0068;

    //main method
    public static void main(String[] args) {
        getDate();
        getTime();
        getPlace();
        dateTimeToJulianDay();
        System.out.println("Latitude/Longitude ------------------------------------------------------------ " + latitude + ", " + (-longitude));
        meanAnomaly();
        equationCentre();
        eclipticalCoordinates();
        equatorialCoordinates();
        observer();
        solarTransit();
        solarTransitTime();
        riseAndSet();
        riseAndSetTime();
    }

    //get calender date from user
    public static void getDate() {
        System.out.println("Enter year: ");
        year = keyboard.nextDouble();
        System.out.println("Enter month: ");
        month = keyboard.nextDouble();
        System.out.println("Enter day: ");
        day = keyboard.nextDouble();
    }

    //get time of day from user
    public static void getTime() {
        System.out.println("Enter time of day:(24hh.mm) ");
        timeDay = keyboard.nextDouble();
    }

    //get latitude and longitude from user
    public static void getPlace() {
        System.out.println("Enter latitude (north positive): ");
        latitude = keyboard.nextDouble();
        System.out.println("Enter longitude (east positive): ");
        longitude = -(keyboard.nextDouble());
    }

    //convert (gregorian) calender date to julian day number - exact method including time of day
    //exact whole number date gives 12.00 hours (noon), date given as e.g. 0.5 corresponds to 00.00 hours (midnight, UTC)
    public static double dateTimeToJulianDay() {
        double hour = Math.floor(timeDay);
        double minutes = timeDay - hour;
        double minDecFraction = (minutes/60)*100;
        double decTime = (hour + minDecFraction)/24;

        //keep year and month as entered for display
        double yearDisplay = year;
        double monthDisplay = month;

        if ((month == 1) || (month == 2)) {
            year = (year - 1);
            month = (month + 12);
        }

        double a = Math.floor(year / 100);
        double b = Math.floor(a / 4);
        double c = 2 - a + b;
        double e = Math.floor(365.25*(year+4716));
        double f = Math.floor(30.6001*(month+1));

        julianDayRaw = c + day + e + f - 1524.5;

        System.out.println("Date/time --------------------------------------------------------------------- " + dfw.format(yearDisplay) + "/"
                + dfw.format(monthDisplay) + "/" + dfw.format(day) + " " + df2.format(timeDay) + " (UTC)");
        System.out.println("Julian Day (date only) is ----------------------------------------------------- " + julianDayRaw);

        julianDay = julianDayRaw + decTime;
        System.out.println("Exact Julian Day Number (including time) is ----------------------------------- " + julianDay);
        return julianDay;
    }

    //calculate the mean anomaly in degrees
    public static double meanAnomaly() {
        final double M0 = 357.5291;//mean anomaly degrees
        final double M1 = 0.98560028;//mean anomaly degrees per day

        double theta= M0 + M1*(julianDay - J2K);
        meanAnomDeg = theta%360;

        System.out.println("Mean anomaly, earth (degrees) ------------------------------------------------- " + meanAnomDeg);
        return meanAnomDeg;
    }

    //calculate the equation of centre - trig functions take arguments in radians
    public static double equationCentre() {
        final double C1 = 1.9148;
        final double C2 = 0.0200;
        final double C3 = 0.0003;

        //convert to radians
        double meanAnomRad = degreesToRadians(meanAnomDeg);
        double equationOfCentre = (C1 * Math.sin(meanAnomRad)) + (C2 * Math.sin(2 * meanAnomRad)) + (C3 * Math.sin(3 * meanAnomRad));
        trueAnomaly = meanAnomDeg + equationOfCentre;

        System.out.println("Equation of centre ------------------------------------------------------------ " + equationOfCentre);
        System.out.println("True anomaly ------------------------------------------------------------------ " + trueAnomaly);
        return trueAnomaly;
    }

    //calculates the ecliptic longitude of the sun
    public static double eclipticalCoordinates() {
        double lambda = trueAnomaly + ELOP + 180;
        eclipticLongSun = lambda%360;

        System.out.println("Ecliptic longitude of sun ----------------------------------------------------- " + eclipticLongSun);
        return eclipticLongSun;
    }

    //calculates the right ascension and declination of the sun
    public static double equatorialCoordinates() {
        //constant approximations for right ascension and declination (earth)
        final double A2 = -2.4657;
        final double A4 = 0.0529;
        final double A6 = -0.0014;
        final double D1 = 22.7908;
        final double D3 = 0.5991;
        final double D5 = 0.0492;

        //convert to radians
        double eclipticLongSunRad = degreesToRadians(eclipticLongSun);
        double sinLambdaSun = Math.sin(eclipticLongSunRad); //sine of the ecliptic longitude
        double cosLambdaSun = Math.cos(eclipticLongSunRad);
        //calculate right ascension
        rASun = eclipticLongSun + (A2*Math.sin(2*eclipticLongSunRad)) + (A4*Math.sin(4*eclipticLongSunRad)) +
                (A6*Math.sin(6*eclipticLongSunRad));
        //alternative method using αsun = arctan(sin λsun cos ε, cos λsun). Returns angle in radians.
        double rASunAlternativeMethod = Math.atan2(sinLambdaSun*Math.cos(degreesToRadians(OOTE)), cosLambdaSun);
        rASunAlternativeMethod = radiansToDegrees(rASunAlternativeMethod);
        //calculate declination
        decSun = (D1*Math.sin(eclipticLongSunRad)) + (D3*Math.pow(sinLambdaSun, 3)) + (D5*Math.pow(sinLambdaSun, 5));
        //alternative method using δsun = arcsin(sin λsun sin ε). Returns angle in radians.
        double decSunAlternativeMethod = Math.asin(sinLambdaSun*Math.sin(degreesToRadians(OOTE)));
        decSunAlternativeMethod = radiansToDegrees(decSunAlternativeMethod);

        System.out.println("Sun's right ascension ------  " + rASun + " or (alternative method) ------ " + rASunAlternativeMethod);
        System.out.println("Sun's declination ---------- " + decSun + " or (alternative method) ------- " + decSunAlternativeMethod);
        return rASun;
    }

    //calculates sun position for earth observer
    public static double observer() {
        //sidereal constants for earth
        final double THETAZERO = 280.1470;
        final double THETAONE = 360.9856235;

        double thetaSidereal = THETAZERO + THETAONE * (julianDay - J2K) - longitude; //raw angle, remove multiples of 360
        double sidereal = thetaSidereal%360; //sidereal time
        double hourAngle = sidereal - rASun;

        System.out.println("Sidereal time (degrees) ------------------------------------------------------- " + sidereal);
        System.out.println("Hour angle -------------------------------------------------------------------- " + hourAngle);

        //convert to radians
        double hourAngleRad = degreesToRadians(hourAngle);
        double latitudeRad = degreesToRadians(latitude);
        double decSunRad = degreesToRadians(decSun);
        double azimuth = Math.atan2(Math.sin(hourAngleRad), Math.cos(hourAngleRad)*Math.sin(latitudeRad) - Math.tan(decSunRad)*Math.cos(latitudeRad));
        azimuth = radiansToDegrees(azimuth) + 180; //convert from radians to degrees, add 180 degrees for conventional azimuth (north=0)

        double altitude = Math.asin(Math.sin(latitudeRad)*Math.sin(decSunRad) + Math.cos(latitudeRad)*Math.cos(decSunRad)* Math.cos(hourAngleRad));
        altitude = radiansToDegrees(altitude); //convert from radians to degrees
        System.out.println("Azimuth ----------------------------------------------------------------------- " + azimuth);
        System.out.println("Altitude ---------------------------------------------------------------------- " + altitude);
        //azimuth 0 = north
        return sidereal;
    }

    //calculate transit time (julian day)
    public static double solarTransit() {
        double n;
        double nNearInt; //nearest whole integer to n
        double jApprox; //approx value of j

        double meanAnomRad = degreesToRadians(meanAnomDeg);
        double lSun = meanAnomDeg + ELOP + 180;
        double lSunRad = degreesToRadians(lSun);

        //n(*) = (J - J2000 - J0)
        n = (julianDay-J2K-JZERO) - longitude/360;
        nNearInt = Math.round(n);//nearest whole number
        jApprox = J2K + JZERO + longitude/360 + nNearInt;
        //Jtransit J(*) + J1 sin M + J2 sin(2 Lsun)
        jTransit = jApprox + (JONE*(Math.sin(meanAnomRad))) + (JTWO*(Math.sin(2*lSunRad)));
        //(repetition method available to increase accuracy) jTransit = jTransit - hourAngleForJTransit/360
        System.out.println("Time of transit (Julian Day) -------------------------------------------------- " + jTransit);
        return jTransit;
    }

    //calculate transit time (hours & minutes)
    public static void solarTransitTime() {
        System.out.println("Time of transit (hours and minutes) ------------------------------------------- " + julianDayToTime(jTransit, julianDayRaw));
    }

    //calculate sunrise & sunset times
    public static void riseAndSet() {
        //constants for daytime length coefficients (earth, degrees)
        //H1 = 22.137; H3 = 0.599; H5 = 0.016;
        //constant for atmospheric refraction effect and apparent solar disc diameter (earth, degrees)
        final double H0 = -0.83;
        final double DSUN = 0.53;

        //H = arccos((sin h₀ − sin φ sin δ)/(cos φ cos δ))
        double H =  Math.acos((Math.sin(degreesToRadians(H0)) - (Math.sin(degreesToRadians(latitude)) *
                Math.sin(degreesToRadians(decSun))))/(Math.cos(degreesToRadians(latitude)) *
                Math.cos(degreesToRadians(decSun))));
        System.out.println("Value of H (degrees) is ------------------------------------------------------- " + radiansToDegrees(H));

        //Jrise = jTransit - H/360 (*J3=1)
        jRise = jTransit - (radiansToDegrees(H)/360);
        System.out.println("Time of sunrise (Julian day) -------------------------------------------------- " + jRise);
        //Jset = jTransit + H/360 (*J3=1)
        jSet = jTransit + (radiansToDegrees(H)/360);
        System.out.println("Time of sunset (Julian day) --------------------------------------------------- " + jSet);
    }

    public static void riseAndSetTime() {
        System.out.println("Time of sunrise (hours and minutes) ------------------------------------------- " + julianDayToTime(jRise, julianDayRaw));
        System.out.println("Time of sunset (hours and minutes) -------------------------------------------- " + julianDayToTime(jSet, julianDayRaw));
    }

    //converts Julian day number to time of day in hours and minutes
    public static String julianDayToTime (double jTime, double julianDayZero) {
        double timeHours = (jTime - julianDayZero)*24;
        double timeMinutes = (timeHours - (Math.floor(timeHours)));
        timeHours = timeHours - timeMinutes;
        if (timeHours < 0) {
            timeHours = timeHours + 24;
        }
        if (timeHours >= 24) {
            timeHours = timeHours - 24;
        }
        timeMinutes = timeMinutes*60;
        //System.out.println("Time of  (hours and minutes) -------------------------------------------- " + dfw.format(timeHours) + "."
        //+ dfm.format(timeMinutes) + " (UTC)");
        String outputString = dfw.format(timeHours) + "." + dfm.format(timeMinutes) + " (UTC)";
        return outputString;
    }

    //convert degrees to radians
    public static double degreesToRadians(double degrees)
    {
        double radians = degrees * 0.0174532925199;
        return radians;
    }

    //convert radians to degrees
    public static double radiansToDegrees(double radians)
    {
        double degrees = radians / 0.0174532925199;
        return degrees;
    }
}
