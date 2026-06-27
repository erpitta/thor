"""
Julian Day and Local Sidereal Time Calculator
Calculates Julian Day Number and Local Sidereal Time from date, time, and geographical position.
"""

import numpy as np
from datetime import datetime
import pytz


def gregorian_to_julian_day(year, month, day, hour, minute, second):
    """
    Convert Gregorian calendar date and time to Julian Day Number (JD).
    
    Parameters:
    -----------
    year : int
        Year (e.g., 2026)
    month : int
        Month (1-12)
    day : int
        Day of month (1-31)
    hour : int
        Hour (0-23)
    minute : int
        Minute (0-59)
    second : float
        Seconds (0-59.999...)
    
    Returns:
    --------
    jd : float
        Julian Day Number (including fractional day)
    jd_integer : float
        Integer part of JD (noon on reference day)
    """
    
    # Adjust for hours/minutes/seconds as fractional day
    day_fraction = (hour + minute / 60.0 + second / 3600.0) / 24.0
    
    # Algorithm from USNO (U.S. Naval Observatory)
    # Adjust month and year for Jan/Feb (astronomical convention)
    if month <= 2:
        year -= 1
        month += 12
    
    # Calculate Julian Day Number
    A = np.floor(year / 100.0)
    B = 2 - A + np.floor(A / 4.0)
    
    JD = (np.floor(365.25 * (year + 4716)) + 
          np.floor(30.6001 * (month + 1)) + 
          day + day_fraction + B - 1524.5)
    
    # Integer JD (at noon UTC)
    JD_integer = np.floor(JD + 0.5)
    
    return JD, JD_integer


def julian_day_to_greenwich_mean_sidereal_time(jd):
    """
    Calculate Greenwich Mean Sidereal Time (GMST) from Julian Day Number.
    
    Parameters:
    -----------
    jd : float
        Julian Day Number
    
    Returns:
    --------
    gmst_hours : float
        Greenwich Mean Sidereal Time in hours (0-24)
    gmst_radians : float
        GMST in radians
    """
    
    # Julian centuries since J2000.0 (2000 January 1.5 TT)
    T = (jd - 2451545.0) / 36525.0
    
    # GMST at 0h UT (in seconds)
    # Equation from USNO Circular 179 (simplified Aoki et al. 1982)
    GMST_0h = (67310.54841 +
               (876600.0 * 3600.0 + 8640184.812866) * T +
               0.093104 * T**2 -
               6.2e-6 * T**3)
    
    # Reduce to 0-86400 range
    GMST_0h = GMST_0h % 86400
    
    # Get fractional day from JD
    jd_fraction = (jd - np.floor(jd + 0.5)) + 0.5
    
    # Universal Time in seconds since 0h UT
    UT_seconds = jd_fraction * 86400.0
    
    # GMST including UT variation
    # Earth rotates 1.00273790935 times per mean solar day
    GMST_seconds = GMST_0h + 1.00273790935 * UT_seconds
    
    # Reduce to 24-hour format
    GMST_seconds = GMST_seconds % 86400
    
    # Convert to hours
    gmst_hours = GMST_seconds / 3600.0
    
    # Convert to radians (360° = 2π rad, so 1 hour = 2π/24 rad)
    gmst_radians = gmst_hours * (2 * np.pi / 24.0)
    
    return gmst_hours, gmst_radians


def local_sidereal_time(gmst_hours, longitude_deg):
    """
    Calculate Local Sidereal Time (LST) from Greenwich Mean Sidereal Time and longitude.
    
    Parameters:
    -----------
    gmst_hours : float
        Greenwich Mean Sidereal Time in hours
    longitude_deg : float
        Geographical longitude in degrees (East is positive, West is negative)
    
    Returns:
    --------
    lst_hours : float
        Local Sidereal Time in hours (0-24)
    lst_radians : float
        LST in radians
    """
    
    # Convert longitude from degrees to hours (15° per hour)
    longitude_hours = longitude_deg / 15.0
    
    # LST = GMST + longitude (in hours)
    lst_hours = gmst_hours + longitude_hours
    
    # Ensure LST is in 0-24 hour range
    lst_hours = lst_hours % 24.0
    
    # Convert to radians
    lst_radians = lst_hours * (2 * np.pi / 24.0)
    
    return lst_hours, lst_radians


def hms_to_decimal(hours, minutes, seconds):
    """
    Convert hours, minutes, seconds to decimal hours.
    """
    return hours + minutes / 60.0 + seconds / 3600.0


def decimal_to_hms(decimal_hours):
    """
    Convert decimal hours to hours, minutes, seconds format.
    
    Returns:
    --------
    (hours, minutes, seconds) : tuple of (int, int, float)
    """
    hours = int(decimal_hours)
    remainder = (decimal_hours - hours) * 60
    minutes = int(remainder)
    seconds = (remainder - minutes) * 60
    
    return hours, minutes, seconds


def decimal_to_dms(decimal_degrees):
    """
    Convert decimal degrees to degrees, arcminutes, arcseconds format.
    
    Returns:
    --------
    (degrees, arcminutes, arcseconds) : tuple of (int, int, float)
    """
    degrees = int(decimal_degrees)
    remainder = abs(decimal_degrees - degrees) * 60
    arcminutes = int(remainder)
    arcseconds = (remainder - arcminutes) * 60
    
    return degrees, arcminutes, arcseconds


# ============================================================================
# MAIN CALCULATION
# ============================================================================

if __name__ == "__main__":
    
    # Input parameters
    print("=" * 70)
    print("Julian Day and Local Sidereal Time Calculator")
    print("=" * 70)
    print()
    
    # Example: June 27, 2026, 12:14:45, 30° E
    year = 2026
    month = 6
    day = 27
    hour = 12
    minute = 14
    second = 45
    longitude_deg = 30.0  # East is positive
    
    print(f"Date/Time: {year}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:{second:02d} UTC")
    print(f"Longitude: {longitude_deg}° (East)")
    print()
    
    # Step 1: Calculate Julian Day Number
    jd, jd_integer = gregorian_to_julian_day(year, month, day, hour, minute, second)
    
    print("-" * 70)
    print("STEP 1: JULIAN DAY NUMBER (JD)")
    print("-" * 70)
    print(f"JD (full):    {jd:.10f}")
    print(f"JD (integer): {jd_integer:.1f}")
    print()
    
    # Step 2: Calculate Greenwich Mean Sidereal Time
    gmst_hours, gmst_radians = julian_day_to_greenwich_mean_sidereal_time(jd)
    gmst_h, gmst_m, gmst_s = decimal_to_hms(gmst_hours)
    
    print("-" * 70)
    print("STEP 2: GREENWICH MEAN SIDEREAL TIME (GMST)")
    print("-" * 70)
    print(f"GMST (decimal): {gmst_hours:.10f} hours")
    print(f"GMST (HMS):     {gmst_h:02d}:{gmst_m:02d}:{gmst_s:06.3f}")
    print(f"GMST (radians): {gmst_radians:.10f} rad")
    print()
    
    # Step 3: Calculate Local Sidereal Time
    lst_hours, lst_radians = local_sidereal_time(gmst_hours, longitude_deg)
    lst_h, lst_m, lst_s = decimal_to_hms(lst_hours)
    
    print("-" * 70)
    print("STEP 3: LOCAL SIDEREAL TIME (LST)")
    print("-" * 70)
    print(f"LST (decimal): {lst_hours:.10f} hours")
    print(f"LST (HMS):     {lst_h:02d}:{lst_m:02d}:{lst_s:06.3f}")
    print(f"LST (radians): {lst_radians:.10f} rad")
    print()
    
    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Julian Day (JD)':<30} = {jd:.10f}")
    print(f"{'Greenwich Mean Sidereal Time':<30} = {gmst_h:02d}:{gmst_m:02d}:{gmst_s:06.3f}")
    print(f"{'Local Sidereal Time (LST)':<30} = {lst_h:02d}:{lst_m:02d}:{lst_s:06.3f}")
    print("=" * 70)
