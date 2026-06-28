"""
Julian Day and Local Sidereal Time Calculator
Calculates Julian Day Number and Local Sidereal Time from date, time, and geographical position.
Also retrieves planetary data from NASA's JPL Horizons system.
"""

import sys
import numpy as np
from astroquery.jplhorizons import Horizons
from datetime import datetime
import requests

def gregorian_to_julian_day(year,month,day,hour,minute,second) :
    """
    Convert gregorian time to julian day
    ————————————————————————————————————
    Input
      > year    (1901-2099)
      > month   (1-12)
      > day     (1-31)
      > hour    (0–23)
      > minute  (0-59)
      > second  (0-59.999...)
    
    Output
      > UT :    Universal time (hours)
      > J0 :    Julian Day at 0h UT (days)
      > JD :    Julian Day at actual UT (days)
      
    """
    
    # Universal Time in hours
    UT = hour + (minute/60) + (second/3600)
    
    # JD Integers :
    a0 = int((month + 9) / 12)
    a  = int(7 * (year + a0) / 4)
    b  = int(275 * month / 9)
    
    # Julian Day number at 0h UT
    J0 = 367 * year - a + b + day + 1721013.5
    
    # Julian Day number at actual UT
    JD = J0 + UT / 24
    
    
    return UT, J0, JD
    
def greenwich_sidereal_time(UT, J0):
    """
    Evaluate greenwich sideral time from julian day number
    ––––––––––––––––––––––––––––––––––––––––––––––––––––––
    Input
      > UT, Universal Time (hours)
      > J0, Julian Day at 0h UT (days)
    
    Output
      > GST, Greenwich Sideral Time (degrees)
      > GST_hours, Greenwich Sideral Time (hours)
    """
    
    # Time in Julian centuries between Julian day @ 0h UT and J2000
    T0 = (J0 - 2451545.0) / 36525.0
    
    # Greenwich sidereal time at 0h UT (in degrees)
    GST0 = 100.4606184 + 36000.77004*T0 + 0.000387933*power_2(T0) - (2.583E-8)*power_3(T0)
    
    # Correction if the value is outside the range [0, 360] degrees
    if GST0 >= 360 :
        GST0 -= int(GST0/360)*360
        
    if GST0 <= 0 :
        GST0 += int(GST0/360)*360
    
    # Greenwich sideral time at current UT
    GST = GST0 + 360.98564724 * UT / 24
    
    # Correction if higher than 360
    if GST >= 360 :
        GST -= int(GST/360)*360
    
    # Convert in hours
    GST_hours = GST * 24 / 360
    
    return GST, GST_hours
    
def local_sidereal_time(GST, longitude) :
    """
    Evaluate Local Sideral Time from Greenwich Sidereal Time
    ––––––––––––––––––––––––––––––––––––––––––––––––––––––––
    Input
      > GST, Greenwich Sideral Time (degrees)
      > longitude_deg, global longitude of the site (degrees)
    
    Output
      > LST, Local Sidereal Time (degrees)
      > LST_hours, Local Sidereal Time (hours)
    """
    
    # Local Sidereal Time @ site
    LST = GST + longitude
    
    # Correction if exceeds 360
    if LST >= 360 :
        LST -= int(LST/360)*360
    
    # Convert in hours
    LST_hours = LST * 24 / 360
    
    return LST, LST_hours
    
def deg_to_hms(a):
    """
    Convert the numer a from degrees in hour-minutes-second units
    –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
    Output
      > h, hours    (0-23)
      > m, minutes  (0-59)
      > s, seconds  (0-59.999...)
    """
    
    # Extract hours
    h = int(a * 24 / 360)
    
    # Extract minutes
    m = int((a * 24 / 360 - h) * 60)
    
    # Extract seconds
    s = ((a * 24 / 360 - h) * 60 - m) * 60
    
    return h, m, s
    
def earth_radii_on_site(latitude, longitude, elevation) :
    """
    Estimation of the radius of the Earth at the side
    –––––––––––––––––––––––––––––––––––––––––––––––––
    Input:
      > latidude  (degrees)
      > elevation (metres)
      
    Output:
      > R, magnitude of earth radii on site (km)
      > R_vector, vector components
    """
    
    # Costant parameters
    Re = 6378       # Equatorial radius of the Earth (km)
    f  = 0.003353   # Oblateness or flattening of the Earth
    
    # Conversion
    rad = np.pi/180
    
    # Normal radius at the zenith direction of the site
    R_phi = 6378 / np.sqrt(1 - (2*f - power_2(f)) * power_2( np.sin(latitude * rad) ))
    
    # Adding the elevation
    Rc = R_phi + elevation/1000
    Rs = elevation/1000 + R_phi * power_2(1 - f)
    
    # Vector components of the Earth radius at site
    x = Rc * np.cos(latitude * rad) * np.cos(longitude * rad)
    y = Rc * np.cos(latitude * rad) * np.sin(longitude * rad)
    z = Rs * np.sin(latitude * rad)
    
    # Earth radius at site
    R = np.sqrt(power_2(x) + power_2(y) + power_2(z))
    
    return R, x, y, z

def parse_arguments():
    """
    Read command line or default
    ––––––––––––––––––––––––––––
    """
    
    # Se è passato "today here", usa data di oggi e posizione rilevata
    if len(sys.argv) >= 3 and sys.argv[1].lower() == "today" and sys.argv[2].lower() == "here":
        adesso = datetime.now()
        year = adesso.year
        month = adesso.month
        day = adesso.day
        hour = adesso.hour
        minute = adesso.minute
        second = adesso.second
        
        # Rileva posizione dal GPS o IP
        latitude, longitude, city, country = get_location_from_ip()
        elevation = 0
        if len(sys.argv) >= 4:
            elevation = float(sys.argv[3])
        else:
            elevation = 0
    
    # Se è passato solo "today"
    elif len(sys.argv) >= 2 and sys.argv[1].lower() == "today":
        adesso = datetime.now()
        year = adesso.year
        month = adesso.month
        day = adesso.day
        hour = adesso.hour
        minute = adesso.minute
        second = adesso.second
        
        # Longitudine (se non passata, usa default)
        if len(sys.argv) >= 4:
            site_str = sys.argv[2]
            elevation = float(sys.argv[3])
            latitude, longitude = map(float, site_str.split('-'))
            city, country = "Custom", "Location"
        else:
            latitude, longitude = 45.7195, 11.3428
            city, country = "Schio", "Italy"
            elevation = 0
    
    # Altrimenti, leggi i parametri normali
    elif len(sys.argv) >= 5:
        data_str = sys.argv[1]
        ora_str = sys.argv[2]
        site_str = sys.argv[3]
        elevation = float(sys.argv[4])
        
        year, month, day = map(int, data_str.split('-'))
        hour, minute, second = map(float, ora_str.split(':'))
        hour = int(hour)
        minute = int(minute)
        second = int(second)
        latitude, longitude = map(float, site_str.split('-'))
        city, country = "Custom", "Location"
    
    else:
        # Default values
        year, month, day = 2026, 6, 27
        hour, minute, second = 14, 54, 29
        latitude, longitude = 45.7195, 11.3428
        city, country = "Schio", "Italy"
        elevation = 200
    
    return year, month, day, hour, minute, second, latitude, longitude, elevation, city, country
    
def power_2(a):
    """
    Elevate the number a to elevation 2
    –––––––––––––––––––––––––––––––––––
    """
    b = a*a
    return b
    
def power_3(a):
    """
    Elevate the number a to elevation 3
    –––––––––––––––––––––––––––––––––––
    """
    b = a*a*a
    return b
    
def get_location_from_ip():
    """
    Rileva latitudine e longitudine dal tuo IP
    ––––––––––––––––––––––––––––––––––––––––––
    """
    try:
        response = requests.get('https://ipapi.co/json/', timeout=5)
        data = response.json()
        
        latitude = float(data['latitude'])
        longitude = float(data['longitude'])
        city = data.get('city', 'Unknown')
        country = data.get('country_name', 'Unknown')
        
        print(f"📍 Location detected from IP: {city}, {country}")
        return latitude, longitude, city, country
    
    except Exception as e:
        print(f"⚠️  Could not detect location: {e}")
        print("Using default location")
        return 45.7195, 11.3428, "Schio", "Italy"

def get_planetary_data_from_nasa(latitude, longitude, date_time):
    """
    Get planetary data from NASA's JPL Horizons system
    ––––––––––––––––––––––––––––––––––––––––––––––––––
    
    Input
      > latitude, longitude (degrees)
      > date_time (datetime object)
    
    Output
      > planetary_data (dictionary with position/velocity vectors)
    """
    try:
        # Coordinate del sito
        location = f"@{latitude},{longitude}"
        
        # Pianeti (ID Horizons)
        planets = {
            'Mercury': '199',
            'Venus': '299',
            'Mars': '499',
            'Jupiter': '599',
            'Saturn': '699',
            'Uranus': '799',
            'Neptune': '899'
        }
        
        planetary_data = {}
        
        print(f"\n🌐 Fetching planetary data from NASA JPL Horizons...")
        print(f"   (This may take a moment...)")
        
        for planet_name, planet_id in planets.items():
            print(f"   • {planet_name:10}...", end="", flush=True)
            
            try:
                # Convert datetime to JD for astroquery
                jd = date_time.toordinal() + 1721425.5 + (date_time.hour + date_time.minute/60 + date_time.second/3600) / 24
                
                obj = Horizons(id=planet_id, location=location, epochs=jd)
                v = obj.vectors()
                
                planetary_data[planet_name] = {
                    'x': float(v['x'][0]),
                    'y': float(v['y'][0]),
                    'z': float(v['z'][0]),
                    'vx': float(v['vx'][0]),
                    'vy': float(v['vy'][0]),
                    'vz': float(v['vz'][0]),
                    'distance': float(np.sqrt(v['x'][0]**2 + v['y'][0]**2 + v['z'][0]**2))
                }
                print(" ✓")
                
            except Exception as e:
                print(f" ✗ ({str(e)[:30]}...)")
                continue
        
        return planetary_data if planetary_data else None
    
    except Exception as e:
        print(f"\n❌ Errore nel caricamento dati planetari: {e}")
        return None

def display_planetary_data(planetary_data):
    """
    Display planetary data in a formatted table
    """
    if not planetary_data:
        print("\n⚠️  No planetary data available")
        return
    
    print()
    print("=" * 100)
    print("PLANETARY DATA (NASA JPL Horizons)")
    print("=" * 100)
    print()
    print(f"{'Planet':<12} | {'X (AU)':<12} | {'Y (AU)':<12} | {'Z (AU)':<12} | {'Distance (AU)':<14}")
    print("-" * 100)
    
    for planet, data in planetary_data.items():
        print(f"{planet:<12} | {data['x']:>11.6f} | {data['y']:>11.6f} | {data['z']:>11.6f} | {data['distance']:>12.6f}")
    
    print()
    print(f"{'Planet':<12} | {'VX (AU/day)':<14} | {'VY (AU/day)':<14} | {'VZ (AU/day)':<14}")
    print("-" * 100)
    
    for planet, data in planetary_data.items():
        print(f"{planet:<12} | {data['vx']:>13.8f} | {data['vy']:>13.8f} | {data['vz']:>13.8f}")
    
    print()
    print("=" * 100)
    print()
    
# ============================================================================
# MAIN CALCULATION
# ============================================================================

# Title
print()
print("=" * 70)
print("Julian Day and Local Sidereal Time Calculator")
print("=" * 70)
print()

# Recall input parameters
year, month, day, hour, minute, second, latitude, longitude, elevation, city, country = parse_arguments()

print(f"Date/Time: {year}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:{second:02d} UTC")
print(f"Location:  {latitude:.4f}° N, {longitude:.4f}° E –– {city}, {country}")
print(f"Elevation: {elevation:.1f} m")
    
# Step 1: Calculate Universal Time and Julian Day Numbers
UT, J0, JD = gregorian_to_julian_day(year, month, day, hour, minute, second)
    
# Step 2: Calculate Greenwick Sidereal Time
GST, GST_hours = greenwich_sidereal_time(UT, J0)
GST_h, GST_m, GST_s = deg_to_hms(GST)

# Step 3: Calculate Local Sidereal Time
LST, LST_hours = local_sidereal_time(GST, longitude)
LST_h, LST_m, LST_s = deg_to_hms(LST)

# Step 4: Calculate Earth radius at site
R, x, y, z = earth_radii_on_site(latitude, longitude, elevation)

print()
print("=" * 70)
print()
print(f"Universal Time (UT) ——————> {UT:.4f} hours")
print(f"Julian Day @ 0h UT –––––––> {J0:.1f} days")
print(f"Julian Day @ current UT ––> {JD:.10f} days")
print(f"Greenwich Sidereal Time ––> {GST:.4f}° / {GST_h:02d}h {GST_m:02d}m {GST_s:06.3f}s")
print(f"Local Sidereal Time ––––––> {LST:.4f}° / {LST_h:02d}h {LST_m:02d}m {LST_s:06.3f}s")
print()
print("=" * 70)
print()
print("Earth Radius at site")
print()
print(f"Module ————————————————––> {R:.4f} km")
print(f"Vector ————————————————––> ({x:.4f}, {y:.4f}, {z:.4f}) [km]")

# Step 5: Get planetary data from NASA
date_time = datetime(year, month, day, hour, minute, second)
planetary_data = get_planetary_data_from_nasa(latitude, longitude, date_time)

# Step 6: Display planetary data
if planetary_data:
    display_planetary_data(planetary_data)
else:
    print("\n⚠️  Could not retrieve planetary data")

print()
