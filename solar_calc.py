import numpy as np
import pandas as pd

def declination_angle(n):
    return 23.45 * np.sin(np.deg2rad(360*((284+n)/365)))

def celcius2kelvin(c):
    return c + 273

def kelvin2celcius(k):
    return k - 273

def celcius2fahrenheit(c):
    return c * 1.8 + 32

def t_sky_berdahl(t_db, t_dp, cloud_cov):
    t_db = celcius2kelvin(c=t_db)
    cloud_cov = float(cloud_cov)/10.0
    eff_cloud_cov = 1 + (0.0224 * cloud_cov) + (0.0035 * (cloud_cov**2)) + (0.00028 * (cloud_cov**3))
    clear_sky_emm = 0.711 + (0.56 * (t_dp/100)) + (0.73 * (t_dp/100)**2)
    t_clear_sky = t_db * clear_sky_emm**0.25
    t_sky = (eff_cloud_cov**0.25) * t_clear_sky
    return kelvin2celcius(k=t_sky)

def hour_angle(month_loc, hour_loc, long_loc, long_std):
    et_list = {
        1:-10.6,
        2:-14.0,
        3:-7.9,
        4:1.2,
        5:3.7,
        6:-1.3,
        7:-6.4,
        8:-3.6,
        9:6.9,
        10:15.5,
        11:13.8,
        12:3.2
    }
    # Not very accurate

    ast = hour_loc + (et_list[month_loc]/60) + (long_loc-long_std)/15
    return 15 * (ast - 12)

def zenith_angle(declination, latitude, hour_angle):
    declination = np.deg2rad(declination)
    latitude = np.deg2rad(latitude)
    hour_angle = np.deg2rad(hour_angle)
    zenith_angle = np.arccos((np.cos(latitude) * np.cos(declination) * np.cos(hour_angle)) + (np.sin(latitude) * np.sin(declination)))
    return np.rad2deg(zenith_angle)

def diff_sky_rad(glob_hor_rad, dir_nor_rad, zenith_angle, slope):
    zenith_angle = np.deg2rad(zenith_angle)
    slope = np.deg2rad(slope)
    return (glob_hor_rad - (dir_nor_rad * np.cos(zenith_angle))) * (1 + np.cos(slope))/2

def diff_refl_rad(glob_hor_rad, gr_refl, slope):
    slope = np.deg2rad(slope)
    return glob_hor_rad * gr_refl * (1 - np.cos(slope))/2

def get_surf_azimuth(orientation):
    az_dict = {'north':180, 'south':0, 'east':-90, 'west':90, 'horizontal':0}
    return az_dict[orientation]

def incidence_angle(orientation, declination, latitude, slope, hour_angle):
    declination = np.deg2rad(declination)
    latitude = np.deg2rad(latitude)
    slope = np.deg2rad(slope)
    surf_azimuth = np.deg2rad(get_surf_azimuth(orientation=orientation))
    hour_angle = np.deg2rad(hour_angle)
    incidence_angle = (np.sin(declination) * np.sin(latitude) * np.cos(slope)) - (np.sin(declination) * np.cos(latitude) * np.sin(slope) * np.cos(surf_azimuth)) + (np.cos(declination) * np.cos(latitude) * np.cos(slope) * np.cos(hour_angle)) + (np.cos(declination) * np.sin(latitude) * np.sin(slope) * np.cos(surf_azimuth) * np.cos(hour_angle)) + (np.cos(declination) * np.sin(slope) * np.sin(surf_azimuth) * np.sin(hour_angle))
    return np.rad2deg(np.arccos(incidence_angle))

def loc_vel(w_spd,w_dir,orientation):
    
    if orientation == 'south':
        deg = 180
        windward = w_dir <= deg+90 and w_dir >= deg-90
    elif orientation == 'east':
        deg = 90
        windward = w_dir <= deg+90 and w_dir >= deg-90
    elif orientation == 'west':
        deg = 270
        windward = w_dir <= deg+90 and w_dir >= deg-90
    elif orientation == 'north':
        deg = 0
        if w_dir == 0 or w_dir == 360:
            windward = True
        else:
            windward = w_dir <= 90 or w_dir >= 270
        #north
    elif orientation == 'horizontal':
        windward = True
    else:
        return None
    
    if windward:
        if w_spd < 2:
            return 0.5
        else:
            return w_spd * 0.25
    else:
        return w_spd * 0.25

def conv_coef(w_spd,w_dir,orientation):
    return 3.5 + (5.6 * loc_vel(w_spd,w_dir,orientation))

def dir_rad(dir_nor_rad, incidence_angle):
    dir_rad = dir_nor_rad * np.cos(np.deg2rad(incidence_angle))
    if dir_rad > 0:
        return dir_rad
    else:
        return 0

def get_df_time():
    day_dict = {
        1:31,
        2:28,
        3:31,
        4:30,
        5:31,
        6:30,
        7:31,
        8:31,
        9:30,
        10:31,
        11:30,
        12:31,
    }
    month_list = []
    tot_days_list = []
    day_list = []
    hour_list = []
    day_n = 1
    for i in range(1,13):
        count = day_dict[i]
        for j in range(1,count+1):
            for k in range(1,25):
                month_list.append(i)
                tot_days_list.append(day_n)
                day_list.append(j)
                hour_list.append(k)
            day_n += 1       
    df_time = pd.DataFrame({
        'DAY_N':tot_days_list,
        'MONTH': month_list,
        'DAY': day_list,
        'HOUR': hour_list
    })
    return df_time

def get_rad_df(tmy_file, latitude, longitude, std_meridian):
    weather_data = pd.read_csv(tmy_file, skiprows=1)
    df_weather = pd.DataFrame(weather_data)
    cols = [
           'GHI (W/m^2)', 'DNI (W/m^2)', 'DHI (W/m^2)',
            'Dry-bulb (C)', 'Dew-point (C)', 'TotCld (tenths)','Wdir (degrees)','Wspd (m/s)'
    ]
    df_weather = df_weather[cols]
    df_weather.columns = ['GHI','DNI','DHI','DBT','DPT','CC','WDIR','WSPD']
    df_time = get_df_time()
    
    df_list = []
    orientation_list = ['horizontal','north','south','east','west']
    slope_list = [0,90,90,90,90]
    
    for i in range(0,len(orientation_list)):
        df = None
        df = df_time.join(df_weather)
        orientation = orientation_list[i]
        slope = slope_list[i]

        df['SKY_TEMP'] = df.apply(lambda x: t_sky_berdahl(x['DBT'], x['DPT'], x['CC']), axis=1)
        df['CONV_COEF'] = df.apply(lambda x: conv_coef(x['WSPD'],x['WDIR'],orientation), axis=1)
        df['DECLINATION'] = df.apply(
        lambda x: 
        declination_angle(n=x['DAY_N']),
        axis=1
        )
        df['HOUR_ANGLE'] = df.apply(
            lambda x:
            hour_angle(
                month_loc=x['MONTH'],
                hour_loc=x['HOUR'],
                long_loc=longitude,
                long_std=std_meridian
            ),
            axis=1
        )
        df['ZENITH_ANGLE'] = df.apply(
            lambda x: 
            zenith_angle
            (declination=x['DECLINATION'],
             latitude=latitude,
             hour_angle=x['HOUR_ANGLE']
            ),
            axis=1
        )
        df['ICD_ANG'] = df.apply(
            lambda x:
            incidence_angle(
                orientation=orientation,
                declination=x['DECLINATION'],
                latitude=latitude,
                slope=slope,
                hour_angle=x['HOUR_ANGLE']
            ),
            axis=1
        )
        df['DSI'] = df.apply(
            lambda x:
            diff_sky_rad(
                glob_hor_rad = x['GHI'],
                dir_nor_rad = x['DNI'],
                zenith_angle = x['ZENITH_ANGLE'],
                slope = slope
            ),
            axis=1
        )
        df['DRI'] = df.apply(
            lambda x:
            diff_refl_rad(
                glob_hor_rad=x['GHI'],
                gr_refl=0.2,
                slope=slope
            ),
            axis=1
        )
        df['DI'] = df.apply(
            lambda x:
            dir_rad(
                dir_nor_rad=x['DNI'],
                incidence_angle=x['ICD_ANG']
            ),
            axis=1
        )
        df['TOT_I'] = df.apply(lambda x: x['DSI'] + x['DRI'] + x['DI'], axis=1)
        # df['SOL_AIR_T_LT'] = df.apply(
        #     lambda x: sol_air_t(t_db = x['DBT'], tot_rad=x['TOT_I'], orientation=orientation, color='light'), axis=1
        # )
        # df['SOL_AIR_T_DK'] = df.apply(
        #     lambda x: sol_air_t(t_db = x['DBT'], tot_rad=x['TOT_I'], orientation=orientation, color='dark'), axis=1
        # )
        
        df_list.append(df)
    
    return df_list