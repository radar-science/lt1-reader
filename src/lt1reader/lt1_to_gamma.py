"""
Core functions for converting LT1 data to Gamma format
"""

import xml.etree.ElementTree as ET
import numpy as np
import rasterio
import math
from datetime import datetime

C = 299792458.0

WGS84_SEMI_MAJOR_AXIS = 6378137.0 
WGS84_ECCENTRICITY_SQUARED = 0.0066943799901
WGS84_SEMI_MINOR_AXIS = WGS84_SEMI_MAJOR_AXIS * math.sqrt(1 - WGS84_ECCENTRICITY_SQUARED)


def parse_iso_to_sod(t):
    """parse ISO time format to seconds of day"""
    _, time_part = t.split('T', 1)
    hh, mm, ss = time_part.split(':')
    sec = float(ss)
    return int(hh)*3600 + int(mm)*60 + sec


def get_text(root, path):
    el = root.find(path)
    return el.text.strip() if el is not None and el.text else None


def parse_meta(meta_xml):
    tree = ET.parse(meta_xml)
    root = tree.getroot()

    rows = int(get_text(root, 'productInfo/imageDataInfo/imageRaster/numberOfRows'))
    cols = int(get_text(root, 'productInfo/imageDataInfo/imageRaster/numberOfColumns'))

    rng_sp = float(get_text(root, 'productSpecific/complexImageInfo/projectedSpacingRange/slantRange'))
    azi_sp = float(get_text(root, 'productSpecific/complexImageInfo/projectedSpacingAzimuth'))

    t_start = get_text(root, 'productInfo/sceneInfo/start/timeUTC')
    t_stop  = get_text(root, 'productInfo/sceneInfo/stop/timeUTC')
    sod_start = parse_iso_to_sod(t_start)
    sod_stop  = parse_iso_to_sod(t_stop)
    sod_center = 0.5*(sod_start + sod_stop)

    inc = float(get_text(root, 'productInfo/sceneInfo/sceneCenterCoord/incidenceAngle'))
    heading = float(get_text(root, 'productInfo/sceneInfo/headingAngle'))

    center_lat = float(get_text(root, 'productInfo/sceneInfo/sceneCenterCoord/lat'))
    center_lon = float(get_text(root, 'productInfo/sceneInfo/sceneCenterCoord/lon'))

    rt_first = float(get_text(root, 'productInfo/sceneInfo/rangeTime/firstPixel'))
    rt_center = float(get_text(root, 'productInfo/sceneInfo/sceneCenterCoord/rangeTime'))
    rt_last  = float(get_text(root, 'productInfo/sceneInfo/rangeTime/lastPixel'))
    near_rg = 0.5*C*rt_first
    center_rg = 0.5*C*rt_center
    far_rg = 0.5*C*rt_last

    freq = float(get_text(root, 'instrument/radarParameters/centerFrequency'))
    rsf  = float(get_text(root, 'instrument/settings/RSF'))
    prf  = float(get_text(root, 'instrument/settings/settingRecord/PRF'))
    bw   = float(get_text(root, 'processing/processingParameter/rangeCompression/chirps/referenceChirp/pulseBandwidth'))
    az_bw= float(get_text(root, 'processing/processingParameter/totalProcessedAzimuthBandwidth'))

    orbit_header = root.find('platform/orbit/orbitHeader')
    svs = root.findall('platform/orbit/stateVec')
    state_vecs = []
    for sv in svs:
        t = get_text(sv, 'timeUTC')
        px = float(get_text(sv, 'posX')); py = float(get_text(sv, 'posY')); pz = float(get_text(sv, 'posZ'))
        vx = float(get_text(sv, 'velX')); vy = float(get_text(sv, 'velY')); vz = float(get_text(sv, 'velZ'))
        state_vecs.append((t, px, py, pz, vx, vy, vz))

    # calculate sar_to_earth_center
    sar_to_earth_center = WGS84_SEMI_MAJOR_AXIS
    if state_vecs:
        _, px, py, pz, _, _, _ = state_vecs[0]
        sar_to_earth_center = math.sqrt(px*px + py*py + pz*pz)
        print(f'[DEBUG] orbit height: {sar_to_earth_center - WGS84_SEMI_MAJOR_AXIS:.1f} meters', flush=True)

    return {
        'rows': rows, 'cols': cols,
        'rng_sp': rng_sp, 'azi_sp': azi_sp,
        'sod_start': sod_start, 'sod_center': sod_center, 'sod_stop': sod_stop,
        'inc': inc, 'heading': heading,
        'center_lat': center_lat, 'center_lon': center_lon,
        'near_rg': near_rg, 'center_rg': center_rg, 'far_rg': far_rg,
        'freq': freq, 'rsf': rsf, 'prf': prf, 'bw': bw, 'az_bw': az_bw,
        'state_vecs': state_vecs,
        'sar_to_earth_center': sar_to_earth_center,
    }


def write_par(meta, out_par):
    with open(out_par, 'w') as f:
        f.write('Gamma Interferometric SAR Processor (ISP) - Image Parameter File\n\n')
        
        # basic information
        f.write('title: LT1 SLC Image\n')
        if meta['state_vecs']:
            first_time = meta['state_vecs'][0][0]
            date_part = first_time.split('T')[0]
            f.write(f'date: {date_part}\n')
        
        f.write('sensor:    LT1 STRIP1 HH\n')
        
        # time information
        f.write(f'start_time: {meta["sod_start"]:26.6f}   s\n')
        f.write(f'center_time:{meta["sod_center"]:26.6f}   s\n')
        f.write(f'end_time:  {meta["sod_stop"]:26.6f}   s\n')
        
        # image basic information
        f.write('line_header_size:                           0\n')
        f.write(f'range_samples:{meta["cols"]:27d}\n')
        f.write(f'azimuth_lines:{meta["rows"]:28d}\n')
        f.write('range_looks:                                1\n')
        f.write('azimuth_looks:                              1\n')
        f.write('image_format:               FCOMPLEX\n')
        f.write('image_geometry:             SLANT_RANGE\n')
        
        # pixel spacing and scale factor
        f.write(f'range_pixel_spacing:{meta["rng_sp"]:20.6f}   m\n')
        f.write(f'azimuth_pixel_spacing:{meta["azi_sp"]:18.6f}   m\n')
        f.write('range_scale_factor:                         1.0\n')
        f.write('azimuth_scale_factor:                       1.0\n')
        
        # slant range information
        f.write(f'near_range_slc:{meta["near_rg"]:23.4f}   m\n')
        f.write(f'center_range_slc:{meta["center_rg"]:21.4f}   m\n')
        f.write(f'far_range_slc:{meta["far_rg"]:24.4f}   m\n')
        
        # geometric parameters
        f.write(f'incidence_angle:{meta["inc"]:24.4f}  degrees\n')
        f.write(f'center_latitude: {meta["center_lat"]:24.6f}  degrees\n')
        f.write(f'center_longitude: {meta["center_lon"]:23.6f}  degrees\n')
        f.write(f'azimuth_angle: {meta["heading"]:30.7f}  degrees\n')
        
        # radar parameters
        f.write('azimuth_deskew:          OFF\n')
        f.write(f'azimuth_line_time:{(1.0/meta["prf"]):18.7e}   s\n')
        f.write(f'radar_frequency:{meta["freq"]:20.7e}  Hz\n')
        f.write(f'adc_sampling_rate:{meta["rsf"]:16.7e}  Hz\n')
        f.write(f'chirp_bandwidth:{meta["bw"]:19.7e}  Hz\n')
        f.write(f'prf: {meta["prf"]:30.7f}  Hz\n')
        f.write(f'azimuth_proc_bandwidth:{meta["az_bw"]:11.5f}  Hz\n')
        f.write(f'heading:{meta["heading"]:30.7f}  degrees\n')
        
        # gain parameters
        f.write('receiver_gain:                             0.0  dB\n')
        f.write('calibration_gain:                          0.0  dB\n')
        
        # earth parameters
        f.write(f'sar_to_earth_center:                   {meta["sar_to_earth_center"]:20.1f}  m\n')
        f.write(f'earth_radius_below_sensor:             {WGS84_SEMI_MAJOR_AXIS:20.1f}  m\n')
        f.write(f'earth_semi_major_axis:                 {WGS84_SEMI_MAJOR_AXIS:20.1f}  m\n')
        f.write(f'earth_semi_minor_axis:                 {WGS84_SEMI_MINOR_AXIS:20.1f}  m\n')
        
        # doppler parameters
        f.write('doppler_polynomial:                       0.0  0.0  0.0  0.0  Hz Hz/m Hz/m^2 Hz/m^3\n')
        f.write('doppler_poly_dot:                         0.0  0.0  0.0  0.0  Hz Hz/m Hz/m^2 Hz/m^3\n')
        f.write('doppler_poly_ddot:                        0.0  0.0  0.0  0.0  Hz Hz/m Hz/m^2 Hz/m^3\n')
        
        # slant range polynomial
        f.write('first_slant_range_polynomial:             0.0  0.0  0.0  0.0  0.0  0.0  s m 1 m^-1 m^-2 m^-3\n')
        f.write('center_slant_range_polynomial:            0.0  0.0  0.0  0.0  0.0  0.0  s m 1 m^-1 m^-2 m^-3\n')
        f.write('last_slant_range_polynomial:              0.0  0.0  0.0  0.0  0.0  0.0  s m 1 m^-1 m^-2 m^-3\n')
        
        # orbit information
        if meta['state_vecs']:
            ref_time = parse_iso_to_sod(meta['state_vecs'][0][0])
            f.write(f'time_of_first_state_vector: {ref_time:20.6f}   s\n')
            if len(meta['state_vecs']) > 1:
                time_interval = parse_iso_to_sod(meta['state_vecs'][1][0]) - ref_time
                f.write(f'state_vector_interval: {time_interval:20.6f}   s\n')
            else:
                f.write('state_vector_interval:                       1.0   s\n')
            
            f.write(f'number_of_state_vectors:{len(meta["state_vecs"]):16d}\n')
            
            for i, sv in enumerate(meta['state_vecs'], 1):
                _, px, py, pz, vx, vy, vz = sv
                f.write(f'state_vector_position_{i}:  {px:11.4f}  {py:11.4f}  {pz:11.4f}   m   m   m\n')
                f.write(f'state_vector_velocity_{i}:    {vx:9.5f}    {vy:9.5f}    {vz:9.5f}   m/s m/s m/s\n')


def get_orbit_times_from_meta(meta_file):
    """get all orbit times from meta.xml file"""
    tree = ET.parse(meta_file)
    root = tree.getroot()
    
    orbit_times = []
    
    # get all orbit state vectors
    state_vecs = root.findall('platform/orbit/stateVec')
    if not state_vecs:
        raise ValueError("no orbit state vectors found")
    
    print(f"[ORBIT] found {len(state_vecs)} orbit state vectors")
    
    for i, sv in enumerate(state_vecs):
        time_elem = sv.find('timeUTC')
        if time_elem is not None:
            time_str = time_elem.text.strip()
            try:
                try:
                    dt = datetime.fromisoformat(time_str.replace('T', ' '))
                except ValueError:
                    date_part, time_part = time_str.split('T')
                    year, month, day = map(int, date_part.split('-'))
                    hh, mm, ss = time_part.split(':')
                    second = float(ss)
                    dt = datetime(year, month, day, int(hh), int(mm), int(second))
                
                target_date = f"{dt.year:04d}{dt.month:02d}{dt.day:02d}"
                target_time_sod = dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond / 1e6
                
                orbit_times.append({
                    'index': i + 1,
                    'date': target_date,
                    'time_sod': target_time_sod,
                    'time_str': time_str
                })
                
                print(f"[ORBIT] orbit vector {i+1}: {time_str} -> {target_date} {target_time_sod:.3f}s")
                
            except Exception as e:
                print(f"[ORBIT] warning: failed to parse orbit vector {i+1} time: {time_str} (error: {e})")
                continue
    
    return orbit_times


def find_matching_orbit(orbit_file, target_date, target_time_sod):
    """find the orbit vector with the closest time to the target time in the orbit file"""
    
    min_time_diff = float('inf')
    best_sv = None
    
    with open(orbit_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            try:
                # orbit file format: year month day hour minute second X Y Z VX VY VZ ...
                parts = line.split()
                if len(parts) >= 12:
                    # parse date and time
                    year = int(parts[0])
                    month = int(parts[1])
                    day = int(parts[2])
                    hour = int(parts[3])
                    minute = int(parts[4])
                    second = float(parts[5])
                    
                    current_date = f"{year:04d}{month:02d}{day:02d}"
                    
                    if current_date == target_date:
                        time_sod = hour * 3600 + minute * 60 + second
                        
                        time_diff = abs(time_sod - target_time_sod)
                        if time_diff < min_time_diff:
                            min_time_diff = time_diff
                            
                            pos_x = float(parts[6])
                            pos_y = float(parts[7]) 
                            pos_z = float(parts[8])

                            vel_x = float(parts[9])
                            vel_y = float(parts[10])
                            vel_z = float(parts[11])
                            
                            best_sv = {
                                'time': time_sod,
                                'pos_x': pos_x,
                                'pos_y': pos_y,
                                'pos_z': pos_z,
                                'vel_x': vel_x,
                                'vel_y': vel_y,
                                'vel_z': vel_z,
                                'time_diff': time_diff
                            }
                    
            except (ValueError, IndexError) as e:
                continue
    
    return best_sv


def update_par_file_with_precise_orbit(par_file, orbit_file, meta_file):
    """update orbit parameters in .par file, replace coarse orbit data with precise orbit data for the corresponding time"""
    
    print(f"[ORBIT] reading meta file: {meta_file}")
    orbit_times = get_orbit_times_from_meta(meta_file)
    
    print(f"[ORBIT] parsing precise orbit file: {orbit_file}")
    orbit_mappings = []
    
    for orbit_info in orbit_times:
        matching_orbit = find_matching_orbit(orbit_file, orbit_info['date'], orbit_info['time_sod'])
        if matching_orbit:
            orbit_mappings.append({
                'par_index': orbit_info['index'],
                'meta_time': orbit_info['time_str'],
                'orbit_data': matching_orbit
            })
            print(f"[ORBIT] matching successful: orbit vector {orbit_info['index']} -> precise orbit time {matching_orbit['time']:.3f}s")
        else:
            print(f"[ORBIT] warning: orbit vector {orbit_info['index']} no matching precise orbit data")
    
    if not orbit_mappings:
        raise ValueError("no matching precise orbit data found")
    
    with open(par_file, 'r') as f:
        lines = f.readlines()
    
    orbit_start = None
    orbit_end = None
    
    for i, line in enumerate(lines):
        if 'number_of_state_vectors:' in line:
            orbit_start = i
            break
    
    if orbit_start is None:
        raise ValueError("no orbit parameters found")
    
    # find the end position of the orbit parameters
    for i in range(orbit_start + 1, len(lines)):
        if not any(keyword in lines[i] for keyword in ['state_vector_position_', 'state_vector_velocity_']):
            orbit_end = i
            break
    else:
        orbit_end = len(lines)
    
    # calculate the number of original orbit vectors
    original_count = 0
    for i in range(orbit_start + 1, orbit_end):
        if 'state_vector_position_' in lines[i]:
            original_count += 1
    
    print(f"[ORBIT] original file contains {original_count} orbit vectors")
    
    # build new orbit parameters
    new_orbit_lines = []
    new_orbit_lines.append(f'number_of_state_vectors:{original_count:16d}\n')
    
    # find the precise orbit data for each orbit vector position
    for i in range(1, original_count + 1):
        # find the precise orbit data
        matching_orbit = None
        for mapping in orbit_mappings:
            if mapping['par_index'] == i:
                matching_orbit = mapping['orbit_data']
                break
        
        if matching_orbit:
            new_orbit_lines.append(f'state_vector_position_{i}:  {matching_orbit["pos_x"]:11.4f}  {matching_orbit["pos_y"]:11.4f}  {matching_orbit["pos_z"]:11.4f}   m   m   m\n')
            new_orbit_lines.append(f'state_vector_velocity_{i}:    {matching_orbit["vel_x"]:9.5f}    {matching_orbit["vel_y"]:9.5f}    {matching_orbit["vel_z"]:9.5f}   m/s m/s m/s\n')
            print(f"[ORBIT]   orbit vector {i}: using precise orbit data (time difference: {matching_orbit['time_diff']:.3f}s)")
        else:
            print(f"[ORBIT] orbit vector {i}: no matching precise orbit data, keep original value")
            pos_line = None
            vel_line = None
            for j in range(orbit_start + 1, orbit_end):
                if f'state_vector_position_{i}:' in lines[j]:
                    pos_line = lines[j]
                elif f'state_vector_velocity_{i}:' in lines[j]:
                    vel_line = lines[j]
            if pos_line and vel_line:
                new_orbit_lines.append(pos_line)
                new_orbit_lines.append(vel_line)
    
    # replace orbit parameters
    new_lines = lines[:orbit_start] + new_orbit_lines + lines[orbit_end:]
    
    # write new file (overwrite original)
    with open(par_file, 'w') as f:
        f.writelines(new_lines)
    
    print(f"[ORBIT] updated orbit parameters: {len(orbit_mappings)} precise orbit data replaced")


def write_slc(tiff_path, out_slc):
    # read GeoTIFF: support (1 complex band) or (2 bands: real/imag)
    with rasterio.open(tiff_path) as src:
        print(f'[TIFF] bands={src.count} size={src.height}x{src.width} dtype={src.dtypes[0]}', flush=True)
        
        if src.count == 2:
            real = src.read(1).astype(np.float32, copy=False)
            imag = src.read(2).astype(np.float32, copy=False)
            
            arr = np.empty((src.height, src.width), dtype=np.complex64)
            arr.real = real
            arr.imag = imag
            
            print(f'[DEBUG] Real range: {real.min():.2f} to {real.max():.2f}')
            print(f'[DEBUG] Imag range: {imag.min():.2f} to {imag.max():.2f}')
            print(f'[DEBUG] Amplitude range: {np.min(np.abs(arr)):.2f} to {np.max(np.abs(arr)):.2f}')
            
        elif src.count == 1:
            arr = src.read(1)
            if not np.iscomplexobj(arr):
                raise RuntimeError('detected 1 band but not complex, cannot determine if I/Q or amplitude/phase.')
        else:
            raise RuntimeError('expected 1 complex band or 2 I/Q bands, actual number of bands is %d' % src.count)
    
    # ensure data layout is correct (row-major)
    arr = np.ascontiguousarray(arr)
    
    # convert to Gamma FCOMPLEX (float32 real + float32 imag)
    arr = arr.astype(np.complex64, copy=False)
    
    # write using little-endian byte order (Gamma usually expects little-endian)
    arr.byteswap().tofile(out_slc)
    
    print(f'[SLC] write: shape={arr.shape} dtype={arr.dtype}', flush=True)


def lt1_to_gamma(meta_xml, tiff_path, out_par, out_slc, orbit_file=None):
    """
    Convert LT1 data to Gamma format
    
    Parameters:
    -----------
    meta_xml : str
        Path to meta.xml file
    tiff_path : str
        Path to input GeoTIFF file
    out_par : str
        Path to output .par file
    out_slc : str
        Path to output .slc file
    orbit_file : str, optional
        Path to precise orbit file (if provided, will update orbit parameters)
    """
    print('[1/3] parse meta.xml', flush=True)
    meta = parse_meta(meta_xml)
    print('[2/3] generate SLC', flush=True)
    write_slc(tiff_path, out_slc)
    print('[3/3] write .par', flush=True)
    write_par(meta, out_par)
    
    # if precise orbit file is provided, update orbit parameters
    if orbit_file:
        print('[4/4] update orbit with precise orbit data', flush=True)
        try:
            update_par_file_with_precise_orbit(out_par, orbit_file, meta_xml)
            print('[OK] orbit parameters updated', flush=True)
        except Exception as e:
            print(f'[WARNING] failed to update orbit parameters: {e}', flush=True)
            print('[OK] done (with orbit update warning)', flush=True)
    else:
        print('[OK] done', flush=True)

