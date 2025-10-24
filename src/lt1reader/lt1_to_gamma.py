import xml.etree.ElementTree as ET
import numpy as np
import rasterio
from datetime import datetime
import math

C = 299792458.0

WGS84_SEMI_MAJOR_AXIS = 6378137.0 
WGS84_ECCENTRICITY_SQUARED = 0.0066943799901
WGS84_SEMI_MINOR_AXIS = WGS84_SEMI_MAJOR_AXIS * math.sqrt(1 - WGS84_ECCENTRICITY_SQUARED)

def parse_iso_to_sod(t):
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
        f.write('title: LT1B SLC Image\n')
        if meta['state_vecs']:
            first_time = meta['state_vecs'][0][0]
            date_part = first_time.split('T')[0]
            f.write(f'date: {date_part}\n')
        
        f.write('sensor:    LT1B STRIP1 HH\n')
        
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

def lt1b_to_gamma(meta_xml, tiff_path, out_par, out_slc):
    print('[1/3] parse meta.xml', flush=True)
    meta = parse_meta(meta_xml)
    print('[2/3] generate SLC', flush=True)
    write_slc(tiff_path, out_slc)
    print('[3/3] write .par', flush=True)
    write_par(meta, out_par)
    print('[OK] done', flush=True)

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--meta', required=True)
    ap.add_argument('--tiff', required=True)
    ap.add_argument('--out-par', required=True)
    ap.add_argument('--out-slc', required=True)
    args = ap.parse_args()
    lt1b_to_gamma(args.meta, args.tiff, args.out_par, args.out_slc)