#!/usr/bin/env python3
"""
Update orbit parameters in Gamma SLC parameter file
Usage: python fix_orbit.py --orbit <orbit file> --meta <meta.xml file> --par <original.par file> --out <new.par file>
"""

import argparse
import re
import xml.etree.ElementTree as ET
from datetime import datetime

def parse_iso_to_sod(t):
    """parse ISO time format to seconds of day"""
    _, time_part = t.split('T', 1)
    hh, mm, ss = time_part.split(':')
    sec = float(ss)
    return int(hh)*3600 + int(mm)*60 + sec

def get_orbit_times_from_meta(meta_file):
    """get all orbit times from meta.xml file"""
    tree = ET.parse(meta_file)
    root = tree.getroot()
    
    orbit_times = []
    
    # get all orbit state vectors
    state_vecs = root.findall('platform/orbit/stateVec')
    if not state_vecs:
        raise ValueError("no orbit state vectors found")
    
    print(f"found {len(state_vecs)} orbit state vectors")
    
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
                
                print(f"orbit vector {i+1}: {time_str} -> {target_date} {target_time_sod:.3f}s")
                
            except Exception as e:
                print(f"warning: failed to parse orbit vector {i+1} time: {time_str} (error: {e})")
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

def update_par_file(par_file, new_par_file, orbit_mappings):
    """update orbit parameters in .par file, replace coarse orbit data with precise orbit data for the corresponding time"""
    
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
    
    print(f"original file contains {original_count} orbit vectors")
    
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
            print(f"  orbit vector {i}: using precise orbit data (time difference: {matching_orbit['time_diff']:.3f}s)")
        else:
            print(f"orbit vector {i}: no matching precise orbit data, keep original value")
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
    
    # write new file
    with open(new_par_file, 'w') as f:
        f.writelines(new_lines)
    
    print(f"updated orbit parameters: {len(orbit_mappings)} precise orbit data replaced")
    print(f"new file saved as: {new_par_file}")

def main():
    parser = argparse.ArgumentParser(description='update orbit parameters in Gamma SLC parameter file')
    parser.add_argument('--orbit', required=True, help='precise orbit file path')
    parser.add_argument('--meta', required=True, help='meta.xml file path')
    parser.add_argument('--par', required=True, help='original.par file path')
    parser.add_argument('--out', required=True, help='new.par file path')
    
    args = parser.parse_args()
    
    try:
        print(f"reading meta file: {args.meta}")
        orbit_times = get_orbit_times_from_meta(args.meta)
        
        print(f"parsing precise orbit file: {args.orbit}")
        orbit_mappings = []
        
        for orbit_info in orbit_times:
            matching_orbit = find_matching_orbit(args.orbit, orbit_info['date'], orbit_info['time_sod'])
            if matching_orbit:
                orbit_mappings.append({
                    'par_index': orbit_info['index'],
                    'meta_time': orbit_info['time_str'],
                    'orbit_data': matching_orbit
                })
                print(f"matching successful: orbit vector {orbit_info['index']} -> precise orbit time {matching_orbit['time']:.3f}s")
            else:
                print(f"warning: orbit vector {orbit_info['index']} no matching precise orbit data")
        
        if not orbit_mappings:
            raise ValueError("no matching precise orbit data found")
        
        print(f"updating parameter file: {args.par}")
        update_par_file(args.par, args.out, orbit_mappings)
        
        print("orbit parameters updated")
        
    except Exception as e:
        print(f"error: {e}")
        return 1
    
    return 0

if __name__ == '__main__':
    exit(main())