#!/usr/bin/env python3
"""
Command-line interface for LT1 to Gamma conversion

test data are in the folder /home/hucy/penguin/lt1togamma/test_data

"""

import argparse
import sys
import os

# Add parent directory to path to import src modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.lt1_to_gamma import lt1_to_gamma


def main():
    """Main entry point for lt1_to_gamma command"""
    parser = argparse.ArgumentParser(
        description='Convert LT1 data to Gamma format, optionally update with precise orbit',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic conversion without precise orbit
  lt1_to_gamma --meta meta.xml --tiff image.tiff --out-par output.par --out-slc output.slc
  
  # Conversion with precise orbit update
  lt1_to_gamma --meta meta.xml --tiff image.tiff --out-par output.par --out-slc output.slc --orbit precise_orbit.txt
        """
    )
    
    parser.add_argument('--meta', required=True,
        help='Path to meta.xml file')
    parser.add_argument('--tiff', required=True,
        help='Path to input GeoTIFF file')
    parser.add_argument('--out-par', required=True,
        help='Path to output .par file')
    parser.add_argument('--out-slc', required=True,
        help='Path to output .slc file')
    parser.add_argument('--orbit', required=False, default=None,
        help='Path to precise orbit file (optional, if provided will update orbit parameters)'
    )
    
    args = parser.parse_args()
    
    # Validate input files exist
    if not os.path.exists(args.meta):
        print(f"Error: meta.xml file not found: {args.meta}", file=sys.stderr)
        return 1
    
    if not os.path.exists(args.tiff):
        print(f"Error: GeoTIFF file not found: {args.tiff}", file=sys.stderr)
        return 1
    
    if args.orbit and not os.path.exists(args.orbit):
        print(f"Error: orbit file not found: {args.orbit}", file=sys.stderr)
        return 1
    
    try:
        lt1_to_gamma(
            args.meta,
            args.tiff,
            args.out_par,
            args.out_slc,
            orbit_file=args.orbit
        )
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

