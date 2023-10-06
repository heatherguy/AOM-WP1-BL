#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH -o /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.out 
#SBATCH -e /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.err
#SBATCH -t 24:00:00



# Check we're in the right directory
cd /gws/nopw/j04/ncas_radar_vol1/heather/AoM2023/AOM-WP1-BL/

# activate python environment
source ~/miniconda3_b/bin/activate
conda activate icecaps

# Generate netcdf file
netcdf_out='/gws/nopw/j04/ncas_radar_vol1/heather/AoM2023/ice-station-data/final_nc/'
in_loc='/gws/nopw/j04/ncas_radar_vol1/heather/AoM2023/ice-station-data/'

start_dat='202305160000'
stop_dat='202306120000'

#python parse_surface-met.py $in_loc $netcdf_out $start_dat $stop_dat

python parse_flux_estimates.py $in_loc $netcdf_out $start_dat $stop_dat 30
python parse_flux_estimates.py $in_loc $netcdf_out $start_dat $stop_dat 15


