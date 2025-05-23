#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH -o /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.out 
#SBATCH -e /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.err
#SBATCH -t 24:00:00
#SBATCH --mem=10000

# Check we're in the right directory
cd /gws/nopw/j04/ncas_radar_vol1/heather/AOM2023/AOM-WP1-BL/

# activate python environment
source ~/miniconda3_b/bin/activate
conda activate icecaps

# Generate netcdf file
netcdf_out='/gws/nopw/j04/ncas_radar_vol1/heather/AOM2023/ice-station-data/final_nc/'
in_loc='/gws/nopw/j04/ncas_radar_vol1/heather/AOM2023/ice-station-data/'

# Ice station 1
start_dat='202305160000'
stop_dat='202305220000'
#python parse_surface-met.py $in_loc $netcdf_out $start_dat $stop_dat

# Ice station 1
start_dat='202305290000'
stop_dat='202306120000'
#python parse_surface-met.py $in_loc $netcdf_out $start_dat $stop_dat

start_dat='202305160000'
stop_dat='202306120000'

python parse_flux_estimates.py $in_loc $netcdf_out $start_dat $stop_dat 30
python parse_flux_estimates.py $in_loc $netcdf_out $start_dat $stop_dat 20
python parse_flux_estimates.py $in_loc $netcdf_out $start_dat $stop_dat 10


