#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH -o /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.out 
#SBATCH -e /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.err
#SBATCH -t 24:00:00
#SBATCH --mem=10000

# Check we're in the right directory
cd /gws/nopw/j04/ncas_radar_vol1/heather/AoM2023/AOM-WP1-BL/

# activate python environment
source ~/miniconda3_b/bin/activate
conda activate icecaps

in_loc='/gws/nopw/j04/ncas_radar_vol1/heather/AoM2023/ice-station-data/raw/'
metek_save='/gws/nopw/j04/ncas_radar_vol1/heather/AoM2023/ice-station-data/metek_processed/'
licor_save='/gws/nopw/j04/ncas_radar_vol1/heather/AoM2023/ice-station-data/licor_processed/'
all_start='20230510'
all_stop='20230612'

# Loop through each day (i.e. one file per day)
start_date=`date -d $all_start '+%Y%m%d'`
end_date=`date -d $all_stop '+%Y%m%d'`

while [ "$start_date" != "$end_date" ]; do
	echo $start_date
    start=$start_date
    stop=$(date -d "$start_date + 1 day" '+%Y%m%d')
    python process_metek.py $in_loc $start $stop $metek_save
    #python process_licor.py $in_loc $start $stop $licor_save
	start_date=$(date -d "$start_date + 1 day" '+%Y%m%d')
done
