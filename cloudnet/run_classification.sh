#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH -o /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.out 
#SBATCH -e /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.err
#SBATCH -t 24:00:00

# activate python environment
cd ~/cloudnetpy
source venv/bin/activate

# Define start and end dates
#start_date='230510'
#end_date='230612'
start_date='230611'
end_date='230612'

# Convert the start and end dates to date objects for comparison
#start_date_obj=$(date -d "$start_date" "+%y%m%d")
#end_date_obj=$(date -d "$end_date" "+%y%m%d")

# on mac
start_date_seconds=$(date -jf "%y%m%d" "$start_date" "+%s")
end_date_seconds=$(date -jf "%y%m%d" "$end_date" "+%s")

# Loop through date strings
current_date_seconds=$start_date_seconds
while [ "$current_date_seconds" -le "$end_date_seconds" ]; do
    # Convert current date object back to date string
    # current_date_string=$(date -d "$current_date_obj" "+%y%m%d")
    current_date_string=$(date -r "$current_date_seconds" "+%y%m%d")

    # Call Python function with the date string as an argument
    python generate_classification.py "$current_date_string" "/Users/heather/cloudnetpy_1.46.4/data/"
    echo $current_date_string

    # Increment the date for the next iteration
    current_date_seconds=$((current_date_seconds + 86400))  # 86400 seconds in a day
    #current_date_obj=$(date -d "$current_date_obj + 1 day" "+%y%m%d")
done