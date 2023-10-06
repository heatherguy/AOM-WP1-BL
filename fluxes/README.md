# AOM-WP1-BL

Data processing scripts for ARTofMELT 2023 WP1: Boundary Layer

## Flux processing code for near surface heat and momentum fluxes from instrumentation on the sea ice. 

Turbulent fluxes of sensible heat, latent heat, and momentum heat are calculated from high resolution measurements of 3D wind, temperature, and humidity flucutations using Eddy Covariance.

### The following corrections are applied: 

- De-spiking: Replace outliers with and 11-data point median filter. Where outliers are defined as data that are greater than 3 standard deviations away from the median. 
- Licor data are corrected for a known delay in data transmission of 200 ms. 
- Cross wind temperature correction: Correct contamination of sonic temperature measurement for lengthening of sonic path by sidewind.
- Tilt correction: Correct tilt by alignment with horizontal streamline over a single averaging period. 
- Sonic acoustic temperature water vapor correction. (T = Ts / (1 + 0.51 q))
- Linear interpolation over data gaps that are smaller than 60% of the averaging period.  
- Absolute wind direction correction for ice floe rotation. 

### The following QC tests are implemented: 

- Basic instrumentation qc (unrealistic values and instrument errors).
- QC correction for licor window contamination
- Stationarity testing
- Calculations of skew and kurtosis
- Calculation of flux development

### The following CF compliant netcdf files are generated:

flux-estimates: 
Flux calculation results and associated qc-flags. 

flux-components: 
- Raw high resolution measurements used for the eddy covariance calculations.
- Results of intermediate calculations. 
- Results of qc-calculations (skew, kurtosis, ect.). 
- Derived Monin-Obukhov parameters.

## Example usage: 

### Preprocess data

    in_loc='/path_to_data/'
    save='/path_to_output_directory/'
    start='202305100000'
    stop='202306120000'

    python process_metek.py $in_loc $start $stop $save
    python process_licor.py $in_loc $start $stop $save

## Calculate fluxes

    in_loc='/path_to_data/'
    save='/path_to_output_directory/'
    start='202305100000'
    stop='202306120000'
    avp=30 # Averaging period for eddy covariance flux calculation (minutes)

    python parse_flux_estimates.py $in_loc $out_loc $start $start $avp
