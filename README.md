# AOM-WP1-BL

Data processing scripts for ARTofMELT 2023 WP1: Boundary Layer

### 1. Ice station surface meteorology data. 

Quality controlled, 1-minute averages of air temperature, relative humidity, wind speed and direction, SW and LW radiative fluxes, snow-to-ice heat flux, snow surface temperature, and 64 cm snow temperature profile from the two on-ice stations during ARTofMELT. 

#### Notes: 
- The radiometer stand (all hemispheric radiosometers and the KT15-2) was located approximately 40 m away from the met mast (all other instruments).
- True wind directions are calculated using the heading between the two GPS units on the Met Mast and the radiometer stand. 
- Vector wind averages are calculated before derving the True wind direction in degrees.
- No tilt correction is applied to the 3D sonic anemometer data. Horizontal wind speed is calculated from all 3D componenets, assuming that the contribution of vertical velocity is negligible.  
- The latitude and longitude data from the GPS units experienced several unexplained and unrealistic step changes. A 'best guess' correction has been applied to remove these step changes based on information from the Sip GPS data, however the latitude and longitude in this file should be considered to have an uncertainty of +/- 0.4 degrees. This impacted both GPS units and therefore did not impact the True wind direction calculation. 
- Air pressure in this file is derived from the internal sensor in the Licor control box, no corrections have been applied. 
- The Met mast KT15 (skin_temperature_1) had an window heater activated. The radiometer stand KT15 (skin_temperature_2) did not.
- The field log-book for both ice stations is attached: ARTofMELT23_icestation_logbook.xls
- The factory calibration is used for the heat-flux plate data. 
- 1-minute averages are labelled for the LHS of the averaging window. 

### 2. Flux processing code for near surface heat and momentum fluxes from instrumentation on the sea ice. 

Turbulent fluxes of sensible heat, latent heat, and momentum heat are calculated from 20 Hz resolution measurements of 3D wind, virtual temperature, and humidity flucutations using Eddy Covariance. Winds are rotated into the streamline prior to flux calculation. 

#### The following corrections are applied: 

- De-spiking: Replace outliers with and 11-data point median filter. Where outliers are defined as data that are greater than 3 standard deviations away from the median. 
- Licor data are corrected for a known delay in data transmission of 200 ms. 
- Cross wind temperature correction: Correct contamination of sonic temperature measurement for lengthening of sonic path by sidewind.
- Tilt correction: Correct tilt by alignment with horizontal streamline over a single averaging period. 
- Sonic acoustic temperature water vapor correction. (T = Ts / (1 + 0.51 q))
- Linear interpolation over data gaps that are smaller than 60% of the averaging period.  
- Absolute wind direction correction for ice floe rotation. 

#### The following QC tests are implemented: 

- Basic instrumentation qc (unrealistic values and instrument errors).
- QC correction for licor window contamination
- Stationarity testing
- Calculations of skew and kurtosis
- Calculation of flux development

#### The following CF-compliant netcdf files are generated:

flux-estimates: 
Flux calculation results and associated qc-flags. 

flux-components: 
- Raw high resolution measurements used for the eddy covariance calculations.
- Results of intermediate calculations. 
- Results of qc-calculations (skew, kurtosis, ect.). 
- Derived Monin-Obukhov parameters.

### Example usage: 

#### Preprocess data

    in_loc='/path_to_data/'
    save='/path_to_output_directory/'
    start='202305100000'
    stop='202306120000'

    python process_metek.py $in_loc $start $stop $save
    python process_licor.py $in_loc $start $stop $save

#### Generate surface met files

    python parse_surface-met.py $in_loc $out_loc $start $stop

#### Calculate fluxes

    start='202305100000'
    stop='202306120000'
    avp=30 # Averaging period for eddy covariance flux calculation (minutes)
    python parse_flux_estimates.py $in_loc $out_loc $start $start $avp
