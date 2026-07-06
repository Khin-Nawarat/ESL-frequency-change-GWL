## Author
Khin Nawarat (March 2026, IHE Delft)
## Code origin

Parts of this workflow adapt functions from the analysis of Tebaldi et al. (2021a), which examined the global warming levels (GWLs) at which present-day 100-year extreme sea level (ESL) events become annual events: https://github.com/DOE-ICoM/tebaldi-etal_2021_natclimchange

## Modifications implemented to produce the datasets

- Extension of the analysis to multiple ESL return periods (10-, 20-, 50-, 100-, and 200-year) using automated looping across RPs
- Removal of two GWL scenarios (2°C+ and 5°C+) from Tebaldi et al. (2021a), which represent alternative SEJ-based high-end ice-sheet melt assumptions
- Restructuring of outputs and generation of NetCDF files for dataset publication
- Additional scripts for technical validation of the datasets

## Reproducing the datasets

Create a working directory with any name of your choice. Within this directory, create three subdirectories:'Rcodes','Rdatasets' and 'Data'. The input data used in this workflow are identical to those used in Tebaldi et al. (2021a,c). Download the input dataset from: https://doi.org/10.5281/zenodo.5095675. After downloading and unzipping the dataset, you should see a directory named: 'tebaldi-etal_2021_natclimchange_data'. This directory contains three subdirectories: 'Kirezci', 'Vousdoukas', and 'Rasmussen'. Rename the directory 'tebaldi-etal_2021_natclimchange_data' to: ‘CSV’. Place this directory inside your working directory. Next, download the R scripts from this repository and place them in the subdirectory: 'Rcodes'.

## Configuration

In all scripts except allfunctions.R and convolve_alltheway_fordatapub_allRPs.R, modify the variable *filedir* so that it points to the path of your working directory.

## Running the workflow

Run the R scripts in the following order to reproduce the datasets.

1. To read and restructure input data, run Read_in_data.R. This script reads the CSV files and restructures them into R arrays. The CSV files contain ESL estimates from the three ESL datasets paired with two alternative sea-level rise projections, organized by projection decade and GWL.
2. To generate probabilistic projections, run probabilisticprojections_fordatapub_allRPs.R. This script applies the Fisher Information Matrix approach to the ESL parameter estimates and convolves samples from their distributions with samples from the sea-level rise projections across different RPs and decades.
3. To apply synthesis methods, run votingsystem_fordatapub_allRPs.R. This script computes the results using the voting system (VS) synthesis method. It also automatically calls function from: convolve_alltheway_fordatapub_allRPs.R to compute the full convolution (FC) synthesis results. The outputs from both approaches are collected for later restructuring into the final datasets.
4. To generate GWL datasets, run VS_FC_GWL_as_netcdf.R. This script saves the results for the GWL at which different ESL return periods become annual events, using both VS and FC approaches, as NetCDF files.
5. To generate timinig datasets, run FC_timing.R. This script extracts the first decade at which present-day ESL return levels become annual events. Then run: FC_timing_as_netcdf.R to save the timing results as NetCDF files.
6. To generate the return-level datasets, run return_levels_as_netcdf.R. This script saves extreme sea-level return levels for each return period under present-day conditions and for each global warming level (GWL) as NetCDF files.
7. To reproduce Figure 2 from the paper,illustrating the dataset structure and spatial coverage, run dataset_structure_and_scale.R.
8. For technical validation, run technical validation.R. This script checks the internal consistency of the datasets and reproduces Figure 3 from the paper.
9. To assess the stability of the Monte Carlo sample size used for uncertainty propagation in the 200-year ESL estimates, run mc_sample_size_stability.R. This script evaluates the sensitivity of the estimated 200-year ESL values to the number of Monte Carlo realizations (1,000, 5,000, and 10,000). To reproduce Figure 4 from the paper, which summarizes the stability analysis, run mc_convergence_plot.py.

## References

Tebaldi, C., et al. (2021a). Extreme sea levels at different global warming levels. Nature Climate Change, 11, 746–751. https://doi.org/10.1038/s41558-021-01127-1

Tebaldi, C., et al. (2021b). Supporting data for Tebaldi et al. 2021 – Nature Climate Change 𝐷𝑎𝑡𝑎𝑠𝑒𝑡. Zenodo. https://doi.org/10.5281/zenodo.5095675

Tebaldi, C., et al. (2021c). Supporting code for Tebaldi et al. 2021 – Nature Climate Change 𝐶𝑜𝑑𝑒 Zenodo. https://doi.org/10.5281/zenodo.7551951

