Scripts used to produce the data and figures for the ARC3O Part 1
=================================================================

What's this?
------------

These are the scripts used to compute data and produce figures for the paper:

Burgard, C., Notz, D., Pedersen, L. T., and Tonboe, R. T.: The Arctic Ocean Observation Operator for 6.9 GHz (ARC3O) – Part 1: How to obtain sea ice brightness temperatures at 6.9 GHz from climate model output, The Cryosphere, 14, 2369–2386, https://doi.org/10.5194/tc-14-2369-2020, 2020.

Computing data
--------------

The reference sea-ice evolution was computed with the model SAMSIM, which can be downloaded on github
`here <https://github.com/pgriewank/SAMSIM>`_. The version used for the analysis was downloaded on March 10th, 2017.

The forcing data for SAMSIM was downloaded from the ERA-Interim dataset, with the script: :file:`data/download_ERA_forcing_data.py`.#

For the analysis of the simplification of temperature and salinity (Section 4), the following were used:
    * `run_simplifications.py </data/run_simplifications.py>`_ run_simplifications.py` (different options can be given in the file), using functions from :file:`data/simplification_functions.py`: produces dat-files as input for MEMLS for each timestep
    * :file:`/data/write_input_netcdf.py`: transformes the single .dat-files into one netcdf-file.
    * :file:`data/run_memls.py`, using :file:`data/memls_functions.py`: simulates the brightness temperatures with MEMLS and stores the results in .dat-files
    * :file:`data/write_output_netcdf.py`: writes the MEMLS output to netcdf

For the sensitivity study looking at the influence of the number of layers (Section 5), the following were used:
    * :file:`data/run_simplifications_layers.py` (can change the amount of layers as an option in the beginning), using functions from :file:`data/simplification_functions.py`: produces dat-files as input for MEMLS for each timestep and each interest in layering
    * :file:`data/write_input_netcdf_layers.py`: transformes the single .dat-files into one netcdf-file, additional sorting by layers
    * :file:`data/run_memls_layers.py`, using :file:`data/memls_functions.py`: simulates the brightness temperatures with MEMLS qnd stores the results in .dat-files
    * :file:`data/write_output_netcdf_layers.py`: writes the MEMLS output to netcdf


Producing figures
-----------------

The final processing and visualization was done using the following scripts:
    * Figure 2: :file:`figures/Figure2.ipynb`
    * Figure 3, 5, 6: :file:`Figure3_5_6.py`
    * Figure 4: :file:`Figure4.py`

Signed: Clara Burgard, 14.08.2020
