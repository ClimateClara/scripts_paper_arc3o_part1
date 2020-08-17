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

The forcing data for SAMSIM was downloaded from the ERA-Interim dataset, with the script: `download_ERA_forcing_data.py </scripts_simulation/download_ERA_forcing_data.py>`_.
This data was converted to the .txt-format necessary as input for SAMSIM, with the script: `convert_nctotxt.py </scripts_simulation/convert_nctotxt.py>`_

For the analysis of the simplification of temperature and salinity (Section 4), the following were used:
    * `run_simplifications.py </data/run_simplifications.py>`_ (different options can be given in the file), using functions from `simplification_functions.py </scripts_simulation/simplification_functions.py>`_: produces dat-files as input for MEMLS for each timestep
    * `write_input_netcdf.py </scripts_simulation/write_input_netcdf.py>`_: transformes the single .dat-files into one netcdf-file.
    * `run_memls.py </scripts_simulation/run_memls.py>`_, using `memls_functions.py </scripts_simulation/memls_functions.py>`_: simulates the brightness temperatures with MEMLS and stores the results in .dat-files
    * `write_output_netcdf.py </scripts_simulation/write_output_netcdf.py>`_: writes the MEMLS output to netcdf

For the sensitivity study looking at the influence of the number of layers (Section 5), the following were used:
    * `run_simplifications_layers.py </scripts_simulation/run_simplifications_layers.py>`_ (can change the amount of layers as an option in the beginning), using functions from `simplification_functions.py </scripts_simulation/simplification_functions.py>`_: produces dat-files as input for MEMLS for each timestep and each interest in layering
    * `write_input_netcdf_layers.py </scripts_simulation/write_input_netcdf_layers.py>`_: transformes the single .dat-files into one netcdf-file, additional sorting by layers
    * `run_memls_layers.py </scripts_simulation/run_memls_layers.py>`_, using `memls_functions.py </scripts_simulation/memls_functions.py>`_: simulates the brightness temperatures with MEMLS qnd stores the results in .dat-files
    * `write_output_netcdf_layers.py </scripts_simulation/write_output_netcdf_layers.py>`_: writes the MEMLS output to netcdf


Producing figures
-----------------

The final processing and visualization was done using the following scripts:
    * Figure 2: `Figure2.ipynb </scripts_figures/Figure2.ipynb>`_
    * Figure 3, 5, 6: `Figures_3_4_6_7.ipynb </scripts_figures/Figures_3_4_6_7.ipynb>`_
    * Figure 4: `Figure5.py </scripts_figures/Figure5.py>`_
    * Values in Table 3: `Table_layers.ipynb </scripts_figures/Table_layers.ipynb>`_

Signed: Clara Burgard, 17.08.2020
