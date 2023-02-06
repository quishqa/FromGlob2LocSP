# FromGlob2LocSP
Scripts to reproduce figures and results of ["From Global to local: A multi-scale air quality modeling study over the Metropolitan Area of São Paulo"](https://www.teses.usp.br/teses/disponiveis/14/14133/tde-15042021-115527/fr.php)

* [`01_Ibirapuera_20yrs_data`](https://github.com/quishqa/FromGlob2LocSP/tree/main/01_Ibirapuera_20yrs_data): Example on how to use [`qualR`](https://github.com/quishqa/qualR),
and how to make the concentration trend plot Figure in Introduction section.
* [`02_CETESB_data`](https://github.com/quishqa/FromGlob2LocSP/tree/main/02_CETESB_data): We used `qualar_py` from [`wrf_sp_eval`](https://github.com/quishqa/WRF-Chem_SP) module to download the data for models performance evaluation.
* [`03_CAM-Chem_scripts`](https://github.com/quishqa/FromGlob2LocSP/tree/main/03_CAM-Chem_scripts): Scripts to compare CAM-Chem simulation against CETESB air quality measurements.
* [`04_Vehicular_Emissions`](https://github.com/quishqa/FromGlob2LocSP/tree/main/04_Vehicular_Emissions): We described how we temporally and spatially disaggregate vehicular emissions in WRF-Chem domain by using road length information. There are also scripts to combine local vehicular emissions with `wrfchemi_d0x` files from global emissions inventories produced by ANTHRO_EMISS.
* [`05_WRF-Chem_scripts`](https://github.com/quishqa/FromGlob2LocSP/tree/main/05_WRF-Chem_scripts): Scripts to compare WRF-Chem simulations against CETESB meteorological and air quality measurements.

* [MUNICH_VEIN_SP](https://github.com/quishqa/MUNICH_VEIN_SP): Here are the scripts to work with MUNICH over São Paulo.

Other useful tools:

* [`qualR`](https://github.com/quishqa/qualR): An R package to download São Paulo air pollution data.
* [`PyChEmiss`](https://github.com/quishqa/PyChEmiss): A Python emissions preprocessor for WRF-Chem regional modeling.
* [`wrf_sp_eval`](https://github.com/quishqa/WRF-Chem_SP): Tools to perform WRF-Chem model evaluation in São Paulo State.
