# Downloading CETESB data for model evaluation

To download CETESB data, we used `qualar_py` from [wrf_sp_eval](https://github.com/quishqa/WRF-Chem_SP)
module.
We first download the pollutant information and then the meteorology information from all available air quality stations,
from October 6th to October 13th, 2014.
Each air quality station information is saved in a Pandas DataFrame,
and all DataFrames are saved in a dictionary.
We saved this dictionary as a pickle.

* `cetesb_pol_2014.pkl` : Contains pollutant information from all air quality station
* `cetesb_met_2014.pkl` : Contains meterological information from all air quality station
