# Creating wrfchemi file based on vehicular emissions

We will spatially disaggregate the total vehicular emissions based on the total street length in each domain cell.
For that, first we'll download the street information for our simulation domain using `1_Download_OSM.py`,
then we'll calculate the total road length in each cell of WRF-domain with `2_OSM2WRF.py`,
to finally disaggregate the emissions into the domain using `VeicEmiss2PyChEmiss.py`,
to finally build the wrfchemi netCDF file using [PyChEmiss](https://github.com/quishqa/PyChEmiss).

## Installation

We recommend using [miniconda](https://docs.conda.io/en/latest/miniconda.html) to work with python.

Then you can create an enviroment to work with this scripts by:

```
conda create --name veic_emiss
conda activate veic_emiss
```

Then install the required packages. You need `osmnx` to download openstreetmap data,
and `xarray` to read and create netCDF.

```
conda install -c conda-forge osmnx
conda install -c conda-forge xarray dask netCDF4 bottleneck
```

When installing `osmx`, you will also have other packages needed as `pandas`, `geopandas`, and `numpy`.

## How to run

### Download OSM data to your domain

First, download the openstreetmap data to your domain with `1_Download_OSM.py`.
You'll need the `geo_em.d01.nc` file.

```
python 1_Download_OSM.py
```

Depending on the size of your domain, It could take some time.
For a domain of 52 x 63 points with &Delta;x = 3 km, the data was download in less than 10 minutes.
However for a domain of 100 x 150 points with &Delta;x = 9 km, it took like 5 hours to finish.
In that case, maybe It is easier to run the script like this:

```
nohup python 1_Download_OSM.py &
```

### Calculate total length in each WRF-Chem domain cell

The script `2_OSM2WRF.py`, will calculate the total road length in each cell.

```
python 2_OSM2WRF.py
```

These two step can be done by using `OSM2VeicEmiss.py`.

### Spatial and temporal disaggregation

Then, used the output from the previous step and run  `VeicEmiss2PyChEmiss.py`:

```
python 3_VeicEmiss2PyChEmiss.py
```
This step is very fast, and will create the output to run PyChEmiss.

### Building wrfchemi file

Finally, you will need to install [PyChEmiss](https://github.com/quishqa/PyChEmiss) and run it.
The chemical speciation is based on CBMZ/MOSAIC mechanism.
The species are detailed in `simple_model_d01.yml`, so you just need to used it as template and edit the **nx** and **ny** ans **cell_area** part.

## Other scripts

### Number of vehicles in Domain
In `utils` folder, you can use `veh_in_domain.py` to calculate the number of vehicles in your domain. It requires the Installation of `geobr` packages by:

```
pip install geobr
```

More information of `geobr`, in [this link](https://github.com/ipeaGIT/geobr)

### Completing anthropogenic emissions
You can also complement your anthropogenic emissions using global inventories that separates emissions by sources, like EDGAR.

First you need to create the `wrfchemi_00z_d0X`using [ANTHRO_EMISS](https://www2.acom.ucar.edu/wrf-chem/wrf-chem-tools-community) for each source.
That is, you need a `wrfchemi_00z_d0X` for industrial emissions, other for domestic emissions, etc. Then you can add them to your previously created vehicular emissions.



## Acknowledgments

Special thanks to LAPAT, Angel Vara-Vela and Carlos M. Gonzalez who inspired these scripts.
