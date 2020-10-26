# Creating wrfchemi file based on vehicular emissions

We will spatially disaggregate the total vehicular emissions based on the total street length in each domain cell.
For that, first we'll download the street information for our simulation domain using `OSM2wrfinput.py`,
then we'll disaggregate the emissions into the domain using `VeicEmiss2PyChEmiss.py`,
to finally build the wrfchemi netCDF file using [PyChEmiss](https://github.com/quishqa/PyChEmiss).

## Installation

We recommend using [miniconda](https://docs.conda.io/en/latest/miniconda.html) to work with python.

Then you can create an enviroment to work with this scripts by:

```
conda create --name veic_emiss
conda activate veic_emiss
```

Then install the required packages

```
conda install -c conda-forge osmnx
conda install -c conda-forge xarray dask netCDF4 bottleneck
```

## How to run

### Openstreetmap data to your domain
First, download the openstreetmap data to your domain with `OSM2wrfinput.py`.
You'll need the `geo_em.d01.nc` file.

```
python OSM2wrfinput.py
```

Depending on the size of your domain, It could take some time.
For a domain of 52 x 63 points with &Delta;x = 3 km, the data was download in less than 10 minutes.
However for a domain of 100 x 150 points with &Delta;x = 9 km, it took like 5 hours to finish.
In that case, maybe It is easier to run the script like this:

```
nohup python OSM2wrfinput.py &
```


### Spatial and temporal disaggregation

Then, used the output from the previous step and run  `VeicEmiss2PyChEmiss.py`:

```
python VeicEmiss2PyChEmiss.py
```
This step is very fast, and will create the output to run PyChEmiss.

### Building wrfchemi file

Finally, you will need to install [PyChEmiss](https://github.com/quishqa/PyChEmiss) and run it.
The chemical speciation is based on CBMZ/MOSAIC mechanism.
The species are detailed in `simple_model_d01.yml`, so you just need to used it as template and edit the **nx** and **ny** ans **cell_area** part.

## Acknowledgments

Special thanks to LAPAT, Angel Vara-Vela and Carlos M. Gonzalez who inspired these scripts.
