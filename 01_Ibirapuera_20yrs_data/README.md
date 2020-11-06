# Plotting 20 years of pollutant concentration information in Ibirapuera Air quality station

Here is an example on how to use [`qualR`](https://github.com/quishqa/qualR) to download 20 years of information.
We then used this information to see the trend in concentration, whether they decline in time or not.

* `data_download.R`: Here we download the data, because it's 20 year it could take a couple of hours. In case there is a problem in the connection, we decided to save each pollutant as soon it is downloaded.
* `plot_concentration.R`: We calculate the monthly mean for each pollutant to plot.
