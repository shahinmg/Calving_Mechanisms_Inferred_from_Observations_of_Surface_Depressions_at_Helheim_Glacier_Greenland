## Calving Mechanisms Inferred from Observations of Surface Depressions at Helheim Glacier Greenland

### Recreate Figures
Notebooks to recreate figures are in the `figure_notebooks` directory 

### Figure Data
Data used for the figures are in the `data` directory. If you are looking for the grounding zone and ground line data, please see the `data/ground_lines_and_grounding_zones` directory. The 2018 grounding lines are in the `2018_grounding_lines_merged_v2.gpkg` file and the 2019 grounding lines are in `2018_grounding_lines_merged_v2.gpkg`. 

### Analysis
Code used for some analysis are in the `analysis` directory. The fuctions for the force balance code are in `force_balance.py`. The script to run the force balance is `force_balance_stresses_its_live_arctic_dem.py`. Note this code was ran with ArcticDem and ITS_LIVE image-pair velocities. The force balance data file names are in `data/force_balance_data_table_S3.csv`. 

### ATLAS Data
The Autonomous Terrestrial Laser Scanner (ATLAS) cloud optimized point clouds are currently in these public s3 buckets, `s3://grid-public-ept/atlas/flat/ATLAS-South/` and `s3://grid-public-ept/atlas/flat/ATLAS-North/`. The velocities are here, `s3://grid-public-ept/atlas/velocity/` and the dems here, `s3://grid-public-ept/atlas/dem/`. 

### Example of how to list and copy data from an s3 bucket

Use the [amazon web services command line interface](https://aws.amazon.com/cli/). 

```
conda install -c conda-forge awscli
```

To list the atlas directory

```
aws s3 ls --no-sign-request s3://atlas-lidar-helheim/
```

To copy a file to your local machine 


```
aws s3 cp --no-sign-request s3://atlas-lidar-helheim/dem/ATLAS-North/191215_121054_idw_geoid_rm.cog.tif . 
```

