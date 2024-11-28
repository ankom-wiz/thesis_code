pip install git+https://github.com/ITC-Water-Resources/gnssr4water.git

#general purpose libraires
import math
import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import cartopy
import shapefile
#visualization libraries
%matplotlib inline
import matplotlib.pyplot as plt

# ^^^ SO FAR SO GOOD IT WORKS

#project modules
import gnssr4river.fresnel.fresnelzone
import gnssr4river.fresnel.plotfresnel
import gnssr4river.fresnel.iterfresnel
import gnssr4river.fresnel.getorbits
import gnssr4river.fresnel.geod
import gnssr4river.fresnel.intersect
import gnssr4river.refl.nmea
