#!/home/ankom/gimathesis/thesis_code/pygima/bin/python3
# the Shebang statement above is only needed when this script is called from bash 

#import the SkyMask class from the gnssr4water module
from gnssr4water.sites.skymask import SkyMask
from gnssr4water.core.gnss import GPSL1
import os
import geopandas as gpd
import matplotlib.pyplot as plt

#load overall workflow config
import os
from gnssr4water.sites.arcbuilder import SatArcBuilder
from gnssr4water.sites.skymask import SkyMask
from gnssr4water.io.nmeastream import NMEAFileStream
from gnssr4water.refl.waterlevel import WaterLevelArc
import asyncio
import xarray as xr
import numpy as np
from glob import glob


def save_as_netcdf(h_est,netcdf_file,group):
    nest=len(h_est)
    time=[est['time'] for est in h_est]
    height=np.full([nest,3],np.nan)
    height_sig=np.full([nest,3],np.nan)
    prominence=np.full([nest,3],np.nan)
    for i,est in enumerate(h_est):
        npk=len(est['height'])
        height[i,0:npk]=est['height'] 
        height_sig[i,0:npk]=est['height_sigma']
        prominence[i,0:npk]=est['prominence']

    ds=xr.Dataset(dict(height=(['time','peak'],height),
                      height_sig=(['time','peak'],height_sig),
                      prominence=(['time','peak'],prominence)),
                  coords=dict(time=time))
    ds=ds.sortby('time')
    if os.path.exists(netcdf_file):
        os.remove(netcdf_file)
    ds.to_netcdf(netcdf_file,mode='w',group=group)




#read jinja mask
dataroot='data'
nc_jinjaout=os.path.join(dataroot,'jinja_masks.nc')
nc_jinjaout_ls=os.path.join(dataroot,'jinja_timeseries.nc')
jinjamask01=SkyMask.load(nc_jinjaout,group='jinjamask01')


def get_JinjaStream(searchstr="eoaf*lz4"):
    """
    Produce a nmea stream from given data files
    """

    nmeadir=os.path.join(dataroot,'nmea/jinja')

    #all available nmea files (Sorted chronologically)!
    nmealist=sorted(glob(nmeadir+"/"+searchstr))
    nmeastream=NMEAFileStream(nmealist)
    return nmeastream,nmealist

#initialize an arc builder
def get_arcbuilder(mask):
    #currently only a subset of the dat is read
    jinjastream,jinjalist=get_JinjaStream('*2024-02*')
    arcfilterConfig={"minLengthSec":600,"split":True,"minElevationSpan":6,"mask":jinjamask01}
    return SatArcBuilder(jinjastream,**arcfilterConfig)




async def compute_timeseries():

    arcbuild01=get_arcbuilder(jinjamask01)
    heightbounds=[1.5,7.0]

    preprocess={"npoly":1,"noiseBandwidth":1,"maxpeaks":3,"resolution":0.01}
    
    i=0

    height_est=[]
    async for arc in arcbuild01.arcs():

        i+=1
        wlarc=WaterLevelArc(arc)
        h_est=wlarc.estimateAntennaHeight_multi_LombScargle(heightbounds,**preprocess)

        if i%10 == 0: 
            print(f"Processed {i:05d} arcs {arc.time[0]}",end='\r')
        height_est.append(h_est)
    
    print(f"Writing time series results to {nc_jinjaout_ls}")
    save_as_netcdf(height_est,netcdf_file=nc_jinjaout_ls,group="jinja_ls")


if __name__ == "__main__":

    asyncio.run(compute_timeseries())


