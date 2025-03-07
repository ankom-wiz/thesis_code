#!/home/ankom/gimathesis/thesis_code/pygima/bin/python3
# the Shebang statement above is only needed when this script is called from bash 

#import the SkyMask class from the gnssr4water module
from gnssr4water.sites.skymask import SkyMask
from gnssr4water.core.gnss import GPSL1
import os
import geopandas as gpd
import matplotlib.pyplot as plt


def main():
    """
    Main function to be called
    """
    dataroot='data'
    nc_jinjaout=os.path.join(dataroot,'jinja_masks.nc')

    #setup the path the jinja shapefile
    jinjashape=os.path.join(dataroot,'gis/poly')
    
    JinjaLoc={"lat":0.414459,"lon":33.207464,"ellipse_height":1135,"recv_height":2.6,"GNSSlambda":GPSL1.length,"noisebandwidth":1}

    #load the Jinja shapefile
    jinja_geopoly=gpd.read_file(jinjashape).loc[0].geometry
    
    jinjamask1=SkyMask(geopoly=jinja_geopoly,lon=JinjaLoc['lon'],lat=JinjaLoc['lat'],ellipsHeight=JinjaLoc['ellipse_height'],antennaHeight=JinjaLoc['recv_height'],wavelength=JinjaLoc['GNSSlambda'])
    jinjamask1.title="Description of the mask (change me)"

    #save/overwrite the mask to a netcdf file 
    jinjamask1.save(nc_jinjaout,group="jinjamask01")


    #possibly create a  quick plot for visualization
    ax=jinjamask1.skyplot()
    plt.savefig('media/jinjamask01_test.png')


if  __name__ == "__main__":
    # This part is executed when this script is called as a script (e.g. python thisscript.py)

    # call the main function as defined above
    main()
