import os
#from nmeastream import NMEAFileStream
from gnssr4water.io.nmeastream import NMEAFileStream
#from skymaskmod import SkyMask
from gnssr4water.sites.skymask import SkyMask
#from arcbuilder import SatArcBuilder
from gnssr4water.sites.arcbuilder import SatArcBuilder
#from waterlevelestimator import WaterLevelEstimator
from gnssr4water.refl.waterlevelestimator import WaterLevelEstimator
from gnssr4water.core.gnss import GPSL1
#from waterlevel import WaterLevelArc
from gnssr4water.refl.waterlevel import WaterLevelArc
import asyncio
import lz4.frame
import geopandas as gpd
import matplotlib.pyplot as plt
import xarray as xr

def scan_lz4_files(folder):
    for fname in os.listdir(folder):
        if fname.endswith(".lz4"):
            try:
                with lz4.frame.open(os.path.join(folder, fname), 'rt') as f:
                    f.read()
            except Exception as e:
                print(f"❌ File {fname} is corrupt: {e}")
            else:
                print(f"✅ File {fname} is OK")

### Create the NMEA stream object
nmea_files_path = "data/nmea/jinja/"

# ### Collect all valid NMEA files
nmea_files = [os.path.join(nmea_files_path, f) for f in os.listdir(nmea_files_path) if f.endswith(('.gz', '.lz4', '.txt'))]
nmea_stream = NMEAFileStream(nmea_files)

# poly_gdf = gpd.read_file("data/gis/poly/jinja_poly2.shp")
# #poly = poly_gdf.geometry.iloc[0]  # Get the first (and presumably only) geometry
# poly_2=gpd.read_file("data/gis/poly/jinja_poly2.shp").loc[0].geometry
# skyMask = SkyMask(poly_2)
# Load the water mask polygon
# geopoly_gdf = gpd.read_file("data/gis/waterm_geopoly/geopoly.shp")
# geopoly = geopoly_gdf.geometry.iloc[0]  # Get the first geometry

# # Verification that the polygons loaded correctly
# print("Cylinder polygon type:", type(poly))
# print("Water mask polygon type:", type(geopoly))
# ###
# skymask = SkyMask(poly=geopoly)
# # Make sure skymask is already created (with either poly or geopoly)
# ax = skymask.skyplot()
# #print(skymask.poly)
# #print(skymask.poly.bounds)
# plt.show()
# plt.show(block=True)

dataroot='data'
#nc_jinjaout=os.path.join(dataroot,'jinja_masks.nc')
#setup the path the jinja shapefile
jinjashape=os.path.join(dataroot,'gis/poly')

JinjaLoc={"lat":0.414459,"lon":33.207464,"ellipse_height":1135,"recv_height":2.6,"GNSSlambda":GPSL1.length,"noisebandwidth":1}
#load the Jinja shapefile
jinja_geopoly=gpd.read_file(jinjashape).loc[0].geometry

jinjamask1=SkyMask(geopoly=jinja_geopoly,lon=JinjaLoc['lon'],lat=JinjaLoc['lat'],ellipsHeight=JinjaLoc['ellipse_height'],antennaHeight=JinjaLoc['recv_height'],wavelength=JinjaLoc['GNSSlambda'])
#jinjamask1.title="Description of the mask (change me)"


builder = SatArcBuilder(snrStream=nmea_stream, 
                        mask=jinjamask1)
#builder.minpoints = 0

# Create output directory
os.makedirs("plots", exist_ok=True)

async def process_arcs():     
    print("Inside process arcs")
    async for arc in builder.arcs():
        print(f"📡 Got arc from {arc.prn} with {len(arc)} points")
        try:
            wl_arc = WaterLevelArc(arc)
            
            t, height, err = wl_arc.estimateAntennaHeight(antennaHeightBounds=[1, 10])
            wl_arc.setAntennaHeight(height)
            wl_arc.npoly=3
            print(f"🌊 Estimated reflector height: {height:.2f} ± {err:.2f} m")

            # Save to CSV or log
            with open("reflector_heights.csv", "a") as f:
                f.write(f"{t},{arc.prn},{height:.2f},{err:.2f}\n")

            # Create a new figure for this arc
            fig, ax = plt.subplots()
            wl_arc.plot(ax=ax, showfit=True)  # Plot with fitted curve
            plt.title(f"SNR Fit - {arc.prn} - {arc.centralT.strftime('%H:%M:%S')}")

            # Save the figure as PNG
            filename = f"snr_fit_{arc.prn}_{arc.centralT.strftime('%Y%m%d_%H%M%S')}.png"
            filepath = os.path.join("plots", filename)
            plt.savefig(filepath)
            plt.close()
            print(f"✅ Saved SNR plot to {filepath}")

        except Exception as e:
            print(f"⚠️ Failed processing arc from {arc.prn}: {e}")

#scan_lz4_files("data/nmea/jinja/")
asyncio.run(process_arcs())


#estimator = WaterLevelEstimator(arcbuilder=builder,zarrlog="outputs_zarr")
#asyncio.run(estimator.start())
#ds = xr.open_zarr("outputs_zarr", group="waterlevel_ema")
#ds.waterlevel.plot(label="Water Level")
#ds.err_waterlevel.plot(label="Uncertainty", linestyle="--")
#plt.legend()
#plt.title("GNSS-R Estimated Water Level")
#plt.xlabel("Time")
#plt.ylabel("Height (m)")
#os.makedirs("plots/waterestimator",exist_ok=True)
#plt.savefig("plots/waterestimator/water_level_timeseries.png", dpi=300)
#plt.show()
#plt.close()


# print(ds)
# ds.waterlevel.plot(label='Water Level (m)')
# ds.err_waterlevel.plot(label='Uncertainty', linestyle='--')
# #plt.legend()
# plt.title("GNSS-R Water Level Time Series")
# plt.xlabel("Time")
# plt.ylabel("Height (m)")
# plt.legend()
# plt.savefig("plots/water_level_timeseries.png", dpi=300)
# plt.show()
# plt.close()
