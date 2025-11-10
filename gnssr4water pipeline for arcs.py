from datetime import datetime
import os
from gnssr4water.io.nmeastream import NMEAFileStream
from gnssr4water.sites.skymask import SkyMask
from gnssr4water.sites.arcbuilder import SatArcBuilder
from gnssr4water.refl.waterlevel import WaterLevelArc
from gnssr4water.refl.waterlevelestimator import WaterLevelEstimator
from gnssr4water.core.gnss import GPSL1

import asyncio
import lz4.frame
import geopandas as gpd
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from scipy.io import savemat

# Modify later to incorporate simulator data
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
 
def save_snr_from_arc(arc,wl_arc):
    """
    Plots the water level arc with a fitted curve and saves the plot as a PNG file.
    Args:
        arc: An object containing information about the arc, including `prn` (satellite PRN) 
             and `centralT` (central time of the arc as a datetime object).
        wl_arc: An object with a `plot` method that generates the water level plot. 
                The `plot` method should accept `ax` and `showfit` parameters.
    Behavior:
        - Creates a new figure for the given arc.
        - Plots the water level arc with a fitted curve.
        - Saves the generated plot as a PNG file in the "plots" directory.
        - The filename is formatted as "snr_fit_<prn>_<centralT>.png", where `<prn>` is the 
          satellite PRN and `<centralT>` is the central time of the arc in the format 
          "YYYYMMDD_HHMMSS".
        - Prints a confirmation message with the file path of the saved plot.
    Raises:
        Any exceptions raised by the `plot` method or file-saving operations will propagate.
    """
    # Create a new figure for this arc
    ax = wl_arc.plot(ax=None, showfit=True)  # Plot with fitted curve
    fig = ax.figure
    # Save the figure as PNG
    filename = f"snr_fit_{arc.prn}_{arc.centralT.strftime('%Y%m%d_%H%M%S')}.png"
    filepath = os.path.join("plots", filename)
    fig.savefig(filepath)
    print(f"✅ Saved SNR plot to {filepath}")

def save_arc_to_matlab(arc, output_dir="data/matlab_output"):
    """
    Saves arc data to a MATLAB .mat file
    
    Args:
        arc: Arc object containing satellite data
        output_dir (str): Directory path where to save the .mat file
    
    Returns:
        str: Path to the saved file
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Prepare the data dictionary with reshaped arrays
    mat_data = {
        'elevation': arc.elev.reshape(-1, 1),  # Reshape to Nx1
        'azimuth': arc.az.reshape(-1, 1),      # Reshape to Nx1
        'snr': arc.cn0.reshape(-1, 1),         # Reshape to Nx1
        'time': np.array([t.timestamp() for t in arc.time]).reshape(-1, 1),  # Reshape to Nx1
        'prn': arc.prn,
        'system': arc.system
    }
    
    # Generate filename with timestamp
    filename = f"arc_{arc.prn}_{arc.centralT.strftime('%Y%m%d_%H%M%S')}.mat"
    filepath = os.path.join(output_dir, filename)
    
    # Save to .mat file
    savemat(filepath, mat_data)
    print(f"✅ Saved MATLAB file: {filepath}")
    
    return filepath

def save_merged_matlab(arcs_data, output_dir="data/matlab_output"):
    """
    Saves all arcs data to a single MATLAB .mat file
    
    Args:
        arcs_data (list): List of dictionaries containing arc data
        output_dir (str): Directory path where to save the .mat file
    
    Returns:
        str: Path to the saved file
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Merge all data into a single dictionary
    merged_data = {
        'elevation': [],
        'azimuth': [],
        'snr': [],
        'time': [],
        'prn': []
    }
    
    for data in arcs_data:
        merged_data['elevation'].append(data['elevation'])
        merged_data['azimuth'].append(data['azimuth'])
        merged_data['snr'].append(data['snr'])
        merged_data['time'].append(data['time'])
        merged_data['prn'].append(data['prn'])
        #merged_data['system'].append(data['system'])
    
    # Generate filename with timestamp
    filename = f"merged_arcs_{datetime.now().strftime('%Y%m%d_%H%M%S')}.mat"
    filepath = os.path.join(output_dir, filename)
    
    # Save to .mat file
    savemat(filepath, merged_data)
    print(f"✅ Saved merged MATLAB file: {filepath}")
    
    return filepath
 
def save_LombScargle(arc,wl_arc):
    """
    Saves a Lomb-Scargle periodogram plot as a PNG file.
    This function generates a Lomb-Scargle periodogram plot using the provided
    `wl_arc` object, saves the plot as a PNG file in the "plots" directory, and
    prints a confirmation message with the file path.
    Args:
        arc: An object containing metadata for the plot, including:
            - `prn`: A unique identifier for the arc.
            - `centralT`: A datetime object representing the central time of the arc.
        wl_arc: An object with a `plotLombScargle` method that generates the plot.
                The method should accept a `showFit` parameter.
    Returns:
        None
    """
    ax = wl_arc.plotLombScargle(showFit=True)
    fig = ax.figure
    fileName = f"LombScargle_{arc.prn}_{arc.centralT.strftime('%Y%m%d_%H%M%S')}.png"
    filePath = os.path.join("plots", fileName)
    fig.savefig(filePath)
    # Save the figure as PNG
    print(f"✅ Saved LombScargle plot to {filePath}")
 
def plot_water_level_estimator():
    """
    Plots the estimated water level and its uncertainty over time, saves the plot as a PNG file, 
    and displays it.
 
    The function reads data from a Zarr dataset located in the "outputs_zarr" directory under 
    the group "waterlevel_ema". It plots the water level and its associated uncertainty 
    (with a dashed line) using Matplotlib. The plot is saved in the "plots/waterestimator" 
    directory as "water_level_timeseries.png" with a resolution of 300 DPI. If the directory 
    does not exist, it is created. The plot is displayed and then closed.
 
    Dependencies:
        - xarray (xr)
        - matplotlib.pyplot (plt)
        - os
 
    Raises:
        FileNotFoundError: If the "outputs_zarr" directory or the specified group does not exist.
    """
    dataset = xr.open_zarr("outputs_zarr", group="waterlevel_ema")
    dataset.waterlevel.plot(label="Water Level")
    dataset.err_waterlevel.plot(label="Uncertainty", linestyle="--")
    plt.legend()
    plt.title("GNSS-R Estimated Water Level")
    os.makedirs("plots/waterestimator",exist_ok=True)
    plt.savefig("plots/waterestimator/water_level_timeseries.png", dpi=300)
    plt.show()
    plt.close()
 
def createMaskFromPoly(JinjaLoc):
    # Read properly shapefiles and extract their geometry using geopandas
    #setup the path the jinja shapefile
    jinjashape=os.path.join("data",'gis/poly')
 
    #load the Jinja shapefile
    jinja_poly=gpd.read_file(jinjashape).loc[0].geometry
    return SkyMask(poly=jinja_poly,lon=JinjaLoc['lon'],lat=JinjaLoc['lat'],ellipsHeight=JinjaLoc['ellipse_height'],antennaHeight=JinjaLoc['recv_height'],wavelength=JinjaLoc['GNSSlambda'])
 
def createMaskFromGeoPoly(JinjaLoc):
    # Read properly shapefiles and extract their geometry using geopandas
    jinjashape = os.path.join("data", 'gis/poly')
    #Should have been: jinjashape = os.path.join("data", 'gis/waterm_geopoly')
    #Because we're passing the geopoly below the code should have worked but for some reason the code runs fine with the poly is assigned as geopoly
    jinja_geopoly = gpd.read_file(jinjashape).loc[0].geometry
    return SkyMask(geopoly=jinja_geopoly,lon=JinjaLoc['lon'],lat=JinjaLoc['lat'],ellipsHeight=JinjaLoc['ellipse_height'],antennaHeight=JinjaLoc['recv_height'],wavelength=JinjaLoc['GNSSlambda'])
 
### Create the NMEA stream object
nmea_files_path = "data/nmea/jinja/"
 
# ### Collect all valid NMEA files
nmea_files = [os.path.join(nmea_files_path, f) for f in os.listdir(nmea_files_path) if f.endswith(('.gz', '.lz4', '.txt'))]
for i in range(len(nmea_files)):
    print(f"✅ NMEA file found: {nmea_files[i]}")
#Create the nmea stream with a single file
nmea_stream = NMEAFileStream(nmea_files)

snrStream = nmea_stream
 
JinjaLoc={"lat":0.414459,"lon":33.207464,"ellipse_height":1135,"recv_height":2.6,"GNSSlambda":GPSL1.length,"noisebandwidth":1}
 
#skyMask = createMaskFromPoly(JinjaLoc)
skyMask = createMaskFromGeoPoly(JinjaLoc)
 
#Create skymask plot
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
skyMask.skyplot(ax=ax)

# Save to file
fig.savefig("skyplot.png", dpi=300, bbox_inches='tight')
 
# builder = SatArcBuilder(snrStream=nmea_stream, 
#                         mask=skyMask)
builder = SatArcBuilder(snrStream=snrStream,mask=skyMask) 
# Create output directory
os.makedirs("plots", exist_ok=True)
 
async def save_snr_plots():     
    print("Inside process arcs")
    arcs_data = []  # List to store all arc data
    total_arcs = 0
    
    # First count total arcs
    # async for _ in builder.arcs():
    #     total_arcs += 1
    #     print(f"Total arcs to process: {total_arcs}")

    # Reset the builder
    # builder.arcqueue = asyncio.Queue(builder.maxarcs)
    # builder.arccache = {}

    # Process arcs with progress tracking
    processed = 0
    async for arc in builder.arcs():
        processed += 1
        print(f"📡 Got arc from {arc.prn} with {len(arc)} points")
        try:
            
            # Check if we have enough points
            if len(arc) < 20:  # Minimum points needed for reliable computation
                #print(f"⚠️ Skipping arc from {arc.prn}: Not enough points ({len(arc)})")
                continue

            #save_arc_to_matlab(arc)

            # Store arc data for merged file
            arc_data = {
                'elevation': arc.elev.reshape(-1, 1),
                'azimuth': arc.az.reshape(-1, 1),
                'snr': arc.cn0.reshape(-1, 1),
                'time': np.array([t.timestamp() for t in arc.time]).reshape(-1, 1),
                'prn': arc.prn
            }
            arcs_data.append(arc_data)
            print(f"✅ Processed {processed} arcs")
            # --- start working codebase --- #
            # wl_arc = WaterLevelArc(arc)
            # t, height, err = wl_arc.estimateAntennaHeight(antennaHeightBounds=[1, 10])
            # wl_arc.setAntennaHeight(height)
            # wl_arc.npoly=3
            # print(f"🌊 Estimated reflector height: {height:.2f} ± {err:.2f} m")
 
            # # Save to CSV or log
            # with open("reflector_heights.csv", "a") as f:
            #     f.write(f"{t},{arc.prn},{height:.2f},{err:.2f}\n")
                
            # save_snr_from_arc(arc, wl_arc)

            # --- end working codebase --- #

            # try:
            #     save_LombScargle(arc, wl_arc)
            #                     #save_LombScargle(arc, wl_arc)
            # except np.linalg.LinAlgError:
            #     print(f"⚠️ Skipping arc from {arc.prn}: Singular matrix in height estimation")
            #     continue
        except Exception as e:
            print(f"⚠️ Failed processing arc from {arc.prn}: {e}")
            # Save merged data after all arcs are processed
    if arcs_data:
        save_merged_matlab(arcs_data)
 
 
# Uncomment or comment to turn process each individual arc and create plots for signal-to-noise ratio
asyncio.run(save_snr_plots())
 
# Uncomment to run water level estimator and create plots for water level and uncertainty
# os.makedirs("outputs_zarr", exist_ok=True)
# estimator = WaterLevelEstimator(arcbuilder=builder,zarrlog="outputs_zarr")
 
# # Starts the estimator for all water level arcs
# asyncio.run(estimator.start())
 
# Plotting the water levels and their uncertainties
##plot_water_level_estimator()
