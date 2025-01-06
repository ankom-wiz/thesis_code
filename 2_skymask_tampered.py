# This file is part of gnssr4water
# gnssr4water is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.

# geoslurp is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with Frommle; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# Author Roelof Rietbroek (r.rietbroek@utwente.nl), 2024

from gnssr4water.core.logger import log
from gnssr4water.io.cf import global_attrs
from pymap3d.enu import geodetic2enu,enu2aer,aer2enu,enu2geodetic
from pymap3d import Ellipsoid
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.patches as mplpatch
wgs84=Ellipsoid.from_name('wgs84')
from scipy.stats import  binned_statistic_2d
from shapely.geometry import Polygon,Point
from shapely import to_wkt,from_wkt
import xarray as xr
from datetime import datetime
from gnssr4water.core.gnss import GPSL1
from gnssr4water.fresnel import firstFresnelZone,elev_from_radius
from numba import jit

### Added imports
import geopandas as gpd
import matplotlib.pyplot as plt

### Read properly shapefiles and extract their geometry using geopandas
# Load the cylindrical polygon around receiver
poly_gdf = gpd.read_file(r"C:\Users\Anastasios_Komiotis\Desktop\data_LakeVictoria\poly\jinja_poly.shp")
poly = poly_gdf.geometry.iloc[0]  # Get the first (and presumably only) geometry

# Load the water mask polygon
geopoly_gdf = gpd.read_file(r"C:\Users\Anastasios_Komiotis\Desktop\data_LakeVictoria\waterm_geopoly\water_mask.shp")
geopoly = geopoly_gdf.geometry.iloc[0]  # Get the first geometry

# Verification that the polygons loaded correctly
print("Cylinder polygon type:", type(poly))
print("Water mask polygon type:", type(geopoly))
###


# Converts a polygon into an azimuth-elevation polygon (azelpoly)
def geo2azelpoly(geopoly,lon,lat,ellipsHeight,antennaHeight,wavelength=GPSL1.length):
    ### Debugging check
    print("geo2azelpoly: Input geopoly:", geopoly)
    print("geo2azelpoly: Reference location - lon:", lon, "lat:", lat, "ellipsHeight:", ellipsHeight,
          "Antenna height:", antennaHeight, "Wavelength:", wavelength)
    ###
    
    # Checks if the input polygon has a simple structure (no self-intersections)
    if not geopoly.is_simple:
        log.warning("Cannot (currently) handle polygons with interiors, taking exterior only")

    # Extracts the polygon's exterior boundary longitude and latitude coordinates
    plon,plat=geopoly.exterior.coords.xy
    ph=ellipsHeight*np.ones(len(plon))
    # We need to convert the lon,lat polygons,fixed to the plane  in the local ENU frame
    e,n,u=geodetic2enu(lat=plat, lon=plon, h=ph, lat0=lat, lon0=lon, h0=ellipsHeight, ell=wgs84, deg=True)
    # Converts ENU coordinates into azimuth, elevation and range
    az,e,r=enu2aer(e,n,u)
    
    #compute the actual elevation assuming the reflection point is a specular point 
    # elev=np.rad2deg(np.arctan2(antennaHeight,r))
    

    #compute the elevation correpsonding to the centroids of the First Fresnel zone
    elev=elev_from_radius(r,antennaHeight,wavelength)
    
    #convert into a shapely polygon
    azelpoly=Polygon(zip(az,elev))
    return azelpoly

# Converts an azimuth-elevation polygon into a geodetic polygon
def azel2geopoly(azelpoly,lon,lat,ellipsHeight,antennaHeight,wavelength=GPSL1.length):
    ### Debugging check
    print("azel2geopoly: Input azelpoly:", azelpoly)
    print("azel2geopoly: Reference location - lon:", lon, "lat:", lat, "ellipsHeight:", ellipsHeight,
          "Antenna height:", antennaHeight, "Wavelength:", wavelength)
    ###
    
    if not azelpoly.is_simple:
        log.warning("Cannot (currently) handle polygons with interiors, taking exterior only")
    
    az,el=azelpoly.exterior.coords.xy
    #compute the radius of the location of the specular point
    #radius=antennaHeight/np.tan(np.deg2rad(el))
    
    #compute radius of fresnel central point
    _,_,radius,_ = firstFresnelZone(wavelength, antennaHeight,np.array(el)) 
    # For the reflection points we actually assume the point has a 0 upward component (both up and elevation component are set to 0) 
    el0=np.zeros(len(el))
    u0=np.zeros(len(el))
    
    e,n,_=aer2enu(az,el0,radius,deg=True)
    # import pdb;pdb.set_trace()
    plat,plon,_=enu2geodetic(e,n,u0,lat0=lat,lon0=lon,h0=ellipsHeight,ell=wgs84,deg=True)

    #convert into a shapely polygon
    geopoly=Polygon(zip(plon,plat))
    return geopoly

# Determines if a given azimuth-elevation point lies inside a polygon
@jit(nopython=True)
def masked_fast(polygon,elevation,azimuth) -> bool:
    """Fast polygon test adapted from this discussion here:https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python"""

    length = polygon.shape[0]-1
    dy2 = elevation - polygon[0][1]
    intersections = 0
    ii = 0
    jj = 1

    while ii<length:
        dy  = dy2
        dy2 = elevation - polygon[jj][1]

        # consider only lines which are not completely above/bellow/right from the point
        if dy*dy2 <= 0.0 and (azimuth >= polygon[ii][0] or azimuth >= polygon[jj][0]):

            # non-horizontal line
            if dy<0 or dy2<0:
                F = dy*(polygon[jj][0] - polygon[ii][0])/(dy-dy2) + polygon[ii][0]

                if azimuth > F: # if line is left from the point - the ray moving towards left, will intersect it
                    intersections += 1
                elif azimuth == F: # point on line
                    return 2

            # point on upper peak (dy2=dx2=0) or horizontal line (dy=dy2=0 and dx*dx2<=0)
            elif dy2==0 and (azimuth==polygon[jj][0] or (dy==0 and (azimuth-polygon[ii][0])*(azimuth-polygon[jj][0])<=0)):
                return 2

        ii = jj
        jj += 1

    #print 'intersections =', intersections
    return intersections & 1


# detailed implementation used to manage and manipulate the geometry of sky masks: 
# areas in the sky that are either "visible" or "masked" based on a receiver's location, elevation and azimuth
# tool for managing and processing sky masks
# functionality includes: saving/loading data, checking if points are masked,
# handling residuals, and visualizing the mask through sky plots
class SkyMask:
    group="skymask" 
    # change poly and geopoly to a shp from local (or download its geometry) 
    # poly -> cylindrical polygon from receiver
    # geopoly -> mask on water alone -> length river about 100m
    # lon and lat is the receiver coordinates
    # ellipsHeight -> ortho_height -> 1135
    def __init__(self,poly=poly,geopoly=geopoly,lon=33.207464,lat=0.414459,ellipsHeight=1135,antennaHeight=2.6,wavelength=GPSL1.length,noisebandwidth=1):
        
        self.res_elev=[]
        self.res_az=[]
        self.res_snr=[]
        self.poly=None
        self.geopoly=None

        # Global attributes that provide metadata for the sky mask on an xarray.Dataset
        globattr=global_attrs()
        globattr["title"]="GNSS-R selection skymask"
        globattr['GNSSWavelength']=wavelength
        self._ds=xr.Dataset(attrs=globattr) #xarray structure to store things into
        
        # Note all arguments are optional so an empty mask can be created 
        # Receiver spatial parameters are assigned to instance variables
        self.lon=lon
        self.lat=lat
        self.ellipseHeight=ellipsHeight

        
        self.antennaHeight=antennaHeight
        self.noiseBandwidth=noisebandwidth

        # Checks if both poly and geopoly are provided at the same time
        # Raises error because both formats for the mask geometry cannot be specified simultaneously
        if poly is not None and geopoly is not None:
            raise RuntimeError("Input is ambigious if both geopoly and poly are provided")
            
        # If poly is provided, it sets self.poly to the given polygon 
        # and uses the azel2geopoly function to convert it into geographic coordinates (longitude/latitude)
        if poly is not None: 
            self.poly=poly
            self.geopoly=azel2geopoly(poly,lon=lon,lat=lat,ellipsHeight=ellipsHeight,antennaHeight=antennaHeight)

        # If geopoly is provided, it sets self.geopoly to the given geographic polygon 
        # and converts it to azimuth-elevation coordinates using the geo2azelpoly function
        if geopoly is not None:
            self.geopoly=geopoly
            self.poly=geo2azelpoly(geopoly,lon=lon,lat=lat,ellipsHeight=ellipsHeight,antennaHeight=antennaHeight)
        self._preppoly()
    
    def _preppoly(self):
        #for fast polygon computation
        # Converts the polygon (if it exists) into a numpy array and stores it in self._poly
        if self.poly is not None:
            self._poly=np.array(self.poly.exterior.xy).T

    # Receiver property getter and setter methods
    @property
    def antennaHeight(self):
        return self._ds.attrs['receiver_antennaheight']

    @antennaHeight.setter
    def antennaHeight(self,height):
        if height is not None:
            self._ds.attrs['receiver_antennaheight']=height
    
    @property
    def ellipseHeight(self):
        return self._ds.attrs['receiver_ellipsheight']

    @ellipseHeight.setter
    def ellipseHeight(self,height):
        if height is not None:
            self._ds.attrs['receiver_ellipsheight']=height
    
    @property
    def lon(self):
        return self._ds.attrs['receiver_lon']

    @lon.setter
    def lon(self,lon):
        if lon is not None:
            self._ds.attrs['receiver_lon']=lon
    
    @property
    def lat(self):
        return self._ds.attrs['receiver_lat']

    @lat.setter
    def lat(self,lat):
        if lat is not None:
            self._ds.attrs['receiver_lat']=lat

    @property
    def noiseBandwidth(self):
        return self._ds.attrs['receiver_noisebandwidth']

    @noiseBandwidth.setter
    def noiseBandwidth(self,bw):
        self._ds.attrs['receiver_noisebandwidth']=bw

    # Allows loading a SkyMask from a file (NetCDF or Zarr)
    @staticmethod
    def load(filename):
        #start with an empty skymask
        skmsk=SkyMask()
        engine=None
        if filename.endswith('.zarr'):
            engine='zarr'
        
        # Reads  polygon data from WKT format -> initializes the geometry -> prepares polygon for computation
        with xr.open_dataset(filename,group=SkyMask.group,engine=engine) as ds:
            skmsk._ds=ds.copy()

        skmsk.poly=from_wkt(skmsk._ds.attrs['azelpoly_wkt'])
        skmsk.geopoly=from_wkt(skmsk._ds.attrs['lonlatpoly_wkt'])
        skmsk._preppoly()
        return skmsk


    # Checks if the given elevation and azimuth are inside the mask polygon
    def masked (self,elevation,azimuth)-> bool:
        return not masked_fast(self._poly,elevation,azimuth)
        # val2= not self.poly.contains(Point(azimuth,elevation))
        # if val != val2:
        # import pdb;pdb.set_trace()
        # return val

    # Checks points to see if they are inside the mask, returning a boolean array
    def isMasked(self,elevation,azimuth):
        """
        returns a boolean array
        """
        return np.array([self.masked(el,az) for el,az in zip(elevation,azimuth)])
        
    # Retrieves weights (SNR error) for a given azimuth and elevation using the nearest neighbor interpolation
    def weights(self,azimuth,elevation):

        return self._ds.snr_error.sel(azimuth=xr.DataArray(azimuth,dims='narc'),elevation=xr.DataArray(elevation,dims='narc'),method='nearest')

    @property
    def title(self):
        return self._ds.attrs['title']

    @title.setter
    def title(self,title):
        self._ds.attrs["title"]=title
    
    # Adds the history attribute of the dataset, including the current timestamp
    def add_history(self,action):
        self._ds.attrs["history"].append(datetime.now().isoformat()+f": {action}")

    # Adds SNR residuals (elevation, azimuth, and SNR residuals) to the class instance
    def append_SNRresidual(self,elev,az,snrres):
        """Append SNR residuals points to the mask

            Parameters
            ----------
            elev : array_like [n]
            Elevation of the points in [degrees]
            az : array_like[n]
            Azimuth of the points in degrees
            snrres: array_like[n]
            SNR residuals
        """
        self.res_elev.extend(elev) 
        self.res_az.extend(az)
        self.res_snr.extend(snrres)
    
    # Computes a weight mask using residual SNR 
    # it calculates a median SNR value for each bin and assigns it to the dataset
    def compute_WeightMask(self,fillmethod='median'):
        if len(self.res_elev) < 10:
            log.info("Enough SMR resildulas must have been added before being able to compute a weight mask")
            return
        
        minaz,minel,maxaz,maxel=self.poly.bounds
        #determine bins
        deltad=2
        elevbins=int((maxel-minel)/deltad)+1
        azbins=int((maxaz-minaz)/deltad)+1
        stat,x_edges,y_edges,bino=binned_statistic_2d(self.res_az,self.res_elev,np.abs(self.res_snr),statistic='median',bins=[azbins,elevbins])
        
        counts,_,_,_=binned_statistic_2d(self.res_az,self.res_elev,None,statistic='count',bins=[azbins,elevbins])
        
        weights=np.where(counts > 100,stat,np.nan)
        self._ds["grd_azimuth"]=(("grd_azimuth",),x_edges) 
        self._ds["grd_elevation"]=(("grd_elevation",),y_edges) 

        self._ds["snr_error"]=(("azimuth","elevation"),weights)

        if fillmethod == "median":
            self._ds['snr_error']=self._ds['snr_error'].fillna(np.median(self.res_snr))

        self._ds=self._ds.assign_coords(azimuth=(("azimuth",),x_edges[0:-1]+deltad/2),elevation=(("elevation",),y_edges[0:-1]+deltad/2))
                                    

    # Creates a skyplot showing the visibility of satellite signals based on the defined polygon/mask
    # It can optionally include the weight map (SNR error)
    def skyplot(self,ax=None,**kwargs):
        maskcolor='#e15989'
        #maskcolor='blue'
        
        if ax is None:
            fig,ax=mpl.subplots(1,1,subplot_kw={'projection': 'polar'})
            ax.title='Skyplot mask'

            ax.set_rlim([90,0])
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)

        #create a matplotlib patch from the polygon data
        polyfine=self.poly.segmentize(2)
        x,y=polyfine.exterior.xy
        ppatch=mplpatch.Polygon(np.array(np.array([[x,y] for x,y in zip(np.deg2rad(x),y)])),edgecolor=maskcolor,fill=False,lw=2)

        #plot the weights within the mask
        if "snr_error" in self._ds:
            # self._ds.snr_error.plot.pcolormesh(x=np.deg2rad(self._ds.azimuth),y=self._ds.elevation,ax=ax)
            axob=ax.pcolormesh(np.deg2rad(self._ds.grd_azimuth),self._ds.grd_elevation,self._ds.snr_error.data.T,shading='flat',cmap='Reds',vmin=0,vmax=30)
            fig=ax.get_figure()
            fig.colorbar(axob,label="SNR error [v/v]")
        
        ax.add_patch(ppatch)

        return ax


    # Saves the mask to a file (NetCDF or Zarr)
    # polygon geometries are stored as WKT in the file attributes
    def save(self,arName,mode='a'):
        """ Save the mask to an archive (netcdf or zarr), for later reuse

            Parameters
            ----------
            Arname: str
            File name to write to
            mode: str
            Mode to write ('a' for appending,'w' for overwriting)
        """
        #save the polygon  as a WKT attribute to the netcdf file
        self._ds.attrs["azelpoly_wkt"]=to_wkt(self.poly)
        self._ds.attrs["lonlatpoly_wkt"]=to_wkt(self.geopoly)
        if arName.endswith('.nc'):
            self._ds.to_netcdf(arName,mode=mode,group=SkyMask.group)
        elif arName.endswith(".zarr"):
            self._ds.to_zarr(arName,mode=mode,group=SkyMask.group)
        else:
            raise RuntimeError(f"archive not supported {arName}")
    
        
    # Segments the azimuth-elevation polygon into smaller parts (if necessary), creating a new SkyMask object
    def segmentize(self,max_segment_length=1):
        """
        Segmentize the azimuth-elevation polygon and create a new Skymask 
        """

        lon=self._ds.attrs["receiver_lon"]
        lat=self._ds.attrs["receiver_lat"]
        oh=self._ds.attrs['receiver_ellipsheight']
        ah=self._ds.attrs['receiver_antennaheight']

        skmsk=SkyMask(poly=self.poly.segmentize(max_segment_length=max_segment_length),lon=lon,lat=lat,ellipsHeight=oh,antennaHeight=ah)
        return skmsk



class SimpleMask(SkyMask):
    def __init__(self,lon,lat,ellipsHeight,antennaHeight,elevations=[5,40],azimuths=[0,360],wavelength=GPSL1.length):
        pnts=[(azimuths[0],elevations[0]),(azimuths[1],elevations[0]),(azimuths[1],elevations[1]),(azimuths[0],elevations[1]),(azimuths[0],elevations[0])]
        super().__init__(poly=Polygon(pnts),lon=lon,lat=lat,ellipsHeight=ellipsHeight,antennaHeight=antennaHeight,wavelength=wavelength)
        self.elevBnds=elevations
        self.azBnds=azimuths

    def masked (self,elevation,azimuth)-> bool:

        if elevation < self.elevBnds[0] or elevation > self.elevBnds[1]:
            # import pdb;pdb.set_trace()
            return True
        if azimuth < 0:
            azimuth+=360
        if azimuth < self.azBnds[0] or azimuth > self.azBnds[1]:
            # import pdb;pdb.set_trace()
            return True
        #ok point is not masked
        return False

### Make an object
# Initialize SkyMask with the geopoly
# skymask = SkyMask(poly=poly,geopoly=geopoly, lon=lon, lat=lat, ellipsHeight=ellipsHeight, antennaHeight=antennaHeight)
