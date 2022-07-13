### Calculate Ventilation Index parameters for WRF-AHW simulations ###
### This method uses hourly WRF output ###
### Created by Michael Fischer on 10/25/17 ###

# Import libraries:
import numpy as np
import datetime as DT
import netCDF4 as NC4
import wrf
import glob
import pcmin
from scipy.interpolate import RectBivariateSpline as rbs
from scipy.interpolate import LinearNDInterpolator as ninterp
from scipy import interpolate
import os

# Establish working directory:
os.chdir('/network/rit/lab/tcdynasty/fischer/Data/Modeling/RI/Gonzalo/Ventilation/')
outputfile = 'Gonzalo_ventilation_params.npy' #be sure to delete this file if want to start from scratch

# Establish constants:
num_ens = 60 #Number of ensemble members
Ae = 6.37E6 #radius of Earth (m)
Rd = 287.04 #gas constant for dry air
Rv = 461.50 #gas constant for water vapor
eps = Rd/Rv
cpd = 1005.7 #specific heat at const pressure for dry air
LvB = 2.555E6 #inflated latent heat of vaporization
r1 = 100.E3 #hypothetical inner-core storm radius (m)
r2 = 300.E3 #hypothetical outer radius (m) within which environmental air can easily intrude
sd = 300.E3 #Radius (m) to use in calculation of pressure centroid

# Establish starting and ending times of simulations:
time_start = DT.datetime(2014,10,13,12,0) #Format: year, month, day, hour, minute
loop_h = 72 # This is the time to loop forward for in hours
num_hours = loop_h+1 # Since Python is not inclusive

# Establish directory path of ensemble output:
path = '/network/rit/lab/tcdynasty/fischer/Data/Modeling/RI/Gonzalo/'

# Establish domain of inner nest to calculate centroids with:
sel_domain = 'd02' #This is the 12-km nest

# Read in TC center locations:
cent_lon = np.genfromtxt('/network/rit/lab/tcdynasty/fischer/Data/Modeling/RI/Gonzalo/Tracking/wrf_tc_center_lons.txt',delimiter=',')
cent_lat = np.genfromtxt('/network/rit/lab/tcdynasty/fischer/Data/Modeling/RI/Gonzalo/Tracking/wrf_tc_center_lats.txt',delimiter=',')

# Establish pressure levels for ventilation calculations:
ptop = 10000 #Top tropospheric level (Pa) to read in
phi = 20000 #Top pressure level (Pa)
plo = 85000 #Bottom pressure level (Pa)
pmid = 60000 #Pressure level to calculate midlevel entropy deficit

# Filter radius around TCs/INVESTs:
pifilt = 2. #for potential intensity (degrees)
envfilt = 8 #for shear (degrees)
shearrad = 500.0E3 #radius for inversion (km)

# Set up polar interpolation grid for area-averaged entropy calculations:
th = np.arange(0,2*np.pi,np.pi/4)
radi = np.arange(50.E3,r1+1,50.E3)
rado = np.arange(r1,r2+1,100.E3)
xi = np.zeros(len(th)*len(radi))
yi = np.zeros(len(th)*len(radi))
ri = np.zeros(len(th)*len(radi))
xo = np.zeros(len(th)*len(rado))
yo = np.zeros(len(th)*len(rado))
ro = np.zeros(len(th)*len(rado))
c1 = 0
c2 = 0
for i in xrange(len(th)):
    for j in xrange(len(radi)):
        xi[c1] = radi[j]*np.cos(th[i])/110.E3
        yi[c1] = radi[j]*np.sin(th[i])/110.E3
        ri[c1] = radi[j]
        c1 = c1+1
    for k in xrange(len(rado)):
        xo[c2] = rado[k]*np.cos(th[i])/110.E3
        yo[c2] = rado[k]*np.sin(th[i])/110.E3
        ro[c2] = rado[k]
        c2 = c2+1

# Get sample pressure levels from ERA-Interim:
dummyfile = '/network/daes/eraint/2005/u.2005.nc'

fin = NC4.Dataset(dummyfile,'r')

lev_hpa = np.array(fin.variables['lev'])
lev = lev_hpa*100. #Convert to Pa
fin.close()

# Find indices of vertical domain layers/levels
pmidi = np.where(lev==pmid)[0][0]
ptopi = np.where(lev==ptop)[0][0]+1
phii = np.where(lev==phi)[0][0]
ploi = np.where(lev==plo)[0][0]

# Create array of zeros for vi, shear, xdef, and mpi:
vi     = np.zeros((num_ens,num_hours)) #Ventilation Index
shear  = np.zeros_like(vi) #Vertical wind shear
ushear = np.zeros_like(vi) #U-component of shear
vshear = np.zeros_like(vi) #V-component of shear
chim   = np.zeros_like(vi) #Entropy deficit
mpi    = np.zeros_like(vi) #Maximum potential intensity

# Initialize arrays with NaNs:
vi[:,:]     = np.nan
shear[:,:]  = np.nan
ushear[:,:] = np.nan
vshear[:,:] = np.nan
chim[:,:]   = np.nan
mpi[:,:]    = np.nan

# Begin loop:
for eni in xrange(num_ens): #Current ensemble index
    curr_ens = str(eni+1).zfill(3) #Current ensemble number, adding 1 since first ensemble is 'e001', not 'e000'
    sel_ens = 'e'+curr_ens+'/' #String of current ensemble member's directory
    print 'Current ensemble member is',sel_ens
    for ts in xrange(num_hours): #Current forecast hour
    	print 'Current forecast hour is:', ts

        curr_time = time_start + DT.timedelta(hours=ts) #The current time in the loop
        curr_time_str = curr_time.strftime('%Y-%m-%d_%H') #String with format: year, month, day, hour

        # Note, some of the ensemble outputs aren't to the exact second,
        # So let's use the "glob" package to search for the output file
        # With the same minute... This method probably won't work if
        # The WRF output is with a frequency > 1 file per minute:
        for file in glob.glob(path+sel_ens+'wrfout_'+sel_domain+'_'+curr_time_str+':00:*'):
            datafile = file

        fin = NC4.Dataset(datafile) #Read in the NetCDF file

        latreg = wrf.getvar(fin, "latitude", meta=False) #Latitude
        lonreg = wrf.getvar(fin, "longitude", meta=False) #Longitude
        tmpsfc = np.array(fin.variables['SST'][0,:,:]) #Sea-surface temperature
        pmsl = wrf.getvar(fin, "slp", meta=False) #Mean sea-level pressure

        # Get meteorological winds at shear levels:
        uv_met = wrf.getvar(fin, "uvmet") #Meteorological winds
        p = wrf.getvar(fin, "pressure") #Pressure levels
        uv_top = wrf.interplevel(uv_met, p, phi/100.) #Meteorological winds at top shear level (hPa)
        uv_bot = wrf.interplevel(uv_met, p, plo/100.) #Meteorological winds at bottom shear level (hPa)
        u = np.zeros((2,uv_top.shape[-2],uv_top.shape[-1]))
        u[0,:,:] = uv_bot[0,:,:]
        u[1,:,:] = uv_top[0,:,:]

        v = np.zeros_like(u)
        v[0,:,:] = uv_bot[1,:,:]
        v[1,:,:] = uv_top[1,:,:]
        
        # Get temperature and specific humidity at multiple pressure levels
        tmpprs = np.zeros((np.size(lev),np.size(latreg),np.size(lonreg)))
        rhprs = np.zeros_like(tmpprs)

        for levi in xrange(tmpprs.shape[0]):
        	t_temp = wrf.getvar(fin, "tk", meta=False) #Temperature (K)
        	tmpprs[levi,:,:] = wrf.interplevel(t_temp, p, lev_hpa[levi]) #Temperature (K) at specified level (hPa)

        	rh_temp = (1./100.)*wrf.getvar(fin, "rh", meta=False) #Relative humidity (convert from %)
        	rhprs[levi,:,:] = wrf.interplevel(rh_temp, p, lev_hpa[levi]) #Relative humidity (dimensionless) at specified level (hPa)

        LAM,PHI = lonreg*np.pi/180.0,latreg*np.pi/180.0 #Convert longitude and latitude to radians
        LAM3d = np.tile(LAM,(np.size(lev[0:ptopi]),1,1))
	PHI3d = np.tile(PHI,(np.size(lev[0:ptopi]),1,1))
	
	print np.shape(LAM3d),np.shape(PHI3d)

        # Calculate potential intensity:
        print 'Calculating potential intensity'
        pp = np.transpose(np.tile(lev[0:ptopi],(np.shape(rhprs)[1],np.shape(rhprs)[2],1)),(2,0,1))
        estar = 611.2*np.exp(17.67*(tmpprs-273.15)/(tmpprs-273.15+243.5))
        vappres = rhprs*estar
        rv = eps*vappres/(pp-vappres)
        vmax_nofilt = np.zeros(tmpsfc.shape) #Unfiltered MPI
        airsea_nofilt = np.zeros(tmpsfc.shape) #Unfiltered air-sea disequilibrium
        for i in xrange(len(latreg)):
            for j in xrange(len(lonreg)):
            	(pmin,vmax_nofilt[i,j],airsea_nofilt[i,j],ifl) = pcmin.pcmin(tmpsfc[i,j]-273.15,pmsl[i,j]/100.00,lev[0:ptopi]/100.00,tmpprs[:,i,j]-273.15,rv[:,i,j]*1000.0,ptopi,ptopi)
            
        # Calculate distance of each grid point from the TC center:
        dist = 2*Ae*np.arcsin(np.sqrt((np.sin((PHI-cent_lat[eni,ts]*np.pi/180)/2))**2 + \
	    	   np.cos(cent_lat[eni,ts]*np.pi/180)*np.cos(PHI)*(np.sin((LAM-cent_lon[eni,ts]*np.pi/180)/2))**2))
        
        # Remove area around TC/INVEST within specified distance (sd), and fill rest with NaNs:
        xlonlim = np.where(dist<=sd,np.nan,lonreg) #Longitudes where distance criterion is satisfied
        xlatlim = np.where(dist<=sd,np.nan,latreg) #Latitudes where distance criterion is satisfied
        vmax    = np.where(dist<=sd,np.nan,vmax_nofilt) #MPI where distance criterion is satisfied
        airsea  = np.where(dist<=sd,np.nan,airsea_nofilt) #Air-sea disequilibrium where distance criterion is satisfied
        
        # Interpolate removed area back in
        lonnonan = np.extract(np.isfinite(vmax),lonreg)
        latnonan = np.extract(np.isfinite(vmax),latreg)
        vmaxnonan = np.extract(np.isfinite(vmax),vmax)
        airseanonan = np.extract(np.isfinite(vmax),airsea)
        
        F1 = ninterp(np.column_stack((lonnonan,latnonan)),vmaxnonan)
        F2 = ninterp(np.column_stack((lonnonan,latnonan)),airseanonan)
        
        # Extract PI value at the TC location:
        mpi[eni,ts] = F1([cent_lon[eni,ts],cent_lat[eni,ts]])
        airsea = F2([cent_lon[eni,ts],cent_lat[eni,ts]])
        
        print 'MPI:',mpi[eni,ts]
        
        # Calculate relative vorticity:
	dvdx    = (1./(Ae*np.cos(PHI3d)))*np.gradient(v,LAM,axis=2) #Zonal gradient of meridional winds
	dudy    = (1./Ae)*np.gradient(u,PHI,axis=1) #Meridional gradient of zonal winds
	vort_s  = u*np.tan(PHI3d)/Ae #Additional term which arises from conversion to spherical coordinates
	relvort = dvdx - dudy + vort_s

	# Calculate divergence:
	dudx  = (1./(Ae*np.cos(PHI3d)))*np.gradient(u,LAM,axis=2) #Zonal gradient of zonal winds
	dvdy  = (1./Ae)*np.gradient(v,PHI,axis=1) #Meridional gradient of meridional winds
	divg_s = v*np.tan(PHI3d)/Ae #Additional term which arises from conversion to spherical coordinates
	diverg = dudx + dvdy - divg_s
            
        #Calculate differential divergence and vorticity, and associated differential velocity potential and streamfunction
        diffdivg = np.where(dist<=shearrad,diverg[1,:,:]-diverg[0,:,:],0)
        diffvort = np.where(dist<=shearrad,relvort[1,:,:]-relvort[0,:,:],0)
        
        #Set up inversion operator
        Q = np.zeros((np.size(diffvort),np.size(diffvort)))
        c = 0
        for i in xrange(np.shape(diffvort)[0]):
            for j in xrange(np.shape(diffvort)[1]):
                tmp = np.zeros(np.shape(diffvort))
                if (j>0):
                    dlam = LAM[i,j] - LAM[i,j-1]
                    tmp[i,j-1] = 1.0/(dlam*np.cos(PHI[i,j]))**2
                if (j<(np.shape(diffvort)[1]-1)):
                    dlam = LAM[i,j+1] - LAM[i,j]
                    tmp[i,j+1] = 1.0/(dlam*np.cos(PHI[i,j]))**2
                if (i>0):
                    dphi = PHI[i,j] - PHI[i-1,j]
                    tmp[i-1,j] = 1.0/dphi**2 - np.tan(PHI[i,j])/(2*dphi)
                if (i<(np.shape(diffvort)[0]-1)):
                    dphi = PHI[i+1,j] - PHI[i,j]
                    tmp[i+1,j] = 1.0/dphi**2 + np.tan(PHI[i,j])/(2*dphi)
                tmp[i,j] = tmp[i,j] - 2.0/(dlam*np.cos(PHI[i,j]))**2 - 2.0/dphi**2
                
                Q[c,:] = np.reshape(tmp,(1,-1))
                c = c+1
                
        vzeta = Ae**2*np.reshape(diffvort,(-1,1))
        vchi = Ae**2*np.reshape(diffdivg,(-1,1))
        
        # Perform inversion:
        diffpsi = np.reshape(np.linalg.solve(Q,vzeta),np.shape(diffvort))
        diffchi = np.reshape(np.linalg.solve(Q,vchi),np.shape(diffdivg))
        
        # Differential nondivergent wind:
        deltaupsi = -1.0/Ae*(diffpsi-diffpsi)/(2.0*dphi)
        deltavpsi = 1.0/(Ae*np.cos(PHI))*(diffpsi-diffpsi)/(2.0*dlam)
        
        # Differential rotational wind:
        deltauchi = 1.0/(Ae*np.cos(PHI))*(diffchi-diffchi)/(2.0*dlam)
        deltavchi = 1.0/Ae*(diffchi-diffchi)/(2.0*dphi)
        
        # Environmental vertical wind shear:
        deltauenv = u[1,:,:]-u[0,:,:]-deltaupsi-deltauchi
        deltavenv = v[1,:,:]-v[0,:,:]-deltavpsi-deltavchi
            
        # Extract shear at TC location:
        Fshear = rbs(lonreg,latreg,np.transpose((deltauenv**2+deltavenv**2)**0.5))
        shear[eni,ts] = Fshear.ev(cent_lon[eni,ts],cent_lat[eni,ts])
        Fushear = rbs(lonreg,latreg,np.transpose(deltauenv))
        ushear[eni,ts] = Fushear.ev(cent_lon[eni,ts],cent_lat[eni,ts])
        Fvshear = rbs(lonreg,latreg,np.transpose(deltavenv))
        vshear[eni,ts] = Fvshear.ev(cent_lon[eni,ts],cent_lat[eni,ts])
        
        # Calculate entropy deficit:
        tmid = tmpprs[pmidi,:,:]
        
        # s*_m:
        rvstarm = eps*estar[pmidi,:,:]/(pmid-estar[pmidi,:,:])
        pd = pmid/(1.0+rvstarm/eps)
        sstarm = cpd*np.log(tmid)-Rd*np.log(pd)+LvB*rvstarm/tmid
        
        # Hypothetical inner-core storm average:
        F = rbs(lonreg, latreg, np.transpose(sstarm))
        disclon = cent_lon[eni,ts]+xi
        disclat = cent_lat[eni,ts]+yi
        sstarminner = F.ev(disclon,disclat)
        sstarminnerbar = np.nansum(sstarminner*ri)/np.sum((np.isfinite(sstarminner))*ri)
                
        # s_m:
        rvmid = rv[pmidi,:,:]
        rhmid = rhprs[pmidi,:,:]
        sm = cpd*np.log(tmid)-Rd*np.log(pd)+LvB*rvmid/tmid-Rv*rvmid*np.log(np.maximum(rhmid,0.01))
        
        # Hypothetical environmental average:
        F = rbs(lonreg, latreg, np.transpose(sm))
        disclon = cent_lon[eni,ts]+xo
        disclat = cent_lat[eni,ts]+yo
        smouter = F.ev(disclon,disclat)
        smouterbar = np.nansum(smouter*ro)/np.sum((np.isfinite(smouter))*ro)
        
        # Calculate entropy deficit at TC location:
        if (np.isfinite(airsea)):
            chim[eni,ts] = np.maximum(sstarminnerbar-smouterbar,0.0)/np.maximum(airsea,1.0)
        else:
            chim[eni,ts] = -99
            
        # Calculate ventilation index:
        if (np.isfinite(airsea)):
            vi[eni,ts] = shear[eni,ts]*chim[eni,ts]/np.maximum(mpi[eni,ts],1.0)
        else:
            vi[eni,ts] = -99
            
                
        # Close files to prevent memory leaks"
        fin.close()






### END ###