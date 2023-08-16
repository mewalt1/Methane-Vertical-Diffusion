#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('reset', '')


# In[2]:


# Lori Neary BIRA-IASB 12/03/2022; Load in files
import numpy as np
import scipy
import scipy.interpolate
import matplotlib.pyplot as plt
import h5py
import glob
import matplotlib.colors as colors
from cmcrameri import cm


# In[3]:


# get the list of files
filelist=sorted(glob.glob('./dust/*.h5'))   #take all the .h5 files from folder (GEM outputs)

lsstr=[]                                    #create empty arrays for each parameter
surf_temp=[]
loct=[]
csza=[]
temp=[]
height=[]	
pressure=[]
km=[]
windu=[]
windv=[]
ustar=[]

for file in filelist:                       
    print('read file: ',file)
    lsstr.append(file[-6:-3])

    fhdf = h5py.File(file,'r')
    
#surface temperature, local time, cosine of solar zenith angle (for each Ls file, 
#there are 48 values together representing 1 sol, starting around midnight)
    surf_temp.append(fhdf["/surf_temp"][:])
    loct.append(fhdf["/localtime"][:])   
    csza.append(fhdf["/csza"][:])
    
#profiles temperature, vertical diffusion coefficient, pressure and height 
#(for each Ls file, 48 profiles of 103 levels) note that level 0 is the top. 
    temp.append(fhdf["/temperature"][:])
    height.append(fhdf["/height_t"][:])
    pressure.append(fhdf["/pressure_t"][:])
    km.append(fhdf["/vdiff_coef"][:])
    windu.append(fhdf["/wind_u"][:])      #zonal wind speed
    windv.append(fhdf["/wind_v"][:])      #meridional wind speed
    ustar.append(fhdf["/ustar"][:])       #friction wind speed 




#surf_temp, loct, csza 36x48
#temp, height, pressure, km 36x48x103


# # Make files for fortran

# In[4]:


nrepetitions=1                              #number of model spin-up repetitions 
heights_new=np.arange(1,3002,1)
for n,file in enumerate(filelist):
    myfile = h5py.File(file,"r")
    filetxt = '\n'.join([" Input data for MARS methane model:",
    " Number of vertical levels (tl and in-program variables):",
    " 3001",
    " Levels temp,sigw data constant after level:",
    " 120",
    " Lat./Long./start Yr Mo D/T. Zone (hr from UTC) - Still Earth coords.",
    " 45.83050 121.94190 2014 8 1 8.0",
    " Model levels (m):"])+"\n"
   
    filetxt+= "\n".join(heights_new.astype(str))+"\n"
    filetxt+= "    ch4 initial value[ppb]: \n"
    filetxt += "1.0000 \n"

    for repetition in range(nrepetitions):
        for nstep in range(48):
           
            heights = myfile["height_t"][nstep,:][::-1]
            vdiff = myfile["vdiff_coef"][nstep,:][::-1]
            temperature = myfile["temperature"][nstep,:][::-1]

            #Interpolate diffusion
            vdiff_f=scipy.interpolate.interp1d(heights,vdiff)
            heights_new=np.arange(1,3002,1)
            vdiff_new = vdiff_f(heights_new)
            #Interpolate temperature
            temp_f=scipy.interpolate.interp1d(heights,temperature)
            heights_temp_new=np.arange(1,121,1)
            temp_new_values=temp_f(heights_temp_new)

            filetxt += "    endtim      pres[mbar]:    par[umol/m2/s]: csza[degrees]: \n"
            filetxt += "%f %f %f %f \n"%(myfile["localtime"][nstep],myfile["pressure_t"][nstep,-1]*1.0e-3,330, myfile["csza"][nstep])
            filetxt += "   temp[K] \n"
            filetxt += "\n".join(temp_new_values.astype(str)) + "\n"
            filetxt += "    K[m2/s] \n"
            filetxt += "\n".join(vdiff_new.astype(str)) + "\n"

   
    with open("K-only/input/czafullinput%02d.dat"%n,"w") as tmpfile:
        tmpfile.write(filetxt)
    myfile.close()
    print(str(n)+"\tInterpolated, wrote, and closed for "+file)


# # Output

# In[6]:


output = np.loadtxt("FileName/output.out") #Methane mixing ratio output has shape 48, 3001
output1 = np.loadtxt("FileName/ech4")
output2 = np.loadtxt("FileName/massdiff")

