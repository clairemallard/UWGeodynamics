#!/usr/bin/env python
# coding: utf-8

# ## Rift Model - no Badlands

# ### 20 degrees obliquity and 2cm/yr extension 

# In[1]:


import UWGeodynamics as GEO
#import glucifer
import badlands

# In[2]:


u = GEO.UnitRegistry


# In[3]:


# Characteristic values of the system - Scalling
half_rate = 1. * u.centimeter / u.year
model_length = 512e3 * u.meter #1024
model_width = 256e3 * u.meter
surfaceTemp = 273.15 * u.degK
baseModelTemp = 1523.15 * u.degK
bodyforce = 3300 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2 #3300

KL = model_length
Kt = KL / half_rate
KM = bodyforce * KL**2 * Kt**2
KT = (baseModelTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"] = KM
GEO.scaling_coefficients["[temperature]"] = KT


# In[4]:


#Model Box
resolution = (256, 128, 64)
Model = GEO.Model(elementRes=resolution,
                  minCoord=(0. * u.kilometer, 0. * u.kilometer, -98. * u.kilometer),
                  maxCoord=(512. * u.kilometer, 256. * u.kilometer, 30. * u.kilometer),
                  gravity=(0.0, 0.0, -9.81 * u.meter / u.second**2))


# In[5]:


#Model Global parameters
Model.minViscosity  = 1e18 * u.pascal * u.second
Model.maxViscosity  = 1e23 * u.pascal * u.second
Model.stressLimiter = 300. * u.megapascal
Model.diffusivity = 1e-6 * u.metre**2 / u.second
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)
#Output Folder
Model.outputDir="outputs_Rift_20_008_Bls"


# In[6]:


# ## Define Materials
#
air_shape         = GEO.shapes.Layer3D(top=Model.top, bottom=10.0 * u.kilometer)
air               = Model.add_material(name="Air", shape=air_shape)
stickyAir_shape   = GEO.shapes.Layer3D(top=air.bottom, bottom= 0.0 * u.kilometer)
stickyAir         = Model.add_material(name="StickyAir", shape=stickyAir_shape)

Ucrust = Model.add_material(name="UpperCrust", shape=GEO.shapes.Layer3D(top=stickyAir.bottom, bottom=-20 * u.kilometer))
Lcrust = Model.add_material(name="LowerCrust", shape=GEO.shapes.Layer3D(top=Ucrust.bottom, bottom=-40 * u.kilometer))
mantleLithosphere = Model.add_material(name="MantleLithosphere", shape=GEO.shapes.Layer3D(top=Lcrust.bottom, bottom=-90 * u.kilometer))
mantle_shape   = GEO.shapes.Layer3D(top=mantleLithosphere.bottom, bottom=Model.bottom)
mantle = Model.add_material(name="Mantle", shape=mantle_shape)
sediment   = Model.add_material(name="Sediment")


# In[7]:


# Notche
top    = GEO.shapes.HalfSpace(normal=(0.,0.,1.), origin=(256.*u.kilometer,144.*u.kilometer,-40.*u.kilometer))
right  = GEO.shapes.HalfSpace(normal=(1.,-0.3,-0.7), origin=(262.*u.kilometer,128.*u.kilometer,-43.*u.kilometer))
left   = GEO.shapes.HalfSpace(normal=(-1.,0.3,-0.7), origin=(250.*u.kilometer,128.*u.kilometer,-43.*u.kilometer))
CompositeShape_notche = left & top & right
notche = Model.add_material(name="notche", shape=CompositeShape_notche)



air.capacity       = 1000. * u.joule / (u.kelvin * u.kilogram)
stickyAir.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)


# In[8]:


## Densities

air.density       = 1. * u.kilogram / u.metre**3
stickyAir.density = 1. * u.kilogram / u.metre**3

Ucrust.density    = GEO.LinearDensity(reference_density=2800. * u.kilogram / u.metre**3)
Lcrust.density    = GEO.LinearDensity(reference_density=2900. * u.kilogram / u.metre**3)
notche.density    = Lcrust.density
mantleLithosphere.density  = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3)
mantle.density     = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3)
sediment.density = GEO.LinearDensity(reference_density=2700. * u.kilogram / u.metre**3)


# In[9]:


#Thermal global parameters
Ucrust.radiogenicHeatProd = 1.2 * u.microwatt / u.meter**3 #0.7
Lcrust.radiogenicHeatProd = 0.6 * u.microwatt / u.meter**3 #0.4
notche.radiogenicHeatProd = Lcrust.radiogenicHeatProd
sediment.radiogenicHeatProd = 1.2 * u.microwatt / u.meter**3 #0.7
mantleLithosphere.radiogenicHeatProd = 0.01154 * u.microwatt / u.meter**3
mantle.radiogenicHeatProd = 0.01154 * u.microwatt / u.meter**3


# In[10]:


#Viscosities
rh = GEO.ViscousCreepRegistry()

air.viscosity                = 1e18 * u.pascal * u.second
stickyAir.viscosity          = 1e19 * u.pascal * u.second

Ucrust.viscosity              = rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
Ucrust.minViscosity = 1e19 * u.pascal * u.second
Ucrust.maxViscosity = 1e23 * u.pascal * u.second

Lcrust.viscosity              = rh.Dry_Mafic_Granulite_Dislocation_Wang_et_al_2012
Lcrust.minViscosity = 1e19 * u.pascal * u.second
Lcrust.maxViscosity = 1e23 * u.pascal * u.second

notche.viscosity = rh.Dry_Mafic_Granulite_Dislocation_Wang_et_al_2012
notche.minViscosity = 1e19 * u.pascal * u.second
notche.maxViscosity = 1e23 * u.pascal * u.second

mantleLithosphere.viscosity  = rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003
mantleLithosphere.minViscosity = 5e20 * u.pascal * u.second
mantleLithosphere.maxViscosity = 1e23 * u.pascal * u.second

mantle.viscosity             = rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003
mantle.minViscosity = 1e20 * u.pascal * u.second
mantle.maxViscosity = 1e23 * u.pascal * u.second

sediment.viscosity           = rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
sediment.minViscosity = 1e19 * u.pascal * u.second
sediment.maxViscosity = 1e23 * u.pascal * u.second


# In[11]:


#Plasticities
pl = GEO.PlasticityRegistry()

Ucrust.plasticity = pl.Rey_and_Muller_2010_UpperCrust
Ucrust.stressLimiter  = 300. * u.megapascal
Ucrust.plasticity.frictionCoefficient = 0.55
Ucrust.plasticity.frictionAfterSoftening = 0.055

Lcrust.plasticity = pl.Rey_and_Muller_2010_LowerCrust
Lcrust.plasticity.frictionAfterSoftening = 0.028
Lcrust.stressLimiter  = 300. * u.megapascal

notche.plasticity = pl.Rey_and_Muller_2010_LowerCrust
notche.plasticity.frictionAfterSoftening = 0.028
notche.stressLimiter  = 150. * u.megapascal

mantleLithosphere.plasticity = pl.Rey_and_Muller_2010_LithosphericMantle
mantleLithosphere.stressLimiter  = 300. * u.megapascal

mantle.plasticity = pl.Rey_and_Muller_2010_LithosphericMantle
mantle.stressLimiter  = 350. * u.megapascal

sediment.plasticity = pl.Rey_and_Muller_2010_UpperCrust
sediment.stressLimiter  = 300. * u.megapascal
sediment.plasticity.frictionCoefficient = 0.55
sediment.plasticity.frictionAfterSoftening = 0.055


# In[12]:


# ## Passive Tracers

import numpy as np

xp = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), 100)
yp = np.linspace(GEO.nd(Model.minCoord[1]), GEO.nd(Model.maxCoord[1]), 100)

xp, yp = np.meshgrid(xp, yp)
xp = xp.flatten()
yp = yp.flatten()
zp = np.zeros(xp.shape)

#Merge into a 3D numpy array
surface_tracersNp=np.zeros(shape=(len(xp),3))
moho_tracersNp=np.zeros(shape=(len(xp),3))
#lab_tracersNp=np.zeros(shape=(len(xp),3))

for k in range(0,len(xp)):
    surface_tracersNp[k][0]=xp[k]
    surface_tracersNp[k][1]=yp[k]
    surface_tracersNp[k][2]=zp[k]
                        
    moho_tracersNp[k][0]=xp[k]
    moho_tracersNp[k][1]=yp[k]
    moho_tracersNp[k][2]=zp[k]+GEO.nd(mantleLithosphere.top)

surface_tracers = Model.add_passive_tracers(name="Surface", vertices=surface_tracersNp)
moho_tracers = Model.add_passive_tracers(name="Moho", vertices=moho_tracersNp)
#Lab_tracers = Model.add_passive_tracers(name="Lab", vertices=[xp,yp,-110.*u.kilometer])


# In[13]:


# ## Temperature Boundary Condition
Model.set_temperatureBCs(top=273.15 * u.degK,
                         nodeSets = [(air_shape, 273.15 * u.degK), (stickyAir_shape, 273.15 * u.degK),(mantle_shape, 1523.15 * u.degK)],
                         bottom=1523.15 * u.degK)



#air.compressibility = 1e4  # Not sure what's a good value is

#P, bottomPress = Model.get_lithostatic_pressureField()


# Get the average of the pressure along the bottom, and make it a pressure BC along the bottom
#bottomPress = GEO.Dimensionalize(numpy.average(bottomPress), u.megapascal)


# In[14]:


# ## Velocity Boundary Conditions
#
# We pull on the left and right side. The back and front wall are freeslip. We use a pseudo isostatic support condition at the bottom.

# In[30]:



Model.set_velocityBCs(front=[None, 0., None],
                      back=[None, 0., None],
                      left=[-0.8* u.centimeter / u.year, 0., None],
                      right=[0.8 * u.centimeter / u.year, 0., None],
                      bottom=GEO.LecodeIsostasy(reference_mat=mantle,
                                                average=False),
                      order_wall_conditions=["front", "back", "left", "right", "bottom", "top"])



# In[ ]:


#Model.set_stressBCs(bottom=[0.0, 0.0, 3.635e9 * u.pascal])


# This would avoid that stresses are concentrated in the walls after inicializing the model using the notche

# In[15]:



#
# # ## Initial Plastic Strain
#
# import numpy as np
#
# def gaussian(xx, centre, width):
#     return ( np.exp( -(xx - centre)**2 / width ))
#
# modelLength = 0.768
# modelWidth = 0.384
#
# tanm2 = (0.36*modelWidth)/2 # tangente de l'angle qui est en radian!
# tanScale = tanm2/modelLength
#
# maxDamage = 0.2577
# Model.plasticStrain.data[:] = maxDamage * np.random.rand(*Model.plasticStrain.data.shape[:])
# angle = 0.35 # 20 degree
# centers = angle * Model.swarm.particleCoordinates.data[:,1] + (((GEO.nd(Model.maxCoord[0]) / 2.0)-tanScale))
#
# Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,0], centers , GEO.nd(4. * u.kilometer))
# Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,2], GEO.nd(-35. * u.kilometer) , GEO.nd(4. * u.kilometer))
#


# In[25]:


# # Surface Processes / BADLANDS

Model.surfaceProcesses = GEO.surfaceProcesses.Badlands(airIndex=[stickyAir.index, air.index],sedimentIndex=sediment.index, XML="./ressources_riftM/badlands-20_2.xml", resolution=1. * u.kilometer,checkpoint_interval=0.05 * u.megayears,restartFolder="outbdls",restartStep=498)


Model.init_model()


Model.set_temperatureBCs(top=273.15 * u.degK,
                         nodeSets = [(air_shape, 273.15 * u.degK), (stickyAir_shape, 273.15 * u.degK)],
                         bottom=1523.15 * u.degK)





import underworld.function as fn

def post_hook():
    coords = fn.input()
    zz = coords[0] / (GEO.nd(Model.maxCoord[0]) - GEO.nd(Model.minCoord[0]))
    fact = fn.math.pow(fn.math.tanh(zz*20.0) + fn.math.tanh((1.0-zz)*20.0) - fn.math.tanh(20.0), 4)
    Model.plasticStrain.data[:] = Model.plasticStrain.data[:] * fact.evaluate(Model.swarm)

Model.post_solve_functions["B"] = post_hook


# In[26]:


#Testing Julians suggestion
scr_rtol = 1e-6
# Schur complement solver options
Model.solver.options.scr.ksp_rtol = scr_rtol
# Inner solve (velocity), A11 options
Model.solver.options.A11.ksp_rtol = 1e-1 * scr_rtol

# In[ ]:


# In[ ]:


##Run verification script and got plots


# In[ ]:


Model.run_for(30. * u.megayears,checkpoint_interval=0.05 * u.megayears,restartStep=498)
import sys
sys.exit()

# In[ ]:


#Model.run_for(10. * u.megayears,restartStep=100,checkpoint_interval=0.05 * u.megayears)#restartStep=52,


# In[ ]:




