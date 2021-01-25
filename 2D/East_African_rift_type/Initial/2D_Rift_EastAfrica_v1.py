#!/usr/bin/env python
# coding: utf-8

# # Tutorial 11 : Coupling a thermo-mechanical model with Badlands (2D)
# 
# This Tutorial requires pyBadlands.
# You can run it using the docker images provided
# 
# ![Tutorial11](images/Tutorial11.gif)

# In[55]:




# In[56]:


import UWGeodynamics as GEO


# In[57]:


u = GEO.UnitRegistry


# ## Scaling

# In[58]:


half_rate = 1.8 * u.centimeter / u.year
model_length = 360e3 * u.meter
surfaceTemp = 273.15 * u.degK
baseModelTemp = 1603.15 * u.degK
bodyforce = 3300 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2

#interfaceX = 600 * u.kilometre
interfaceX = 180 * u.kilometre
UpperCrustHL = -30 * u.kilometre
UpperCrustHR = -22 * u.kilometre
CrustHL = -44 * u.kilometre
CrustHR = -34 * u.kilometre
LithL = -150 * u.kilometre
LithR = -130 * u.kilometre

KL = model_length
Kt = KL / half_rate
KM = bodyforce * KL**2 * Kt**2
KT = (baseModelTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT


# # UW Model

# In[59]:


Model = GEO.Model(elementRes=(192, 96), 
                  minCoord=(0. * u.kilometer, -170. * u.kilometer), 
                  maxCoord=(360. * u.kilometer, 10. * u.kilometer), 
                  gravity=(0.0, -9.81 * u.meter / u.second**2))


# In[60]:


Model.outputDir="output_tutorialCoupling2D"


# ## Materials

# In[61]:


Model.diffusivity = 1e-6 * u.metre**2 / u.second 
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)


# The Model we are building is essentially a layered cake. The geodynamics module provide and easy way to define a layer by defining shape as *layer* and specifying its *top* and *bottom*. The order is important: when 2 shapes overlay each other, only the second is used.

# In[62]:


air = Model.add_material(name="Air", shape=GEO.shapes.Layer(top=Model.top, bottom=2.0 * u.kilometer))
stickyAir = Model.add_material(name="StickyAir", shape=GEO.shapes.Layer(top=air.bottom, bottom= 0.0 * u.kilometer))

#uppercrust = Model.add_material(name="UppperCrust", shape=GEO.shapes.Layer(top=stickyAir.bottom, bottom=-35.0 * u.kilometer))
#mantleLithosphere = Model.add_material(name="MantleLithosphere", shape=GEO.shapes.Layer(top=uppercrust.bottom, bottom=-100.0 * u.kilometer))
#mantle = Model.add_material(name="Mantle", shape=GEO.shapes.Layer(top=mantleLithosphere.bottom, bottom=Model.bottom))
#sediment = Model.add_material(name="Sediment")


# In[63]:


#Left Box Definitions
uppercrust_Geo_L = GEO.shapes.Box(top=0.* u.kilometre, bottom=UpperCrustHL, minX=0.*u.kilometre, maxX=interfaceX)
uppercrust_L = Model.add_material(name="Uppercrust_L", shape=uppercrust_Geo_L)

lowercrust_Geo_L = GEO.shapes.Box(top=uppercrust_L.bottom, bottom=CrustHL, minX=0.*u.kilometre, maxX=interfaceX)
lowercrust_L = Model.add_material(name="Lowercrust_L", shape=lowercrust_Geo_L)

mantleLithosphere_Geo_L = GEO.shapes.Box(top=lowercrust_L.bottom, bottom=LithL, minX=0.*u.kilometre, maxX=interfaceX)
mantleLithosphere_L = Model.add_material(name="MantleLithosphere_L", shape=mantleLithosphere_Geo_L)

sediment = Model.add_material(name="Sediment")


# In[64]:


#Right Box Definitions
uppercrust_Geo_R = GEO.shapes.Box(top=stickyAir.bottom, bottom=UpperCrustHR, minX=interfaceX, maxX=model_length)
uppercrust_R = Model.add_material(name="Uppercrust_R", shape=uppercrust_Geo_R)

lowercrust_Geo_R = GEO.shapes.Box(top=uppercrust_R.bottom, bottom=CrustHR, minX=interfaceX, maxX=model_length)
lowercrust_R = Model.add_material(name="Lowercrust_R", shape=lowercrust_Geo_R)

mantleLithosphere_Geo_R = GEO.shapes.Box(top=lowercrust_R.bottom, bottom=LithR, minX=interfaceX, maxX=model_length)
mantleLithosphere_R = Model.add_material(name="MantleLithosphere_R", shape=mantleLithosphere_Geo_R)

mantleAsthenosphere_Geo_R = GEO.shapes.Box(top=mantleLithosphere_R.bottom, bottom=mantleLithosphere_L.bottom, minX=interfaceX, maxX=model_length)
mantleAsthenosphere_R = Model.add_material(name="MantleAsthenosphere_R", shape=mantleAsthenosphere_Geo_R)


# In[65]:


mantleAsthenosphere_Geo_comb = GEO.shapes.Layer(top=mantleAsthenosphere_R.bottom, bottom=Model.bottom)
mantleAsthenosphere_comb = Model.add_material(name="MantleAsthenosphere_comb", shape=mantleAsthenosphere_Geo_comb)


# In[66]:


air.diffusivity = 1.0e-6 * u.metre**2 / u.second
stickyAir.diffusivity = 1.0e-6 * u.metre**2 / u.second

air.capacity = 100. * u.joule / (u.kelvin * u.kilogram)
stickyAir.capacity = 100. * u.joule / (u.kelvin * u.kilogram)


# In[67]:


air.density = 1. * u.kilogram / u.metre**3
stickyAir.density = 1. * u.kilogram / u.metre**3

uppercrust_L.density = GEO.LinearDensity(2620. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
uppercrust_R.density = GEO.LinearDensity(2620. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

lowercrust_L.density = GEO.LinearDensity(2620. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
lowercrust_R.density = GEO.LinearDensity(2620. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

mantleLithosphere_L.density = GEO.LinearDensity(3370. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
mantleLithosphere_R.density = GEO.LinearDensity(3370. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

mantleAsthenosphere_R.density = GEO.LinearDensity(3395. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
mantleAsthenosphere_comb.density = GEO.LinearDensity(3395. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

sediment.density = 2300. * u.kilogram / u.metre**3


# In[ ]:





# In[68]:


uppercrust_L.radiogenicHeatProd = 0.7 * u.microwatt / u.meter**3
uppercrust_R.radiogenicHeatProd = 0.7 * u.microwatt / u.meter**3

lowercrust_L.radiogenicHeatProd = 0.7 * u.microwatt / u.meter**3
lowercrust_R.radiogenicHeatProd = 0.7 * u.microwatt / u.meter**3

mantleLithosphere_L.radiogenicHeatProd = 0.02e-6 * u.microwatt / u.meter**3
mantleLithosphere_R.radiogenicHeatProd = 0.02e-6 * u.microwatt / u.meter**3

sediment.radiogenicHeatProd = 0.6 * u.microwatt / u.meter**3


# In[69]:


rh = GEO.ViscousCreepRegistry()


# In[70]:


#air.viscosity                = 1e19 * u.pascal * u.second
#stickyAir.viscosity          = 1e20 * u.pascal * u.second
#uppercrust.viscosity         = 5.0 * rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995
#mantleLithosphere.viscosity  = 5.0 * rh.Dry_Olivine_Dislocation_Karato_and_Wu_1990
#mantle.viscosity             = 1.0 * rh.Dry_Olivine_Dislocation_Karato_and_Wu_1990
#sediment.viscosity           = rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995


# In[71]:


air.viscosity                = 1e19 * u.pascal * u.second
stickyAir.viscosity          = 1e20 * u.pascal * u.second

uppercrust_L.viscosity         = 10.0 * rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995
uppercrust_R.viscosity         = 1.0 * rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995

lowercrust_L.viscosity         = 5.0 * rh.Wet_Anorthite_Dislocation_Ribacki_et_al_2000
lowercrust_R.viscosity         = 1.0 * rh.Wet_Anorthite_Dislocation_Ribacki_et_al_2000

mantleLithosphere_L.viscosity  = 1.0 * rh.Dry_Olivine_Dislocation_Karato_and_Wu_1990
mantleLithosphere_R.viscosity  = 1.0 * rh.Dry_Olivine_Dislocation_Karato_and_Wu_1990

mantleAsthenosphere_R.viscosity             = 1.0 * rh.Dry_Olivine_Dislocation_Karato_and_Wu_1990
mantleAsthenosphere_comb.viscosity             = 1.0 * rh.Dry_Olivine_Dislocation_Karato_and_Wu_1990

sediment.viscosity           = rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995


# In[72]:


#air.viscosity                = 1e19 * u.pascal * u.second
#stickyAir.viscosity          = 1e20 * u.pascal * u.second

#uppercrust_L.viscosity          = 1e22 * u.pascal * u.second
#lowercrust_L.viscosity         = 5e21 * u.pascal * u.second
#mantleLithosphere_L.viscosity  = 1e22 * u.pascal * u.second

#uppercrust_R.viscosity          = 1e22 * u.pascal * u.second
#lowercrust_R.viscosity         = 5e21 * u.pascal * u.second
#mantleLithosphere_R.viscosity  = 1e22 * u.pascal * u.second
#mantleAsthenosphere_R.viscosity             = 5e20 * u.pascal * u.second

#mantleAsthenosphere_comb.viscosity             = 5e20 * u.pascal * u.second

#sediment.viscosity         = 1e22 * u.pascal * u.second


# In[73]:


pl = GEO.PlasticityRegistry()


# In[74]:


uppercrust_L.plasticity         = pl.Huismans_et_al_2011_Crust
uppercrust_R.plasticity         = pl.Huismans_et_al_2011_Crust

lowercrust_L.plasticity         = pl.Huismans_et_al_2011_Crust
lowercrust_R.plasticity         = pl.Huismans_et_al_2011_Crust

mantleLithosphere_L.plasticity  = pl.Huismans_et_al_2011_Crust
mantleLithosphere_R.plasticity  = pl.Huismans_et_al_2011_Crust

mantleAsthenosphere_R.plasticity             = pl.Huismans_et_al_2011_Crust
mantleAsthenosphere_comb.plasticity             = pl.Huismans_et_al_2011_Crust

sediment.plasticity           = pl.Huismans_et_al_2011_Crust


# In[75]:


sediment


# ## Boundary Conditions

# In[76]:


Model.set_temperatureBCs(top=293.15 * u.degK, 
                         bottom=1603.15 * u.degK, 
                         materials=[(air, 293.15 * u.degK), (stickyAir, 293.15 * u.degK)])

Model.set_velocityBCs(left=[-2.5 * u.centimeter / u.year, None],
                       right=[2.5 * u.centimeter / u.year, None],
                       bottom=GEO.LecodeIsostasy(reference_mat=mantleAsthenosphere_comb, average=True))


# ## Initial Damage

# In[77]:


import numpy as np

def gaussian(xx, centre, width):
    return ( np.exp( -(xx - centre)**2 / width ))

maxDamage = 0.7
Model.plasticStrain.data[:] = maxDamage * np.random.rand(*Model.plasticStrain.data.shape[:])
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,0], (GEO.nd(Model.maxCoord[0] - Model.minCoord[0])) / 2.0, GEO.nd(5.0 * u.kilometer))
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,1], GEO.nd(-35. * u.kilometer) , GEO.nd(5.0 * u.kilometer))



# In[80]:


Model.init_model()



# # Badlands

# In[23]:


Model.surfaceProcesses = GEO.surfaceProcesses.Badlands(airIndex=[stickyAir.index, air.index], sedimentIndex=sediment.index,
                                          XML="ressources/badlands.xml", resolution=1. * u.kilometer, 
                                          checkpoint_interval=0.01 * u.megayears)


# ## Run Model

# In[24]:


Model.solver.set_inner_method("mumps")
Model.solver.set_penalty(1e6)

# In[48]:


Model.run_for(5.0 * u.megayear, checkpoint_interval=0.01*u.megayear)


# In[ ]:




