#!/usr/bin/env python
# coding: utf-8

# ## Rift Model - no Badlands

# ### 20 degrees obliquity and 2cm/yr extension

# In[1]:


import UWGeodynamics as GEO

# In[2]:
import badlands

u = GEO.UnitRegistry


# In[3]:


# Characteristic values of the system - Scalling
half_rate = 0.5 * u.centimeter / u.year
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
Model.outputDir="outputs_Rift_20_050_Bls"


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

air.viscosity    = 1e18 * u.pascal * u.second
stickyAir.viscosity = 1e19 * u.pascal * u.second

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

surface_tracersNp = np.zeros(shape=(len(xp),3))
moho_tracersNp = np.zeros(shape=(len(xp),3))

surface_tracersNp[:, 0] = xp
surface_tracersNp[:, 1] = yp
surface_tracersNp[:, 2] = zp

moho_tracersNp[:, 0] = xp
moho_tracersNp[:, 1] = yp
moho_tracersNp[:, 2] = zp + GEO.nd(mantleLithosphere.top)

surface_tracers = Model.add_passive_tracers(name="Surface", vertices=surface_tracersNp)
moho_tracers = Model.add_passive_tracers(name="Moho", vertices=moho_tracersNp)


# In[13]:


## Temperature Boundary Condition
#Model.set_temperatureBCs(top=273.15 * u.degK, bottom=1523.15 * u.degK)
Model.set_temperatureBCs(top=273.15 * u.degK,
nodeSets = [(air_shape, 273.15 * u.degK), (stickyAir_shape, 273.15 * u.degK),(mantle_shape, 1523.15 * u.degK)],
bottom=1523.15 * u.degK)


Model.set_velocityBCs(front=[None, 0., None],
                      back=[None, 0., None],
                      left=[-0.5* u.centimeter / u.year, 0., None],
                      right=[0.5 * u.centimeter / u.year, 0., None],
                      bottom=GEO.LecodeIsostasy(reference_mat=mantle,
                                                average=True),
                      order_wall_conditions=["front", "back", "left", "right", "bottom", "top"])

Model.surfaceProcesses = GEO.surfaceProcesses.Badlands(airIndex=[stickyAir.index, air.index], sedimentIndex=sediment.index,
XML="./ressources_riftM/badlands_20_05.xml", resolution=1. * u.kilometer,checkpoint_interval=0.05 * u.megayears,restartFolder="outbdls",restartStep=1)

Model.init_model()

#second Tempt BCs to fix temp in the mantle
## Temperature Boundary Condition
Model.set_temperatureBCs(top=273.15 * u.degK,
nodeSets = [(air_shape, 273.15 * u.degK), (stickyAir_shape, 273.15 * u.degK)],
bottom=1523.15 * u.degK)




import underworld.function as fn

def post_hook():
    coords = fn.input()
    zz = coords[0] / (GEO.nd(Model.maxCoord[0]) - GEO.nd(Model.minCoord[0]))
    fact = fn.math.pow(fn.math.tanh(zz*20.0) + fn.math.tanh((1.0-zz)*20.0) - fn.math.tanh(20.0), 4)
    Model.plasticStrain.data[:] = Model.plasticStrain.data[:] * fact.evaluate(Model.swarm)

#Model.postSolveHook = post_hook
Model.post_solve_functions["B"] = post_hook

#Testing for SOlver options - fragment from Calire original script
#Model.solver.options.set_inner_method = "mg"
#Model.solver.options.A11.mg_coarse_pc_factor_mat_solver_type = "mumps"

#Testing Julians suggestion
scr_rtol = 1e-6
# Schur complement solver options
Model.solver.options.scr.ksp_rtol = scr_rtol
# Inner solve (velocity), A11 options
Model.solver.options.A11.ksp_rtol = 1e-1 * scr_rtol

#Model.run_for(nstep=5, checkpoint_interval=1)

Model.run_for(40. * u.megayears,checkpoint_interval=0.05 * u.megayears,restartStep=1,restartDir="outputs_Rift_20_050_Bls")

