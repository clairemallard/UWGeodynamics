
# coding: utf-8

# # 3D passive margin

# - configuration A: (128 x 128 x 64 nodes) 128 CPUS
# - configuration B: (256 x 256 x 96 nodes) 256 CPUS
# The user should remove visualisation from the python script before running the model on raijin.
# In [2]:
# 

# In[1]:


import matplotlib.pyplot as plt
import UWGeodynamics as GEO
u = GEO.UnitRegistry
import glucifer


# In[2]:


##### Characteristic values of the system

half_rate = 1.8 * u.centimeter / u.year # What is it?
model_length = 200e3 * u.meter
model_width = 100e3 * u.meter
surfaceTemp = 273.15 * u.degK
baseModelTemp = 1603.15 * u.degK
bodyforce = 3300 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2 # ckeck publi schellart

KL = model_length # Characteristic length
Kt = KL / half_rate # Characteristic time 
KM = bodyforce * KL**2 * Kt**2 # Characteristic mass
KT = (baseModelTemp - surfaceTemp) # Characteristic temperature

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT


# In[3]:


Model = GEO.Model(elementRes=(128, 64, 64), 
                  minCoord=(0. * u.kilometer, 0. * u.kilometer, -110. * u.kilometer), 
                  maxCoord=(200. * u.kilometer, 100. * u.kilometer, 10. * u.kilometer), 
                  gravity=(0.0 ,0.0, -9.81 * u.meter / u.second**2))


# ## Output

# In[4]:

Model.mesh_advector(axis=0)
Model.outputDir="3D_passive_margin"

# ## Limiters
# 

# In[6]:


Model.stressLimiter = 300.0 * u.megapascal
Model.maxViscosity = 5.0e23 * u.pascal * u.second
Model.minViscosity = 1.0e19 * u.pascal * u.second
Model.diffusivity = 1.1e-6 * u.metre**2 / u.second 
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)


# ## Define Materials
# The model has initially 4 materials (air, crust, mantle lithosphere and mantle). We add a fifth material for the sediment. Sediment will only appear if surface processes are turned on...(and if there is sedimentation of course)

# In[5]:


## Definition of materials

air               = Model.add_material(name="Air", shape=GEO.shapes.Layer3D(top=Model.top, bottom=4.0 * u.kilometer))
stickyAir         = Model.add_material(name="StickyAir", shape=GEO.shapes.Layer3D(top=air.bottom, bottom= 0.0 * u.kilometer))
sedimentMarge     = Model.add_material(name="SedimentMarge", shape=GEO.shapes.Layer3D(top=stickyAir.bottom, bottom= -9.0 * u.kilometer))
uppercrust        = Model.add_material(name="UpperCrust", shape=GEO.shapes.Layer3D(top=sedimentMarge.bottom, bottom=-15.0 * u.kilometer))
lowercrust        = Model.add_material(name="LowerCrust", shape=GEO.shapes.Layer3D(top=-15.0 * u.kilometer, bottom=-30.0 * u.kilometer))
mantleLithosphere = Model.add_material(name="MantleLithosphere", shape=GEO.shapes.Layer3D(top=-30.0 * u.kilometer, bottom=-40.0 * u.kilometer))
mantle            = Model.add_material(name="Mantle", shape=GEO.shapes.Layer3D(top=mantleLithosphere.bottom, bottom=Model.bottom))

# Definition of the Fault between continent and sediments

wz_slope = GEO.shapes.HalfSpace(normal=(-0.8,0.2,0.5), origin=(158.*u.kilometer,50.*u.kilometer, 0.*u.kilometer))
wz_bottom = GEO.shapes.HalfSpace(normal=(0.8,-0.2,-0.5), origin=(160.*u.kilometer,50.*u.kilometer,0.*u.kilometer))
wz_top = GEO.shapes.HalfSpace(normal=(0.,0.,1.), origin=(200.*u.kilometer,100.*u.kilometer,0.*u.kilometer))
wz_fond = GEO.shapes.HalfSpace(normal=(0.,0.,-1.), origin=(200.*u.kilometer,100.*u.kilometer,-33.*u.kilometer))
compositeShape = wz_slope & wz_top & wz_bottom & wz_fond
Fault = Model.add_material(name="Fault", shape=compositeShape)

# Continent part
top = GEO.shapes.HalfSpace(normal=(0.,0.,1.), origin=(200.*u.kilometer,100.*u.kilometer,2.*u.kilometer))
right = GEO.shapes.HalfSpace(normal=(1.,0.,0.), origin=(200.*u.kilometer,100.*u.kilometer,-0.*u.kilometer))
bottom = GEO.shapes.HalfSpace(normal=(0.,0.,-1.), origin=(200.*u.kilometer,100.*u.kilometer,-30.*u.kilometer))
slope = GEO.shapes.HalfSpace(normal=(-0.8,0.2,0.5), origin=(160.*u.kilometer,50.*u.kilometer,0.*u.kilometer))
CompositeShape_conti = slope & top & right & bottom
continent = Model.add_material(name="continent", shape=CompositeShape_conti)

# Fig = glucifer.Figure()
# Fig.Points(Model.swarm, Model.materialField, cullface=False, opacity=1.)
# viewer = Fig.viewer(resolution=(1200,600))
# viewer = Fig.viewer(axis=True)
# viewer.window()




# ### Material specific definitions
# 

# In[7]:


# I Have changes the values here. We can use more realistic values later on.

#air.capacity = 100. * u.joule / (u.kelvin * u.kilogram)
#stickyAir.capacity = 100. * u.joule / (u.kelvin * u.kilogram)


# In[8]:


# Density

air.density = 1. * u.kilogram / u.metre**3
stickyAir.density = 1. * u.kilogram / u.metre**3
Fault.density = GEO.LinearDensity(2000. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
# There was an error here. You had sedimentMarge = GEO.Lin... instead of sedimentMarge.density
sedimentMarge.density = GEO.LinearDensity(2200. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
uppercrust.density = GEO.LinearDensity(2400. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
continent.density = GEO.LinearDensity(2720. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
lowercrust.density = GEO.LinearDensity(2720. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
mantleLithosphere.density = GEO.LinearDensity(3370. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
mantle.density = GEO.LinearDensity(3395. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)


# In[9]:


# Radiogenic Heat Production

sedimentMarge.radiogenicHeatProd   = 0.70 * u.microwatt / u.meter**3
uppercrust.radiogenicHeatProd = 0.70 * u.microwatt / u.meter**3
lowercrust.radiogenicHeatProd = 0.30 * u.microwatt / u.meter**3
continent.radiogenicHeatProd = 0.30 * u.microwatt / u.meter**3
mantleLithosphere.radiogenicHeatProd = 0.01358 * u.microwatt / u.meter**3
Fault.radiogenicHeatProd = 0.25 * u.microwatt / u.meter**3 # faccenda 2008: 1-5
mantle.radiogenicHeatProd = 0.02e-6 * u.microwatt / u.meter**3


# ### Viscous Rheologies

# In[10]:


rh = GEO.ViscousCreepRegistry() 

Model.maxViscosity = 5e23 * u.pascal * u.second
Model.minViscosity = 1e18 * u.pascal * u.second

air.viscosity           = 1e18 * u.pascal * u.second
stickyAir.viscosity     = 1e19 * u.pascal * u.second
sedimentMarge.viscosity = rh.Kronenberg_et_al_1984
uppercrust.viscosity    = rh.Goetze_et_al_1978
lowercrust.viscosity    = rh.Wang_et_al_2012
continent.viscosity     = 30 * rh.Goetze_et_al_1978
mantleLithosphere.viscosity = 5 * rh.Brace_and_Kohlstedt_1980
mantle.viscosity        = rh.Brace_and_Kohlstedt_1980
Fault.viscosity         = rh.Goetze_et_al_1978


# ## Plasticities

# In[11]:


pl = GEO.PlasticityRegistry()


# In[12]:


# I have cleaned up the following

uppercrust.plasticity         = pl.Rey_et_al_2010_UpperCrust
lowercrust.plasticity         = pl.Rey_et_al_2010_UpperCrust
sedimentMarge.plasticity      = pl.Rey_et_al_2010_UpperCrust
continent.plasticity          = pl.Rey_et_al_2010_UpperCrust
mantleLithosphere.plasticity  = pl.Rey_et_al_2010_Mantle
mantle.plasticity             = pl.Rey_et_al_2010_Mantle

# Let's make the fault weak
Fault.plasticity              = GEO.DruckerPrager(cohesion=2.0 * u.megapascal, frictionCoefficient=0.1154)

uppercrust.stressLimiter         = 300. * u.megapascal
lowercrust.stressLimiter         = 300. * u.megapascal
sedimentMarge.stressLimiter      = 300. * u.megapascal
continent.stressLimiter          = 300. * u.megapascal
mantleLithosphere.stressLimiter  = 380. * u.megapascal
mantle.stressLimiter             = 380. * u.megapascal
Fault.stressLimiter              = 300. * u.megapascal


# ## Add Temperature Boundary Conditions
# 
# Temperature is 293.15K at the top and 3570K at the bottom. Temperature is constant in the mantle and the air layers. 

# In[13]:


Model.set_temperatureBCs(top=293.15 * u.degK, 
                         bottom=1603.15 * u.degK, 
                         indexSets=[(stickyAir.indices, 293.15 * u.degK ),
                                    (air.indices, 293.15 * u.degK )]) 


# ## Add Velocity Boundary Conditions
# 
# We push on the right side. The back and front wall are freeslip. We use a pseudo isostatic support condition at the bottom.

# In[19]:


# The Model was free slip everywhere, which is OK for gravitationally unstable set up but not here...

Model.set_velocityBCs(left=[-2.5 * u.centimeter / u.year, None, 0.],
                      right=[0., None, 0.],
                      back=[None, 0., None],
                      front=[None, 0., None],
                      bottom=GEO.LecodeIsostasy(reference_mat=mantle.index,
                                                average=False))


# # Add Passive Tracers 

# In[15]:


import numpy as np

xp = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), 100)
yp = np.linspace(GEO.nd(Model.minCoord[1]), GEO.nd(Model.maxCoord[1]), 100)

xp, yp = np.meshgrid(xp, yp)
xp = xp.flatten()
yp = yp.flatten()
zp = np.zeros(xp.shape)

surface_tracers = Model.add_passive_tracers(name="Surface2", vertices=[xp, yp, zp])
moho_tracers = Model.add_passive_tracers(name="Moho", vertices=[xp, yp, zp+GEO.nd(mantle.top)])


# In[16]:


# revert Back to normal values
GEO.rcParams["initial.nonlinear.tolerance"]= 5e-2
GEO.rcParams["nonlinear.tolerance"]= 1e-2
#GEO.rcParams["advection.diffusion.method"] = "SLCN"

# In[17]:


Model.init_model()


# In[18]:




# In[ ]:

Model.restart(step=8)
Model.run_for(5.0* u.megayears, checkpoint_interval=0.1*u.megayears)

