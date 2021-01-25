

import UWGeodynamics as GEO
#import glucifer

import badlands
# In[2]:

u = GEO.UnitRegistry

##### Characteristic values of the system

half_rate = 1.0 * u.centimeter / u.year # What is it?
model_length = 512e3 * u.meter
model_width = 256e3 * u.meter
surfaceTemp = 273.15 * u.degK
baseModelTemp = 1623.15 * u.degK
bodyforce = 3300. * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2 # ckeck publi schellart

KL = model_length # Characteristic length
Kt = KL / half_rate # Characteristic time
KM = bodyforce * KL**2 * Kt**2 # Characteristic mass
KT = (baseModelTemp - surfaceTemp) # Characteristic temperature

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT


# In[3]:

resolution = (256, 128, 128)
Model = GEO.Model(elementRes=resolution,
                  minCoord=(0. * u.kilometer, 0. * u.kilometer, -120. * u.kilometer),
                  maxCoord=(512. * u.kilometer, 256. * u.kilometer, 40. * u.kilometer),
                  gravity=(0.0 ,0.0, -9.81 * u.meter / u.second**2))


# ## Output
Model.outputDir="CM_5cmyr"
Model.minViscosity = 1e18 * u.pascal * u.second
Model.maxViscosity = 1e23 * u.pascal * u.second
Model.stressLimiter = 300. * u.megapascal
Model.diffusivity = 1e-6 * u.metre**2 / u.second
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)


# ## Define Materials
air_shape         = GEO.shapes.Layer3D(top=Model.top, bottom=10.0 * u.kilometer)


air               = Model.add_material(name="air", shape=air_shape)

stickyAir_shape   = GEO.shapes.Layer3D(top=air.bottom, bottom= 0.0 * u.kilometer)

stickyAir         = Model.add_material(name="StickyAir", shape=stickyAir_shape)
sedimentMarge     = Model.add_material(name="SedimentMarge", shape=GEO.shapes.Layer3D(top=stickyAir.bottom, bottom= -9.0 * u.kilometer))
uppercrust        = Model.add_material(name="UppperCrust", shape=GEO.shapes.Layer3D(top=sedimentMarge.bottom, bottom=-15.0 * u.kilometer))
midcrust          = Model.add_material(name="MidCrust", shape=GEO.shapes.Layer3D(top=-15.0 * u.kilometer, bottom=-25.0 * u.kilometer))
lowercrust        = Model.add_material(name="LowerCrust", shape=GEO.shapes.Layer3D(top=-25.0 * u.kilometer, bottom=-35.0 * u.kilometer))
mantleLithosphere = Model.add_material(name="MantleLithosphere", shape=GEO.shapes.Layer3D(top=-35.0 * u.kilometer, bottom=-110.0 * u.kilometer))

mantle_shape      = GEO.shapes.Layer3D(top=mantleLithosphere.bottom, bottom=Model.bottom)

mantle            = Model.add_material(name="Mantle", shape=mantle_shape)
sediment          = Model.add_material(name="Sediment")

wz_slope = GEO.shapes.HalfSpace(normal=(-0.8,0.0,0.5), origin=(328.*u.kilometer,50.*u.kilometer, 0.*u.kilometer))
wz_bottom = GEO.shapes.HalfSpace(normal=(0.8,-0.0,-0.5), origin=(330.*u.kilometer,50.*u.kilometer,0.*u.kilometer))
wz_top = GEO.shapes.HalfSpace(normal=(0.,0.,1.), origin=(200.*u.kilometer,100.*u.kilometer,0.*u.kilometer))
wz_fond = GEO.shapes.HalfSpace(normal=(0.,0.,-1.), origin=(200.*u.kilometer,100.*u.kilometer,-38.*u.kilometer))
compositeShape = wz_slope & wz_top & wz_bottom & wz_fond
Fault = Model.add_material(name="Fault", shape=compositeShape)

# Continent part
top = GEO.shapes.HalfSpace(normal=(0.,0.,1.), origin=(200.*u.kilometer,100.*u.kilometer,5.*u.kilometer))
right = GEO.shapes.HalfSpace(normal=(1.,0.,0.), origin=(512.*u.kilometer,256.*u.kilometer,-0.*u.kilometer))
bottom = GEO.shapes.HalfSpace(normal=(0.,0.,-1.), origin=(200.*u.kilometer,100.*u.kilometer,-35.*u.kilometer))
slope = GEO.shapes.HalfSpace(normal=(-0.8,0.0,0.5), origin=(330.*u.kilometer,50.*u.kilometer,0.*u.kilometer))
CompositeShape_conti = slope & top & right & bottom
continent = Model.add_material(name="continent", shape=CompositeShape_conti)



# ### Material specific definitions
#

#air.diffusivity = 1.0e-6 * u.metre**2 / u.second
air.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)
# stickyAir.diffusivity = 1.0e-6 * u.metre**2 / u.second
stickyAir.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)


#air.density = 0.1 * u.kilogram / u.metre**3 #1.

# In[8]:

# Density

air.density         = 1. * u.kilogram / u.metre**3
stickyAir.density   = 1. * u.kilogram / u.metre**3
Fault.density = GEO.LinearDensity(2000. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
#sedimentMarge.density = GEO.LinearDensity(reference_density=2500. * u.kilogram / u.metre**3)
sedimentMarge.density = GEO.LinearDensity(2500. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

#uppercrust.density  = GEO.LinearDensity(reference_density=2620. * u.kilogram / u.metre**3)
uppercrust.density  = GEO.LinearDensity(2620. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
#midcrust.density          = GEO.LinearDensity(reference_density=2750. * u.kilogram / u.metre**3)
midcrust.density          = GEO.LinearDensity(2750. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
#lowercrust.density  = GEO.LinearDensity(reference_density=2850. * u.kilogram / u.metre**3)
lowercrust.density  = GEO.LinearDensity(2850. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
#continent.density   = GEO.LinearDensity(reference_density=2800. * u.kilogram / u.metre**3)
continent.density = GEO.LinearDensity(2750. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
#mantleLithosphere.density = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3)
mantleLithosphere.density = GEO.LinearDensity(3370. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
#mantle.density      = GEO.LinearDensity(reference_density=3395. * u.kilogram / u.metre**3)
mantle.density      = GEO.LinearDensity(3395. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
#sediment.density          = GEO.LinearDensity(reference_density=2300. * u.kilogram / u.metre**3)
sediment.density       = GEO.LinearDensity(2400. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

# Radiogenic Heat Production
Fault.radiogenicHeatProd = 0.4 * u.microwatt / u.meter**3 # faccenda 2008: 1-5

sedimentMarge.radiogenicHeatProd   = 0.95  * u.microwatt / u.meter**3

uppercrust.radiogenicHeatProd = 0.75 * u.microwatt / u.meter**3
midcrust.radiogenicHeatProd   = 0.60 * u.microwatt / u.meter**3
lowercrust.radiogenicHeatProd = 0.45 * u.microwatt / u.meter**3
continent.radiogenicHeatProd = 0.50 * u.microwatt / u.meter**3

mantleLithosphere.radiogenicHeatProd = 0.01154 * u.microwatt / u.meter**3
mantle.radiogenicHeatProd = 0.01154  * u.microwatt / u.meter**3

sediment.radiogenicHeatProd   = 0.85 * u.microwatt / u.meter**3



rh = GEO.ViscousCreepRegistry()

air.viscosity           = 1e18 * u.pascal * u.second
stickyAir.viscosity     = 1e19 * u.pascal * u.second

sedimentMarge.viscosity = rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
sedimentMarge.minViscosity = 1e19 * u.pascal * u.second
#sedimentMarge.maxViscosity = 5e22 * u.pascal * u.second

uppercrust.viscosity    = rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
uppercrust.minViscosity = 1e19 * u.pascal * u.second
#uppercrust.maxViscosity = 5e22 * u.pascal * u.second

midcrust.viscosity    = rh.Wet_Quartz_Dislocation_Paterson_and_Luan_1990
midcrust.minViscosity = 1e19 * u.pascal * u.second
#midcrust.maxViscosity = 5e22* u.pascal * u.second

lowercrust.viscosity    = rh.Wet_Olivine_Dislocation_Goetze_et_al_1978#Wet_Quartz_Dislocation_Goetze_et_al_1978
lowercrust.minViscosity = 1e19 * u.pascal * u.second
#lowercrust.maxViscosity = 5e22 * u.pascal * u.second

continent.viscosity    = 5e23 * u.pascal * u.second
continent.minViscosity = 5e19 * u.pascal * u.second
#continent.maxViscosity = 5e22 * u.pascal * u.second

mantleLithosphere.viscosity  = rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003
mantleLithosphere.minViscosity = 5e20 * u.pascal * u.second
mantleLithosphere.maxViscosity = 5e22 * u.pascal * u.second

mantle.viscosity             = rh.Wet_Olivine_Dislocation_Hirth_and_Kohlstedt_2003
mantle.minViscosity          = 5e20 * u.pascal * u.second

sediment.viscosity           = rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995
sediment.minViscosity        = 1e19 * u.pascal * u.second

Fault.viscosity    = 1e20 * u.pascal * u.second # * rh.Goetze_et_al_1978
Fault.minViscosity = 1e18 * u.pascal * u.second
#Fault.maxViscosity = 1e20 * u.pascal * u.second




# ## Plasticities

pl = GEO.PlasticityRegistry()

uppercrust.plasticity = 0
uppercrust.stressLimiter           = 300. * u.megapascal
#uppercrust.plasticity.frictionCoefficient = 0.55
#uppercrust.plasticity.frictionAfterSoftening = 0.055

midcrust.plasticity = 0
midcrust.stressLimiter           = 300. * u.megapascal
#midcrust.plasticity.frictionCoefficient = 0.55
#midcrust.plasticity.frictionAfterSoftening = 0.055

lowercrust.plasticity = 0
#lowercrust.plasticity.frictionCoefficient = 0.55
#lowercrust.plasticity.frictionAfterSoftening = 0.028
lowercrust.stressLimiter         = 300. * u.megapascal

continent.plasticity = 0
continent.stressLimiter         = 300. * u.megapascal

mantle.plasticity = pl.Rey_and_Muller_2010_LithosphericMantle
mantle.stressLimiter  = 350. * u.megapascal

mantleLithosphere.plasticity = pl.Rey_and_Muller_2010_LithosphericMantle
mantleLithosphere.stressLimiter  = 300. * u.megapascal


sedimentMarge.plasticity = 0
sedimentMarge.stressLimiter           = 300. * u.megapascal
#sedimentMarge.plasticity.frictionCoefficient = 0.55
#sedimentMarge.plasticity.frictionAfterSoftening = 0.055

sediment.plasticity = 0
sediment.stressLimiter           = 300. * u.megapascal
#sediment.plasticity.frictionCoefficient = 0.55
#sediment.plasticity.frictionAfterSoftening = 0.055


Fault.plasticity              = GEO.DruckerPrager(cohesion=2.0 * u.megapascal, frictionCoefficient=0.1154)
Fault.stressLimiter              = 300. * u.megapascal


# ## Passive Tracers

import numpy as np

xp = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), 256)
yp = np.linspace(GEO.nd(Model.minCoord[1]), GEO.nd(Model.maxCoord[1]), 256)

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


#surface_tracers = Model.add_passive_tracers(name="Surface", vertices=[xp, yp, zp])
surface2_tracers = Model.add_passive_tracers(name="Surface2", vertices=surface_tracersNp)
moho_tracers = Model.add_passive_tracers(name="Moho", vertices=moho_tracersNp)
#Lab_tracers = Model.add_passive_tracers(name="Lab", vertices=[xp,yp,-90.*u.kilometer])


# ## Add Temperature Boundary Conditions
#
# Temperature is 293.15K at the top and 3570K at the bottom. Temperature is constant in the mantle and the air layers.




Model.set_temperatureBCs(top=273.15 * u.degK,
                         nodeSets = [(air_shape, 273.15 * u.degK), (stickyAir_shape, 273.15 * u.degK),(mantle_shape, 1623.15 * u.degK)],
                         bottom=1623.15 * u.degK)


Model.set_velocityBCs(front=[None, 0., None],
                      back=[None, 0., None],
                      left=[-5.0 * u.centimeter / u.year, 0., None],
                      right=[0., 0. , None],
                      bottom=GEO.LecodeIsostasy(reference_mat=mantle, average=False),
                      top= [None, None, None],
                      order_wall_conditions=["front", "back", "left", "right", "bottom", "top"])



#
Model.surfaceProcesses = GEO.surfaceProcesses.Badlands(airIndex=[stickyAir.index, air.index], sedimentIndex=sediment.index,XML="ressources/badlands_0.xml", resolution=1. * u.kilometer,restartFolder="outbdls",
            checkpoint_interval=0.05 * u.megayears,restartStep=47)



Model.init_model()



Model.set_temperatureBCs(top=273.15 * u.degK,
                         nodeSets = [(air_shape, 273.15 * u.degK), (stickyAir_shape, 273.15 * u.degK)],
                         bottom=1623.15 * u.degK)


import underworld.function as fn
def post_hook():
    coords = fn.input()
    zz = coords[0] / (GEO.nd(Model.maxCoord[0]) - GEO.nd(Model.minCoord[0]))
    fact = fn.math.pow(fn.math.tanh(zz*50.0) + fn.math.tanh((1.0-zz)*50.0) - fn.math.tanh(50.0), 10.0)
    Model.plasticStrain.data[:] = Model.plasticStrain.data[:] * fact.evaluate(Model.swarm) # if symetric
    #Model.plasticStrain.data[250:] = 1
    #uw.barrier()

Model.post_solve_functions["B"] = post_hook

##Testing Julians suggestion
#scr_rtol = 1e-6
## Schur complement solver options
#Model.solver.options.scr.ksp_rtol = scr_rtol
## Inner solve (velocity), A11 options
#Model.solver.options.A11.ksp_rtol = 1e-1 * scr_rtol

GEO.rcParams["popcontrol.particles.per.cell.3D"] = 60
GEO.rcParams["swarm.particles.per.cell.3D"] = 60



GEO.rcParams["initial.nonlinear.tolerance"] = 5e-2
GEO.rcParams["nonlinear.tolerance"] = 5e-2

Model.run_for(20.0* u.megayears,checkpoint_interval=0.05*u.megayears,restartStep=47)

