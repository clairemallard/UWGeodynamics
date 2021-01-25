

import UWGeodynamics as GEO
import glucifer


# In[2]:

u = GEO.UnitRegistry


# In[3]:

# Characteristic values of the system
half_rate = 1.8 * u.centimeter / u.year
model_length = 500e3 * u.meter
model_width = 500e3 * u.meter
surfaceTemp = 273.15 * u.degK
baseModelTemp = 1623.15 * u.degK
bodyforce = 3370 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2 #3300

KL = model_length
Kt = KL / half_rate
KM = bodyforce * KL**2 * Kt**2
KT = (baseModelTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"] = KM
GEO.scaling_coefficients["[temperature]"] = KT


# In[4]:

Model = GEO.Model(elementRes=(128, 128, 64),
                  minCoord=(0. * u.kilometer, 0. * u.kilometer, -180. * u.kilometer),
                  maxCoord=(500. * u.kilometer, 500. * u.kilometer, 20. * u.kilometer),
                  gravity=(0.0, 0.0, -9.81 * u.meter / u.second**2))


# ## Global definitions

# In[5]:

Model.outputDir="outputs"


# In[6]:

Model.maxViscosity  = 5e23 * u.pascal * u.second
Model.minViscosity  = 1e19 * u.pascal * u.second
Model.stressLimiter = 300. * u.megapascal
Model.diffusivity = 1e-6 * u.metre**2 / u.second
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)


# ## Define Materials
#


air = Model.add_material(name="Air", shape=GEO.shapes.Layer3D(top=Model.top, bottom=0 * u.kilometer))
Ucrust = Model.add_material(name="UpperCrust", shape=GEO.shapes.Layer3D(top=air.bottom, bottom=-20 * u.kilometer))
Lcrust = Model.add_material(name="LowerCrust", shape=GEO.shapes.Layer3D(top=Ucrust.bottom, bottom=-35 * u.kilometer))
mantleLithosphere = Model.add_material(name="MantleLithosphere", shape=GEO.shapes.Layer3D(top=Lcrust.bottom, bottom=-160 * u.kilometer))
mantle = Model.add_material(name="Mantle", shape=GEO.shapes.Layer3D(top=mantleLithosphere.bottom, bottom=Model.bottom))
sediment   = Model.add_material(name="Sediment")



# ### Material specific definitions
#

# In[10]:

air.diffusivity = 1.0e-5 * u.metre**2 / u.second
air.capacity = 1000. * u.joule / (u.kelvin * u.kilogram)

## Density

air.density = 1. * u.kilogram / u.metre**3
Ucrust.density  = GEO.LinearDensity(reference_density=2800. * u.kilogram / u.metre**3) # 2800.
Lcrust.density  = GEO.LinearDensity(reference_density=2900. * u.kilogram / u.metre**3) # 28900.
mantleLithosphere.density  = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3)
mantle.density  = GEO.LinearDensity(reference_density=3370. * u.kilogram / u.metre**3)
sediment.density = 2700. * u.kilogram / u.metre**3 # 2700.

# In[11]:

Ucrust.radiogenicHeatProd = 1.3 * u.microwatt / u.meter**3 #0.7
Lcrust.radiogenicHeatProd = 0.2 * u.microwatt / u.meter**3 #0.4
sediment.radiogenicHeatProd = 0.7 * u.microwatt / u.meter**3 #0.7
mantleLithosphere.radiogenicHeatProd = 0 * u.microwatt / u.meter**3
mantle.radiogenicHeatProd = 0 * u.microwatt / u.meter**3


# ### Viscous Rheologies


rh = GEO.ViscousCreepRegistry()


air.viscosity                = 1e19 * u.pascal * u.second
Ucrust.viscosity              = rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995
Lcrust.viscosity              = rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995
mantleLithosphere.viscosity  = rh.Dry_Olivine_Diffusion_Hirth_and_Kohlstedt_2003
mantle.viscosity             = rh.Wet_Olivine_Diffusion_Hirth_and_Kohlstedt_2003
sediment.viscosity           = rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995

# ### Plasticities


pl = GEO.PlasticityRegistry()


Ucrust.plasticity = GEO.DruckerPrager(
    cohesion=10.*u.megapascal,
    cohesionAfterSoftening=2.0*u.megapascal,
    frictionCoefficient=0.12,
    frictionAfterSoftening=0.02,
    epsilon1=0.0, epsilon2=0.2
    )
Lcrust.plasticity = GEO.DruckerPrager(
    cohesion=20.*u.megapascal,
    cohesionAfterSoftening=4.*u.megapascal,
    frictionCoefficient=0.577,
    frictionAfterSoftening=0.1154,
    epsilon1=0.0, epsilon2=0.2
    )
Lcrust.stressLimiter  = 150. * u.megapascal

mantleLithosphere.plasticity = GEO.DruckerPrager(
    cohesion=10.*u.megapascal,
    cohesionAfterSoftening=2.*u.megapascal,
    frictionCoefficient=0.577,
    frictionAfterSoftening=0.1154,
    epsilon1=0.0, epsilon2=0.2
    )
mantle.plasticity = GEO.DruckerPrager(
    cohesion=10.*u.megapascal,
    cohesionAfterSoftening=2.*u.megapascal,
    frictionCoefficient=0.577,
    frictionAfterSoftening=0.1154,
    epsilon1=0.0, epsilon2=0.2
    )
sediment.plasticity = GEO.DruckerPrager(
    cohesion=10.*u.megapascal,
    cohesionAfterSoftening=2.0*u.megapascal,
    frictionCoefficient=0.12,
    frictionAfterSoftening=0.02,
    epsilon1=0.0, epsilon2=0.2
    )

# ## Passive Tracers
#


import numpy as np

xp = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), 100)
yp = np.linspace(GEO.nd(Model.minCoord[1]), GEO.nd(Model.maxCoord[1]), 100)

xp, yp = np.meshgrid(xp, yp)
xp = xp.flatten()
yp = yp.flatten()
zp = np.zeros(xp.shape)

surface_tracers = Model.add_passive_tracers(name="Surface", vertices=[xp, yp, zp])
moho_tracers = Model.add_passive_tracers(name="Moho", vertices=[xp, yp, zp+GEO.nd(mantleLithosphere.top)])


# ## Temperature Boundary Condition


Model.set_temperatureBCs(top=273.15 * u.degK,
                         bottom=1623.15 * u.degK,
                         materials = [(air, 273.15 * u.degK), (mantle, 1623.15 * u.degK)])

Model.init_model()

Model.set_temperatureBCs(top=273.15 * u.degK,
                         bottom=1623.15 * u.degK)


#air.compressibility = 1e4  # Not sure what's a good value is

#P, bottomPress = Model.get_lithostatic_pressureField()


# Get the average of the pressure along the bottom, and make it a pressure BC along the bottom
#bottomPress = GEO.Dimensionalize(numpy.average(bottomPress), u.megapascal)



# ## Velocity Boundary Conditions
#
# We pull on the left and right side. The back and front wall are freeslip. We use a pseudo isostatic support condition at the bottom.

# In[30]:

Model.set_velocityBCs(left=[-0.35* u.centimeter / u.year, 0., None],
                      right=[0.35 * u.centimeter / u.year, 0., None],
                      back=[None, 0., None],
                      front=[None, 0., None],
                      bottom=GEO.LecodeIsostasy(reference_mat=mantle,
                                                average=True))





# ## Initial Plastic Strain

import numpy as np

def gaussian(xx, centre, width):
    return ( np.exp( -(xx - centre)**2 / width ))

maxDamage = 0.25
Model.plasticStrain.data[:] = maxDamage * np.random.rand(*Model.plasticStrain.data.shape[:])
centers = 0.3 * Model.swarm.particleCoordinates.data[:,1] + 0.34
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,0], centers, GEO.nd(5.0 * u.kilometer))
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,2], GEO.nd(-35. * u.kilometer) , GEO.nd(5.0 * u.kilometer))




# # Surface Processes / BADLANDS

Model.surfaceProcesses = GEO.surfaceProcesses.Badlands(airIndex=[air.index], sedimentIndex=sediment.index,
                                            XML="resources/badlands.xml", resolution=1. * u.kilometer,
                                            checkpoint_interval=0.005 * u.megayears)



GEO.rcParams["initial.nonlinear.tolerance"] = 1e-2 #1e-3
GEO.rcParams["nonlinear.tolerance"] = 1e-3 #1e-3
GEO.rcParams["swarm.particles.per.cell.3D"] = 60
GEO.rcParams["popcontrol.particles.per.cell.3D"] = 60



Model.run_for(50.0 * u.megayears, checkpoint_interval=0.005 * u.megayears, restart_checkpoint=5)
