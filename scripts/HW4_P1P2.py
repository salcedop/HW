import openmc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import pandas as pd
import argparse

#material 1

def create_material(boron_conc,pitch,lattice,enrichment,moderator,temp1,temp2):
  

            uo2 = openmc.Material(1,"fuel",temperature=temp1)
            uo2.add_element('U', 1.0, enrichment=enrichment)
            uo2.add_element('O', 2.0)
            uo2.set_density('g/cc', 10.0)
            zirconium = openmc.Material(2,"zirconium", temperature=temp2)
            zirconium.add_element('Zr', 1.0),
            zirconium.set_density('g/cm3', 6.6),
            #this function creates borated light water, for other moderators you will need to replace\n",
            if (moderator == 'light_water'):
                mod = openmc.model.borated_water(boron_ppm=boron_conc, temperature=temp2,pressure=15)
            elif (moderator == 'heavy_water'):
                mod = openmc.Material(name='heavy water')
                mod.add_nuclide('H2',2.0)
                mod.add_nuclide('O16',1.0)
                mod.add_s_alpha_beta('c_D_in_D2O')
                mod.set_density('g/cm3',1.1)
            elif(moderator == 'graphite'):
               mod = openmc.Material(name='Buffer')
               mod.set_density('g/cm3', 1.0)
               mod.add_element('C', 1.0)
               mod.add_s_alpha_beta('c_Graphite')
            else:
               print('no valid moderator was provided')


            mats = openmc.Materials([uo2, zirconium, mod])
            mats.export_to_xml()
            fuel_or = openmc.ZCylinder(R=0.50)
            clad_ir = openmc.ZCylinder(R=0.51)
            clad_or = openmc.ZCylinder(R=0.57)
            fuel_region = -fuel_or
            gap_region = +fuel_or & -clad_ir
            clad_region = +clad_ir & -clad_or
            fuel = openmc.Cell(1, 'fuel')
            fuel.fill = uo2
            fuel.region = fuel_region
            gap = openmc.Cell(2, 'air gap')
            gap.region = gap_region
            clad = openmc.Cell(3, 'clad')
            clad.fill = zirconium
            clad.region = clad_region
            if (lattice == 'rectangular'):
               box = openmc.get_rectangular_prism(width=pitch, height=pitch,boundary_type='reflective')
            elif (lattice == 'triangular'):
                box = openmc.get_hexagonal_prism(edge_length=pitch,boundary_type='reflective')
            else:

               print("wrong lattice specified! Can only handle rectangular/triangular lattice")
            water_region = box & +clad_or
            moderator = openmc.Cell(4, 'moderator')
            moderator.fill = mod
            moderator.region = water_region
            root = openmc.Universe(cells=(fuel, gap, clad, moderator))
            geom = openmc.Geometry(root)
            geom.export_to_xml()
            
            cell_filter = openmc.CellFilter([fuel,gap, clad, moderator])
            
            #tallies over 2 energy groups with 4 eV being thermal bound\n",
            #energy_filter = openmc.EnergyFilter([0., 4.0, 20.0e6])
            #t = openmc.Tally(1)
            tallies_file = openmc.Tallies()
            
  
            #t.filters = [cell_filter, energy_filter]
            
            #thermal, fast and total leakage tallies

            leak = openmc.Tally(name = 'leakage')
            leak.filters = [cell_filter]
            leak.scores = ['current']
            tallies_file.append(leak)

            thermal_leak = openmc.Tally(name='thermal leakage')
            thermal_leak.filters = [cell_filter,openmc.EnergyFilter([0.,0.625])]
            thermal_leak.scores = ['current']
            tallies_file.append(thermal_leak) 
          
            fast_leak = openmc.Tally(name='fast leakage')
            fast_leak.filters = [cell_filter, openmc.EnergyFilter([0.625,20e6])]
            fast_leak.scores = ['current']
            tallies_file.append(fast_leak)
             

            #K-eigenvalue (infinity)
            fiss_rate = openmc.Tally(name='fiss.rate')
            abs_rate = openmc.Tally(name ='abs.rate')
            fiss_rate.scores = ['nu-fission']
            abs_rate.scores = ['absorption']
            tallies_file += (fiss_rate, abs_rate)

            #Resonance Escape Prob
 
            therm_abs_rate = openmc.Tally(name='therm.abs.rate')
            therm_abs_rate.scores = ['absorption']
            therm_abs_rate.filters = [openmc.EnergyFilter([0.,0.625])]
            tallies_file.append(therm_abs_rate)

            #Fast Fission Factor

            therm_fiss_rate = openmc.Tally(name='therm.fiss.rate')
            therm_fiss_rate.scores = ['nu-fission']
            therm_fiss_rate.filters = [openmc.EnergyFilter([0.,0.625])]
            tallies_file.append(therm_fiss_rate)


            #Thermal Flux Utilization

            fuel_therm_abs_rate = openmc.Tally(name='fuel.therm.abs.rate')
            fuel_therm_abs_rate.scores = ['absorption']
            fuel_therm_abs_rate.filters = [openmc.EnergyFilter([0.,0.625]), openmc.CellFilter([fuel])]
            tallies_file.append(fuel_therm_abs_rate)
         
            # these are the main reaction rates you should need\n",
            #t.scores = ['absorption','nu-fission','fission']
            #tallies = openmc.Tallies([t])
            tallies_file.export_to_xml()



def post_process_statepoint(statepoint):
       
       fiss_rate = statepoint.get_tally(name='fiss.rate')
       abs_rate = statepoint.get_tally(name='abs.rate')
       leak = statepoint.get_tally(name='leakage')
       
       #get k-infinity
       keff = fiss_rate/(abs_rate + leak)
       df1 = keff.get_pandas_dataframe()
       #resonance escape prob
       therm_abs_rate = statepoint.get_tally(name='therm.abs.rate')
       thermal_leak = statepoint.get_tally(name='thermal leakage')
       res_esc = (therm_abs_rate + thermal_leak) / (abs_rate + thermal_leak)
       df2 = res_esc.get_pandas_dataframe()
       #fast fission factor
       therm_fiss_rate = statepoint.get_tally(name='therm.fiss.rate')
       fast_fiss = fiss_rate / therm_fiss_rate
       df3 = fast_fiss.get_pandas_dataframe()
       #thermal utilization factor
       fuel_therm_abs_rate = statepoint.get_tally(name='fuel.therm.abs.rate')
       therm_util = fuel_therm_abs_rate / therm_abs_rate
       df4 = therm_util.get_pandas_dataframe()
       #eta
       eta = therm_fiss_rate / fuel_therm_abs_rate
       df5 = eta.get_pandas_dataframe()

       frames = [df1,df2,df3,df4,df5]
       result = pd.concat(frames)
       result.to_csv('five_factor.csv')
       
       

parser = argparse.ArgumentParser()

parser.add_argument('-p','--problem',choices=('p1','p2','p3'), default='p1',
                   help='problem number')
parser.add_argument('-subp','--subproblem',choices=('a','b'),default='a',
                   help='subproblem number')
parser.add_argument('-n','--node',action='store_true',
                   help='Are you using Skylake node?')

args = parser.parse_args()
#set MPI argument
if args.node:

   mpi_args=['mpiexec','-bind-to','numa','-n','2']
else:
   mpi_args=['mpiexec','-n','2']

if (args.problem == 'p1'):

  moderator = 'heavy_water'
  enrichment = 0.007
  lattice = ['rectangular']

  if (args.subproblem == 'a'):
     pitch = [1.40,1.50,1.60]
     proc_tallies = False
  else:
     pitch = [1.40]
     proc_tallies = True

elif(args.problem == 'p2'):
  
  moderator = 'graphite'
  enrichment = 0.007
  lattice = ['rectangular']

  if (args.subproblem == 'a'):
     pitch = [1.40,1.50,1.60]
     proc_tallies = False
  else:
     pitch = [1.40]
     proc_tallies = True
elif(args.problem == 'p3'):

   moderator = 'light_water'
   enrichment = 4
   lattice = ['rectangular','triangular']

   if (args.subproblem == 'a'):
     pitch = [1.40,1.50,1.60]
     proc_tallies = False
   else:
     pitch  = "Error! Only processing subpart==a for p==3"

#Settings

settings = openmc.Settings()
total_batches = 10
settings.batches = total_batches
settings.inactive = 5
settings.particles = 100000

bounds = [-0.4, -0.4, -0.4, 0.4, 0.4, 0.4]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],only_fissionable=True)
settings.source = openmc.source.Source(space=uniform_dist)
settings.temperature={'tolerance':10000,'multipole':True}

settings.export_to_xml()

#pitch = [0.93,1.0,1.05,1.25,1.45]
#pitch = [1.4]
boron = [0]
nt = len(pitch)
nb = len(boron)
nl = len(lattice)
temp1 = 900
temp2 = 600
k = np.zeros([nt,nb])


for m in range(nb):
   for j in range(nt):
       for n in range(nl):

            create_material(boron[m],pitch[j],lattice[n],enrichment,moderator,temp1,temp2)
            
            openmc.run(mpi_args=mpi_args)
            
            #make sure the number in this file name matches the number of cycles you are running!!!\n",
            sp = openmc.StatePoint('statepoint.'+str(total_batches)+'.h5')
            # this reads the tally with all reaction rates, not just absorption\n",
            #tally = sp.get_tally(scores=['absorption'])
            # this is the final k-effective of the simulation\n
            k[j,m] = str(sp.k_combined).split('+')[0]
            if proc_tallies:
               post_process_statepoint(sp)
            #os.remove('statepoint.'+str(total_batches)+'.h5')
            os.remove('summary.h5')
            del sp
 
print(k)

plt.figure(figsize=(20,10))
plt.plot(pitch,k[:,0], linewidth=10)
plt.legend(['0ppm Boron'], fontsize=30)
plt.xlabel('Pitch', fontsize=30)
plt.ylabel('k', fontsize=30)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.savefig('keff vs pin pitch.png')
