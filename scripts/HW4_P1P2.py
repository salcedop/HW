import openmc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import pandas as pd
import argparse
import math


#material 1

def create_material(boron_conc,pitch,lattice,enrichment,moderator,temp1,temp2):
  
            rf = 0.50
            r2 = rf + 0.01
            r3 = rf + 0.06
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
            fuel_or = openmc.ZCylinder(R=rf)
            clad_ir = openmc.ZCylinder(R=r2)
            clad_or = openmc.ZCylinder(R=r3)
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
       factors = np.zeros([4])
       fiss_rate = statepoint.get_tally(name='fiss.rate')
       abs_rate = statepoint.get_tally(name='abs.rate')
       leak = statepoint.get_tally(name='leakage')
       leak = leak.summation(filter_type=openmc.CellFilter,remove_filter = True)
       #get k-infinity
       keff = fiss_rate/(abs_rate + leak)
       df1 = keff.get_pandas_dataframe()
       #resonance escape prob
       therm_abs_rate = statepoint.get_tally(name='therm.abs.rate')
       thermal_leak = statepoint.get_tally(name='thermal leakage')
       thermal_leak = thermal_leak.summation(filter_type=openmc.CellFilter,remove_filter = True)
       res_esc = (therm_abs_rate + thermal_leak) / (abs_rate + thermal_leak)
       df2 = res_esc.get_pandas_dataframe()
       
       #factors[0] = (df2['mean'])
       
       #fast fission factor
       therm_fiss_rate = statepoint.get_tally(name='therm.fiss.rate')
       fast_fiss = fiss_rate / therm_fiss_rate
       df3 = fast_fiss.get_pandas_dataframe()
       #factors[1] = (df3['mean'])
       #thermal utilization factor
       fuel_therm_abs_rate = statepoint.get_tally(name='fuel.therm.abs.rate')
       therm_util = fuel_therm_abs_rate / therm_abs_rate
       df4 = therm_util.get_pandas_dataframe()
       #factors[2] = (df4['mean'])
       #eta
       eta = therm_fiss_rate / fuel_therm_abs_rate
       df5 = eta.get_pandas_dataframe()
       #factors[3] = (df5['mean'])
       frames = [df1,df2,df3,df4,df5]
       #result = pd.concat(frames)
       #result.to_csv('five_factor.csv')

       for i in range(len(factors)):
          factors[i] = frames[i+1]['mean']

       return factors
       

parser = argparse.ArgumentParser()

parser.add_argument('-p','--problem',choices=('p1','p2','p3'), default='p1',
                   help='problem number')
parser.add_argument('-subp','--subproblem',choices=('a','b'),default='a',
                   help='subproblem number')
parser.add_argument('-n','--node',action='store_true',
                   help='Are you using Skylake node?')

args = parser.parse_args()

#pitch = [8.5,9.,9.5,10.,10.5,11.]
#pitch = [8,8.5,9.,9.5,10.]

pitch  = [1.15,1.2,1.3,1.8,2.,2.2,2.4,2.6,2.8,3.] #light-water
#pitch = [3.,6.,9.,12.,15.] #graphite 
#[1.,7.,10.,15.,20.] #heavy water
x = pitch
#set MPI argument
if args.node:

   mpi_args=['mpiexec','-bind-to','numa','-n','2']
else:
   mpi_args=['mpiexec','-n','2']

if (args.problem == 'p1'):

  moderator = 'heavy_water'
  enrichment = None #fresh
  lattice = ['rectangular']
  xtitle = 'pitch (cm)'

  if (args.subproblem == 'a'):
     
     ran = 1
     offset = 0
     proc_tallies = False
     title = 'k_versus_rod_pitch_heavy-water'
     ytitle = 'k'
     
     plot_legend = ['Rectangular lattice']
  else:
     
     ran = 4
     offset = 1
     title = 'Four_factors_versus_rod_pitch_heavy-water'
     ytitle = 'Four factors'
     
     plot_legend = ['esc.prob','fast.fiss','thermal.uti','eta']
     proc_tallies = True

elif(args.problem == 'p2'):
  
  moderator = 'graphite'
  enrichment = None #fresh
  lattice = ['rectangular']
  xtitle = 'pitch (cm)'
  
  if (args.subproblem == 'a'):
     ran = 1
     offset = 0
     title = 'k_versus_rod_pitch_graphite'
     ytitle = 'k'

     plot_legend = ['Rectangular lattice']
     proc_tallies = False
  else:
      
     ran = 4
     offset = 1
     title = 'Four_factors_versus_rod_pitch_graphite'
     ytitle = 'Four factors'
     
     plot_legend = ['esc.prob','fast.fiss','thermal.uti','eta']
     proc_tallies = True
elif(args.problem == 'p3'):

   moderator = 'light_water'
   enrichment = 4
   lattice = ['rectangular','triangular']

   plot_legend = ['Rectangular lattice','Hexagonal lattice']
   ytitle = 'k'
   ran = 1
   offset = 0
   if (args.subproblem == 'a'):
     
     title = 'k_versus_rod_pitch_light-water'
     xtitle = 'pitch (cm)'
     
     proc_tallies = False
   else:
     #raise ValueError("cannot specify combination of p3 and subp b")
     title = 'k_versus_volume_ratio_mod_to_fuel_light-water'
     xtitle = 'Moderator to fuel volumne ratio'
     proc_tallies = False    
     l_pitch = len(pitch)
     x = np.zeros([l_pitch])
     mod_to_fuel_hex = np.zeros([l_pitch])
     for ii in range(l_pitch):
         x[ii] = ((pitch[ii]**2)/(math.pi*(0.5**2)))-1
         mod_to_fuel_hex[ii] = ((3*math.sqrt(3)/2)*pitch[ii]**2)/(math.pi*0.5**2)-1
     
#Settings

settings = openmc.Settings()
total_batches = 10
settings.batches = total_batches
settings.inactive = 5
settings.particles = 100000

bounds = [-0.5, -0.5, -0.5, 0.5, 0.5, 0.5]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],only_fissionable=True)
settings.source = openmc.source.Source(space=uniform_dist)
settings.temperature={'tolerance':10000,'multipole':True}

settings.export_to_xml()

#pitch = [0.93,1.0,1.05,1.25,1.45]
#pitch = [1.4]
boron = [0]
nt = len(pitch)
nl = len(lattice)
temp1 = 900
temp2 = 600
#k = np.zeros([nt,nl])
factors = np.zeros([nt,nl,5])

for j in range(nt):
    for n in range(nl):

           create_material(0,pitch[j],lattice[n],enrichment,moderator,temp1,temp2)
            
           openmc.run(mpi_args=mpi_args)
            
           sp = openmc.StatePoint('statepoint.'+str(total_batches)+'.h5')
           factors[j,n,0] = str(sp.k_combined).split('+')[0]
           if proc_tallies:
               factors[j,n,1:5] = post_process_statepoint(sp)
               
            #os.remove('statepoint.'+str(total_batches)+'.h5')
           os.remove('summary.h5')
           del sp
 
#print(k)
print(factors)

plt.figure(figsize=(20,10))

for i in range(ran):
  for j in range(nl):
      if ((j != 0) and (args.subproblem == 'b')):
         x = mod_to_fuel_hex
      plt.plot(x,factors[:,j,i+offset], linewidth=5)

plt.legend(plot_legend, fontsize=12)

plt.title(title,fontsize=30)
plt.xlabel(xtitle, fontsize=30)
plt.ylabel(ytitle, fontsize=30)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig(title+'_'+args.problem+'.png') 
