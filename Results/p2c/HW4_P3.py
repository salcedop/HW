#!/usr/bin/env python

import openmc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import pandas as pd
import argparse

#material 1

def build_first_model(ppm):
            rfuel=0.5
            r2 = rfuel + 0.01
            r3 = rfuel + 0.06
            lattice="rectangular"
            pitch = 1.25
            uo2 = openmc.Material(1,"fuel",temperature=900)
            uo2.add_element('U', 1.0, enrichment=0.7)
            uo2.add_element('O', 2.0)
            uo2.set_density('g/cc', 10.0)
            zirconium = openmc.Material(2,"zirconium", temperature=600)
            zirconium.add_element('Zr', 1.0),
            zirconium.set_density('g/cm3', 6.6),
            
            mod = openmc.Material(3,'heavy water',temperature=600)
            mod.add_nuclide('H2',2.0)
            mod.add_nuclide('O16',1.0)
            mod.add_s_alpha_beta('c_D_in_D2O')
            mod.set_density('g/cm3',1.1)

            #adding Dysprosium Oxide Dy2O3
            #dy2o3 = openmc.Material(8,"crit_mat",temperature=900)
            #dy2o3.add_element('Dy',2.0*ppm*1E-6)
            #dy2o3.add_element('O',3.0*ppm*1E-6)
            uo2.add_element('Dy',2.0*ppm*1E-6)
            uo2.add_element('O',3.0*ppm*1E-6)
         
            materials = openmc.Materials([uo2,zirconium,mod])
            
            fuel_or = openmc.ZCylinder(R=rfuel)
            clad_ir = openmc.ZCylinder(R=r2)
            clad_or = openmc.ZCylinder(R=r3)
            fuel_region = -fuel_or
            gap_region = +fuel_or & -clad_ir
            clad_region = +clad_ir & -clad_or
            fuel = openmc.Cell(1, 'fuel')
            fuel.fill = uo2 
            #fuel.fill = dy2o3
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
            geometry = openmc.Geometry(root)
                      
            #Settings

            settings = openmc.Settings()
            total_batches = 100
            settings.batches = total_batches
            settings.inactive = 25
            settings.particles = 100000

            bounds = [-rfuel, -rfuel, -rfuel, rfuel, rfuel, rfuel]
            uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],only_fissionable=True)
            settings.source = openmc.source.Source(space=uniform_dist)
            settings.temperature={'tolerance':10000,'multipole':True}

            settings.output = {'tallies': False}

            model = openmc.model.Model(geometry,materials,settings)

            return model

def build_second_model(enrichment):
            rfuel=0.5
            r2 = rfuel + 0.01
            r3 = rfuel + 0.06
            lattice="rectangular"
            pitch = 9
            uo2 = openmc.Material(1,"fuel",temperature=900)
            uo2.add_element('U', 1.0, enrichment=enrichment)
            uo2.add_element('O', 2.0)
            uo2.set_density('g/cc', 10.0)
            zirconium = openmc.Material(2,"zirconium", temperature=600)
            zirconium.add_element('Zr', 1.0),
            zirconium.set_density('g/cm3', 6.6),
            #this function creates borated light water, for other moderators you will need to replace\n",
            mod = openmc.Material(name='Buffer',temperature=600)
            mod.set_density('g/cm3', 1.0)
            mod.add_element('C', 1.0)
            mod.add_s_alpha_beta('c_Graphite')
          
            materials = openmc.Materials([uo2, zirconium, mod])
            
            fuel_or = openmc.ZCylinder(R=rfuel)
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
            geometry = openmc.Geometry(root)
                      
            #Settings

            settings = openmc.Settings()
            total_batches = 100
            settings.batches = total_batches
            settings.inactive = 25
            settings.particles = 100000

            bounds = [-rfuel, -rfuel, -rfuel, rfuel, rfuel, rfuel]
            uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],only_fissionable=True)
            settings.source = openmc.source.Source(space=uniform_dist)
            settings.temperature={'tolerance':10000,'multipole':True}

            settings.output = {'tallies': False}

            model = openmc.model.Model(geometry,materials,settings)

            return model
            
      
       

parser = argparse.ArgumentParser()

parser.add_argument('-p','--problem',choices=('p1','p2'), default='p1',
                   help='problem number')
parser.add_argument('-n','--node',action='store_true',
                   help='Are you using Skylake node?')

args = parser.parse_args()
#set MPI argument
if args.node:

   mpi_args=['mpiexec','-bind-to','numa','-n','2']
else:
   mpi_args=['mpiexec','-n','2']

if (args.problem == 'p1'):

  crit_var,guess,keffs = openmc.search_for_keff(build_first_model,initial_guess=1500,
                                              tol=1.E-2,print_iterations=True)

else:

   crit_var,guess,keffs = openmc.search_for_keff(build_second_model,initial_guess=1.,
                                              tol=1.E-2,print_iterations=True)


#crit_var,guess,keffs = openmc.search_for_keff(build_first_model,bracket=[1900,2800],
#                                              tol=1.E-2,bracketed_method='bisect',
#                                              print_iterations=True)

#crit_var,guess,keffs = openmc.search_for_keff(build_second_model,initial_guess=8.,
                                              #tol=1.E-2,print_iterations=True)


print('Critical Variable Concentration: '+str(crit_var))
openmc.run(mpi_args=mpi_args)
            
os.remove('summary.h5')
