import sys
import os
import openmc

import json
import pandas as pd

# MATERIALS

breeder_material = openmc.Material()  # lithium lead chemical equation is Pb84.2Li15.8
breeder_material.add_element('Pb', 84.2, percent_type='ao')
# 50% enriched lithium 6, natural percentage of lithium 6 is just 7% 
breeder_material.add_element('Li', 15.8, percent_type='ao', enrichment={{{:Enrich}}}, enrichment_target='Li6', enrichment_type='ao')
breeder_material.set_density('atom/b-cm', 3.2720171e-2)  # around 11 g/cm3


steel = openmc.Material()
steel.set_density('g/cm3', 7.75)
steel.add_element('Fe', 0.95, percent_type='wo')
steel.add_element('C', 0.05, percent_type='wo')

mats = openmc.Materials([breeder_material, steel])


# GEOMETRY

# surfaces
vessel_inner = openmc.Sphere(r=500)
first_wall_outer_surface = openmc.Sphere(r=510)
breeder_blanket_outer_surface = openmc.Sphere(r={{{:OuterWall}}}, boundary_type='vacuum')


# cells
inner_vessel_region = -vessel_inner
inner_vessel_cell = openmc.Cell(region=inner_vessel_region)

first_wall_region = -first_wall_outer_surface & +vessel_inner
first_wall_cell = openmc.Cell(region=first_wall_region)
first_wall_cell.fill = steel

breeder_blanket_region = +first_wall_outer_surface & -breeder_blanket_outer_surface
breeder_blanket_cell = openmc.Cell(region=breeder_blanket_region)
breeder_blanket_cell.fill = breeder_material

universe = openmc.Universe(cells=[inner_vessel_cell, first_wall_cell, breeder_blanket_cell])
geom = openmc.Geometry(universe)


# SIMULATION SETTINGS

# Instantiate a Settings object
sett = openmc.Settings()
sett.batches = 10
sett.inactive = 0
sett.particles = 500
sett.run_mode = 'fixed source'

# Create a DT point source
source = openmc.Source()
source.space = openmc.stats.Point((0, 0, 0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14e6], [1])
sett.source = source


tallies = openmc.Tallies()

# added a cell tally for tritium production
cell_filter = openmc.CellFilter(breeder_blanket_cell)
tbr_tally = openmc.Tally(name='TBR')
tbr_tally.filters = [cell_filter]
tbr_tally.scores = ['(n,Xt)']  # Where X is a wildcard character, this catches any tritium production
# this allows the tally to be recorded per nuclide so we can see which one contributes to tritium production more
tbr_tally.nuclides = [openmc.Nuclide('Li6'), openmc.Nuclide('Li7')] 
tallies.append(tbr_tally)


# Run OpenMC!
model = openmc.model.Model(geom, mats, sett, tallies)

model.export_to_xml()

os.system('openmc')

# open the results file
sp = openmc.StatePoint("statepoint.10.h5")

# access the tally using pandas dataframes
tbr_tally = sp.get_tally(name='TBR')
df = tbr_tally.get_pandas_dataframe()

# prints the contents of the dataframe
df

# sums up all the values in the mean column
tbr_tally_result = df['mean'].sum()

# sums up all the values in the std. dev. column
tbr_tally_std_dev = df['std. dev.'].sum()

TBR_dict = {'TBR': tbr_tally_result,
            'TBR_std': tbr_tally_std_dev
            }

# json_object = json.dumps(TBR_dict)

with open("openmc.out", "w") as outfile:
    outfile.write(tbr_tally_result)


#df_out = pd.DataFrame([tbr_tally_result], columns=['TBR'])
#df_out.to_csv("TBR.csv")
#df_out.to_json("TBR.json")

