import sys
import os
import openmc

import json
import pandas as pd
   #MATERIALS#

min = 1
max = 100

R1 = {{{:R1}}}
R2 = 96.723

mats = openmc.Materials()

tungsten = openmc.Material(name='Tungsten', material_id=1)
tungsten.set_density('g/cm3', 19.0)
tungsten.add_element('W', 1.0)
mats.append(tungsten)

breeder_material = openmc.Material(name='PbLi', material_id=2) #Pb84.2Li15.8 with natural enrichment of Li6
# enrichment_fraction = 0.50
enrichment_fraction = {{{:E}}}
breeder_material.add_element('Pb', 84.2,'ao')
breeder_material.add_nuclide('Li6', enrichment_fraction*15.8, 'ao')
breeder_material.add_nuclide('Li7', (1.0-enrichment_fraction)*15.8, 'ao')
breeder_material.set_density('atom/b-cm',3.2720171e-2)
mats.append(breeder_material)

eurofer = openmc.Material(name='EUROFER97',material_id=3)
eurofer.set_density('g/cm3', 7.75)
eurofer.add_element('Fe', 89.067, percent_type='wo')
eurofer.add_element('C', 0.11, percent_type='wo')
eurofer.add_element('Mn', 0.4, percent_type='wo')
eurofer.add_element('Cr', 9.0, percent_type='wo')
eurofer.add_element('Ta', 0.12, percent_type='wo')
eurofer.add_element('W', 1.1, percent_type='wo')
eurofer.add_element('N', 0.003, percent_type='wo')
eurofer.add_element('V', 0.2, percent_type='wo')
mats.append(eurofer)
mats.export_to_xml()
#GEOMETRY#

sphere1 = openmc.Sphere(r=min)
sphere2 = openmc.Sphere(r=R1)
sphere3 = openmc.Sphere(r=R2)
sphere4 = openmc.Sphere(r=max, boundary_type='vacuum')

vac1 = -sphere1
mat1 = +sphere1 & -sphere2
mat2 = +sphere2 & -sphere3
mat3 = +sphere3 & -sphere4
vac2 = +sphere4

vacuum1 = openmc.Cell(region=vac1, cell_id=1)
first = openmc.Cell(region=mat1, cell_id=2)
first.fill = tungsten
second = openmc.Cell(region=mat2, cell_id=3)
second.fill = breeder_material
third = openmc.Cell(region=mat3, cell_id=4)
third.fill = eurofer
vacuum2 = openmc.Cell(region=vac2, cell_id=5)

root = openmc.Universe(cells=(vacuum1, first, second, third, vacuum2))
geom = openmc.Geometry(root)
geom.export_to_xml()

#SETTINGS#

batches = 10
inactive = 0
particles = 1000

source = openmc.Source()
source.space = openmc.stats.Point((0,0,0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14e6], [1])

sett = openmc.Settings()
sett.batches = batches
sett.inactive = inactive
sett.particles = particles
sett.output = {'tallies': False}
sett.run_mode = 'fixed source'
sett.source = source
sett.export_to_xml()

#TALLIES#

tallies = openmc.Tallies()

cell_filter = openmc.CellFilter(second)
tbr_tally = openmc.Tally(name='TBR')
tbr_tally.filters = [cell_filter]
tbr_tally.scores = ['(n,Xt)'] # MT 205 is the (n,Xt) reaction where X is a wildcard, if MT 105 or (n,t) then some tritium production will be missed, for example (n,nt) which happens in Li7 would be missed
tallies.append(tbr_tally)

filter = openmc.SurfaceFilter(sphere4)
leakage_tally = openmc.Tally(name='leakage')
leakage_tally.filters = [filter]
leakage_tally.scores = ['current']
tallies.append(leakage_tally)
tallies.export_to_xml()

model = openmc.model.Model(geom, mats, sett, tallies)

# model.export_to_xml()

os.system('openmc')

# open the results file
sp = openmc.StatePoint("statepoint.10.h5")

# access the tally using pandas dataframes
tbr_tally = sp.get_tally(name='TBR')
df = tbr_tally.get_pandas_dataframe()

# sums up all the values in the mean column
tbr_tally_result = df['mean'].sum()

# sums up all the values in the std. dev. column
tbr_tally_std_dev = df['std. dev.'].sum()


with open("openmc.out", "w") as outfile:
    outfile.write(f"{tbr_tally_result} {tbr_tally_std_dev}")
