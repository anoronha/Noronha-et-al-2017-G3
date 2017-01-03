# -*- coding: utf-8 -*-
"""
Created on Tuesday, January 3, 2017 10:31 CST

@author: anoronha

Implemntation of stalagmite calcite deposition rate model described by
Dreybrodt and others.  See Noronha et al., 2017, G3
dx.doi.org/10.1002/2016GC006644 for discussion of this
implementation and relevant references.  Deposition model applied to
measured values from Jinpasan Cave, Guam.
"""

from pylab import *
import pandas as pd
import numpy as np
import os
import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
from math import *

delta = 0.01 #thin film thickness in cm
drop_volume = 0.07
molar_mass_Ca = 40.08
molar_mass_calcite = 100.09
sec_per_min = 60.0
min_per_hour = 60.0
hour_per_day = 24.0
seconds_per_day = sec_per_min * min_per_hour * hour_per_day
datpath = "/Users/alexandramagana/Data Scripting/Python/GWA Wells"

JinapsanDripwater = pd.read_pickle('JinapsanDripwater.pickle')
JinapsanPlateCalcite = pd.read_pickle('JinapsanPlateCalcite.pickle')

#generate input files to pass to phreeqc
def generate_input_files(drip_data):
    input_files = []
    index_list = []
    selected_output = 'SELECTED_OUTPUT\n-molalities Ca+2 HCO3-\n-activities H+ Ca+2 CO2 HCO3- CO3-2'
    drip_data = drip_data[(np.isnan(drip_data.Temp) == False) & \
                          (np.isnan(drip_data.CO2) == False)]
    for i in range(0,len(drip_data)):
        index_rec = drip_data['index'].iloc[i]
        solution = 'SOLUTION %s\ntemp %5.2f\nunits ppm\npH 8.0' \
                    %(i, drip_data.Temp.iloc[i])
        eql_phase = 'EQUILIBRIUM_PHASES\nCO2(g) %5.2f\nCalcite 0.0' \
                    %(log10(drip_data.CO2.iloc[i]/(1e6)))
        ind_file = '%s\n%s\n%s\nEND' \
                    %(solution, eql_phase, selected_output)
        input_files.append(ind_file)
        index_list.append(index_rec)
    return input_files, index_list


#calculate drip water [Ca] in equilbirium with cave pCO2 using phreeqc via phreeqcpy
def calc_co2(input_files, index_list, data):
    def selected_array(db_path, input_string):
        dbase = phreeqc_mod.IPhreeqc("/usr/local/lib/libiphreeqc.dylib")
        dbase.load_database(db_path)
        dbase.run_string(input_string)
        return dbase.get_selected_output_array()
    results_table = pd.DataFrame(columns = ['index', 'Ca_eq'])
    for i in range(0,len(input_files)):
        result = selected_array(os.path.join(datpath, 'phreeqc.dat'),
                                input_files[i])
        dicts = {'index': index_list[i],
                 'Ca_eq': result[2][8]}
        results_table = results_table.append([dicts])
    return results_table

#calcluate calcite deposition rate
def calc_deposition(water_data):
    GR_holder = []
    delta_Ca_holder = []
    Ca_initial_holder = []
    for i in range(0,len(water_data)):
        #Equation from Hansen2013 pg 244, units of alpha in cm/s
        alpha_p = (0.52 + 0.04 * water_data.Temp.iloc[i] + \
                0.004 * water_data.Temp.iloc[i]**2) * 1e-5
        #convert from ppb to mol/L
        Ca_initial = water_data.Ca.iloc[i]/(1e6 * molar_mass_Ca)
        Ca_initial_holder.append(Ca_initial)
        #calcuate c - c_eq, mmol/cm^3 = mol/L
        delta_Ca = (Ca_initial - water_data.Ca_eq.iloc[i])
        delta_Ca_holder.append(delta_Ca)
        if delta_Ca < 0:
            delta_Ca = 0
        #drip interval units of s, W_0 has units of mmol/cm^2*s
        W_0 = delta_Ca * (delta/water_data.DripInterval.iloc[i]) * \
              (1 - exp((-alpha_p/delta) * water_data.DripInterval.iloc[i]))
        #cm^3/cm = cm^2
        area = water_data.DripVol.iloc[i] / delta
        #convert from mmol/cm^2*s to mg/day
        GR = W_0 * molar_mass_calcite * seconds_per_day * area
        GR_holder.append(GR)
    water_data['GrowthRate_Modeled'] = GR_holder
    water_data['Delta_Ca'] = delta_Ca_holder
    water_data['Ca_initial'] = Ca_initial_holder
    return water_data

JinapsanDripwater = JinapsanDripwater.reset_index(drop = False)
input_files, index_list = generate_input_files(JinapsanDripwater)
results_table = calc_co2(input_files, index_list, JinapsanDripwater)
JinapsanDripwater = pd.merge(JinapsanDripwater,
                              results_table,
                              how='left',
                              on='index')
GrowthModeled = calc_deposition(JinapsanDripwater)
GrowthModeled.to_pickle('GrowthModeled.pickle')
