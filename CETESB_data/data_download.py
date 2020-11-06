#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 22:36:07 2020

@author: mgavidia
"""


import wrf_sp_eval.qualar_py as qr
import pandas as pd
import pickle


# Loading names of CETESB air quality station 
aqs = pd.read_csv("cetesb2017_latlon.dat")

cetesb_login = "XXXXXXXXX"
cetesb_pass = "XXXXXXXXX"

start_date = "06/10/2014"
end_date = "14/10/2014"

# Downloading pollutant data
cetesb_data = {}

for code in aqs.code:
    print(aqs.name[aqs.code == code].values[0], ":", code)
    cetesb_data[aqs.name[aqs.code == code].values[0]] = qr.all_pols(cetesb_login,
                                                                    cetesb_pass,
                                                                    start_date,
                                                                    end_date,
                                                                    code)


a_dict = open("cetesb_pol_2014.pkl", "wb")
pickle.dump(cetesb_data, a_dict)
a_dict.close()


a_dict = open("cetesb_pols_2014.pkl", "rb")
cetesb = pickle.load(a_dict)
a_dict.close()

# Downloading metorological data

cetesb_met = {}
for code in aqs.code:
    print(aqs.name[aqs.code == code].values[0], ":", code)
    cetesb_met[aqs.name[aqs.code == code].values[0]] = qr.all_met(cetesb_login,
                                                                   cetesb_pass,
                                                                   start_date,
                                                                   end_date,
                                                                   code, 
                                                                   in_k=True)
    

a_dict = open("cetesb_met_2014.pkl", "wb")
pickle.dump(cetesb_met, a_dict)
a_dict.close()
