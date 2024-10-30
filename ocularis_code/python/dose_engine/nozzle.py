#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 15:43:56 2021

@author: bjoerk_c

This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""


import numpy as np
import pandas as pd


class Nozzle:
    def __init__(self, data, data_path): #, config

        self.data = data
        # self.config = config
        self.data_path = data_path
        self.tune_foil = pd.read_csv(
            data_path / "BaseData" / (data["tune_foil"] + ".csv"))
        self.tune = []
        self.target_range = []
        self.foil = []
        self.modulation = []
        self.modulator_wheel = []
        self.MWwatereq = []
        self.MWweights = []

        self.title = []
        self.inflection_shift_ratio = []
        self.inflection_shift_mm = []
        
        # Linear fit parameters for collimation_modifier
        self.beta_k, self.beta_m = [], []
        
        self.max_range = None

    def set_title(self):
        #self.title = "Tune: {}, Target Range: {} mm, Modulation Range: {} mm, Scatter foil: {}, Modulator wheel: {}".format(self.tune, self.target_range, self.modulation , self.foil, self.modulator_wheel )
        #self.title = "Tune: {}, Target Range: {} mm, Modulation Range: {} mm, Scatter foil: {}, Modulator wheel: {}".format(
        #   "ESS12", self.target_range, self.modulation, self.foil, self.modulator_wheel)
        self.title = "Tune: {}, Proton Range: {} mm, Modulation Range: {} mm, Scatter foil: {}, Modulator wheel: {}".format(self.tune,
            round(self.target_range, 3), round(self.modulation, 3), self.foil, self.modulator_wheel)

    def set_basedata(self):

        # Factor to shift the theoretical lateral infleciton line
        self.inflection_shift_ratio = np.asarray(pd.read_csv(
            self.data_path / "BaseData" / "Inflection_lines_ratio.csv")[self.foil])

        self.inflection_shift_mm = np.asarray(pd.read_csv(
            self.data_path / "BaseData" / "Inflection_shift_mm.csv")[self.foil])

        # k and m describing the linear relationship of the beta parameter over depth z for each scattering foil
        # Beta = k*z = m, where beta = m at isocenter
        # self.beta_k, self.beta_m = np.asarray(pd.read_csv(
        #     self.data_path / "BaseData" / "Beta_Parametrized.csv")[self.foil])
        
        
        # if "linear_betas" in self.config.keys():
        #     print("k and m from linear_betas")
        # Loading data from CCD measurements taken on the 26th of September 2023
        df = pd.read_csv(
            self.data_path / "BaseData" / "linear_betas.csv", index_col=0) #[self.foil]

        self.beta_k, self.beta_m = float(df.loc[df.index == self.foil]["k"]), float(df.loc[df.index == self.foil]["m"])
    
    
        self.max_range = self.tune_foil.loc[self.tune_foil['Foil'] == self.foil]["MaxRange"]

    def config_nozzle(self, config):
        
        # Rounding to 1 decimal place like practiced in the clinic
        self.target_range = round(config["Target_range"], 1)
        self.modulation = round(config["Modulation_range"], 1)

        MWTable_filename = []
        

        for i in range(len(self.tune_foil)):
            if self.target_range <= self.tune_foil.loc[i, :]["MaxRange"] and self.target_range >= self.tune_foil.loc[i, :]["MinRange"]:
                self.tune = self.tune_foil.loc[i, :]["Tune"]
                self.foil = self.tune_foil.loc[i, :]["Foil"]
                MWTable_filename = self.tune_foil.loc[i, :]["MWTable"]
                break

        assert MWTable_filename != []

        self.MWTable = pd.read_csv(
            self.data_path / "BaseData" / (MWTable_filename + ".csv"))

        for i in range(len(self.MWTable) - 1, -1, -1):
            if self.MWTable.loc[i, "MaxMod"] >= self.modulation:
                self.modulator_wheel = self.MWTable.loc[i, "WheelNo"]
                _mw_is_allowed(self.modulator_wheel)
                break
            
        if self.modulator_wheel == []:
            print("Note: Requested modulation range > expected range. Chose MW associated with longest SOBP profile")
            self.modulator_wheel = self.MWTable.loc[0, "WheelNo"]
            _mw_is_allowed(self.modulator_wheel)

        self.MWwatereq = np.asarray(pd.read_csv(
            self.data_path / "BaseData" / (self.data["MWwatereq"] + ".csv"))[str(self.modulator_wheel)])
        self.MWweights = np.asarray(pd.read_csv(
            self.data_path / "BaseData" / (self.data["MWweights"] + ".csv"))[str(self.modulator_wheel)])

        self.set_title()
        self.set_basedata()

    def config_nozzle_by_mw(self, config, wheel):

        self.modulator_wheel = int(wheel)

        _mw_is_allowed(self.modulator_wheel)

        found = False

        for item in self.data["MW"]:
            MWTable = pd.read_csv(
                self.data_path / "BaseData" / (item + ".csv"))
            for i in range(len(MWTable)):
                if self.modulator_wheel == MWTable["WheelNo"][i]:
                    found = True
                    self.modulation = MWTable["MaxMod"][i]
                    MWTable_filename = item
                    break

        assert found, f"Data for wheel {wheel} not found"

        for i in range(len(self.tune_foil)):
            if self.tune_foil.loc[i, :]["MWTable"] == MWTable_filename:
                self.target_range = self.tune_foil.loc[i, :]["MaxRange"]

        for i in range(len(self.tune_foil)):

            if self.target_range <= self.tune_foil.loc[i, :]["MaxRange"] and self.target_range >= self.tune_foil.loc[i, :]["MinRange"]:
                self.tune = self.tune_foil.loc[i, :]["Tune"]
                self.foil = self.tune_foil.loc[i, :]["Foil"]
                break

        self.MWwatereq_all = pd.read_csv(
            self.data_path / "BaseData" / (self.data["MWwatereq"] + ".csv"))
        self.MWweights_all = pd.read_csv(
            self.data_path / "BaseData" / (self.data["MWweights"] + ".csv"))
        self.MWwatereq = np.asarray(
            self.MWwatereq_all[str(self.modulator_wheel)])
        self.MWweights = np.asarray(
            self.MWweights_all[str(self.modulator_wheel)])

        self.set_title()
        self.set_basedata()

    def __repr__(self):
        return self.title


def _mw_is_allowed(mw):

    # allowed = [901, 811, 701, 611, 511, 401, 301, 201, 101] #101	102	103	201	202	301	302	303	401	402	403	404	501	502	503	504	505	511	601	611	602	603	604	605	701	702	703	704	705	706	801	802	803	804	805	806	807	811	901	902	903	904	905	906	907	908


    if mw in [801, 601, 501]:
        print("You specified mw 801, 601 or 501. Did you mean 811, 611, 511 ?")
        raise ValueError

    # if mw not in allowed:
    #     raise ValueError




# mydict = {}
# mylist = []

# for key in keys:
#     mydict[key] = a[key][0]*mag - mag
#     mylist.append(a[key][0])


# my_df = pd.DataFrame([mydict])


# my_df.to_csv(algo.data_path / "BaseData" / "Inflection_shift_mm.csv")
