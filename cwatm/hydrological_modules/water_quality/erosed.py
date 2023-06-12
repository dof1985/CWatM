# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Name:        Water quality - erosion-sediment module (EROSED)
# Purpose:     simulate total suspended solids (TSS) in rivers
#
# Author:      TT, FS, PB, MS, DF
#
# Created:     20/01/2022
# Copyright:   (c) TT, FS, PB, MS, DF 2022
# -------------------------------------------------------------------------
import numpy as np

from cwatm.management_modules.data_handling import *
from cwatm.management_modules.globals import *
from cwatm.hydrological_modules.water_quality.waterquality_vars import waterquality_vars


class waterquality_erosed(object):
    """
        WATER QUALITY - EROSION-SEDIMENT MODULE TREATMENT

        Note:
        ------

        How to use:
        ------

        Optional input:
        ------

        **Global variables**

        =====================================  ======================================================================  =====
        Variable [self.var]                    Description                                                             Unit
        =====================================  ======================================================================  =====
        Here                                   Here                                                                    --
        =====================================  ======================================================================  =====

        **Functions**
        """
    def __init__(self, model):
        self.var = model.var
        self.model = model
        self.waterquality_vars = waterquality_vars(model)

    def initial(self):
        """
                INITIAL PART OF THE EROSED MODULE

        Sediment yield per grid cell is calculated with the Modified Universal Soil Loss Equation (MUSLE)
        Williams (1995)
        """
        # load initial MUSLE maps
        # exponential function for fraction of rock (exp-0.053 * fcr)
        self.var.fcfr = 1

        # K_usle: USLE soil erodibility factor
        self.var.kFactor = loadmap('kFactor')

        # C_usle: USLE cover and management factor
        self.var.cFactor = loadmap('cFactor')

        # ls_usle: USLE topographic factor (slope-length)
        self.var.lsFactor = loadmap('lsFactor')

        # manning overland rougness: values for landcoverclasses from settingsfile
        # do not forget to add reference for chosen values
        overlandManningVars = ['manForest', 'manGrassland', 'manirrPaddy', 'manirrNonPaddy']
        self.var.manOverland = np.tile(globals.inZero, (4, 1))
        i = 0
        for variable in overlandManningVars:
            self.var.manOverland[i] += loadmap(variable)
            i += 1
        # manningsN channel
        self.var.manNChan = loadmap('chanMan')
        # grid slope
        tanslope = loadmap('tanslope')

        # setting slope >= 0.00001 to prevent 0 value
        # underlying datasets for tanslope and slopelength are derived from different DEMS, to keep in mind
        self.var.tanslope = np.maximum(tanslope, 0.00001)

        # channel flow time of concentration: unrealistic values. substituted wth. self.var.travelTime
        #tch = divideArrays(0.62 * self.var.chanLength * np.power(self.var.manNChan, 0.75), np.power(self.var.cellArea, 0.125) * np.power(self.var.chanGrad, 0.375))
        
        # channel sediment [kg]
        self.var.channel_sed = globals.inZero.copy()
        self.var.channel_sedConc = globals.inZero.copy()
        self.var.outlet_sed = globals.inZero.copy()
        
        # lake reservoirs [kg]
        self.var.resLakeInflow_sed = globals.inZero.copy()
        self.var.resLakeOutflow_sed = globals.inZero.copy()
        self.var.resLake_sed = globals.inZero.copy()
        self.var.resLake_sedConc = globals.inZero.copy()
        #### Is there anyway to check for initial balance - i.e. so all soil_P in kg at time step = 0 == self.var.soil_PConc_total

        # abstraction [kg]
        self.var.channel_sed_Abstracted = globals.inZero.copy()
        self.var.resLake_sed_Abstracted = globals.inZero.copy()
        self.var.groundwater_sed_Abstracted = globals.inZero.copy()
        self.var.domestic_sed_Abstracted = globals.inZero.copy()
        self.var.livestock_sed_Abstracted = globals.inZero.copy()
        self.var.industry_sed_Abstracted = globals.inZero.copy()
        self.var.irrigation_sed_Abstracted = globals.inZero.copy()
        self.var.returnflowIrr_sed = globals.inZero.copy()
        
    def dynamic(self):
        """
        Dynamic part of EROSED module
        """
        '''
        # Modified Universal Soil Erosion (MUSLE) for sediment yield
        # M_(in_land)  = 11.8×〖(Q_surf*q_peak*A_grid)〗^0.56×K×C×P×LS*f_cfr
        # 11.8 & 0.56 -> calibration parameters, call from settingsfile
        # Q_surf surface runoff volume in mm from cwatm
        # q_peak...peak runoff rate (m3/s); a_tc*Q_surf*Agrid/3.6*t_conc
            #a_tc...frac. of daily rain falling in time of concentration
            # t_conc...time of concentration for grid (model variable) hour
        #if self.var.a
        #self.var.runoffEnergyFactor = self.var.sum_directRunoff * 2
        '''



        self.waterquality_vars.dynamic() # TO FIX
        directRunoff_m3sec = self.var.directRunoff[0:4] * self.var.cellArea / self.var.DtSec
        
        #tov = divideArrays(np.power(self.var.lsFactor, 0.6) * np.power(self.var.manOverland, 0.6),
        #                   18 * np.power(self.var.tanslope, 0.3))
        
        tov = divideArrays(self.var.lsFactor * np.power(self.var.manOverland, 0.6), 3600 * np.power(directRunoff_m3sec, 0.4) * np.power(self.var.tanslope, 0.3))
        tch = self.var.travelTime / 3600  # converted from seconds to hours

        self.var.tconc = tov + tch  # [hours]
        # a05 load dummy value: fraction of daily rain falling in the half-hour highest intensity
        # if time series read netcdf2
        a05 = np.maximum(np.minimum(loadmap('a05'), 1.0), 0.0) #must be a fraction between 0 and 1

        # atc : fraction of rain falling in the time of concentration
        self.var.atc = 1 - np.exp(2 * self.var.tconc * np.log(1 - a05))

        # qpeak: peak runoffrate m3/s
        self.var.qpeak = divideArrays(self.var.atc * self.var.directRunoff[0:4] * self.var.cellArea, 3600. * self.var.tconc) #[m3s-1]

        # MUSLE: sediment yield per day and grid in [1000 kg]

        self.var.sedYieldLand = loadmap('a') * np.power(self.var.directRunoff[0:4] * self.var.qpeak * self.var.cellArea, loadmap('b')) *\
            self.var.kFactor * self.var.cFactor * self.var.lsFactor * self.var.fcfr
        #self.var.sedYieldLand_sum = np.nansum(self.var.fracVegCover[0:4]*self.var.sedYieldLand, axis=0)
        erosedVarsSum = ['sedYieldLand']
        for variable in erosedVarsSum:
            vars(self.var)["sum_" + variable] = np.nansum(vars(self.var)[variable] * self.var.fracVegCover[0:4], axis=0)
        
        
        
        # channel
        self.var.channel_sed_Abstracted = np.maximum(np.minimum(self.var.act_channelAbst * self.var.cellArea *  self.var.channel_sedConc, self.var.channel_sed), 0.)
        