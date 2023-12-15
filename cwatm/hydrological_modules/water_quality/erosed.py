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

    def sediments_in_channel(self, channel_sed, channel_sedConc, prf, Q, A, csp, spexp, V, Kch, Cch):
        # this function is used for sediment routing sub-steps in the channel
        '''
        prf - peak rate factor, atm. to be defined in settingsfile
        Q - discharge [m3/s]
        A - channel crossectional area [m2]
        qChanPeak - peak channel flow rate [m3s-1]
        vChanPeak - peak channel flow velocity [ms-1]
        concSedMax - maximum sediment transport capacity [kg.m-3]
        csp - user defined coefficient
        spexp - user defined exponent, usually between 1-2, set to 1.5 according to original Bagnold stream power equation (Arnold et. al., 1995)
        sedDep - deposition of sediments in channel [kg]
        sedDeg - degradation of sediments in channel [kg]
        channel_sed - sediment in channel [kg]
        channel_sedConc - sediment concentration in channel [kg.m-3]
        Kch - Channel erodibility factor
        Cch - channel cover factor
        '''
        qChanPeak = prf * Q
        vChanPeak = qChanPeak / A #substituted totalCrossArea wth. crossArea calculated in waterquality_vars.py
        #print('vcahnpeak', np.nanmean(vChanPeak), 'vchanmin' , np.nanmin(vChanPeak), 'vchanmax' ,np.nanmax(vChanPeak))
        dummyvelocity = divideValues(self.var.travelTime, self.var.chanLength)
        #print('dVmean', np.nanmean(dummyvelocity), 'dvmin', np.nanmin(dummyvelocity), 'dvmax', np.nanmax(dummyvelocity))
        #print('Qcahnpeak', np.nanmean(Q))
        concSedMax = csp * np.power(vChanPeak, spexp)
        #print('concsedmaxmean', np.nanmean(concSedMax), 'concsedmaxmin', np.nanmin(concSedMax), 'concsedmaxmax', np.nanmax(concSedMax))
        sedDep = np.where(channel_sedConc > concSedMax, (channel_sedConc - concSedMax) * V, 0.)
        #print('Seddep', np.nanmean(sedDep))
        sedDeg = np.where(channel_sedConc <= concSedMax, (concSedMax - channel_sedConc) * Kch * Cch * V, 0.)
        #print('SedDeg',np.nanmean(sedDeg))

        #dChanSedConc = divideValues(sedDeg - sedDep, V)
        #channel_sedConc += dChanSedConc
        dchannelSed = sedDeg - sedDep
        channel_sed += dchannelSed
        channel_sedConc = divideValues(channel_sed,V)
        return channel_sed, channel_sedConc

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

        # slope length = 50 (Malago et al., 2018 and Vigiak et al., 2015)
        self.var.slopelength = globals.inZero.copy() + 50.

        # manning overland roughness: values for landcoverclasses from settingsfile
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
        # tch = divideArrays(0.62 * self.var.chanLength * np.power(self.var.manNChan, 0.75), np.power(self.var.cellArea, 0.125) * np.power(self.var.chanGrad, 0.375))

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

        # instream routing
        # channel erodibility factor
        self.var.Kch = globals.inZero.copy() + 0.8
        # channel cover factor
        self.var.Cch = globals.inZero.copy() + 0.9
        # peak runoff factor
        self.var.prf = globals.inZero.copy() + loadmap('prf')
        # csp
        self.var.csp = globals.inZero.copy() + loadmap('csp')
        # spexp
        self.var.spexp = globals.inZero.copy() + loadmap('spexp')


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
        
        
        self.waterquality_vars.dynamic()  # TO FIX
        directRunoff_m3sec = self.var.directRunoff[0:4] * self.var.cellArea / self.var.DtSec
        runoffm3s = self.var.runoff * self.var.cellArea / self.var.DtSec

        self.var.directRunoff_mm = self.var.directRunoff[0:4]/1000
        


        #tov = divideArrays(np.power(self.var.lsFactor, 0.6) * np.power(self.var.manOverland, 0.6),
        #                  18 * np.power(self.var.tanslope, 0.3))
        vov = divideArrays(np.power(runoffm3s, 0.4)*np.power(self.var.tanslope, 0.3), np.power(self.var.manOverland, 0.6))

        tov = divideArrays(self.var.lsFactor * np.power(self.var.manOverland, 0.6),
                           3600 * np.power(self.var.directRunoff_mm, 0.4) * np.power(self.var.tanslope, 0.3))
        tov2 = divideArrays(self.var.slopelength, 3600 * vov)
        
        
        #tov2 = divideArrays(np.power(self.var.slopelength, 0.6))
        
        tch = self.var.travelTime / 3600  # converted from seconds to hours
        self.var.tconc = tov + tch  # [hours]
        
        # a05 load dummy value: fraction of daily rain falling in the half-hour highest intensity
        # if time series read netcdf2
        a05 = np.maximum(np.minimum(loadmap('a05'), 1.0), 0.0)  # must be a fraction between 0 and 1

        # atc : fraction of rain falling in the time of concentration
        self.var.atc = 1 - np.exp(2 * self.var.tconc * np.log(1 - a05))

        # qpeak: peak runoffrate m3/s
        self.var.qpeak = divideArrays(self.var.atc * self.var.directRunoff_mm[0:4] * self.var.cellArea,
                                      3600. * self.var.tconc)  # [m3s-1]

        # MUSLE: sediment yield per day and grid in [1000 kg]
        self.var.sedYieldLand = loadmap('a') * np.power(self.var.directRunoff_mm[0:4] * self.var.qpeak * self.var.cellArea,loadmap('b')) * self.var.kFactor * self.var.cFactor * self.var.lsFactor * self.var.fcfr

        # self.var.sedYieldLand_sum = np.nansum(self.var.fracVegCover[0:4]*self.var.sedYieldLand, axis=0)
        erosedVarsSum = ['sedYieldLand', 'channel_sed', 'channel_sedConc']
        for variable in erosedVarsSum:
            vars(self.var)["sum_" + variable] = np.nansum(vars(self.var)[variable] * self.var.fracVegCover[0:4], axis=0)
        
        # CHANNEL
        if checkOption('includeWaterDemand'):
            self.var.channel_sed_Abstracted = np.maximum(
            np.minimum(self.var.act_channelAbst * self.var.cellArea * self.var.channel_sedConc, self.var.channel_sed),
            0.)
        
        # as an output variable
        self.var.sum_sedYieldLand_tonha =  divideValues(self.var.sum_sedYieldLand, self.var.cellArea * 0.0001)
        
        erosed_debug = False
        if erosed_debug:
            print('runoffm3s', np.nansum(runoffm3s), 'runoffmin', np.nanmin(runoffm3s) , 'runoffmax', np.nanmax(runoffm3s))
            #print('runoffm3s', np.nanmean(directRunoff_m3sec), 'runoffmin', np.nanmin(directRunoff_m3sec) , 'runoffmax', np.nanmax(directRunoff_m3sec))
            #print('directRunoffm', self.var.directRunoff)
            #print('directRunoffmm', directRunoff_mm)
            print('vov', np.nanmean(vov), 'vovmin', np.nanmin(vov) , 'vovmax', np.nanmax(vov))
            print('Tov', np.nanmean(tov), 'tovmin', np.nanmin(tov) , 'tovmax', np.nanmax(tov))
            print('Tov2', np.nanmean(tov2), 'tov2min', np.nanmin(tov2) , 'tov2max', np.nanmax(tov))
            #print('maxtimerunoff', np.nanmean(self.var.maxtime_runoff_conc))
            #print('tchmean', np.nanmean(tch), 'tchmin', np.nanmin(tch), 'tchmax',np.nanmax(tch))
            #print('tconcmean', np.nanmean(self.var.tconc), 'tconcmin', np.nanmin(self.var.tconc), 'tconcmax', np.nanmax(self.var.tconc))
            
            ### Stop run at step 282 - for debugging purposes
            if (dateVar['curr'] == 282):
               ii=1

        

