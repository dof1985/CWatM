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
    
    def convert_kc_to_c(self, kc, max_kc, beta):
        fc_est = globals.inZero.copy()
        
        # estimate frac_cover based on kc from Allen, 1998 ().
        fc_est = np.where(kc <= 0.15, 0.01 + 0.09 * (0.15 - (1 - kc)), fc_est)
        fc_est = np.where(np.logical_and(kc > 0.15, kc <= 0.5), 0.1 + 0.25 * (1 - (0.5 - kc)), fc_est)
        fc_est = np.where(np.logical_and(kc > 0.5, kc <= 0.7), 0.35 + 0.15 * (1 - (0.7 - kc)), fc_est)
        fc_est = np.where(kc > 0.7, 0.5 + 0.2 * (1 - (max_kc - kc)), fc_est)
        
        # convert fc to c_factor - Gyssels et al., 2005 (doi: 10.1191/0309133305pp443ra)
        
        c_fact_out = np.exp(beta * 100 * fc_est)
        return(c_fact_out)
        
    
    def sediments_in_channel(self, channel_sed, channel_sedConc, prf, Q, A, csp, spexp, Kch, Cch, Dt):
        # this function is used for sediment routing sub-steps in the channel

        '''
        prf - peak rate factor, atm. to be defined in settingsfile
        Q - discharge [m3/s]
        A - channel crossectional area [m2]
        qChanPeak - peak channel flow rate [m3/s]
        vChanPeak - peak channel flow velocity [m/s]
        concSedMax - maximum sediment transport capacity [kg/m3]
        csp - user defined coefficient
        spexp - user defined exponent, usually between 1-2, set to 1.5 according to original Bagnold stream power equation (Arnold et. al., 1995)
        sedDep - deposition of sediments in channel [kg/sub time step]
        sedDeg - degradation of sediments in channel [kg/sub time step]
        channel_sed - sediment flow in channel [kg/sub time step]
        channel_sedConc - sediment concentration in channel [kg/m3]
        Kch - Channel erodibility factor
        Cch - channel cover factor
        Dt - Number of seconds per routing step
        '''
        #pre_channel_sed = self.var.channel_sed.copy()
        qChanPeak = prf * Q
        vChanPeak = qChanPeak / A #substituted totalCrossArea wth. crossArea calculated in waterquality_vars.py
        #dummyvelocity = divideValues(self.var.travelTime, self.var.chanLength)
        concSedMax = csp * np.power(vChanPeak, spexp)
    

        # Deposition and degreadtion is divided by number of routing steps 
        # sedDt = sedConc * second in time step --> second in day / number of routing time step
        sedDep = Dt * Q * np.where(channel_sedConc > concSedMax, (channel_sedConc - concSedMax), 0.) # kg/m3 -> kg/s -> kg/subtime step
        sedDeg = Dt * Q * np.where(channel_sedConc <= concSedMax, (concSedMax - channel_sedConc) * Kch * Cch, 0.) # kg/subtime step
       
        dchannelSed = sedDeg - sedDep # kg/subtime step
        channel_sed += dchannelSed  # kg/subtime step
        channel_sedConc = divideValues(channel_sed, Dt * Q)
       

        return channel_sed, channel_sedConc, sedDep, sedDeg
    '''    
    def sediments_in_lakes_reservoirs(self, conc_i, conc_eq, ks, d_50, V, t):
        """
        conc_eq ... equilibrium conc. of suspended solids in waterbody (kg/m3)
        conc_i ... initial conc. of suspended solids in waterbody (kg/m3)
        conc_f ... final conc. of suspended solids in waterbody (kg/m3)
        ks ... decay constant (m3/day) ; was (l/day) -> DF
        t ... days of timestep (day)
        d_50 ... median particle size of inflow sediment (um - mikrometer)
        V ... lake/res volume (m3)
        sed_stl ... amount of sediments settled in a day (kg)
        """
        #print(((conc_i - conc_eq) * np.exp(-ks * t * d_50))[conc_i > conc_eq])

        conc_f = np.where(conc_i > conc_eq, (conc_i - conc_eq) * np.exp(-ks * t * d_50) + conc_eq, conc_i)

        sed_stl = (conc_i - conc_f) * V
        mass_f = conc_f * V
       
        return mass_f, sed_stl
    '''
    
    def sediments_in_lakes_reservoirs(self, q0, s_area, v_setl, res_sed, t):
        """
        CALCULATE SEDIMENTATION USING A MODIFIED OVERFLOW RATE MODEL FOLLOWING SWAT
        Returns up-dated res/lake sediments, and sediments settled per sub-time step, and trapping efficiency
        q0 - outflow [m3 / subtimestep]
        s_area - water body surface area in m2
        vsetl - users defined settling velocity [m day-1]
        res_sed - sediments loadings in res/lake [kg/sub timestep]
        t - number of sub-daily timesteps
        conc_eq ... equilibrium conc. of suspended solids in waterbody (kg/m3)

        sed_setl ... amount of sediments settled in a subtime step (kg)
        """
        
        # TRAPPING EFFICIENCY(TE) = v_setl/v_ovflw [day-1]
        
        v_ovflw = divideValues(q0 * t, s_area) # [m]
        # calculate te and restrict between 0 to 1
        te = divideValues(v_setl, v_ovflw) # [day-1]
        te = np.minimum(np.maximum(te, 0), 1)
        
        # Sediment trapping per sub-timestep
        sed_setl = (res_sed * te) / t
        sed_setl = np.where(sed_setl > res_sed, res_sed, sed_setl)

        res_sed -= sed_setl

       
        return res_sed, sed_setl, te

    def initial(self):
        """
                INITIAL PART OF THE EROSED MODULE

        Sediment yield per grid cell is calculated with the Modified Universal Soil Loss Equation (MUSLE)
        Williams (1995)
        """
        # load initial MUSLE maps
        # map with percentage of rock in first soil layer (%) must be provided, e.g., SoilGrids Parameter cfvo for depth 0-5cm
        self.var.CFRG = np.exp(-0.053 * loadmap('rockFrac'))

        i=1

        # K_usle: USLE soil erodibility factor
        self.var.kFactor = loadmap('kFactor')

        # C_usle: USLE cover and management factor
        self.var.cFactor = loadmap('cFactor')
        
        # C_usle  from Kc-factor - monthly
        if 'cfactor_from_kc' in binding  and returnBool('cfactor_from_kc') == True:
            self.var.cfactor_arr = np.tile(globals.inZero,(4,1))
            self.var.c_factor_beta =  -0.048 # https://doi.org/10.1016/j.rse.2018.04.008
            if 'c_factor_beta' in binding:
                self.var.c_factor_beta = loadmap('c_factor_beta')
            
            # max kcmaps
            self.var.max_kcGrassland = readnetcdf2('grassland_cropCoefficientNC', 0, useDaily = "max")
            self.var.max_kcPaddy = readnetcdf2('irrPaddy_cropCoefficientNC', None, "max")
            self.var.max_kcNonPaddy = readnetcdf2('irrNonPaddy_cropCoefficientNC', None, "max")
        
        # p_usle: USLE land conservation factor
        self.var.pFactor = globals.inZero.copy() + 1.
        if 'pFactor' in binding:
            self.var.pFactor = globals.inZero.copy() + loadmap('pFactor')

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
        # grid slope length
        tanslope = loadmap('tanslope')

        # setting slope >= 0.00001 to prevent 0 value
        # underlying datasets for tanslope and slopelength are derived from different DEMS, to keep in mind
        self.var.tanslope = np.maximum(tanslope, 0.00001)

        # channel flow time of concentration: unrealistic values. substituted wth. self.var.travelTime
        # tch = divideArrays(0.62 * self.var.chanLength * np.power(self.var.manNChan, 0.75), np.power(self.var.cellArea, 0.125) * np.power(self.var.chanGrad, 0.375))        
        self.var.sedStor_gridcell = self.var.load_initial('sedStor_gridcell', default = globals.inZero.copy())
        
        if checkOption('includeRunoffConcentration'):
            self.var.sedRunoff_conc = np.tile(globals.inZero,(self.var.maxtime_runoff_conc, 1))
            for i in range(self.var.maxtime_runoff_conc):
                self.var.sedRunoff_conc[i] = self.var.load_initial("sedRunoff_conc", number = i+1)

        self.var.sedToChannel = globals.inZero.copy()
        
        # channel sediment [kg]
        self.var.channel_sed = self.var.load_initial('channel_sed', default = globals.inZero.copy())
        self.var.channel_sed_Dt = self.var.load_initial('channel_sed_Dt', default = globals.inZero.copy()) # channel sed timestep
        self.var.channel_sedConc = self.var.load_initial('channel_sedConc', default = globals.inZero.copy())
        self.var.outlet_sed = globals.inZero.copy()
        self.var.channel_sedDep = globals.inZero.copy()
        # channelbed degradation - input into channels
        self.var.channel_sedDeg = globals.inZero.copy() 
        
        # lake reservoirs [kg]
        self.var.resLakeInflow_sed = globals.inZero.copy()
        self.var.resLakeOutflow_sed = globals.inZero.copy()
        self.var.resLake_sed = self.var.load_initial('resLake_sed', default = globals.inZero.copy())
        self.var.resLake_sedConc = self.var.load_initial('resLake_sedConc', default = globals.inZero.copy())

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
        
        # Sediment loss depth (mm)         
        self.var.sedimentLossDepth_mm = globals.inZero.copy()

        # instream routing
        # channel erodibility factor
        self.var.Kch = globals.inZero.copy() + 0.003
        if 'channel_erodibility' in binding:
            self.var.Kch = globals.inZero.copy() + loadmap('channel_erodibility')
            
        # channel cover factor
        self.var.Cch = globals.inZero.copy() + 0.0015
        if 'channel_cover' in binding:
            self.var.Cch = globals.inZero.copy() + loadmap('channel_cover')
        
        # peak runoff factor
        self.var.prf = globals.inZero.copy() + loadmap('prf')
        # csp
        self.var.csp = globals.inZero.copy() + loadmap('csp')
        # spexp
        self.var.spexp = globals.inZero.copy() + loadmap('spexp')
        
        # sediment delivery ratio
        self.var.sdr_coeff = globals.inZero.copy() + 0.05
        if 'sdr_coeff' in binding:
            self.var.sdr_coeff = globals.inZero.copy() + loadmap('sdr_coeff')
            
        ### Dummy variables for lakes and reservoir function
        if checkOption('includeWaterBodies'):
            self.var.resLakeSedSetlVelocity = np.compress(self.var.compress_LR, globals.inZero.copy() + 0.01)
            
            if 'res_setlVelocity' in binding:
                self.var.resLakeSedSetlVelocity  = np.compress(self.var.compress_LR, globals.inZero.copy() + loadmap('res_setlVelocity'))
            '''    
            if 'ks_sediment' in binding:
                self.var.ks_sed = np.compress(self.var.compress_LR, globals.inZero.copy() + loadmap('ks_sediment')) # day-1 (decay constant)
            else:
                self.var.ks_sed = np.compress(self.var.compress_LR, globals.inZero.copy() + .184)  # day-1 (decay constant)

            if 'd50_sediment' in binding:
                self.var.d50_sed = np.compress(self.var.compress_LR, globals.inZero.copy() + loadmap('d50_sediment'))
            else:
                self.var.d50_sed = np.compress(self.var.compress_LR, globals.inZero.copy() + 32.)

            if 'eq_conc_sediment' in binding:
                self.var.conc_sed_eq = np.compress(self.var.compress_LR, globals.inZero.copy() + loadmap('eq_conc_sediment') / 1000) # mg per l to kg per m3
            else:
                self.var.conc_sed_eq = np.compress(self.var.compress_LR, globals.inZero.copy() + 15) / 1000  # mg per l to kg per m3
            '''






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
        self.var.runoffm3s = self.var.directRunoff[0:4] * self.var.cellArea / self.var.DtSec
        #runoffm3s = self.var.runoff * self.var.cellArea / self.var.DtSec
        
        #self.var.sum_runoffm3s = self.var.sum_directRunoff * self.var.cellArea / self.var.DtSec
        self.var.directRunoff_mm = self.var.directRunoff[0:4] * 1000
        

        # overland flow time of concentration without vov
        #self.var.tov = divideArrays(np.power(self.var.lsFactor, 0.6) * np.power(self.var.manOverland, 0.6), 18 * np.power(self.var.tanslope, 0.3))
        self.var.vov = divideArrays(np.power(self.var.runoffm3s, 0.4) * np.power(self.var.tanslope, 0.3), np.power(self.var.manOverland, 0.6))

        #tov = divideArrays(self.var.lsFactor * np.power(self.var.manOverland, 0.6),
        #                   3600 * np.power(self.var.directRunoff_mm, 0.4) / self.var.DtSec * np.power(self.var.tanslope, 0.3))

        self.var.tov = divideArrays(self.var.slopelength, 3600 * self.var.vov)
        
        self.var.tov = np.where(self.var.runoffm3s > 1, self.var.tov, 24)
        
        #tov2 = divideArrays(np.power(self.var.slopelength, 0.6))
        
        self.var.tch = np.where(self.var.runoffm3s > 1, self.var.travelTime / 3600, 24)  # converted from seconds to hours
        
        #print('mean tch: ', np.nanmean(tch), ' max tch: ', np.nanmax(tch), ' median tch: ', np.median(tch))
        #print('mean tov: ', np.nanmean(tov), ' max tov: ', np.nanmax(tov), ' median tov: ', np.median(tov))
        self.var.tconc = self.var.tov + self.var.tch  # [hours]
        #self.var.tconc = self.var.tch  # [hours]
        
        # a05 load dummy value: fraction of daily rain falling in the half-hour highest intensity
        # if time series read netcdf2
        a05 = np.maximum(np.minimum(loadmap('a05'), 1.0), 0.0)  # must be a fraction between 0 and 1

        # atc : fraction of rain falling in the time of concentration
        self.var.atc = 1 - np.exp(2 * self.var.tconc * np.log(1 - a05))

        # qpeak: peak runoffrate m3/s
        self.var.qpeak = divideArrays(self.var.atc * self.var.directRunoff_mm[0:4] * (self.var.cellArea/10**6), 3.6 * self.var.tconc)  # [m3s-1]
        
        ### CALCULATE C-FACTOR BASED ON KC-valuesof cropland

        if dateVar['newStart'] or (dateVar['currDate'].day in [1,11,21]):
            if 'cfactor_from_kc' in binding  and returnBool('cfactor_from_kc') == True:
                # convert kc of paddy and nonpaddy irr to c factor
                # if irrigation is included No ranges between 0-3; else 0-1; forest = 0 , grassland = 1, paddyIrr = 2, nonPaddyIrr = 3;
                self.var.cfactor_arr[0] = 0.00155 # mean of range from Panagos et al., 2015 (dx.doi.org/10.1016/j.landusepol.2015.05.021)
                
                # grassland combines natural grasslands 0.01 - 0.08 & pasture land 0.05 - 0.15 (Panagos et al., 2015) -- > 0.1  for the unmanaged share; managed share use kc-to-c fucntion
                grassland_cFactor =  self.convert_kc_to_c(kc = self.var.cropKC[1], max_kc = self.var.max_kcGrassland ,beta = self.var.c_factor_beta)
                self.var.cfactor_arr[1] = 0.1 * (1. - self.var.fracManagedGrassland) + self.var.fracManagedGrassland * grassland_cFactor
                
                # irrPaddy
                self.var.cfactor_arr[2] = self.convert_kc_to_c(kc = self.var.cropKC[2], max_kc = self.var.max_kcPaddy ,beta = self.var.c_factor_beta)
                
                # irrNonPaddy
                self.var.cfactor_arr[3] = self.convert_kc_to_c(kc = self.var.cropKC[3], max_kc = self.var.max_kcNonPaddy ,beta = self.var.c_factor_beta)
      
        # MUSLE: sediment yield per day and grid in [1000 kg]
        self.var.sedYieldLand = loadmap('a') * np.power(self.var.directRunoff_mm[0:4] * self.var.qpeak * self.var.cellArea, loadmap('b')) * self.var.kFactor * self.var.cFactor * self.var.lsFactor * self.var.pFactor * self.var.CFRG
        

        # MUSLE: sediment yield per day and grid in [1000 kg]             
        if 'cfactor_from_kc' in binding  and returnBool('cfactor_from_kc') == True:
            self.var.sedYieldLand = loadmap('a') * np.power(self.var.directRunoff_mm[0:4] * self.var.qpeak * self.var.cellArea, loadmap('b')) * self.var.kFactor * self.var.cfactor_arr * self.var.lsFactor * self.var.pFactor * self.var.CFRG
        
        # correct for snow
        self.var.sedYieldLand = divideArrays(self.var.sedYieldLand, np.exp(3 * self.var.SnowCover /  25.4))
        
        # stop sediment yield if frost index > threshold
        self.var.sedYieldLand = np.where(self.var.FrostIndex > self.var.FrostIndexThreshold, 0., self.var.sedYieldLand)
        
        # cap sedimentYield based on volumetric ratio - 0.4 is the top ratio - use bulk density of top soil
        self.var.sedVolumeRatio = globals.inZero.copy() + 0.4
        sedVol =  divideArrays(self.var.sedYieldLand * 1000, self.var.rho1 * self.var.gCm3TokgM3) # ton soil to m3 soil
        sedCap = self.var.directRunoff[0:4]  * self.var.cellArea * self.var.sedVolumeRatio
        sedVolAdj = np.where(sedVol > sedCap, sedCap, sedVol)
        self.var.sedYieldLand = self.var.rho1 * self.var.gCm3TokgM3 * (sedVolAdj / 1000)
        
        # calculate depth of soil loss (mm)
        self.var.sedimentLossDepth_mm = divideValues(self.var.sedYieldLand * np.tile(self.var.soildepth[0], (4, 1)), np.tile(self.var.cellArea,  (4, 1)))
        
        
        # self.var.sedYieldLand_sum = np.nansum(self.var.fracVegCover[0:4]*self.var.sedYieldLand, axis=0)
        #erosedVarsSum = ['sedYieldLand', 'channel_sed', 'channel_sedConc']
        erosedVarsSum = ['sedYieldLand', 'qpeak', 'tconc', 'sedimentLossDepth_mm', 'runoffm3s', 'tov', 'tch', 'directRunoff_mm']
        for variable in erosedVarsSum:
            vars(self.var)["sum_" + variable] = np.nansum(vars(self.var)[variable] * self.var.fracVegCover[0:4], axis=0)
        
        self.var.sedToChannel = (self.var.sum_sedYieldLand * 1000).copy() 
        
        # Calculate sed to channel and lag
        if checkOption('includeRunoffConcentration'):
            runoffConcShare = divideArrays(self.var.runoff_conc, np.nansum(self.var.runoff_conc, axis = 0))
            
            sedToStor = self.var.sedToChannel * runoffConcShare
            
            self.var.sedRunoff_conc = np.roll(self.var.sedRunoff_conc, -1, axis=0)
            self.var.sedRunoff_conc[self.var.maxtime_runoff_conc - 1] = globals.inZero
            self.var.sedRunoff_conc = self.var.sedRunoff_conc + sedToStor
            
            self.var.sedStor_gridcell = self.var.sedStor_gridcell - self.var.sedRunoff_conc[0] + self.var.sedToChannel
            self.var.sedToChannel = self.var.sedRunoff_conc[0, :].copy()
        
        # Sediment delivery ratio

        # apply reduction factor (e.g., deposition or sediment trapped by plants, grasses, etc.)
        # self.var.sdr_coeff  is a calibration factor recommend values are 0.01 - 0.5 (needed to be tested - DF)
        self.var.sdr =  np.exp(-self.var.sdr_coeff * self.var.sum_tconc)
        
        # apply delivery rate
        self.var.sedToChannel = self.var.sdr * self.var.sedToChannel
        
        # as an output variable
        self.var.sum_sedYieldLand_tonha = divideValues(self.var.sum_sedYieldLand, self.var.cellArea * 0.0001)

        #LAKES AND RESERVOIRS
        # detention time (storage/outflow) # is it used ? DF
        if checkOption('includeWaterBodies'):
            self.var.lakeResOutflowM3s = self.var.lakeResOutflowM * self.var.cellArea / 86400
            self.var.detentionTime = self.var.lakeResStorage / self.var.lakeResOutflowM3s