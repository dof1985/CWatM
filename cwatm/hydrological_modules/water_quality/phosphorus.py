# -------------------------------------------------------------------------
# Name:        Water quality - phosphrous module
# Purpose:
#
# Author:      TT, FS, PB, MS, DF
#
# Created:     10/10/2022
# Copyright:   (c) TT, FS, PB, MS, DF 2022
# -------------------------------------------------------------------------

from cwatm.management_modules.data_handling import *
from cwatm.management_modules.globals import *


class waterquality_phosphorus(object):
    """
    WATER QUALITY - PHOSPHORUS MODULE TREATMENT
        
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
    
    def discretizeSoilP(self, Plab, TDP, EPC0, P_in, Qr, Qi, Qp, Kf, soilmass, Vs, runoff_adj):
        
        # Variables:
        # Labile p stocks, Dissolved P stocks, Labile P inputs, P outputs in runoff, interflow, and percolation/groundwater recharge,
        # soil moisture content
        
        
        
        # calculate a and b terms
        
        #a = Plab + P_in
        a =  Kf * soilmass * EPC0 + P_in
        b = divideArrays(Kf * soilmass + Qr + Qi + Qp, Vs)
        
        # calculate new TDP value
        preTDP = TDP.copy()
        prePlab = Plab.copy()
        TDP = (divideArrays(a, b) + (TDP - divideArrays(a, b))* np.exp(-1 * b)) #* self.var.naturalLandFrac# in kg
        
        # Discretized soil labile P
        b0 = b * Vs 
        
        # calculate sorption/absorption
        sorp = Kf * soilmass * (divideArrays(a, b0) -\
            EPC0 + divideArrays(1, b) * (divideArrays(TDP, Vs) - divideArrays(a, b0)) *  (1 -  np.exp(-1 * b)))
        sorp =  np.where(Vs == 0, 0, sorp)
        
        # update Plab
        Plab = Plab + sorp# * self.var.naturalLandFrac
        
        
        # calculate fluxes (in kg)
        P_Qr = (Qr * (divideArrays(preTDP , Vs))) * runoff_adj # calibration parameter runoff_adj > 0 
        P_Qi = Qi * divideArrays(preTDP, Vs)
        P_Qp = Qp * divideArrays(preTDP, Vs)
        
        # Add balance term to Plab
        b_in = P_in
        b_out = P_Qr + P_Qi + P_Qp
        b_sto = Plab - prePlab + TDP - preTDP
        
        bal = b_in - b_out - b_sto
        Plab = Plab + bal
        
        
        
        # calculate dynamic EPC
        EPC0 = divideArrays(Plab, Kf * soilmass)
        
        return Plab, TDP, EPC0, P_Qr, P_Qi, P_Qp
        
    def __init__(self, model):
        self.var = model.var
        self.model = model
	
    def initial(self):
        
        ''' 
        All soil P compartment are numpy arrays nested in lists. 
        The list has 4 elements correspondeing with the following landcover classes:
            * Forest No.0   (Natural)
            * Grasland/managed grassland| non irrigated land No.1 (Semi natural)
            * Paddy irrigation No.2 (cropland)
            * non-Paddy irrigation No.3 (cropland)
        '''
        
        # load initial inactive soil p as a fraction of TP; default - 0.85
        soil_P_fracInactive_managed = globals.inZero.copy() + 0.75
        soil_P_fracInactive_natural = globals.inZero.copy() + 0.9
        if 'fractionInactiveManaged_P' in binding:
            soil_P_fracInactive_managed = loadmap('fractionInactiveManaged_P') # [fraction of TP that is inactive]
        if 'fractionInactiveNatural_P' in binding:
            soil_P_fracInactive_natural = loadmap('fractionInactiveNatural_P') # [fraction of TP that is inactive]
        
        # load soil absorption coefficient Kf | from mm kgsoil-1 to m kgsoil-1
        self.var.Kf = loadmap('Kf') / 1000 
        
        if not self.var.loadInit:
            # load initial total soil p concentration [% of weight]
            # multiply by soil mass [kg] and get mass of initial P per layer in [kg]
            PSoil_init1 = self.var.soilM1 * loadmap('soilTP1_init') / 100
            PSoil_init2 = self.var.soilM2 * loadmap('soilTP2_init') / 100
            PSoil_init3 = self.var.soilM3 * loadmap('soilTP3_init') / 100

            # inactive soil P [kg P] - check if 1e-06 or 1e+06 and also check if it results in per m^2 or per grid cell.
            # natural land - ignores managed grassland for now
            self.var.soil_P_inactive1[0:2] += soil_P_fracInactive_natural * PSoil_init1
            self.var.soil_P_inactive2[0:2] += soil_P_fracInactive_natural * PSoil_init2
            self.var.soil_P_inactive3[0:2] += soil_P_fracInactive_natural * PSoil_init3
            
            # managed land - ignores managed grassland for now
            self.var.soil_P_inactive1[2:4] += soil_P_fracInactive_managed * PSoil_init1
            self.var.soil_P_inactive2[2:4] += soil_P_fracInactive_managed * PSoil_init2
            self.var.soil_P_inactive3[2:4] += soil_P_fracInactive_managed * PSoil_init3

            # labile soil P [kg P]
            # Applied to land cover [1, 2, 3] - natural landcover (forests) initial labile is assumed to be 0. 
            self.var.soil_P_labile1 += (PSoil_init1 - self.var.soil_P_inactive1) * self.var.naturalLandFrac
            self.var.soil_P_labile2 += (PSoil_init2 - self.var.soil_P_inactive2) * self.var.naturalLandFrac
            self.var.soil_P_labile3 += (PSoil_init3 - self.var.soil_P_inactive3) * self.var.naturalLandFrac
            
            # calculate dissolved soil P [kg P /m2]
            '''
                For each landcover i and WQ soil layer j, the soil equilibrium TDP concentration of zero sorption is
                EPC_ij = soil_P_labile_ij / (Kf * soilM_j)
                
                Whereas Kf is Soil P adsorption coefficient, default – 1.1*10-4

                later the TDP is
                self.var.soil_P_dissovled_ij = EPC_ij * wq_soilMoisture_ij 
                
                Where wq_soilMoisture_ij is the mositure content of land cover i and soil layer j  in meters
            '''
            
            # calculate  soil equilibrium TDP concentration of zero sorption  [kg / m]
            self.var.EPC1 = divideArrays(self.var.soil_P_labile1, self.var.Kf * self.var.soilM1)
            self.var.EPC2 = divideArrays(self.var.soil_P_labile2, self.var.Kf * self.var.soilM2)
            self.var.EPC3 = divideArrays(self.var.soil_P_labile3, self.var.Kf * self.var.soilM3)
            
            # calculate TDP [kg P]
            self.var.soil_P_dissolved1 = np.minimum(self.var.EPC1 * self.var.w1, self.var.soil_P_labile1) * self.var.naturalLandFrac
            self.var.soil_P_dissolved2 = np.minimum(self.var.EPC2 * self.var.w2, self.var.soil_P_labile2) * self.var.naturalLandFrac
            self.var.soil_P_dissolved3 = np.minimum(self.var.EPC3 * self.var.w3, self.var.soil_P_labile3) * self.var.naturalLandFrac

            # update labile P to keep P balance
            self.var.soil_P_labile1 -= self.var.soil_P_dissolved1
            self.var.soil_P_labile2 -= self.var.soil_P_dissolved2
            self.var.soil_P_labile3 -= self.var.soil_P_dissolved3
        
        else:
            # if self.var.initLoadFile is set to True
            soil_vars = ["soil_P_inactive1", "soil_P_inactive2", "soil_P_inactive3", "soil_P_labile1", "soil_P_labile2",\
                    "soil_P_labile3", "soil_P_dissolved1", "soil_P_dissolved2", "soil_P_dissolved3"]
            for var in soil_vars:
                for lc_idx in range(4):
                    vars(self.var)[var][lc_idx, :] = self.var.load_initial(var + '_' + str(lc_idx), default = globals.inZero.copy())
                       
        # load data for calculating mineral P weathering and supply to rivers - based on https://doi.org/10.1016/j.chemgeo.2013.10.025
        
        self.var.background_P_mineral = globals.inZero.copy()
        if 'background_P_mineral' in binding:
            self.var.background_P_mineral = loadmap('background_P_mineral')
               
        self.var.soil_shielding = globals.inZero.copy() + 1.
        if 'soil_shielding' in binding:
            self.var.soil_shielding = loadmap('soil_shielding')
        
        self.var.activation_energy = globals.inZero.copy()
        if 'activation_energy' in binding:
            self.var.activation_energy = loadmap('activation_energy')
        
        
        # load groundwater P concentration [kg / m3]
        self.var.GW_P_Conc = globals.inZero.copy()
        if 'GW_P_Conc' in binding:
            self.var.GW_P_Conc = loadmap('GW_P_Conc') * 10 ** 3 
            
        # load runoff P adjustment coefficient
        self.var.runoff_Padj = globals.inZero.copy() + 1.
        if 'runoff_Padj' in binding:
            self.var.runoff_Padj = loadmap('runoff_Padj')
        
        # create empty maps for P source point loading (domestic sewers) and non source point loadings (e.g., irrigation return flow).
        self.var.returnflowNonIrr_P = globals.inZero.copy()
        self.var.returnflowIrr_P = globals.inZero.copy()
        
        # variables for daily Point Source loadings
        self.var.PntSource_NetPload_topsoil = globals.inZero.copy()
        self.var.PntSource_NetPload_channel = globals.inZero.copy()
        self.var.PntSource_NetPload_soil23 = globals.inZero.copy()
        
        # runoff, interflow, baseflow P [kg]
        self.var.directRunoff_P = globals.inZero.copy()
        self.var.interflow_P = globals.inZero.copy()
        self.var.baseflow_P = globals.inZero.copy()
        self.var.sedYieldLand_PP = globals.inZero.copy()
        self.var.sedYieldLand_inactiveP = globals.inZero.copy()
        self.var.mineralWeat_P = globals.inZero.copy()
        
        # channel phsophorus [kg]
        self.var.channel_P = self.var.load_initial('channel_P', default = globals.inZero.copy())
        self.var.channel_PP = self.var.load_initial('channel_PP', default = globals.inZero.copy())
        self.var.channel_inactiveP = self.var.load_initial('channel_inactiveP', default = globals.inZero.copy())
        self.var.channel_PConc = self.var.load_initial('channel_PConc', default = globals.inZero.copy())
        self.var.channel_PPConc = self.var.load_initial('channel_PPConc', default = globals.inZero.copy())
        self.var.channel_inactivePConc = self.var.load_initial('channel_inactivePConc', default = globals.inZero.copy())
        self.var.outlet_P = globals.inZero.copy()
        self.var.outlet_PP = globals.inZero.copy()
        self.var.outlet_inactiveP = globals.inZero.copy()
        
        # lake reservoirs [kg]
        self.var.resLakeInflow_P = globals.inZero.copy()
        self.var.resLakeInflow_PP = globals.inZero.copy()
        self.var.resLakeInflow_inactiveP = globals.inZero.copy()
        self.var.resLakeOutflow_P = globals.inZero.copy()
        self.var.resLakeOutflow_PP = globals.inZero.copy()
        self.var.resLakeOutflow_inactiveP = globals.inZero.copy()
        self.var.resLake_P = self.var.load_initial('resLake_P', default = globals.inZero.copy())
        self.var.resLake_PP = self.var.load_initial('resLake_PP', default = globals.inZero.copy())
        self.var.resLake_inactiveP = self.var.load_initial('resLake_inactiveP', default = globals.inZero.copy())
        self.var.resLake_PConc = self.var.load_initial('resLake_PConc', default = globals.inZero.copy())
        self.var.resLake_PPConc = self.var.load_initial('resLake_PPConc', default = globals.inZero.copy())
        self.var.resLake_inactivePConc = self.var.load_initial('resLake_inactivePConc', default = globals.inZero.copy())
        self.var.resLake_TPConc = globals.inZero.copy()
        self.var.resLake_TDPConc = globals.inZero.copy()
        #### Is there anyway to check for initial balance - i.e. so all soil_P in kg at time step = 0 == self.var.soil_PConc_total

        # abstraction [kg]
        self.var.channel_P_Abstracted = globals.inZero.copy()
        self.var.channel_PP_Abstracted = globals.inZero.copy()
        self.var.channel_inactiveP_Abstracted = globals.inZero.copy()
        self.var.resLake_P_Abstracted = globals.inZero.copy()
        self.var.resLake_PP_Abstracted = globals.inZero.copy()
        self.var.resLake_inactiveP_Abstracted = globals.inZero.copy()
        self.var.groundwater_P_Abstracted = globals.inZero.copy()
        self.var.domestic_P_Abstracted = globals.inZero.copy()
        self.var.livestock_P_Abstracted = globals.inZero.copy()
        self.var.industry_P_Abstracted = globals.inZero.copy()
        self.var.irrigation_P_Abstracted = globals.inZero.copy()
        self.var.irrigation_inactiveP_Abstracted = globals.inZero.copy()
        self.var.returnflowIrr_P = globals.inZero.copy()
        
        # P retention [fraction]
        self.var.channelLake_P_retention = globals.inZero.copy()
        self.var.avg_channelLake_P_retention = globals.inZero.copy()
        
        # In stream/lake sorption/de-sorption
        self.var.n_water = globals.inZero.copy() + 1.
        if 'n_water' in binding:
            self.var.n_water = loadmap('n_water')
        
        # 1.612 - m3/kg <--> 1612 l/kg
        self.var.kf_water = globals.inZero.copy() + 1.612 # [l/kg]
        if 'kf_water' in binding:
            self.var.kf_water = loadmap('kf_water')
        ## initiate all phosphrous stocks -> soil, channel, lakes/reservoirs, groundwater
        ## calculate all conversion factors
        
        # dynamic -> load inputs and ground cover adjustment
        # dynamic_soil, dynamic_...
        
     
        ###############################
        # soil properties
        # load bulkdensity for each soil layer
        
        #create for each soil layer soilmass1 =  bulkdensity1 * cellarea * (thickness of the layer)
        ################################
    
    def dynamic_P_retention(self):
        '''
            Retention is applied as a fraction proportionally to PP, TDP and inactive P in channels, reservoirs and lakes/reservoirs
            
            R = 1 - exp(-(Vf/Hl))
            
            where R is retention fraction, Vf is nutrient uptake velocity  and Hl is hydrological loading
        '''
        
        # calculate water bodies volume (live storage) and depth
        if checkOption('includeWaterBodies'):
            wb_volume = np.where(self.var.waterBodyTypTemp > 0, self.var.lakeResStorage, self.var.substepChannelStorage)
            wb_depth = np.where(self.var.waterBodyTypTemp > 0, divideValues(self.var.lakeResStorage , self.var.lakeArea), self.var.waterLevel)
        else:
            wb_volume = self.var.substepChannelStorage.copy()
            wb_depth = self.var.waterLevel.copy()
            
        # calculate residence time
        r_t = divideValues(wb_volume, self.var.discharge)
        hl = divideValues(wb_depth, r_t)
        vf = 1.411E-06 * 1.06 ** (self.var.waterTemperature - 20)
        
        r_f = 1 - np.exp(-1 * (divideValues(vf, hl)))
        r_f = np.where(self.var.discharge < 0.01, 0., r_f)
        return(r_f)
    
    def dynamic_channel_sorption(self, TDP, PP, Mss, Kf_w, n_w, v, t):
        # function goes here - to be used in routing sub-steps
        '''
        dPP - Sorption/de-soprtion flux [kg / subtimestep]
        TDP - Total dissolved P in channel/lake [kg]
        PP  -  Total particulate P in channel/lake [kg]
        Msss -  Total suspended solids in channel/lake [kg]
        
        Kf_w - Water column sorption coefficient [m3 / kg]
        n_w - Water column Freundlich isotherm constant [unitless]
        v - water volume in channel/lake [m3]
        t - number of sub-timesteps
        
        # check units in doc. and when loading inputs
        '''
        TDPc = divideValues(TDP, v)
        # EPC0_w - Water column equilibrium TDP concentration of zero sorption [kg/m3]
        EPC0_w = divideValues(PP, Kf_w * Mss) ** n_w
        EPC0_w = np.where(Mss <= 1, TDPc, np.where(PP <= 1 , TDPc, EPC0_w))
        EPC0_w = np.where(Kf_w <= 0, TDPc, EPC0_w)
       
        dPP = Kf_w * (TDPc ** (1 / n_w) - EPC0_w ** (1 / n_w)) * v
        
        # restrict by availability of TDP/PP and calculate per sub-timestep
        dPPt = np.where(dPP < 0, -1 * np.minimum(np.abs(dPP / t), PP / t), np.minimum(dPP / t, TDP / t))
        

        # update TDP and PP & return TDP, PP & EPC0_w
        PP += dPPt
        TDP -= dPPt

        return TDP, PP, EPC0_w

    def dynamic(self):
        # phosphrous dynamic part ###
        day_of_year = globals.dateVar['currDate'].timetuple().tm_yday
        wd_date = globals.dateVar['currDate']
        # soil depth ratio - layer 0 out of 0 + 1 : to split P inputs
        soil_depthRatio1 = divideValues(self.var.soildepth[0], self.var.soildepth[0] + self.var.soildepth[1])
    
        # load daily net inputs 
        # should be at [microgram P / m2]
        croplandInputNet = readnetcdf2('P_Cropland_Input', wd_date, useDaily='monthly', value='P_Cropland') / globals.dateVar['daysInMonth']
        croplandInputNet1 = croplandInputNet * soil_depthRatio1 * self.var.cellArea / 10**9 # from microgram to kg
        croplandInputNet2 = croplandInputNet * (1 - soil_depthRatio1) * self.var.cellArea / 10**9
        
        shrManureGrassland = readnetcdf2('shareManure_Grassland', wd_date, useDaily='yearly', value='ManureShare_Grassland')
        
        grasslandInputNet = readnetcdf2('P_Grassland_Input', wd_date, useDaily='yearly', value='P_Grassland') / globals.dateVar['daysInYear']
        
        manure_grasslandInput = grasslandInputNet * shrManureGrassland * self.var.cellArea
        grasslandInputNet1 = grasslandInputNet * (1 - shrManureGrassland) * soil_depthRatio1 * self.var.cellArea
        grasslandInputNet2 = grasslandInputNet * (1 - shrManureGrassland) * (1 - soil_depthRatio1) * self.var.cellArea

        '''
        later change to:
        self.var.soil_PConc_total = readnetcdf2('P_netInput', day_of_year, useDaily='DOY', value='Net_P_Input')
        '''
        
        # read Point source net P loadings - daily
        
        self.var.PntSource_NetPload_topsoil = readnetcdf2('P_sourcePoint', wd_date, useDaily='yearly', value='opendef')
        self.var.PntSource_NetPload_channel = readnetcdf2('P_sourcePoint', wd_date, useDaily='yearly', value='wwtp')
        self.var.PntSource_NetPload_soil23 = readnetcdf2('P_sourcePoint', wd_date, useDaily='yearly', value='latrines')

        pitLatrinesDepth = 3
        if 'pitLatrinesDepth' in binding:
            pitLatrinesDepth = loadmap('pitLatrinesDepth')
        PntSource_toSoilLyr2 = self.var.soildepth[1] + 0.05 >= pitLatrinesDepth
        PntSource_toSoilLyr3 = self.var.soildepth[1] + 0.05 < pitLatrinesDepth
        
        #self.var.PntSource_toSoil =  self.var.PntSource_NetPload_topsoil + self.var.PntSource_NetPload_soil23
        self.var.PntSource_toSoil =   self.var.PntSource_NetPload_soil23
      
        #manureFrac = loadmap('f_Manure')
        
        ## Place holder for P weathering
        self.var.soil_P_weathering = globals.inZero.copy()
        
        # soil dynamic part ###
        
        
        # calculate water movement based on prior timestep concentrations - resulting with mass of layer P per unit area
        '''
        POut_j = (sum of water outflows* from j / soilMoisture_j) * Mass P in join
        * Outflows excluding evaporation
        
        '''
        # calculate relative moisture of soil layer 1 out of 12 - to split interflow of P
        relMoisture2 = divideArrays(self.var.pre_w2, self.var.soildepth[1]) 
        relMoisture3 = divideArrays(self.var.pre_w3, self.var.soildepth[2]) 
        interflowDivider = divideArrays(relMoisture2, relMoisture2 + relMoisture3) #* 0
        
        # Calculate Plab Net Input - Currently only apply on managed grasslands and on irrigated agriculture
        
        self.var.soil_P_input1[2:4] = croplandInputNet1# * self.var.naturalLandFrac
        self.var.soil_P_input2[2:4] = croplandInputNet2# * self.var.naturalLandFrac
        
        # Temporary
        self.var.soil_P_input1[1] = grasslandInputNet1 *  self.var.managedGrassland
        self.var.soil_P_input2[1] = grasslandInputNet2 *  self.var.managedGrassland
        
        # add irrigation to TDP & inactive to inactive
        self.var.soil_P_input1[3] += self.var.sum_irrigation_P_Applied
        self.var.soil_P_inactive1[3] += self.var.sum_irrigation_inactiveP_Applied
        
        # temporary
        pre_lab1 = self.var.soil_P_labile1.copy()
        pre_lab2 = self.var.soil_P_labile2.copy()
        pre_lab3 = self.var.soil_P_labile3.copy()
        
        # pre TDP
        pre_TDP1 = self.var.soil_P_dissolved1.copy()
        pre_TDP2 = self.var.soil_P_dissolved2.copy()
        pre_TDP3 = self.var.soil_P_dissolved3.copy()
        
        #discretizeSoilP(Plab, TDP, EPC0, P_in, Qr, Qi, Qp, Kf, soilmass, Vs)
        '''
        # Only applied for non-natural land covers
        self.var.soil_P_labile1 += self.var.soil_P_input1
        self.var.soil_P_labile2 += self.var.soil_P_input2
        self.var.soil_P_labile3 += 0
        '''
     
        #https://github.com/LeahJB/SimplyP/
        # run dynamic soil P - layer 1 ######
        
        # outputs =  [Plab, TDP, EPC0, P_Qr, P_Qi, P_Qp]
        # self.var.directRunoff[0:4]
        nonNaturalDirectRunoff = self.var.directRunoff[0:4].copy()
        nonNaturalDirectRunoff[0:2] = globals.inZero.copy()
        
        outputs = self.discretizeSoilP(Plab = self.var.soil_P_labile1, TDP = self.var.soil_P_dissolved1,\
            EPC0 =  self.var.EPC1, P_in = self.var.soil_P_input1,\
            # + self.var.PntSource_NetPload_topsoil,\
            Qr = nonNaturalDirectRunoff, Qi = globals.inZero, Qp = self.var.perc1to2,\
            Kf = self.var.Kf, soilmass = self.var.soilM1, Vs = self.var.w1, runoff_adj = self.var.runoff_Padj)

        # update variables
        self.var.soil_P_labile1 =  outputs[0].copy()
        self.var.soil_P_dissolved1 =  outputs[1].copy()
        self.var.EPC1 =  outputs[2].copy()
        
        # update fluxes - no interflow from layer 1
        directRunoff_P =  outputs[3].copy()
        perc1to2_P = outputs[5].copy()
        
        # Add Manure from Grassland to runoff
        directRunoff_P[1] += manure_grasslandInput
        
        # run dynamic soil P - layer 2 ######
        outputs = self.discretizeSoilP(Plab = self.var.soil_P_labile2, TDP = self.var.soil_P_dissolved2,\
            EPC0 =  self.var.EPC2, P_in = self.var.soil_P_input2 + perc1to2_P + np.tile(self.var.PntSource_NetPload_soil23 * PntSource_toSoilLyr2, (4,1)),\
            Qr = globals.inZero, Qi = self.var.interflow[0:4] * interflowDivider, Qp = self.var.perc2to3,\
            Kf = self.var.Kf, soilmass = self.var.soilM2, Vs = self.var.w2, runoff_adj = self.var.runoff_Padj)
        
        # update variables
        self.var.soil_P_labile2 =  outputs[0].copy()
        self.var.soil_P_dissolved2 =  outputs[1].copy()
        self.var.EPC2 =  outputs[2].copy()
        
        # update fluxes - no runoff from layer 2
        interflow2_P =  outputs[4].copy()
        perc2to3_P = outputs[5].copy()
        
        # run dynamic soil P - layer 3 ######
        
        
        outputs = self.discretizeSoilP(Plab = self.var.soil_P_labile3, TDP = self.var.soil_P_dissolved3,\
            EPC0 =  self.var.EPC3, P_in = perc2to3_P + np.tile(self.var.PntSource_NetPload_soil23 * PntSource_toSoilLyr3, (4, 1)),\
            Qr = globals.inZero, Qi = self.var.interflow[0:4] * (1- interflowDivider), Qp = self.var.grossGWrechargeFromSoil[0:4],\
            Kf = self.var.Kf, soilmass = self.var.soilM3, Vs = self.var.w3, runoff_adj = self.var.runoff_Padj)
        
        # update variables
        self.var.soil_P_labile3 =  outputs[0].copy()
        self.var.soil_P_dissolved3 =  outputs[1].copy()
        self.var.EPC3 =  outputs[2].copy()
        
        # update fluxes - no runoff from layer 2
        interflow3_P =  outputs[4].copy()
        toGW = outputs[5].copy()
      
        # PP Delivery to channel: soil flux from EroSed + attched labile
        
        # Enrichment factors based on finer soil praticles

        E_pp = np.where(self.var.sedYieldLand * 1000 > 0.1, np.exp(2.00 - 0.16 * np.log(self.var.sedYieldLand * 0.1)), 1.) # sedYieldLand * 1000 (to kg) / 10000 (per ha)
        
        sedYieldLand_PP = np.maximum(np.minimum(E_pp * self.var.sedYieldLand * 1000 * divideArrays(self.var.soil_P_labile1, self.var.soilM1), self.var.soil_P_labile1) , 0.)
        sedYieldLand_inactiveP = np.maximum(np.minimum(E_pp * self.var.sedYieldLand * 1000 * divideArrays(self.var.soil_P_inactive1, self.var.soilM1), self.var.soil_P_inactive1) , 0.)
        # update soil layer 1 - erosion is only allowed from the top soil
        self.var.soil_P_labile1 -= sedYieldLand_PP
        self.var.soil_P_inactive1 -= sedYieldLand_inactiveP
        # sum outflows
        
        # PP into channel [kg]
        self.var.sedYieldLand_PP = np.nansum(sedYieldLand_PP * self.var.fracVegCover[0:4], axis = 0)
        self.var.sedYieldLand_inactiveP = np.nansum(sedYieldLand_inactiveP * self.var.fracVegCover[0:4], axis = 0)
        
        # runoff [kg]
        self.var.directRunoff_P = np.nansum(directRunoff_P * self.var.fracVegCover[0:4], axis = 0)
        self.var.interflow_P = np.nansum((interflow2_P + interflow3_P) * self.var.fracVegCover[0:4], axis = 0)
        self.var.baseflow_P = self.var.baseflow * self.var.cellArea * self.var.GW_P_Conc 
        
        #self.var.returnflowNonIrr_P = np.where(self.var.returnflowNonIrr > 0, self.var.PntSource_NetPload_channel,0.)
        self.var.returnflowNonIrr_P = np.where(self.var.returnflowNonIrr > 0, self.var.PntSource_NetPload_channel + self.var.PntSource_NetPload_topsoil, 0.)
        
        # to runoff [kg] 
        self.var.runoff_P = self.var.directRunoff_P +  self.var.interflow_P + self.var.baseflow_P + self.var.returnflowNonIrr_P

        # to groundwater [kg]
        self.var.toGroundwater_P =  np.nansum(toGW * self.var.fracVegCover[0:4], axis = 0)
        
        # calculating mineral P weathering and supply to rivers
        # based on https://doi.org/10.1016/j.chemgeo.2013.10.025
        
        self.var.mineralWeat_P = self.var.background_P_mineral * self.var.sum_directRunoff * self.var.soil_shielding *\
            np.exp(-1 * (self.var.activation_energy / 8.3144) * (divideValues(globals.inZero.copy() + 1., self.var.waterTemperature + 273.15) - 1 / 284.15))
        # Update dissolved P mass in soils  - runoff TDP is allowed only when directRunoff > 0
        # to be improved - self.var.runoff_Padj should be also accompanied by max runoff concentration, e.g., how much phosphrous can be removed by the surfacerunoff as a function of runoff volume
      
    ## END SOIL P DYNAMIC
        # calculate soil moisture content for WQ soil layers
        
        # exported phsophorus by outlet [kg]
        self.var.outlet_P = globals.inZero.copy()

        # calculate inputs to dissloved P
        # balance labile/active using EPC
        
                # Water demand ### lift, reservoir type4 is currently excluded
        # Only with self.var.sectorSourceAbstractionFractions = True
        # channel
        if checkOption('includeWaterDemand'):
            self.var.channel_P_Abstracted = np.maximum(np.minimum(self.var.act_channelAbst * self.var.cellArea *  (self.var.channel_PConc / 10**3), self.var.channel_P), 0.)
            self.var.channel_PP_Abstracted = np.maximum(np.minimum(self.var.act_channelAbst * self.var.cellArea * (self.var.channel_PPConc / 10**3), self.var.channel_PP), 0.)
            self.var.channel_inactiveP_Abstracted = np.maximum(np.minimum(self.var.act_channelAbst * self.var.cellArea * (self.var.channel_inactivePConc / 10**3), self.var.channel_inactiveP), 0.)
            self.var.channel_sed_Abstracted = np.maximum(
            np.minimum(self.var.act_channelAbst * self.var.cellArea * self.var.channel_sedConc, self.var.channel_sed),
            0.)
        '''
        # lake/reservoir
        if checkOption('includeWaterBodies'):   
            resLake_P_compress = np.compress(self.var.compress_LR, self.var.resLake_P)
            resLake_P_Abstracted_small = np.minimum(divideValues(resLake_P_compress, self.var.lakeResStorageC + self.var.act_bigLakeAbstC) * (self.var.act_bigLakeAbstC), resLake_P_compress)               
            self.var.resLake_P_Abstracted = globals.inZero.copy()
            np.put(self.var.resLake_P_Abstracted, self.var.decompress_LR, resLake_P_Abstracted_small)
        '''
        
        # groundwater
        self.var.groundwater_P_Abstracted = self.var.nonFossilGroundwaterAbs * self.var.cellArea * self.var.GW_P_Conc 
        
        # domestic P 
        self.var.domestic_P_Abstracted = (self.var.channel_P_Abstracted + self.var.channel_PP_Abstracted) * divideValues(self.var.Channel_Domestic, self.var.domesticDemand) +\
            self.var.resLake_P_Abstracted * divideValues(self.var.Lake_Domestic + self.var.Res_Domestic, self.var.domesticDemand) +\
            self.var.groundwater_P_Abstracted * divideValues(self.var.GW_Domestic, self.var.domesticDemand)
        
        # industry P 
        self.var.industry_P_Abstracted = (self.var.channel_P_Abstracted + self.var.channel_PP_Abstracted) * divideValues(self.var.Channel_Industry, self.var.industryDemand) +\
            self.var.resLake_P_Abstracted * divideValues(self.var.Lake_Industry + self.var.Res_Industry, self.var.industryDemand) +\
            self.var.groundwater_P_Abstracted * divideValues(self.var.GW_Industry, self.var.industryDemand)
            
        # livestock P 
        self.var.livestock_P_Abstracted = (self.var.channel_P_Abstracted + self.var.channel_PP_Abstracted) * divideValues(self.var.Channel_Livestock, self.var.livestockDemand) +\
            self.var.resLake_P_Abstracted * divideValues(self.var.Lake_Livestock + self.var.Res_Livestock, self.var.livestockDemand) +\
            self.var.groundwater_P_Abstracted * divideValues(self.var.GW_Livestock, self.var.livestockDemand)
        
        # irrigation P
        self.var.irrigation_P_Abstracted = (self.var.channel_P_Abstracted + self.var.channel_PP_Abstracted) * divideValues(self.var.Channel_Irrigation, self.var.totalIrrDemand) +\
            self.var.resLake_P_Abstracted * divideValues(self.var.Lake_Irrigation + self.var.Res_Irrigation, self.var.totalIrrDemand) +\
            self.var.groundwater_P_Abstracted * divideValues(self.var.GW_Irrigation, self.var.totalIrrDemand)
        
        self.var.irrigation_inactiveP_Abstracted = (self.var.channel_inactiveP_Abstracted) * divideValues(self.var.Channel_Irrigation, self.var.totalIrrDemand) +\
            self.var.resLake_inactiveP_Abstracted * divideValues(self.var.Lake_Irrigation + self.var.Res_Irrigation, self.var.totalIrrDemand)
            
        # irrigation returnFlows
        self.var.returnflowIrr_P = np.minimum(self.var.irrigation_P_Abstracted  * divideValues(self.var.returnflowIrr, self.var.act_irrWithdrawal), self.var.irrigation_P_Abstracted)
        self.var.returnflowIrr_inactiveP = np.minimum(self.var.irrigation_inactiveP_Abstracted  * divideValues(self.var.returnflowIrr, self.var.act_irrWithdrawal), self.var.irrigation_inactiveP_Abstracted)
        
        self.var.sum_irrigation_P_Applied = self.var.irrigation_P_Abstracted - self.var.returnflowIrr_P
        self.var.sum_irrigation_inactiveP_Applied = self.var.irrigation_inactiveP_Abstracted - self.var.returnflowIrr_inactiveP

        