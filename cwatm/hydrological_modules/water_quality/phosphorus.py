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
    
    def soilP(self, Plab, TDP, EPC0, Plab_in, TDP_in, Qr, Qi, Qp, Kf, soilmass, Vs):
        
        # Variables:
        # Labile p stocks, Dissolved P stocks, Labile P inputs, TDP inputs, P outputs in runoff, interflow, and percolation/groundwater recharge,
        # soil moisture content (t-1)
        
        
        
        # calculate out fluxes [kg]
        
        totalFlux = np.minimum((Qr + Qi + Qp) * divideArrays(TDP, Vs), TDP)
        
        P_Qr = divideArrays(Qr, Qr + Qi + Qp) * totalFlux
        P_Qi = divideArrays(Qi, Qr + Qi + Qp) * totalFlux
        P_Qp = divideArrays(Qp, Qr + Qi + Qp) * totalFlux
        
        # calculate intermediate TDP and Plab
        TDP = TDP + TDP_in - P_Qr - P_Qi - P_Qp
        Plab = Plab + Plab_in
        
        # calculate EPC0
        EPC0 = divideArrays(Plab, Kf * soilmass)
        
        # calculate sorption
        sorp =  Kf * (TDP -  EPC0)
        
        # Update variables
        Plab = Plab + sorp
        TDP = TDP - sorp

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
        # load initial total soil p concentration [% of weight]
        # multiply by soil mass [kg] and get mass of initial P per layer in [kg]
        PSoil_init1 = self.var.soilM1 * loadmap('soilTP1_init') / 100
        PSoil_init2 = self.var.soilM2 * loadmap('soilTP2_init') / 100
        PSoil_init3 = self.var.soilM3 * loadmap('soilTP3_init') / 100
        
        # load initial inactive soil p as a fraction of TP; default - 0.85
        soil_P_fracInactive = globals.inZero.copy() + 0.85
        if 'fractionInactive_P' in binding:
            soil_PConc_inactive = loadmap('fractionInactive_P') # [fraction of TP that is inactive]
        
        # load soil absorption coefficient Kf | from mm kgsoil-1 to m kgsoil-1
        self.var.Kf = loadmap('Kf') / 1000 
        
        # inactive soil P [kg P] - check if 1e-06 or 1e+06 and also check if it results in per m^2 or per grid cell.
        self.var.soil_P_inactive1 += soil_PConc_inactive * PSoil_init1
        self.var.soil_P_inactive2 += soil_PConc_inactive * PSoil_init2
        self.var.soil_P_inactive3 += soil_PConc_inactive * PSoil_init3

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
        
        # runoff, interflow, baseflow P [kg]
        self.var.directRunoff_P = globals.inZero.copy()
        self.var.interflow_P = globals.inZero.copy()
        self.var.baseflow_P = globals.inZero.copy()

        # channel phsophorus [kg]
        self.var.channel_P = globals.inZero.copy()
        self.var.channel_PConc = globals.inZero.copy()
        
        # lake reservoirs [kg]
        self.var.resLakeInflow_P = globals.inZero.copy()
        self.var.resLakeOutflow_P = globals.inZero.copy()
        self.var.resLake_P = globals.inZero.copy()
        self.var.resLake_PConc = globals.inZero.copy()
        #### Is there anyway to check for initial balance - i.e. so all soil_P in kg at time step = 0 == self.var.soil_PConc_total

        # abstraction [kg]
        self.var.channel_P_Abstracted = globals.inZero.copy()
        self.var.resLake_P_Abstracted = globals.inZero.copy()
        self.var.groundwater_P_Abstracted = globals.inZero.copy()
        self.var.domestic_P_Abstracted = globals.inZero.copy()
        self.var.livestock_P_Abstracted = globals.inZero.copy()
        self.var.industry_P_Abstracted = globals.inZero.copy()
        self.var.irrigation_P_Abstracted = globals.inZero.copy()
        self.var.returnflowIrr_P = globals.inZero.copy()
        
        ## initiate all phosphrous stocks -> soil, channel, lakes/reservoirs, groundwater
        ## calculate all conversion factors
        
        # dynamic -> load inputs and ground cover adjustment
        # dynamic_soil, dynamic_...
        
     
        ###############################
        # soil properties
        # load bulkdensity for each soil layer
        
        #create for each soil layer soilmass1 =  bulkdensity1 * cellarea * (thickness of the layer)
        ################################
       
    def dynamic(self):
        # phosphrous dynamic part ###
        day_of_year = globals.dateVar['currDate'].timetuple().tm_yday
        
        # soil depth ratio - layer 0 out of 0 + 1 : to split P inputs
        soil_depthRatio1 = divideValues(self.var.soildepth[0], self.var.soildepth[0] + self.var.soildepth[1])
    
        # load daily net inputs 
        # should be at [kg P / m2]
        PSoil_Input1 = loadmap('P_netInput') * soil_depthRatio1 * self.var.cellArea
        PSoil_Input2 = loadmap('P_netInput') * (1 - soil_depthRatio1) * self.var.cellArea
        '''
        later change to:
        self.var.soil_PConc_total = readnetcdf2('P_netInput', day_of_year, useDaily='DOY', value='Net_P_Input')
        '''
        manureFrac = loadmap('f_Manure')
        
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
        
        self.var.soil_P_input1[2:4] = PSoil_Input1# * self.var.naturalLandFrac
        self.var.soil_P_input2[2:4] = PSoil_Input2# * self.var.naturalLandFrac
        
        # Temporary
        self.var.soil_P_input1[1] = PSoil_Input1 *  self.var.managedGrassland
        self.var.soil_P_input2[1] = PSoil_Input2 *  self.var.managedGrassland
        
        # temporary
        pre_lab1 = self.var.soil_P_labile1.copy()
        pre_lab2 = self.var.soil_P_labile2.copy()
        pre_lab3 = self.var.soil_P_labile3.copy()
        
        # pre TDP
        pre_TDP1 = self.var.soil_P_dissolved1.copy()
        pre_TDP2 = self.var.soil_P_dissolved2.copy()
        pre_TDP3 = self.var.soil_P_dissolved3.copy()
        
        
        # Apply irrigation 
        self.var.irrigation_P_Applied = np.maximum(self.var.irrigation_P_Abstracted - self.var.returnflowIrr_P, 0.) * divideArrays(self.var.infiltration, self.var.availWaterInfiltration)[0:4] * (self.var.onlyIrrPaddy + self.var.onlyIrr)
        shareInfiltrationTopSoil = divideArrays(self.var.infiltration -  self.var.infiltration2, self.var.infiltration)[0:4]
        
        #SoilP(Plab, TDP, EPC0,  Plab_in, TDP_in, Qr, Qi, Qp, Kf, soilmass, Vs)
        # run dynamic soil P - layer 1 ######
        
        # outputs =  [Plab, TDP, EPC0, P_Qr, P_Qi, P_Qp]
        # self.var.directRunoff[0:4]
        nonNaturalDirectRunoff = self.var.directRunoff[0:4].copy()
        nonNaturalDirectRunoff[0:2] = globals.inZero.copy()
        outputs = self.soilP(Plab = self.var.soil_P_labile1, TDP = self.var.soil_P_dissolved1,\
            EPC0 =  self.var.EPC1, Plab_in = self.var.soil_P_input1, TDP_in = self.var.irrigation_P_Applied * shareInfiltrationTopSoil,\
            Qr = nonNaturalDirectRunoff, Qi = globals.inZero, Qp = self.var.perc1to2,\
            Kf = self.var.Kf, soilmass = self.var.soilM1, Vs = self.var.pre_w1)

        # update variables
        self.var.soil_P_labile1 =  outputs[0].copy()
        self.var.soil_P_dissolved1 =  outputs[1].copy()
        self.var.EPC1 =  outputs[2].copy()
        
        # update fluxes - no interflow from layer 1
        directRunoff_P =  outputs[3].copy()
        perc1to2_P = outputs[5].copy()
        
        # run dynamic soil P - layer 2 ######
        outputs = self.soilP(Plab = self.var.soil_P_labile2, TDP = self.var.soil_P_dissolved2,\
            EPC0 =  self.var.EPC2, Plab_in = self.var.soil_P_input2, TDP_in = perc1to2_P + self.var.irrigation_P_Applied * (1 - shareInfiltrationTopSoil),\
            Qr = globals.inZero, Qi = self.var.interflow[0:4] * interflowDivider, Qp = self.var.perc2to3,\
            Kf = self.var.Kf, soilmass = self.var.soilM2, Vs = self.var.pre_w2)
        
        # update variables
        self.var.soil_P_labile2 =  outputs[0].copy()
        self.var.soil_P_dissolved2 =  outputs[1].copy()
        self.var.EPC2 =  outputs[2].copy()
        
        # update fluxes - no runoff from layer 2
        interflow2_P =  outputs[4].copy()
        perc2to3_P = outputs[5].copy()
        
        # run dynamic soil P - layer 3 ######
         #
        outputs = self.soilP(Plab = self.var.soil_P_labile3, TDP = self.var.soil_P_dissolved3,\
            EPC0 =  self.var.EPC3, Plab_in = globals.inZero, TDP_in = perc2to3_P,\
            Qr = globals.inZero, Qi = self.var.interflow[0:4] * (1- interflowDivider), Qp = self.var.grossGWrechargeFromSoil[0:4],\
            Kf = self.var.Kf, soilmass = self.var.soilM3, Vs = self.var.pre_w3)
        
        # update variables
        self.var.soil_P_labile3 =  outputs[0].copy()
        self.var.soil_P_dissolved3 =  outputs[1].copy()
        self.var.EPC3 =  outputs[2].copy()
        
        # update fluxes - no runoff from layer 2
        interflow3_P =  outputs[4].copy()
        toGW = outputs[5].copy()
       
      
        # sum outflows
        
        # runoff [kg]
        self.var.directRunoff_P = np.nansum(directRunoff_P * self.var.fracVegCover[0:4], axis = 0)
        self.var.interflow_P = np.nansum((interflow2_P + interflow3_P) * self.var.fracVegCover[0:4], axis = 0)
        self.var.baseflow_P = self.var.baseflow * self.var.cellArea * self.var.GW_P_Conc 
        
        self.var.returnflowNonIrr_P = np.where(self.var.returnflowNonIrr > 0, loadmap('P_sourcePoint'),0.)
        
        # to runoff [kg] 
        self.var.runoff_P = self.var.directRunoff_P +  self.var.interflow_P + self.var.baseflow_P + self.var.returnflowNonIrr_P

        # to groundwater [kg]
        self.var.toGroundwater_P =  np.nansum(toGW * self.var.fracVegCover[0:4], axis = 0)
        
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
        self.var.channel_P_Abstracted = np.maximum(np.minimum(self.var.act_channelAbst * self.var.cellArea *  self.var.channel_PConc, self.var.channel_P), 0.)
        
        # lake/reservoir
        resLake_P_compress = np.compress(self.var.compress_LR, self.var.resLake_P)
        resLake_P_Abstracted_small = np.minimum(divideValues(resLake_P_compress, self.var.lakeResStorageC + self.var.act_bigLakeAbstC) * (self.var.act_bigLakeAbstC), resLake_P_compress)               
        self.var.resLake_P_Abstracted = globals.inZero.copy()
        np.put(self.var.resLake_P_Abstracted, self.var.decompress_LR, resLake_P_Abstracted_small)
                            
        # groundwater
        self.var.groundwater_P_Abstracted = self.var.nonFossilGroundwaterAbs * self.var.cellArea * self.var.GW_P_Conc 
        
        # domestic P 
        self.var.domestic_P_Abstracted = self.var.channel_P_Abstracted * divideValues(self.var.Channel_Domestic, self.var.domesticDemand) +\
            self.var.resLake_P_Abstracted * divideValues(self.var.Lake_Domestic + self.var.Res_Domestic, self.var.domesticDemand) +\
            self.var.groundwater_P_Abstracted * divideValues(self.var.GW_Domestic, self.var.domesticDemand)
        
        # industry P 
        self.var.industry_P_Abstracted = self.var.channel_P_Abstracted * divideValues(self.var.Channel_Industry, self.var.industryDemand) +\
            self.var.resLake_P_Abstracted * divideValues(self.var.Lake_Industry + self.var.Res_Industry, self.var.industryDemand) +\
            self.var.groundwater_P_Abstracted * divideValues(self.var.GW_Industry, self.var.industryDemand)
            
        # livestock P 
        self.var.livestock_P_Abstracted = self.var.channel_P_Abstracted * divideValues(self.var.Channel_Livestock, self.var.livestockDemand) +\
            self.var.resLake_P_Abstracted * divideValues(self.var.Lake_Livestock + self.var.Res_Livestock, self.var.livestockDemand) +\
            self.var.groundwater_P_Abstracted * divideValues(self.var.GW_Livestock, self.var.livestockDemand)
        
        # irrigation P
        self.var.irrigation_P_Abstracted = self.var.channel_P_Abstracted * divideValues(self.var.Channel_Irrigation, self.var.totalIrrDemand) +\
            self.var.resLake_P_Abstracted * divideValues(self.var.Lake_Irrigation + self.var.Res_Irrigation, self.var.totalIrrDemand) +\
            self.var.groundwater_P_Abstracted * divideValues(self.var.GW_Irrigation, self.var.totalIrrDemand)
        
        # irrigation returnFlows
        self.var.returnflowIrr_P = divideValues(self.var.returnflowIrr, self.var.act_irrWithdrawal)
        
