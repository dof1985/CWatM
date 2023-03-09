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

        # Division between landcovers is assumed proportional to their relative area fraction
        #soil_depthRatio1 = divideValues(self.var.wq_soilDepth1, self.var.wq_soilDepth1 + self.var.wq_soilDepth2)
    
        # labile soil P [kg P]
        # Applied to land cover [1, 2, 3] - natural landcover (forests) initial labile is assumed to be 0. 
        self.var.soil_P_labile1[1:4] += PSoil_init1 - self.var.soil_P_inactive1[1:4]
        self.var.soil_P_labile2[1:4] += PSoil_init2 - self.var.soil_P_inactive2[1:4]
        self.var.soil_P_labile3[1:4] += PSoil_init3 - self.var.soil_P_inactive3[1:4]

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
        self.var.soil_P_dissolved1 = np.minimum(self.var.EPC1 * self.var.w1, self.var.soil_P_labile1)
        self.var.soil_P_dissolved2 = np.minimum(self.var.EPC2 * self.var.w2, self.var.soil_P_labile2)
        self.var.soil_P_dissolved3 = np.minimum(self.var.EPC3 * self.var.w3, self.var.soil_P_labile3)

        # update labile P to keep P balance
        self.var.soil_P_labile1 -= self.var.soil_P_dissolved1
        self.var.soil_P_labile2 -= self.var.soil_P_dissolved2
        self.var.soil_P_labile3 -= self.var.soil_P_dissolved3

        # load groundwater P concentration [kg / M]
        self.var.GW_P_Conc = globals.inZero.copy()
        if 'GW_P_Conc' in binding:
            self.var.GW_P_Conc += loadmap('GW_P_Conc')
            
        # load runoff P adjustment coefficient
        self.var.runoff_Padj = globals.inZero.copy() + 1.
        if 'runoff_Padj' in binding:
            self.var.runoff_Padj = loadmap('runoff_Padj')
        
        
        # channel phsophorus [kg]
        self.var.channel_P = globals.inZero.copy()
        self.var.channel_PConc = globals.inZero.copy()
        
        self.var.ChanDiffuse_P = globals.inZero.copy()
        if 'ChannelDiff_P' in binding:
            self.var.ChanDiffuse_P = loadmap('ChannelDiff_P')
        #### Is there anyway to check for initial balance - i.e. so all soil_P in kg at time step = 0 == self.var.soil_PConc_total

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
        '''
        # dynamic p inactive ###
        # Add only to managed grasslands and irrigated crops
        # rewrite - make sure that weathering is per soil layer, and cap it to the available p inactive
        delta_soil_P_inactive =  -self.var.soil_P_weathering
        
        delta_fromManure1 = self.var.soil_P_inactive1.copy()
        delta_fromManure2 = self.var.soil_P_inactive1.copy()
        
        # managed grasslands - add to grassland * frac_managed_grassland - currently it is assumed that soil layer 3 has no inactive P - I guess this is not correct
        self.var.soil_P_inactive1[1] += delta_soil_P_inactive * soil_depthRatio1 * self.var.fracManagedGrassland
        self.var.soil_P_inactive2[1] += delta_soil_P_inactive * (1 - soil_depthRatio1) * self.var.fracManagedGrassland
        #self.var.soil_P_inactive3[1] += delta_soil_P_inactive * (1 - soil_depthRatio1) * self.var.fracManagedGrassland
        
        # irrigated crops (non Paddy)
        self.var.soil_P_inactive1[3] += delta_soil_P_inactive * soil_depthRatio1
        self.var.soil_P_inactive2[3] += delta_soil_P_inactive * (1 - soil_depthRatio1)
        
        # calculate added inactive soil P from manure
        delta_fromManure1 -= self.var.soil_P_inactive1
        delta_fromManure2 -= self.var.soil_P_inactive2
        '''
        
        ''' 
        Soil dissolved and labile dynamic section 
        
      0 1. Infiltration & Leak to surface runoff  (soil.py line 470 -473, 490 -492)
      0 1.** Irrigation ** soil.py (line 256) ; also see Paddy irrigation in (soil.pyj lines 262 -281), inclue also openwaterEvaporation
      1 2. capillary rise from groundwater soil.py (line 292 - 303, Modflow) - CANCEL
      1 2.** capillary rise from soil layers soil.py (line 337 - 342, Modflow; lines 539 -553 non modflow) - CANCEL
      1 3. Evapotranspiration ?! - assuming plant uptake takes a solution in equilibria (soil.py lines 408 -410) - CANCEL
      - 4. Baresoil evaporation - only water (soil.py  lines 415 -422) 
      1 5. Soil percolation  - (soil.py lines 662 -669); include self.var.perc3toGW
        6. Interflow -- includes percolation + preferential flows (the latter has P = 0) soil.py lines 744 -750
        7. Groundwater recharge ?!
        
        
        ## crops extnetion soil.py -> lines 684 -737
        
        '''
        
        '''
            Soil layers 1 -2 exchange water with their environment - resulting in a different 
            soil water moisture in time  t+1, relative to time 1. Some processes exchange water-P solution, and
            other only exchange water.
            We calculate mass P changes (TDP, self.var.soil_P_dissolve_j) onlt for the first type. The other type is not \
            accounted for in this module, but is captured in the change of soil moisture.
            
            # water P solution - leak to Surface runoff (only from soil layer 1), Irrigation (Infiltration from..),
                Evapotranspiration, Interflow, Groundwater recharge, Soil percolation.
            # Only water - Infiltration, Bare soil evaporation
        '''
        
        
        # calculate water movement based on prior timestep concentrations - resulting with mass of layer P per unit area
        '''
        POut_j = (sum of water outflows* from j / soilMoisture_j) * Mass P in join
        * Outflows excluding evaporation
        
        '''
        # calculate relative moisture of soil layer 1 out of 12 - to split interflow of P
        relMoisture2 = divideArrays(self.var.pre_w2, self.var.soildepth[1]) 
        relMoisture3 = divideArrays(self.var.pre_w3, self.var.soildepth[2]) 
        interflowDivider = divideArrays(relMoisture2, relMoisture2 + relMoisture3)
  
        interflow2 = self.var.interflow[0:4] * interflowDivider * divideArrays(self.var.soil_P_dissolved2, self.var.pre_w2)
        interflow3 = self.var.interflow[0:4] * (1 - interflowDivider) * divideArrays(self.var.soil_P_dissolved3, self.var.pre_w3)
        toGW = self.var.perc3toGW[0:4] * divideArrays(self.var.soil_P_dissolved3, self.var.pre_w3)
        POut1 = self.var.perc1to2 * divideArrays(self.var.soil_P_dissolved1, self.var.pre_w1) # REDESIGN ALL TO BE WRITTEN AS IN THIS LINE
        POut2 = self.var.perc2to3 * divideArrays(self.var.soil_P_dissolved2, self.var.pre_w2) + interflow2
        POut3 = toGW + interflow3
        
        # INPUT TDPs 
        PIn1 = 0
        PIn2 = self.var.perc1to2  * divideArrays(self.var.soil_P_dissolved1, self.var.pre_w1)
        PIn3 = self.var.perc2to3  * divideArrays(self.var.soil_P_dissolved2, self.var.pre_w2)
           
        # Calculate P Net Input
           
        P_netInput1 =  PIn1 - POut1# + self.var.soil_P_weathering  # do we only have weathering in the topsoil?
        P_netInput2 =  PIn2 - POut2
        P_netInput3 =  PIn3 - POut3
        
        # Calculate Plab Net Input - Currently only apply on managed grasslands and on irrigated agriculture

        # managed grasslands 
        self.var.soil_P_labile1[1] += PSoil_Input1 * self.var.fracManagedGrassland
        self.var.soil_P_labile2[1] += PSoil_Input2 * self.var.fracManagedGrassland
        self.var.soil_P_labile3[1] += 0 * self.var.fracManagedGrassland
        
        # irrigated crops (non Paddy)
        self.var.soil_P_labile1[3] += PSoil_Input1 
        self.var.soil_P_labile2[3] += PSoil_Input2 
        self.var.soil_P_labile3[3] += 0
        
        dissolvedtoRunoff1 = (self.var.directRunoff[0:4] > 0) * self.var.runoff_Padj * \
            self.var.wq_SoilDepthRunoff * np.maximum(self.var.soil_P_dissolved1 + P_netInput1, 0.)
        # limit runoff P concentration to soil P concentration
        soil_P_dissolved1Conc = divideArrays(self.var.soil_P_dissolved1, self.var.w1)
        runoff_P_conc = divideArrays(dissolvedtoRunoff1, self.var.directRunoff[0:4])
        dissolvedtoRunoff1 = np.minimum(runoff_P_conc, soil_P_dissolved1Conc) * self.var.directRunoff[0:4]

        # to runoff [kg]
        self.var.runoff_P = np.nansum((dissolvedtoRunoff1  + interflow2 + interflow3) * self.var.fracVegCover[0:4], axis = 0) 

        # to groundwater [kg]
        self.var.toGroundwater_P =  np.nansum(toGW * self.var.fracVegCover[0:4], axis = 0)
        
        # Update dissolved P mass in soils  - runoff TDP is allowed only when directRunoff > 0
        # to be improved - self.var.runoff_Padj should be also accompanied by max runoff concentration, e.g., how much phosphrous can be removed by the surfacerunoff as a function of runoff volume
        # now we get negative TDP in layer 1
        
        self.var.soil_P_dissolved1 = np.maximum(self.var.soil_P_dissolved1 + P_netInput1 - dissolvedtoRunoff1, 0.)
        self.var.soil_P_dissolved2 = np.maximum(self.var.soil_P_dissolved2 + P_netInput2, 0.)
        self.var.soil_P_dissolved3 = np.maximum(self.var.soil_P_dissolved3 + P_netInput3, 0.)
        
        # Advective P fluxes END ### 
        
        
        # Calculate New TDP and Plabile
        # Discretizied soil water TDP - follow simplyP model.py lines 42 - 47 ; https://github.com/LeahJB/SimplyP/
        b1 = divideArrays(self.var.Kf * self.var.soilM1 + self.var.w1 - self.var.pre_w1, self.var.w1)
        a_b1 = divideArrays(self.var.soil_P_labile1, b1) # a/b
        self.var.soil_P_dissolved1 = a_b1 + (self.var.soil_P_dissolved1 - a_b1) * np.exp(-1 * b1)
        
        b2 = divideArrays(self.var.Kf * self.var.soilM2 + self.var.w2 - self.var.pre_w2, self.var.w2)
        a_b2 = divideArrays(self.var.soil_P_labile2, b2) # a/b
        self.var.soil_P_dissolved2 = a_b2 + (self.var.soil_P_dissolved2 - a_b2) * np.exp(-1 * b2)
        
        b3 = divideArrays(self.var.Kf * self.var.soilM3 + self.var.w3 - self.var.pre_w3, self.var.w3)
        a_b3 = divideArrays(self.var.soil_P_labile3, b3) # a/b
        self.var.soil_P_dissolved3 = a_b3 + (self.var.soil_P_dissolved3 - a_b3) * np.exp(-1 * b3)
        
    
        # Discretizied Plabile
        b0_1 = b1 * self.var.w1
        a_b0_1 = divideArrays(self.var.soil_P_labile1, b0_1) # a/b
        dPLabile1 = self.var.Kf * self.var.soilM1 * (a_b1 - self.var.EPC1 + \
            (divideArrays(globals.inZero.copy() + 1., b1) * \
            (divideArrays(self.var.soil_P_dissolved1, self.var.w1) - a_b0_1) * \
            (1 - np.exp(-1 * b1))))
            
        dPLabile1 = np.where(self.var.w1 == 0, 0., dPLabile1)

            
        b0_2 = b2 * self.var.w2
        a_b0_2 = divideArrays(self.var.soil_P_labile2, b0_2) # a/b
        dPLabile2 = self.var.Kf * self.var.soilM2 * (a_b2 - self.var.EPC2 + \
            (divideArrays(globals.inZero.copy() + 1., b2) * \
            (divideArrays(self.var.soil_P_dissolved2, self.var.w2) - a_b0_2) * \
            (1 - np.exp(-1 * b2))))
            
        dPLabile2 = np.where(self.var.w2 == 0, 0., dPLabile2)
        
            
        b0_3 = b3 * self.var.w3
        a_b0_3 = divideArrays(self.var.soil_P_labile3, b0_3) # a/b
        dPLabile3 = self.var.Kf * self.var.soilM3 * (a_b3 - self.var.EPC3 + \
            (divideArrays(globals.inZero.copy() + 1., b3) * \
            (divideArrays(self.var.soil_P_dissolved3, self.var.w3) - a_b0_3) * \
            (1 - np.exp(-1 * b3))))
            
        dPLabile3 = np.where(self.var.w3 == 0, 0., dPLabile3)

        # limit TDP to be >= 0
        self.var.soil_P_dissolved1 = np.maximum(self.var.soil_P_dissolved1, 0.)
        self.var.soil_P_dissolved2 = np.maximum(self.var.soil_P_dissolved2, 0.)
        self.var.soil_P_dissolved3 = np.maximum(self.var.soil_P_dissolved3, 0.)

        # Update labile P mass in soils
        self.var.soil_P_labile1 += dPLabile1
        self.var.soil_P_labile2 += dPLabile2
        self.var.soil_P_labile3 += dPLabile3
        
        
        # update  soil equilibrium TDP concentration of zero sorption  
        self.var.EPC1 = divideArrays(self.var.soil_P_labile1, self.var.Kf * self.var.soilM1)
        self.var.EPC2 = divideArrays(self.var.soil_P_labile2, self.var.Kf * self.var.soilM2)
        self.var.EPC3 = divideArrays(self.var.soil_P_labile3, self.var.Kf * self.var.soilM3)

        
      
    ## END SOIL P DYNAMIC
        # calculate soil moisture content for WQ soil layers
        
        # exported phsophorus by outlet [kg]
        self.var.outlet_P = globals.inZero.copy()

        # calculate inputs to dissloved P
        # balance labile/active using EPC
        
        
        
       