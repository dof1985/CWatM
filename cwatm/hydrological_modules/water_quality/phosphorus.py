# -------------------------------------------------------------------------
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
        # load initial total soil p concentration
        self.var.soil_PConc_total = loadmap('total_Soil_PConc')
        
        # load initial inactive soil p concentration
        self.var.soil_PConc_inactive = globals.inZero.copy()
        if 'inactive_Soil_PConc' in binding:
           self.var.soil_PConc_inactive += loadmap('inactive_Soil_PConc') #[mg P/kg soil]
        
        # load soil absorption coefficient Kf 
        self.var.Kf = globals.inZero.copy() + loadmap('Kf')
        
        
        # inactive soil P [kg P / m2]
        self.var.soil_P_inactive1 += self.var.soil_PConc_inactive * 1e-06 * self.var.soilM1
        self.var.soil_P_inactive2 += self.var.soil_PConc_inactive * 1e-06 * self.var.soilM2
        
        # labile soil P [kg P / m2]
        # natural landcover (forests) initial labile is assumed to be 0.
        self.var.soil_P_labile1[1:4] += self.var.soil_PConc_total * 1e-06 * self.var.soilM1 - self.var.soil_P_inactive1[1:4]
        self.var.soil_P_labile2[1:4] += self.var.soil_PConc_total * 1e-06 * self.var.soilM2 - self.var.soil_P_inactive2[1:4]
        # calculate dissolved soil P [kg P /m2]
        '''
            For each landcover i and WQ soil layer j, the soil equilibrium TDP concentration of zero sorption is
            EPC_ij = soil_P_labile_ij / (Kf * soilM_j)
            
            Whereas Kf is Soil P adsorption coefficient, default – 1.1*10-4

            later the TDP is
            self.var.soil_P_dissovled_ij = EPC_ij * wq_soilMoisture_ij * 1000
            
            Where wq_soilMoisture_ij is the mositure content of land cover i and soil layer j  in meters
        '''
        
        # calculate  soil equilibrium TDP concentration of zero sorption 
        self.var.EPC1 = divideArrays(self.var.soil_P_labile1, self.var.Kf * self.var.soilM1)
        self.var.EPC2 = divideArrays(self.var.soil_P_labile2, self.var.Kf * self.var.soilM2)
        
        # calculate TDP [kg P /m2]
        self.var.soil_P_dissolved1 = self.var.EPC1 * self.var.wq_soilMoisture1 * 1000
        self.var.soil_P_dissolved2 = self.var.EPC2 * self.var.wq_soilMoisture2 * 1000
        
        
        # load groundwater P concentration [kg / M]
        self.var.GW_P_Conc = globals.inZero.copy()
        if 'GW_P_Conc' in binding:
            self.var.GW_P_Conc += loadmap('GW_P_Conc')
            
        # load runoff P adjustment coefficient
        self.var.runoff_Padj = globals.inZero.copy() + 1.
        if 'runoff_Padj' in binding:
            self.var.runoff_Padj = loadmap('runoff_Padj')
        
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
        
        # load daily net inputs 
        # should be at [kg P / m2]
        self.var.soil_PConc_total = loadmap('P_netInput')
        '''
        later change to:
        self.var.soil_PConc_total = readnetcdf2('P_netInput', day_of_year, useDaily='DOY', value='Net_P_Input')
        '''
        manureFrac = loadmap('f_Manure')
        manureInactiveFrac = loadmap('f_InactiveManure')
        
        ## Place holder for P weathering
        self.var.soil_P_weathering = globals.inZero.copy()
        
        # soil dynamic part ###
        
        # dynamic p inactive ###
        # Add only to managed grasslands and irrigated crops
        soil_depthRatio1 = divideValues(self.var.wq_soilDepth1, self.var.wq_soilDepth1 + self.var.wq_soilDepth2)
        delta_soil_P_inactive = self.var.soil_PConc_total * manureFrac * manureInactiveFrac - self.var.soil_P_weathering
        
        delta_fromManure1 = self.var.soil_P_inactive1.copy()
        delta_fromManure2 = self.var.soil_P_inactive1.copy()
        
        # managed grasslands - add to grassland * frac_managed_grassland
        self.var.soil_P_inactive1[1] += delta_soil_P_inactive * soil_depthRatio1 * self.var.fracManagedGrassland
        self.var.soil_P_inactive2[1] += delta_soil_P_inactive * (1 - soil_depthRatio1) * self.var.fracManagedGrassland
        
        # irrigated crops (non Paddy)
        self.var.soil_P_inactive1[3] += delta_soil_P_inactive * soil_depthRatio1
        self.var.soil_P_inactive2[3] += delta_soil_P_inactive * (1 - soil_depthRatio1)
        
        # calculate added inactive soil P from manure
        delta_fromManure1 -= self.var.soil_P_inactive1
        delta_fromManure2 -= self.var.soil_P_inactive2
        
        
        ''' 
        Soil dissolved and labile dynamic section 
        
      0 1. Infiltration & Leak to surface runoff  (soil.py line 470 -473, 490 -492)
      0 1.** Irrigation ** soil.py (line 256) ; also see Paddy irrigation in (soil.pyj lines 262 -281), inclue also openwaterEvaporation
      1 2. capillary rise from groundwater soil.py (line 292 - 303, Modflow)
      1 2.** capillary rise from soil layers soil.py (line 337 - 342, Modflow; lines 539 -553 non modflow)
      1 3. Evapotranspiration ?! - assuming plant uptake takes a solution in equilibria (soil.py lines 408 -410) 
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
        POut1 = divideArrays(self.var.wq_Transpiration1 + self.var.wq_Percolation1to2, self.var.wq_soilMoisture1) * self.var.soil_P_dissolved1
        POut2 = divideArrays(self.var.wq_Transpiration2 + self.var.wq_capRise1 + self.var.wq_Percolation2toGW + \
                self.var.wq_Interflow2, self.var.wq_soilMoisture2) * self.var.soil_P_dissolved2
        
        # PIn1 = In_From_2 +  In_From_GW
        PIn1 = divideArrays(self.var.wq_capRise1, self.var.wq_soilMoisture2)  * self.var.soil_P_dissolved2 + \
                self.var.wq_capRiseFromGW1 * self.var.GW_P_Conc
        # PIn2 = In_From_1 + In_From_GW
        PIn2 = divideArrays(self.var.wq_Percolation1to2, self.var.wq_soilMoisture1)  * self.var.soil_P_dissolved1 + \
                self.var.wq_capRiseFromGW2 * self.var.GW_P_Conc
        
        # Calculate P Net Input
     
        P_netInput1 =  PIn1 + self.var.soil_PConc_total * soil_depthRatio1 - delta_fromManure1 + \
                self.var.soil_P_weathering - POut1
        P_netInput2 =  PIn2 + self.var.soil_PConc_total * (1 - soil_depthRatio1) - delta_fromManure2 - POut2
        
        # Update dissolved P mass in soils 
        dissolvedtoRunoff1 = self.var.runoff_Padj * np.maximum(self.var.soil_P_dissolved1 + P_netInput1, 0.)
        self.var.soil_P_dissolved1 = (1 - self.var.runoff_Padj) * np.maximum(self.var.soil_P_dissolved1 + P_netInput1, 0.)
        self.var.soil_P_dissolved2 = np.maximum(self.var.soil_P_dissolved2 + P_netInput2, 0.)
        
        # Update soil Moisture
        self.var.wq_soilMoisture1 = self.var.w1 * self.var.wq_relSoilDepth1
        self.var.wq_soilMoisture2 = self.var.w1 * (1 - self.var.wq_relSoilDepth1) + self.var.w2 + self.var.w3
        
        '''
        # update  soil equilibrium TDP concentration of zero sorption  - is this step required?   
        self.var.EPC1 = divideArrays(self.var.soil_P_labile1, self.var.Kf * self.var.soilM1)
        self.var.EPC2 = divideArrays(self.var.soil_P_labile2, self.var.Kf * self.var.soilM2)
        '''
        
        # to runoff [kg]
        self.var.runoff_P = np.nansum((dissolvedtoRunoff1 + self.var.wq_Interflow2 * self.var.soil_P_dissolved2) * self.var.fracVegCover[0:4], axis = 0) * \
                self.var.cellArea
        
        # to groundwater
        self.var.toGroundwater_P =  np.nansum(self.var.wq_Percolation2toGW * self.var.soil_P_dissolved2 * self.var.fracVegCover[0:4], axis = 0) * \
                self.var.cellArea
        
        # Calculate delta_P_labile
        delta_P_labile1 = self.var.Kf * self.var.soilM1 * \
                (divideArrays(self.var.soil_P_dissolved1, self.var.wq_soilMoisture1) - self.var.EPC1)     
        delta_P_labile2 = self.var.Kf * self.var.soilM2 * \
                (divideArrays(self.var.soil_P_dissolved2, self.var.wq_soilMoisture2) - self.var.EPC2)
        
        # do not allow negative labile
        delta_P_labile1 = np.where(delta_P_labile1 < 0, np.maximum(self.var.soil_P_labile1, delta_P_labile1), delta_P_labile1)
        delta_P_labile2 = np.where(delta_P_labile2 < 0, np.maximum(self.var.soil_P_labile2, delta_P_labile2), delta_P_labile2)
        
        # do not allow negative TDP
        delta_P_labile1 = np.where(delta_P_labile1 > 0, np.minimum(self.var.soil_P_dissolved1, delta_P_labile1), delta_P_labile1)
        delta_P_labile2 = np.where(delta_P_labile2 > 0, np.minimum(self.var.soil_P_dissolved2, delta_P_labile2), delta_P_labile2)
        
        # Update dissolved P mass in soils - remove change in labile
        self.var.soil_P_dissolved1 -= delta_P_labile1
        self.var.soil_P_dissolved2 -= delta_P_labile2
        
        # Update labile P mass in soils
        self.var.soil_P_labile1 += delta_P_labile1
        self.var.soil_P_labile2 += delta_P_labile2
        
    
        # calculate soil moisture content for WQ soil layers
        
        print(self.var.runoff_P)
        print(self.var.toGroundwater_P)
        # calculate inputs to dissloved P
        # balance labile/active using EPC
        
        
        
        
        
       