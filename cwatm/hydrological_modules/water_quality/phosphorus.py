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
        self.var.EPC1 = divideArrays(self.var.soil_P_labile1, self.var.Kf)
        self.var.EPC2 = divideArrays(self.var.soil_P_labile2, self.var.Kf)
        
        # calculate TDP [kg P /m2]
        self.var.soil_P_dissolved1 = self.var.EPC1 * self.var.wq_soilMoisture1 * 1000
        self.var.soil_P_dissolved2 = self.var.EPC2 * self.var.wq_soilMoisture2 * 1000
        
        
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
        '''
        # variable names for P balance items inputs - all inputs should be provided as total KG phosphrous input
        self.var.P_fertIn_varName = "PFert"
        self.var.P_manureIn_varName = "Pmanure"
        self.var.P_cropsOut_varName = "Pex"
        
        # load initial P stocks in Soil (concentration)
        self.var.soilPConc_inactiveIni = loadmap('soil_inactive_phosphrous')
        self.var.soilPConc_totalIni = loadmap('soil_total_phosphrous')
        
        # load initial P stocks in Groundwater
        if 'gw_dissolved_phosphrous' in binding:
            self.var.gwPConc_dissolvedIni = loadmap('gw_dissolved_phosphrous')
        else:
            self.var.gwPConc_dissolvedIni = global.inZero.copy()
        
        
        # MAYBE AS a pre process - or contiue later to calculate the p_kf
        ## read in organic matter content to calculate P partitioning coefficient
        # load soil organic matter by layer
        self.var.soilOrganicM0 = loadmap('soil_organic0')
        self.var.soilOrganicM1 = loadmap('soil_organic1')
        self.var.soilOrganicM2 = loadmap('soil_organic2')
        
        # P partitioning constant - check number - can also be dynamic with soil organic matter
        self.var.soilP_kf = 0.2 # this can be loaded from settings file
        #self.var.waterP_kf = 0.2 # this can be loaded from settings file
        self.var.P_manureShareInactive = 0.07 # ratio
        ####
        
        # Read initial
        self.var.fractionGrasslandCropGraze = loadmap('grasslandsShare_crop_and_graze')
      
        # calculate P inactive - p. 5 documenation - note forests, natural grasslands only have inactive - so the divission of active is split proprotionally to relative area of all other landcoves ( except of water and sealed areas )
        # calculate initial EPC0 from labile mass
        
     
        e.g layer1
        epc0 = self.var.soilPConc_labileIni  /  kf 
        
        
    
        
        ## array of arrays - each sub array is all p (before partitioning in any land cover)
        ## OR each array is Pin, Pinorganic,  TDP, inorganic, Labile)
        # for all p types
        pSoilStockVars = ['soilP_inactive', 'soilP_labile','soilP_dissolved', 'soilPConc_inactive', 'soilPConc_labile','soilPConc_dissolved']
        # not sure 6 is the right number - maybe water is out, so only 4 
        for variable in pSoilStockVars: vars(self.var)[variable] = np.tile(globals.inZero,(4,1))
        for variable in pSoilStockVars: vars(self.var)["sum_" + variable] = globals.inZero.copy()

        # Initial calculations
        # Soil
        for No in range(4):
            
            # ['forest', 'grassland', 'irrPaddy', 'irrNonPaddy','sealed']
            # Calcuclate inactive P in soil. Forest and natural grassland only has inactive P.
            self.var.soilPConc_inactive[No] = self.var.soilPConc_inactiveIni * self.var.fracVegCover[No]
            
            if No == 0:
                # assume natural/semi natural system only have inactive phosphrous
                self.var.soilPConc_labile[No] = globals.inZero.copy()
                self.var.soilPConc_dissolved[No] = globals.inZero.copy()
            else:
                if No == 1:
                    soilPConc_active = (self.var.soilPConc_total - self.var.soilPConc_inactiveIni)  * self.var.fracVegCover[No]
                else:
                    soilPConc_active = (self.var.soilPConc_total - self.var.soilPConc_inactiveIni)  * self.var.fracVegCover[No]

                self.var.soilPConc_labile = soilPConc_active / ( 1 + self.var.soilP_kf)
                self.var.soilPConc_dissolved = soilPConc_active * ((self.var.soilP_kf) / ( 1 + self.var.soilP_kf))
        '''   
    def dynamic(self):
        print("p dyn")
        '''
        if dateVar['newStart'] or dateVar['newYear']:
            timediv = globals.dateVar['daysInYear']
            
            ## load inputs  - P balance - [t/ha] -> [t/m2]  
            self.var.P_fertIn = readnetcdf2('PTerBalanceMaps', useDaily = "yearly", value=self.var.P_fertIn_varName)  timediv
            self.var.P_manureIn = readnetcdf2('PTerBalanceMaps', useDaily = "yearly", value=self.var.P_manureIn_varName)  / timediv
            self.var.P_cropsOut = readnetcdf2('PTerBalanceMaps', useDaily = "yearly", value=self.var.P_cropsOut_varName)  / timediv
            # currenlty for param (below) can we get a map?
            self.var.P_atmDeposit =  loadmap('PAtmDeposition')  / timediv
        
            ## load share of managed grassland (out of all grasslands) -> change to readnetcdf2
            self.var.fractionGrasslandCropGraze = loadmap('grasslandsShare_crop_and_graze')
        '''