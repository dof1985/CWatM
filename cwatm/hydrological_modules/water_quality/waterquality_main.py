# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Name:        Water quality - main module
# Purpose:
#
# Author:      TT, FS, PB, MS, DF
#
# Created:     10/10/2022
# Copyright:   (c) TT, FS, PB, MS, DF 2022
# -------------------------------------------------------------------------

from cwatm.management_modules.data_handling import *
from cwatm.management_modules.globals import *

# define water quality sub-modules
from cwatm.hydrological_modules.water_quality.phosphorus import waterquality_phosphorus

class water_quality(object):
    """
    WATER QUALITY - MAIN MODULE
        
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
        
        # define water quality sub-modules
        self.waterquality_p = waterquality_phosphorus(self)
	
    def initial(self):
        
        
        # Create Arrays for modules - adapted from landcoverType.py
        
        phosphorusVars = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_labile1', 'soil_P_labile2', \
                          'soil_P_dissolved1', 'soil_P_dissolved2', 'EPC1', 'EPC2']
        
        # only applied to landcovers with soil underneath (grassland and managed grasslands are treated has one landcover
        for variable in phosphorusVars: vars(self.var)[variable] = np.tile(globals.inZero,(4,1))
        #for variable in phosphorusVars: vars(self.var)["sum_" + variable] = globals.inZero.copy()
        
        # Soil layers in the water quality module ###
        '''
        CWatM has three soil layers, and the water quality modules uses only two layers. 
        Soil depth of the water quality layers are: 
        d_1 = D_1 * 0.2 (where as D is CWatM layer, and d is for the water quality)
        d_2 = D_1 * 0.8 + D_2 + D_3
        
        The mass of the soil layers is calculated as [kg / m2]
        m_1 = d_1 * rho_1  (where rho_1 is bulk density in kg/m3)
        m_2 = d_2 * rho_2' 
        
        rho_2' = (rho_1 * D_1 * 0.8 + D_2 * rho_2 + D_3 * rho_3) / d_2
        m_2 = (rho_1 * D_1 * 0.8 + D_2 * rho_2 + D_3 * rho_3) 

        '''
        
        self.var.gCm3TokgM3 = 1000
        self.var.M2mm = 1000

        # soil depth [m]
        self.var.wq_soilDepth1 = self.var.soildepth[0] * 0.2
        self.var.wq_soilDepth2 = self.var.soildepth[0] * 0.8 + self.var.soildepth[1] + self.var.soildepth[2]
        
        # soil mass [kg/m2]
        self.var.soilM1 = self.var.wq_soilDepth1 * loadmap('rho1') * self.var.gCm3TokgM3
        self.var.soilM2 = (self.var.soildepth[0] * loadmap('rho1')  * 0.8 + self.var.soildepth[1] * loadmap('rho2') + self.var.soildepth[2] * loadmap('rho3')) * self.var.gCm3TokgM3 
        
        # calculate soil moisture content for WQ soil layers
        self.var.wq_soilMoisture1 = self.var.w1 * 0.2
        self.var.wq_soilMoisture2 = self.var.w1 * 0.8 + self.var.w2 + self.var.w3
        
     
        # Load managed grassland fraction ###
        
        self.var.managedGrassland = globals.inZero.copy()
        if 'fracManagedGrassland' in binding:
            self.var.managedGrassland = loadmap('fracManagedGrassland')
        
        
        # soil calculation - e.g., water quality soil layers
        # soil depths for water quality
        # soil masses for water quality
        
        ## initiate all phosphrous stocks -> soil, channel, lakes/reservoirs, groundwater
        ## calculate all conversion factors
        
        # dynamic -> load inputs and ground cover adjustment
        # dynamic_soil, dynamic_...
        
        # check if water quality is included
        self.var.includeWaterQuality =  False
        if 'includeWaterQuality' in option:
            self.var.includeWaterQuality =  checkOption('includeWaterQuality')

        if self.var.includeWaterQuality:
            # load sub-modules
             
            self.var.includePhosphorus = False
            if 'includePhosphorus' in binding:
                self.var.includePhosphorus = returnBool('includePhosphorus')
            
            
            # run initial sub-modules
            
            if self.var.includePhosphorus:
                self.waterquality_p.initial()
            
            
            
            # sum fractions of flows/stocks from different landcovers
            
            # sum total soil P stocks [kg / m2] 
            for variable in phosphorusVars:
                vars(self.var)["sum_" + variable] = np.nansum(vars(self.var)[variable] * self.var.fracVegCover[0:4], axis = 0)
            
            soil_P_total1 = (self.var.sum_soil_P_inactive1 + self.var.sum_soil_P_labile1 + self.var.sum_soil_P_dissolved1) * self.var.cellArea
            soil_P_total2 = (self.var.sum_soil_P_inactive2 + self.var.sum_soil_P_labile2 + self.var.sum_soil_P_dissolved2) * self.var.cellArea
            
            print(soil_P_total1 + soil_P_total2)

            
            
    def dynamic(self):
        print("wq dyn.")
    
       