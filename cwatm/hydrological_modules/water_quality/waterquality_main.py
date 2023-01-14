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
        
            
    def dynamic(self):
        print("wq dyn.")
    
       