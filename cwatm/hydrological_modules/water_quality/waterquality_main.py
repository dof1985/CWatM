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

# import water quality sub-modules
from cwatm.hydrological_modules.water_quality.phosphorus import waterquality_phosphorus
from cwatm.hydrological_modules.water_quality.erosed import waterquality_erosed

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
        self.erosed = waterquality_erosed(model)
	
    def initial(self):
        """
        # soil calculation - e.g., water quality soil layers
        # soil depths for water quality
        # soil masses for water quality
        
        ## initiate all phosphrous stocks -> soil, channel, lakes/reservoirs, groundwater
        ## calculate all conversion factors
        
        # dynamic -> load inputs and ground cover adjustment
        # dynamic_soil, dynamic_...
        """

        # check settings file, if water quality should be included
        self.var.includeWaterQuality =  False
        if 'includeWaterQuality' in option:
            self.var.includeWaterQuality =  checkOption('includeWaterQuality')

        if self.var.includeWaterQuality:
            # load sub-modules
            self.var.includePhosphorus = False
            if 'includePhosphorus' in binding:
                self.var.includePhosphorus = returnBool('includePhosphorus')

            self.var.includeErosed = False
            if 'includeErosed' in binding:
                self.var.includeErosed = returnBool('includeErosed')
            
            
            # run initial sub-modules
            
            if self.var.includePhosphorus:
                self.waterquality_p.initial()

            if self.var.includeErosed:
                self.erosed.initial()
        
            
    def dynamic(self):

        if self.var.includeWaterQuality:
            if self.var.includeErosed:
                self.erosed.dynamic()



