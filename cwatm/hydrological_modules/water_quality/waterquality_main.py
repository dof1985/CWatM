# -------------------------------------------------------------------------
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
            

            # Create sub-modules variables
            
            '''
                Water flows adapated to the water quality soil layer strcuture: transpiration, capilary rise, percolation, GW recharge,
                and inteflow
                
                Transpiration - assumed to 'pump' dissolved P proportionally to its concentration - should be dropped, Net P inputs partially account for P removal by plants,
                    and the mechanism for nutrient uptake by plant is more complicated.
                capRise - soil capillary rise between layer 2 to layer 1 takes nutrients with it.
                capRiseFromGW - capillary rise from GW into soil layers can bring phosphorous from groundwater to soil layers (if accounted for).
                Percolation1to2 - simulate dissolved P moving with soil percolation.
                Percolation2toGW - simulate dissolved P 'losses' to groundwater. At this stage it doesn't affect groundwater P concetrations.
                Interflow2 - simulate lateral flow from soil layer 2 of the water quality (and 3 of cwatm) to the runoff.
            '''

            waterqualityVars = ['wq_Transpiration1', 'wq_Transpiration2', 'wq_capRise1', 'wq_capRiseFromGW1', 'wq_capRiseFromGW2', \
                                'wq_Percolation1to2', 'wq_Percolation2toGW', 'wq_Interflow2']

            for variable in waterqualityVars: vars(self.var)[variable] = np.tile(globals.inZero, (4, 1))

            if self.var.includePhosphorus:
                phosphorusVars = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_dissolved1', 'soil_P_dissolved2', 'EPC1', 'EPC2', 'runoff_P']
                phosphorusVarsSum = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_dissolved1', 'soil_P_dissolved2']
        
                # only applied to landcovers with soil underneath (grassland and managed grasslands are treated has one landcover
                for variable in phosphorusVars: vars(self.var)[variable] = np.tile(globals.inZero,(4,1))
            
            
            # Soil layers in the water quality module ###
            '''
            CWatM has three soil layers, and the water quality modules uses only two layers. 
            Soil depth of the water quality layers are: 
            d_1 = D_1 * self.var.wq_relSoilDepth1 (where as D is CWatM layer, and d is for the water quality)
            d_2 = D_1 * (1  -  self.var.wq_relSoilDepth1) + D_2 + D_3
            
            The mass of the soil layers is calculated as [kg / m2]
            m_1 = d_1 * rho_1  (where rho_1 is bulk density in kg/m3)
            m_2 = d_2 * rho_2' 
            
            rho_2' = (rho_1 * D_1 * (1 - self.var.wq_relSoilDepth1) + D_2 * rho_2 + D_3 * rho_3) / d_2
            m_2 = (rho_1 * D_1 * (1 - self.var.wq_relSoilDepth1) + D_2 * rho_2 + D_3 * rho_3) 
    
            '''

            
        
            self.var.wq_SoilDepth1 = 0.01 # [soil depth 1 for WQ in meters; default to 0.01; 10 mm)
            self.var.wq_relSoilDepth1 = divideValues(self.var.wq_SoilDepth1, self.var.soildepth[0])  # need to be revised
        
            # soil depth [m]
            self.var.wq_soilDepth1 = self.var.soildepth[0] * self.var.wq_relSoilDepth1
            self.var.wq_soilDepth2 = self.var.soildepth[0] * (1 - self.var.wq_relSoilDepth1) + self.var.soildepth[1] + self.var.soildepth[2]

            # soil mass [kg/m2]
            
            '''
            CWatM-WQ global is based on bulk density data from openLandMap (Hengl, 2018; https://zenodo.org/record/2525665#.Y-NjZi_TVD_)
            The input includes three soil depths (0 cm, 30 cm, and 100 cm) corresponding to the CWatM three soil layers.
            
            Raw data available, at 250 meters spatial resolution, was aggregated to 10,000 meters and downloaded from 
            GoogleEarth Engine (GEE). It was later resampled and rectified to the global CWatM extent and resolution (5'), 
            and exported as a layered NetCDF.
            '''
            # loaded as bulkDensity g/cm3
            self.var.gCm3TokgM3 = 1000
            
            rho1 = globals.inZero.copy() + 1.
            rho2 = globals.inZero.copy() + 1.
            rho3 = globals.inZero.copy() + 1.

            if 'bulkDensity' in binding:
                rho1 = readnetcdf2('bulkDensity', value= 'rho1')
                rho2 = readnetcdf2('bulkDensity', value= 'rho2')
                rho3 = readnetcdf2('bulkDensity', value= 'rho3')
               
            self.var.soilM1 = self.var.wq_soilDepth1 * rho1 * self.var.gCm3TokgM3
            self.var.soilM2 = (self.var.soildepth[0] * rho1  * (1 - self.var.wq_relSoilDepth1) + self.var.soildepth[1] * rho2 + self.var.soildepth[2] * rho3) * self.var.gCm3TokgM3


            # calculate soil moisture content for WQ soil layers
            self.var.wq_soilMoisture1 = self.var.w1 * self.var.wq_relSoilDepth1
            self.var.wq_soilMoisture2 = self.var.w1 * (1 - self.var.wq_relSoilDepth1) + self.var.w2 + self.var.w3
        
     
            # Load managed grassland fraction ###
        
            self.var.managedGrassland = globals.inZero.copy()
            if 'fracManagedGrassland' in binding:
                self.var.managedGrassland = loadmap('fracManagedGrassland')


            # Run initial sub-modules
            
            if self.var.includePhosphorus:
                self.waterquality_p.initial()
            
            if self.var.includeErosed:
                self.erosed.initial()

            if self.var.includePhosphorus:
                # sum total soil P stocks [kg / m2]
                for variable in phosphorusVarsSum:
                    vars(self.var)["sum_" + variable] = np.nansum(vars(self.var)[variable] * self.var.fracVegCover[0:4], axis = 0) * self.var.cellArea

    def dynamic(self): 
        if dateVar['newStart'] or dateVar['newYear']:
            current_year = globals.dateVar['currDate']
            '''
            later change to:
            self.var.fracManagedGrassland = readnetcdf2('fractionManagedGrassland', current_year,\
            useDaily='yearly', value='fracManagedGrassland')
            '''
            self.var.fracManagedGrassland = globals.inZero.copy()
            if 'fractionManagedGrassland' in binding:
                self.var.fracManagedGrassland = loadmap('fractionManagedGrassland')
         
        # run initial sub-modules
        # Run the Erosion and Sediment Yield (EroSed) module
        if self.var.includeErosed:
            self.erosed.dynamic()
        
        # Run the  module        
        if self.var.includePhosphorus:
            phosphorusVarsSum = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_labile1', 'soil_P_labile2', \
                          'soil_P_dissolved1', 'soil_P_dissolved2']
                          
            self.waterquality_p.dynamic()
            
            # sum total soil P stocks [kg / m2] 
            for variable in phosphorusVarsSum:
                vars(self.var)["sum_" + variable] = np.nansum(vars(self.var)[variable] * self.var.fracVegCover[0:4], axis = 0) * self.var.cellArea

    
       

