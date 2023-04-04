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
            
            # create water quality variables
            waterQualityVars = ['pre_w1', 'pre_w2', 'pre_w3', 'norm_w1', 'norm_w2', 'norm_w3']
            for variable in waterQualityVars: vars(self.var)[variable] = np.tile(globals.inZero,(4,1))            
            
            # to delete
            self.var.norm_w1 = divideArrays(self.var.w1, self.var.cellArea) 
            self.var.norm_w2 = divideArrays(self.var.w2, self.var.cellArea) 
            self.var.norm_w3 = divideArrays(self.var.w3, self.var.cellArea) 
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

            if self.var.includePhosphorus:
                phosphorusVars = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3',\
                                  'soil_P_labile1', 'soil_P_labile2', 'soil_P_labile3',\
                                  'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3',\
                                  'EPC1', 'EPC2', 'EPC3', 'runoff_P', 'toGroundwater_P']
                phosphorusVarsSum = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_labile3', 'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3', 'EPC1', 'EPC2', 'EPC3']
        
                # only applied to landcovers with soil underneath (grassland and managed grasslands are treated has one landcover
                for variable in phosphorusVars: vars(self.var)[variable] = np.tile(globals.inZero,(4,1))

            
            
            # Define depth of top soil from which nutrients leaks to runoff: 10 mm###


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


            
            wq_SoilDepth1 = 0.01 
            self.var.wq_SoilDepthRunoff = divideValues(wq_SoilDepth1, self.var.soildepth[0]) 
        
            # soil mass [kg/m2]
            
            '''
            CWatM-WQ global uses bulk density data from 'The Global Soil Dataset for Earth System Modeling'
            (Shangguan et al., 2014; https://doi.org/10.1002/2013MS000293; http://globalchange.bnu.edu.cn/research/soilw)
           
            The input includes three soil depths (0 -4.5 cm, 16.6 -28.9 cm, and 49.3 -82.9 cm) to represent CWatM's  three soil layers.
            '''
            # loaded as bulkDensity g/cm3
            self.var.gCm3TokgM3 = 1000
            
            rho1 = globals.inZero.copy() + 1.
            rho2 = globals.inZero.copy() + 1.
            rho3 = globals.inZero.copy() + 1.

            if 'rho1' in binding:
                rho1 = loadmap('rho1')
            if 'rho2' in binding:
                rho2 = loadmap('rho2')
            if 'rho3' in binding:
                rho3 = loadmap('rho3')
    
            # Calculate soil mass as [kg]
            self.var.soilM1 = self.var.soildepth[0] * rho1 * self.var.gCm3TokgM3 * self.var.cellArea
            self.var.soilM2 = self.var.soildepth[1] * rho2 * self.var.gCm3TokgM3 * self.var.cellArea
            self.var.soilM3 = self.var.soildepth[2] * rho3 * self.var.gCm3TokgM3 * self.var.cellArea

     
            # Load managed grassland fraction ###
        
            self.var.managedGrassland = globals.inZero.copy()
            if 'fracManagedGrassland' in binding:
                self.var.managedGrassland = loadmap('fracManagedGrassland')


            # Run initial sub-modules
            if self.var.includePhosphorus:
                self.waterquality_p.initial()

            if self.var.includeErosed:

                #erosedVars = ['runoffEnergyMusle']
                #
                #for variable in erosedVars: vars(self.var)[variable] = np.title(globals.inZero, (4, 1))

                self.erosed.initial()

            if self.var.includePhosphorus:
                # sum total soil P stocks [kg / m2]
                for variable in phosphorusVarsSum:
                    vars(self.var)["sum_" + variable] = np.nansum(vars(self.var)[variable] * self.var.fracVegCover[0:4], axis = 0)

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
         
        # Erosion and Sediment Yield (EroSed) dynamic part
        if self.var.includeErosed:

            # erosedVarsSum = ['runoff_energy_term']
            #
            # for variable in erosedVarsSum:
            #     vars(self.var)[variable] = np.title(globals.inZero, (4, 1))

            self.erosed.dynamic()

            # for variable in erosedVarsSum:
            #     vars(self.var)["sum_" + variable] = np.nansum(vars(self.var)[variable] * self.var.fracVegCover[0:4], axis = 0) * self.var.cellArea
        
        # Run the  module        
        if self.var.includePhosphorus:
 
                          
            self.waterquality_p.dynamic()
            
            # sum total soil P stocks [kg / m2] 
            phosphorusVarsSum = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_labile3', 'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3', 'EPC1', 'EPC2', 'EPC3']
                                
            for variable in phosphorusVarsSum:
                vars(self.var)["sum_" + variable] = np.nansum(vars(self.var)[variable] * self.var.fracVegCover[0:4], axis = 0)
            
            # TDP soil concentration [mg / liter]
            self.var.sum_soil_P_dissolvedConc1 = divideValues(self.var.sum_soil_P_dissolved1, self.var.sum_w1) * 10**3
            self.var.sum_soil_P_dissolvedConc2 = divideValues(self.var.sum_soil_P_dissolved2, self.var.sum_w2) * 10**3
            self.var.sum_soil_P_dissolvedConc3 = divideValues(self.var.sum_soil_P_dissolved3, self.var.sum_w3) * 10**3
            
            # PP soil concentration [mg / kg soil] - miligram
            self.var.sum_soil_P_labileConc1 = divideValues(self.var.sum_soil_P_labile1, self.var.soilM1) * 10**6
            self.var.sum_soil_P_labileConc2 = divideValues(self.var.sum_soil_P_labile2, self.var.soilM2) * 10**6
            self.var.sum_soil_P_labileConc3 = divideValues(self.var.sum_soil_P_labile3, self.var.soilM3) * 10**6

            

