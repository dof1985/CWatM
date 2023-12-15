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
    '''
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
    '''

    def __init__(self, model):
        self.var = model.var
        self.model = model
        
        # define water quality sub-modules
        self.waterquality_p = waterquality_phosphorus(self)
        self.erosed = waterquality_erosed(model)
	
    def initial(self):

        if self.var.includeWaterQuality:
            # number of nutrient fluxes simulated in CWatM-WQ (sediment, phosphorus); self.var.n_fluxes is defined in initcondition.py
            
            # create water quality variables
            waterQualityVars = ['pre_w1', 'pre_w2', 'pre_w3', 'naturalLandFrac', 'onlyIrr', 'onlyIrrPaddy']
            for variable in waterQualityVars: vars(self.var)[variable] = np.tile(globals.inZero,(4,1))            
            
            # Create sub-modules variables
            if self.var.includePhosphorus:
                phosphorusVars = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3',\
                                  'soil_P_labile1', 'soil_P_labile2', 'soil_P_labile3',\
                                  'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3',\
                                  'EPC1', 'EPC2', 'EPC3', 'runoff_P', 'toGroundwater_P',\
                                  'soil_P_input1', 'soil_P_input2', 'irrigation_P_Applied', 'irrigation_inactiveP_Applied']
                                  
                phosphorusVarsSum = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_labile3', 'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3', 'EPC1', 'EPC2', 'EPC3',\
                                'soil_P_input1', 'soil_P_input2', 'irrigation_P_Applied', 'irrigation_inactiveP_Applied']
                                
                phosphorusVarsSumCat = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_labile3', 'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3', 'EPC1', 'EPC2', 'EPC3',\
                                'w1', 'w2', 'w3', 'soilM1_f', 'soilM2_f', 'soilM3_f']    
            
            
        
                # only applied to landcovers with soil underneath (grassland and managed grasslands are treated has one landcover
                for variable in phosphorusVars: vars(self.var)[variable] = np.tile(globals.inZero,(4,1))

                self.var.sum_soil_TP_urbanLoss = globals.inZero.copy()
            
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
            
            self.var.soilM1_f = np.tile(self.var.soilM1, (4, 1))
            self.var.soilM2_f = np.tile(self.var.soilM2, (4, 1))
            self.var.soilM3_f = np.tile(self.var.soilM3, (4, 1))
     
            # Load managed grassland fraction ###
        
            self.var.managedGrassland = globals.inZero.copy()
            if 'fracManagedGrassland' in binding:
                self.var.managedGrassland = loadmap('fracManagedGrassland')
            
            # multidimensional array of ones with  managed grasslands as feractions and forest land as 0.
            self.var.naturalLandFrac[1] += 1.#self.var.managedGrassland 
            self.var.naturalLandFrac[2:4] += 1.
            self.var.naturalLandFrac[0] += 1.
            
            self.var.onlyIrrPaddy[2] += 1.
            self.var.onlyIrr[3] += 1.
            # create lakesRes sub-compartments
            if checkOption('includeWaterBodies'):
                '''
                    resLakeSubcompartments dimensions:
                    0: nutrient flux
                    1: water body
                    2: routing step
                '''
                self.var.resLakeSubcompartments = np.tile(1., (self.var.n_fluxes, self.var.waterBodyOutC.shape[0], self.var.noRoutingSteps)) 
                
                self.var.ResLake_Conc_LRC = np.tile(np.compress(self.var.compress_LR, globals.inZero.copy()), (self.var.n_fluxes, 1))           
                self.var.resLakeSubcompartments_wtr = np.tile(1., (self.var.waterBodyOutC.shape[0], self.var.noRoutingSteps))
                
                # Initial storage
                if self.var.includeErosed:
                    # soil index is 0
                    self.var.resLakeSubcompartments[0, :, :] = self.var.resLakeSubcompartments[0, :, :] * np.transpose(np.tile(self.var.lakeResStorageC * self.var.ResLake_Conc_LRC[0, :] / self.var.noRoutingSteps, (self.var.noRoutingSteps, 1)))
                    
                if self.var.includePhosphorus:
                    # phosphorus index is 1
                    self.var.resLakeSubcompartments[1, :, :] = self.var.resLakeSubcompartments[1, :, :] * np.transpose(np.tile(self.var.lakeResStorageC * self.var.ResLake_Conc_LRC[1, :]/ self.var.noRoutingSteps, (self.var.noRoutingSteps, 1)))

           
                # Initital sub-compartments storage water [m3]
                #self.var.resLakeSubcompartments_wtr = self.var.resLakeSubcompartments_wtr * np.transpose(np.tile(self.var.lakeResStorageC / self.var.noRoutingSteps, (self.var.noRoutingSteps, 1)))
                
                
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
                
                for variable in phosphorusVarsSumCat:
                    vars(self.var)["sum_" + variable + "_natural"] = vars(self.var)[variable][0] * self.var.fracVegCover[0] + vars(self.var)[variable][1] * (1 - self.var.managedGrassland) * self.var.fracVegCover[1]
                    vars(self.var)["sum_" + variable + "_managed"] = vars(self.var)[variable][1] * self.var.managedGrassland * self.var.fracVegCover[1] + np.nansum(vars(self.var)[variable][2:4] * self.var.fracVegCover[2:4], axis = 0)

    def dynamic(self): 
        
        # landcover transitions
        if self.var.includePhosphorus:
                self.var.soil_P_inactive_urbanLoss = globals.inZero.copy()
                self.var.soil_P_dissolved_urbanLoss = globals.inZero.copy()

                
                self.var.sum_soil_TP_urbanLoss = globals.inZero.copy()

        if dateVar['newYear'] and not dateVar['newStart']:
            delta_LC = self.var.fracVegCover - self.var.fracVegCover_former
            
            # 0 : forest
            # 1 : grassland
            # 2 : paddy
            # 3 : irrNonPaddy
            # 4 : sealed
            # 5 : water 
        
            # DEFINE THERSHOLD TO CHANGE
          
            # End guess landcover transitions  ####
            

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

        
        
        # Calculate gross GW recharge from 3rd soil layer for water quality
        self.var.grossGWrechargeFromSoil =  (self.var.gwRecharge + self.var.capRiseFromGW) * divideArrays(self.var.perc3toGW, self.var.perc3toGW + self.var.prefFlow)
        self.var.grossGWrechargeFromPrefFlow = (self.var.gwRecharge + self.var.capRiseFromGW) * divideArrays(self.var.prefFlow, self.var.perc3toGW + self.var.prefFlow)
        
        # run initial sub-modules
        # Run the Erosion and Sediment Yield (EroSed) module

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
            
            # update nutrients soil stocks due to land cover change; note sigma(deltaP_j) = 0; i.e., P mass is balanced 

            if dateVar['newYear'] and not dateVar['newStart']:
                
                # 0 : forest
                # 1 : grassland
                # 2 : paddy
                # 3 : irrNonPaddy
                # 4 : sealed
                # 5 : water 
                
                
                delta_LC = self.var.fracVegCover - self.var.fracVegCover_former
                
                urbanGain = np.where(delta_LC[4] > 0, delta_LC[4], 0.)
                allLoss = np.nansum(np.where(delta_LC < 0, np.abs(delta_LC), 0.), axis = 0) 
                wghts = divideArrays(np.where(delta_LC < 0, np.abs(delta_LC), 0.), allLoss)[0:4]
                # update all landcovers to new fractions
                wghts_add = globals.inZero.copy()
                

                '''
                    Tot =  sigma(fi*xi) : i E 0,4
                    loss - total weighted urban loss
                    
                    sigma(fi'*xi') = Tot - loss
                    xi/xj = xi'/xj' = Constant
                    
                    x0 = (Tot - Loss) /  (f'0 + f'1/a + f'2/b + f'3/c)
                    a = x0/x1; b = x0/x2; c = x0/x3
                    
                '''
                
                for var in ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3']:
                    toUrban = np.nansum(np.minimum(np.tile(urbanGain, (4,1)) * wghts * vars(self.var)[var], vars(self.var)[var]), axis = 0)
                    a = divideValues(vars(self.var)[var][0], vars(self.var)[var][1])
                    b = divideValues(vars(self.var)[var][0], vars(self.var)[var][2])
                    c = divideValues(vars(self.var)[var][0], vars(self.var)[var][3])
                    
                    Tot =  np.nansum(vars(self.var)[var] * self.var.fracVegCover_former[0:4] ,axis = 0)
                    x0 = divideValues(Tot - toUrban, self.var.fracVegCover[0] + divideValues(self.var.fracVegCover[1], a) +\
                        divideValues(self.var.fracVegCover[2], b) + divideValues(self.var.fracVegCover[3], c))
                    x1 = divideValues(x0, a)
                    x2 = divideValues(x0, b)
                    x3 = divideValues(x0, c)
                    
                    vars(self.var)[var][0] = x0.copy()
                    vars(self.var)[var][1] = x1.copy()
                    vars(self.var)[var][2] = x2.copy()
                    vars(self.var)[var][3] = x3.copy()
                    self.var.soil_P_inactive_urbanLoss += toUrban

                for var in ['soil_P_labile1', 'soil_P_labile2', 'soil_P_labile3']:
                    toUrban = np.nansum(np.minimum(np.tile(urbanGain, (4,1)) * wghts * vars(self.var)[var], vars(self.var)[var]), axis = 0)
                    a = divideValues(vars(self.var)[var][0], vars(self.var)[var][1])
                    b = divideValues(vars(self.var)[var][0], vars(self.var)[var][2])
                    c = divideValues(vars(self.var)[var][0], vars(self.var)[var][3])
                    
                    Tot =  np.nansum(vars(self.var)[var] * self.var.fracVegCover_former[0:4] ,axis = 0)
                    x0 = divideValues(Tot - toUrban, self.var.fracVegCover[0] + divideValues(self.var.fracVegCover[1], a) +\
                        divideValues(self.var.fracVegCover[2], b) + divideValues(self.var.fracVegCover[3], c))
                    x1 = divideValues(x0, a)
                    x2 = divideValues(x0, b)
                    x3 = divideValues(x0, c)
                    
                    vars(self.var)[var][0] = x0.copy()
                    vars(self.var)[var][1] = x1.copy()
                    vars(self.var)[var][2] = x2.copy()
                    vars(self.var)[var][3] = x3.copy()
                    self.var.soil_P_labile_urbanLoss += toUrban

                     
                for var in ['soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3']:
                    toUrban = np.nansum(np.minimum(np.tile(urbanGain, (4,1)) * wghts * vars(self.var)[var], vars(self.var)[var]), axis = 0)
                    a = divideValues(vars(self.var)[var][0], vars(self.var)[var][1])
                    b = divideValues(vars(self.var)[var][0], vars(self.var)[var][2])
                    c = divideValues(vars(self.var)[var][0], vars(self.var)[var][3])
                    
                    Tot =  np.nansum(vars(self.var)[var] * self.var.fracVegCover_former[0:4] ,axis = 0)
                    x0 = divideValues(Tot - toUrban, self.var.fracVegCover[0] + divideValues(self.var.fracVegCover[1], a) +\
                        divideValues(self.var.fracVegCover[2], b) + divideValues(self.var.fracVegCover[3], c))
                    x1 = divideValues(x0, a)
                    x2 = divideValues(x0, b)
                    x3 = divideValues(x0, c)
                    
                    vars(self.var)[var][0] = x0.copy()
                    vars(self.var)[var][1] = x1.copy()
                    vars(self.var)[var][2] = x2.copy()
                    vars(self.var)[var][3] = x3.copy()
                    self.var.soil_P_dissolved_urbanLoss += toUrban

                self.var.sum_soil_TP_urbanLoss = self.var.soil_P_inactive_urbanLoss + self.var.soil_P_labile_urbanLoss + self.var.soil_P_dissolved_urbanLoss
                
            ## End phosphorus balance #######
                
                           

            # sum total soil P stocks [kg / m2] 
            phosphorusVarsSum = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_labile3', 'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3', 'EPC1', 'EPC2', 'EPC3',\
                                'soil_P_input1', 'soil_P_input2', 'irrigation_P_Applied', 'irrigation_inactiveP_Applied']
            
            
            
            for variable in phosphorusVarsSum:
                vars(self.var)["sum_" + variable] = np.nansum(vars(self.var)[variable] * self.var.fracVegCover[0:4], axis = 0)
            
   
            phosphorusVarsSumCat = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_labile3', 'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3', 'EPC1', 'EPC2', 'EPC3',\
                                'w1', 'w2', 'w3', 'soilM1_f', 'soilM2_f', 'soilM3_f']
            
            for variable in phosphorusVarsSumCat:
                vars(self.var)["sum_" + variable + "_natural"] = vars(self.var)[variable][0] * self.var.fracVegCover[0] + vars(self.var)[variable][1] * (1 - self.var.managedGrassland) * self.var.fracVegCover[1]
                vars(self.var)["sum_" + variable + "_managed"] = vars(self.var)[variable][1] * self.var.managedGrassland * self.var.fracVegCover[1] + np.nansum(vars(self.var)[variable][2:4] * self.var.fracVegCover[2:4], axis = 0)

            
            # EPC [liter /mg]
            self.var.sum_EPC1 = self.var.sum_EPC1 / 10**3
            self.var.sum_EPC2 = self.var.sum_EPC2 / 10**3
            self.var.sum_EPC3 = self.var.sum_EPC3 / 10**3
            
            # TDP soil concentration [mg / liter]
            self.var.sum_soil_P_dissolvedConc1 = divideValues(self.var.sum_soil_P_dissolved1, self.var.sum_w1 * self.var.cellArea) * 10**3
            self.var.sum_soil_P_dissolvedConc2 = divideValues(self.var.sum_soil_P_dissolved2, self.var.sum_w2 * self.var.cellArea) * 10**3
            self.var.sum_soil_P_dissolvedConc3 = divideValues(self.var.sum_soil_P_dissolved3, self.var.sum_w3 * self.var.cellArea) * 10**3
            
            # natural landcovers
            self.var.sum_soil_P_dissolvedConc1_natural =  divideValues(self.var.sum_soil_P_dissolved1_natural, self.var.sum_w1_natural * self.var.cellArea) * 10**3
            self.var.sum_soil_P_dissolvedConc2_natural =  divideValues(self.var.sum_soil_P_dissolved2_natural, self.var.sum_w2_natural * self.var.cellArea) * 10**3
            self.var.sum_soil_P_dissolvedConc3_natural =  divideValues(self.var.sum_soil_P_dissolved3_natural, self.var.sum_w3_natural * self.var.cellArea) * 10**3
            
            self.var.sum_soil_P_dissolvedConc_natural = divideValues(self.var.sum_soil_P_dissolved1_natural + self.var.sum_soil_P_dissolved2_natural + self.var.sum_soil_P_dissolved3_natural,\
                (self.var.sum_w1_natural + self.var.sum_w2_natural + self.var.sum_w3_natural)  * self.var.cellArea) * 10**3
                
            # managed landcovers
            self.var.sum_soil_P_dissolvedConc1_managed =  divideValues(self.var.sum_soil_P_dissolved1_managed, self.var.sum_w1_managed * self.var.cellArea) * 10**3
            self.var.sum_soil_P_dissolvedConc2_managed =  divideValues(self.var.sum_soil_P_dissolved2_managed, self.var.sum_w2_managed * self.var.cellArea) * 10**3
            self.var.sum_soil_P_dissolvedConc3_managed =  divideValues(self.var.sum_soil_P_dissolved3_managed, self.var.sum_w3_managed * self.var.cellArea) * 10**3
            
            self.var.sum_soil_P_dissolvedConc_managed = divideValues(self.var.sum_soil_P_dissolved1_managed + self.var.sum_soil_P_dissolved2_managed + self.var.sum_soil_P_dissolved3_managed,\
                (self.var.sum_w1_managed + self.var.sum_w2_managed + self.var.sum_w3_managed)  * self.var.cellArea) * 10**3
            
            # PP soil concentration [mg / kg soil] - miligram
            self.var.sum_soil_P_labileConc1 = divideValues(self.var.sum_soil_P_labile1, self.var.soilM1) * 10**6
            self.var.sum_soil_P_labileConc2 = divideValues(self.var.sum_soil_P_labile2, self.var.soilM2) * 10**6
            self.var.sum_soil_P_labileConc3 = divideValues(self.var.sum_soil_P_labile3, self.var.soilM3) * 10**6
            
            # natural landcovers
            self.var.sum_soil_P_labileConc1_natural =  divideValues(self.var.sum_soil_P_labile1_natural, self.var.sum_soilM1_f_natural) * 10**6
            self.var.sum_soil_P_labileConc2_natural =  divideValues(self.var.sum_soil_P_labile2_natural, self.var.sum_soilM2_f_natural) * 10**6
            self.var.sum_soil_P_labileConc3_natural =  divideValues(self.var.sum_soil_P_labile3_natural, self.var.sum_soilM3_f_natural) * 10**6
            
            self.var.sum_soil_P_labileConc_natural = divideValues(self.var.sum_soil_P_labile1_natural + self.var.sum_soil_P_labile2_natural + self.var.sum_soil_P_labile3_natural,\
                self.var.sum_soilM1_f_natural + self.var.sum_soilM2_f_natural + self.var.sum_soilM3_f_natural) * 10**6
            
            # managed landcovers
            self.var.sum_soil_P_labileConc1_managed =  divideValues(self.var.sum_soil_P_labile1_managed, self.var.sum_soilM1_f_managed) * 10**6
            self.var.sum_soil_P_labileConc2_managed =  divideValues(self.var.sum_soil_P_labile2_managed, self.var.sum_soilM2_f_managed) * 10**6
            self.var.sum_soil_P_labileConc3_managed =  divideValues(self.var.sum_soil_P_labile3_managed, self.var.sum_soilM3_f_managed) * 10**6
            
            self.var.sum_soil_P_labileConc_managed = divideValues(self.var.sum_soil_P_labile1_managed + self.var.sum_soil_P_labile2_managed + self.var.sum_soil_P_labile3_managed,\
                (self.var.sum_w1_managed + self.var.sum_w2_managed + self.var.sum_w3_managed)  * self.var.cellArea) * 10**3
                
            # P inactive soil concentration [mg / kg soil] - miligram 
            self.var.sum_soil_P_inactiveConc1 = divideValues(self.var.sum_soil_P_inactive1, self.var.soilM1) * 10**6
            self.var.sum_soil_P_inactiveConc2 = divideValues(self.var.sum_soil_P_inactive2, self.var.soilM2) * 10**6
            self.var.sum_soil_P_inactiveConc3 = divideValues(self.var.sum_soil_P_inactive3, self.var.soilM3) * 10**6
            
            # natural landcovers
            self.var.sum_soil_P_inactiveConc1_natural =  divideValues(self.var.sum_soil_P_inactive1_natural, self.var.sum_soilM1_f_natural) * 10**6
            self.var.sum_soil_P_inactiveConc2_natural =  divideValues(self.var.sum_soil_P_inactive2_natural, self.var.sum_soilM2_f_natural) * 10**6
            self.var.sum_soil_P_inactiveConc3_natural =  divideValues(self.var.sum_soil_P_inactive3_natural, self.var.sum_soilM3_f_natural) * 10**6
            
            self.var.sum_soil_P_inactiveConc_natural = divideValues(self.var.sum_soil_P_inactive1_natural + self.var.sum_soil_P_inactive2_natural + self.var.sum_soil_P_inactive3_natural,\
                self.var.sum_soilM1_f_natural + self.var.sum_soilM2_f_natural + self.var.sum_soilM3_f_natural) * 10**6
            
            # managed landcovers
            self.var.sum_soil_P_inactiveConc1_managed =  divideValues(self.var.sum_soil_P_inactive1_managed, self.var.sum_soilM1_f_managed) * 10**6
            self.var.sum_soil_P_inactiveConc2_managed =  divideValues(self.var.sum_soil_P_inactive2_managed, self.var.sum_soilM2_f_managed) * 10**6
            self.var.sum_soil_P_inactiveConc3_managed =  divideValues(self.var.sum_soil_P_inactive3_managed, self.var.sum_soilM3_f_managed) * 10**6
            
            self.var.sum_soil_P_inactiveConc_managed = divideValues(self.var.sum_soil_P_inactive1_managed + self.var.sum_soil_P_inactive2_managed + self.var.sum_soil_P_inactive3_managed,\
                self.var.sum_soilM1_f_managed + self.var.sum_soilM2_f_managed + self.var.sum_soilM3_f_managed) * 10**6
            
            # sum soil P  [kg]
            self.var.tot_soil_P_dissolved = self.var.sum_soil_P_dissolved1 + self.var.sum_soil_P_dissolved2 + self.var.sum_soil_P_dissolved3
            self.var.tot_soil_P_labile = self.var.sum_soil_P_labile1 + self.var.sum_soil_P_labile2 + self.var.sum_soil_P_labile3
            self.var.tot_soil_P_inactive = self.var.sum_soil_P_inactive1 + self.var.sum_soil_P_inactive2 + self.var.sum_soil_P_inactive3
            
            self.var.tot_soil_P_dissolvedConc = divideValues(self.var.tot_soil_P_dissolved, (self.var.sum_w1 + self.var.sum_w2 + self.var.sum_w3)* self.var.cellArea) * 10**3
            self.var.tot_soil_P_labileConc = divideValues(self.var.tot_soil_P_labile, self.var.soilM1 + self.var.soilM2 + self.var.soilM3) * 10**6
            self.var.tot_soil_P_inactiveConc = divideValues(self.var.tot_soil_P_inactive, self.var.soilM1 + self.var.soilM2 + self.var.soilM3) * 10**6

            
            # sum soil P input [kg]
            self.var.tot_soil_P_input = self.var.sum_soil_P_input1 + self.var.sum_soil_P_input2
