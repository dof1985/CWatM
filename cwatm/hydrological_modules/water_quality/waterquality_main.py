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
                                  'soil_P_input1', 'soil_P_input2', 'irrigation_P_Applied']
                                  
                phosphorusVarsSum = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_labile3', 'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3', 'EPC1', 'EPC2', 'EPC3',\
                                'soil_P_input1', 'soil_P_input2', 'irrigation_P_Applied']
                                
                phosphorusVarsSumCat = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_labile3', 'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3', 'EPC1', 'EPC2', 'EPC3',\
                                'w1', 'w2', 'w3', 'soilM1_f', 'soilM2_f', 'soilM3_f']    
            
            
        
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
        
        ## Calculate landcover change fractions to adjust soil nutrient stocks
        # Save previous year landcover at the last day of each year
        #if dateVar['doy'] == dateVar['daysInYear']:
        #    self.var.pre_fracVegCover = self.var.fracVegCover.copy()
        
        if self.var.includePhosphorus:
            self.var.soil_TP_urbanLoss = globals.inZero.copy()
        
        '''
        # calculate change fractions and landcover to landcover changes
        if dateVar['newYear'] and not dateVar['newStart']:
            
            delta_fracVegCover = (self.var.fracVegCover - self.var.pre_fracVegCover) 
            
            natural = delta_fracVegCover[0] + delta_fracVegCover[1]
            managed = delta_fracVegCover[2] + delta_fracVegCover[3]
            urban = delta_fracVegCover[4] + delta_fracVegCover[5]
            # natural/managed loss, urban gained
            natural = np.where(natural < 0, -1 * natural, 0)
            managed = np.where(managed < 0, -1 * managed, 0)
            totLoss = managed + natural
            urban = np.where(urban > 0, urban, 0)
            
            if self.var.includePhosphorus:
                self.var.soil_TP_urbanLoss = urban * (divideValues(managed, totLoss) *\
                (self.var.sum_soil_P_dissolved1_managed + self.var.sum_soil_P_dissolved2_managed + self.var.sum_soil_P_dissolved3_managed +\
                        self.var.sum_soil_P_labile1_managed + self.var.sum_soil_P_labile2_managed + self.var.sum_soil_P_labile3_managed) +\
                        divideValues(natural, totLoss) * (self.var.sum_soil_P_dissolved1_natural + self.var.sum_soil_P_dissolved2_natural + self.var.sum_soil_P_dissolved3_natural +\
                        self.var.sum_soil_P_labile1_natural + self.var.sum_soil_P_labile2_natural + self.var.sum_soil_P_labile3_natural))
            
            self.var.gain_fracVegCover = divideArrays(np.where(delta_fracVegCover > 0, delta_fracVegCover, 0.), self.var.pre_fracVegCover)
           
            #0 - grassland (natural & managed)
            #1 - forest (natural)
            #2 - irrPaddy (managed)
            #3 - irrNonPaddy (managed)
            #4 - sealed (sealed)
            #5 - water (sealed)
            
          
            
            loss_fracVegCoverManaged = np.nansum(np.abs(np.where(delta_fracVegCover < 0, delta_fracVegCover, 0.))[2:4], axis = 0)
        '''    
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
            '''
            if dateVar['newYear'] and not dateVar['newStart']:
                # calculate last time step total P concentration
                pre_tot_soil_P_dissolved = divideValues(self.var.tot_soil_P_dissolved, (self.var.sum_w1 + self.var.sum_w2 + self.var.sum_w3) * self.var.cellArea) * 10**3
                pre_tot_soil_P_labile = divideValues(self.var.tot_soil_P_labile, (self.var.soilM1 + self.var.soilM2 + self.var.soilM3)) * 10**6
                # calculate up-to-date total P concentration
             
                phosphorusVarsSum = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_labile3', 'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3', 'EPC1', 'EPC2', 'EPC3',\
                                'soil_P_input1', 'soil_P_input2', 'irrigation_P_Applied']
            
                phosphorusVarsSum = ['soil_P_labile1', 'soil_P_labile2', 'soil_P_labile3', 'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3']
                                
                for variable in phosphorusVarsSum:
                    vars(self.var)["sum_" + variable] = np.nansum(vars(self.var)[variable] * self.var.fracVegCover[0:4], axis = 0)
                
                tot_soil_P_dissolved = divideValues(self.var.sum_soil_P_dissolved1 + self.var.sum_soil_P_dissolved2 + self.var.sum_soil_P_dissolved3,\
                    (self.var.sum_w1 + self.var.sum_w2 + self.var.sum_w3) * self.var.cellArea) * 10**3
                tot_soil_P_labile = divideValues(self.var.sum_soil_P_labile1 + self.var.sum_soil_P_labile2 + self.var.sum_soil_P_labile3,\
                    (self.var.soilM1 + self.var.soilM2 + self.var.soilM3)) * 10**6
                
                # calculate multipliers
                multiplier_dissolved = divideValues(tot_soil_P_dissolved, pre_tot_soil_P_dissolved)
                multiplier_labile = divideValues(tot_soil_P_labile, pre_tot_soil_P_labile)
                
                shr_dissolve = divideValues(tot_soil_P_dissolved, tot_soil_P_dissolved + tot_soil_P_labile)
                multiplier_avg = shr_dissolve * multiplier_dissolved + (1 - shr_dissolve) * multiplier_labile
                
                print("multipliuer TDP [" + str(np.nanmin(multiplier_dissolved)) + ", " +  str(np.nanmax(multiplier_dissolved)) + "]" + \
                    str(np.nanmean(multiplier_dissolved)))
                print("multipliuer labile [" + str(np.nanmin(multiplier_labile)) + ", " +  str(np.nanmax(multiplier_labile)) + "]" + \
                    str(np.nanmean(multiplier_labile)))
                print("multipliuer avg [" + str(np.nanmin(multiplier_avg)) + ", " +  str(np.nanmax(multiplier_avg)) + "]" + \
                    str(np.nanmean(multiplier_avg)))
               
                # update all P stocks
                self.var.soil_P_dissolved1 = self.var.soil_P_dissolved1 * multiplier_avg
                self.var.soil_P_dissolved2 = self.var.soil_P_dissolved2 * multiplier_avg
                self.var.soil_P_dissolved3 = self.var.soil_P_dissolved3 * multiplier_avg
                self.var.soil_P_labile1 = self.var.soil_P_labile1 * multiplier_avg
                self.var.soil_P_labile2 = self.var.soil_P_labile2 * multiplier_avg
                self.var.soil_P_labile3 = self.var.soil_P_labile3 * multiplier_avg
                
              
                
                # Assume P concentrations are unchanged - TDP
                
                TP_pre = np.nansum(self.var.sum_soil_P_dissolved1 + self.var.sum_soil_P_dissolved2 + self.var.sum_soil_P_dissolved3 +\
                                self.var.sum_soil_P_labile1 + self.var.sum_soil_P_labile2 + self.var.sum_soil_P_labile3, axis = 0)
                
                sum_soil_P_dissolved1 =  self.var.sum_soil_P_dissolvedConc1 * 10**-3 * self.var.sum_w1 * self.var.cellArea
                sum_soil_P_dissolved2 =  self.var.sum_soil_P_dissolvedConc2 * 10**-3 * self.var.sum_w2 * self.var.cellArea
                sum_soil_P_dissolved3 =  self.var.sum_soil_P_dissolvedConc3 * 10**-3 * self.var.sum_w3 * self.var.cellArea
                
                # Assume P concentrations are unchanged - P labile
                sum_soil_P_labile1 =  self.var.sum_soil_P_labileConc1 * 10**-6 * self.var.soilM1
                sum_soil_P_labile2 =  self.var.sum_soil_P_labileConc2 * 10**-6 * self.var.soilM2
                sum_soil_P_labile3 =  self.var.sum_soil_P_labileConc3 * 10**-6 * self.var.soilM3
                
                TP_post = np.nansum(sum_soil_P_dissolved1 + sum_soil_P_dissolved2 + sum_soil_P_dissolved3 +\
                                sum_soil_P_labile1 + sum_soil_P_labile2 + sum_soil_P_labile3, axis = 0)
                
                self.var.soil_TP_urbanLoss = np.maximum(TP_pre- TP_post, 0.)
                print(np.nansum(self.var.soil_TP_urbanLoss))
                
                def splitSumStock(var, var_target):
                    varOut = np.tile(var, (4, 1))
                    varOut = var * divideArrays(var_target, self.var.fracVegCover[0:4])
                    return(varOut)
                print(np.nanmean(np.nansum(self.var.soil_P_dissolved1 * self.var.pre_fracVegCover[0:4], axis = 0)))    
                print(np.nanmean(sum_soil_P_dissolved1, axis = 0))
                # back-calculate TDP mass
                self.var.soil_P_dissolved1 = splitSumStock(sum_soil_P_dissolved1, self.var.soil_P_dissolved1)
                self.var.soil_P_dissolved2 = splitSumStock(sum_soil_P_dissolved2, self.var.soil_P_dissolved2)
                self.var.soil_P_dissolved3 = splitSumStock(sum_soil_P_dissolved3, self.var.soil_P_dissolved3)

                # back-calculate P labile mass
                self.var.soil_P_labile1 = splitSumStock(sum_soil_P_labile1, self.var.soil_P_labile1)
                self.var.soil_P_labile2 = splitSumStock(sum_soil_P_labile2, self.var.soil_P_labile2)
                self.var.soil_P_labile3 = splitSumStock(sum_soil_P_labile3, self.var.soil_P_labile3)
                
                print(np.nanmean(np.nansum(self.var.soil_P_dissolved1 * self.var.fracVegCover[0:4], axis = 0)))
            '''
            # sum total soil P stocks [kg / m2] 
            phosphorusVarsSum = ['soil_P_inactive1', 'soil_P_inactive2', 'soil_P_inactive3', 'soil_P_labile1', 'soil_P_labile2', \
                                'soil_P_labile3', 'soil_P_dissolved1', 'soil_P_dissolved2', 'soil_P_dissolved3', 'EPC1', 'EPC2', 'EPC3',\
                                'soil_P_input1', 'soil_P_input2', 'irrigation_P_Applied']
                                
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
                self.var.sum_soilM1_f_managed + self.var.sum_soilM2_f_managed + self.var.sum_soilM3_f_managed) * 10**6
            
            # sum soil P  [kg]
            self.var.tot_soil_P_dissolved = self.var.sum_soil_P_dissolved1 + self.var.sum_soil_P_dissolved2 + self.var.sum_soil_P_dissolved3
            self.var.tot_soil_P_labile = self.var.sum_soil_P_labile1 + self.var.sum_soil_P_labile2 + self.var.sum_soil_P_labile3
            
            self.var.tot_soil_P_dissolvedConc = divideValues(self.var.tot_soil_P_dissolved, (self.var.sum_w1 + self.var.sum_w2 + self.var.sum_w3)* self.var.cellArea) * 10**3
            self.var.tot_soil_P_labileConc = divideValues(self.var.tot_soil_P_labile, self.var.soilM1 + self.var.soilM2 + self.var.soilM3) * 10**6

            
            # sum soil P input [kg]
            self.var.tot_soil_P_input = self.var.sum_soil_P_input1 + self.var.sum_soil_P_input2
            
            '''
            if dateVar['newYear'] and not dateVar['newStart']:
                pre_TP =  (pre_tot_soil_P_labile / 10**6) *  (self.var.soilM1 + self.var.soilM2 + self.var.soilM3) +\
                (pre_tot_soil_P_dissolved / 10**3) * (self.var.sum_w1 + self.var.sum_w2 + self.var.sum_w3)
                
                self.var.soil_TP_urbanLoss = pre_TP - (self.var.tot_soil_P_dissolved + self.var.tot_soil_P_labile)
                
            '''
           