# import libraries       
import numpy as np
import os

from cwatm.management_modules.data_handling import *
from cwatm.hydrological_modules.groundwater_modflow.modflow6 import ModFlowSimulation
# impotlib to install libraries on the fly if needed e.g. rasterio
import importlib


"""

**Global variables**

=====================================  ======================================================================  =====
Variable [self.var]                    Description                                                             Unit 
=====================================  ======================================================================  =====
modflow                                Flag: True if modflow_coupling = True in settings file                  --   
sum_gwRecharge                         groundwater recharge                                                    m    
cellArea                               Area of cell                                                            m2   
gwdepth_observations                   Input, gw_depth_observations, groundwater depth observations            m    
gwdepth_adjuster                       Groundwater depth adjuster                                              m    
baseflow                               simulated baseflow (= groundwater discharge to river)                   m    
capillar                               Flow from groundwater to the third CWATM soil layer. Used with MODFLOW  m    
capriseindex                                                                                                   --   
soildepth12                            Total thickness of layer 2 and 3                                        m    
leakageriver_factor                                                                                            --   
leakagelake_factor                                                                                             --   
modflow_timestep                       Chosen ModFlow model timestep (1day, 7days, 30days, etc.)               --   
head                                   Simulated ModFlow water level [masl]                                    m    
gwdepth_adjusted                       Adjusted depth to groundwater table                                     m    
gwdepth                                Depth to groundwater table                                              m    
modflow_cell_area                                                                                              --   
modflowsteady                          True if modflow_steadystate = True in settings file                     --   
Ndays_steady                           Number of steady state run before the transient simulation              --   
channel_ratio                                                                                                  --   
modflowtotalSoilThickness              Array (nrows, ncol) used to compute water table depth in post-processi  m    
load_init_water_table                                                                                          --   
GW_pumping                             Input, True if Groundwater_pumping=True                                 bool 
use_complex_solver_for_modflow                                                                                 --   
use_super_complex_solver_for_modflow                                                                           --   
availableGWStorageFraction                                                                                     --   
wells_index                                                                                                    --   
depth                                                                                                          --   
sumed_sum_gwRecharge                                                                                           --   
modflow_compteur                       Counts each day relatively to the chosen ModFlow timestep, allow to ru  --   
modflow_watertable                                                                                             --   
writeerror                                                                                                     --   
modflowdiscrepancy                                                                                             --   
groundwater_storage_top_layer                                                                                  --   
groundwater_storage_available          Groundwater storage. Used with MODFLOW.                                 m    
gwstorage_full                         Groundwater storage at full capacity                                    m    
permeability                                                                                                   --   
modfPumpingM_actual                    Actual groundwater pumping. Used with MODFLOW.                          m    
gwdepth_difference_sim_obs             Difference between simulated and observed groundwater table             m    
modflow_head_adjusted                                                                                          --   
waterdemandFixed                                                                                               --   
modfPumpingM                                                                                                   --   
=====================================  ======================================================================  =====

**Functions**
"""
    
def parseArray(arr, splt = ";"):
    return arr.split(splt)

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def decompress(map, nanvalue=None):
    """
    Decompressing CWatM maps from 1D to 2D with missing values

    :param map: compressed map
    :return: decompressed 2D map
    """

    dmap = maskinfo['maskall'].copy()
    dmap[~maskinfo['maskflat']] = map[:]
    if nanvalue is not None:
        dmap.data[np.isnan(dmap.data)] = nanvalue

    return dmap.data

def load_modflow_basin(self, var):
    rasterio = importlib.import_module("rasterio", package=None)
    with rasterio.open(cbinding(var), 'r') as src:
        mbas = src.read(1).astype(bool)  # read in as 2-dimensional array (nrows, ncols).
        self.domain = {
            'rowsize': abs(src.profile['transform'].e),
            'colsize': abs(src.profile['transform'].a),
            'nrow': int(src.profile['height']),
            'ncol': int(src.profile['width']),
            'west': src.profile['transform'].c,
            'east': src.profile['transform'].c + (src.profile['width']-1) * abs(src.profile['transform'].a),
            'north': src.profile['transform'].f,
            'south': src.profile['transform'].f - (src.profile['height']-1) * abs(src.profile['transform'].e)
        }
        domain['rowsize'] = abs(src.profile['transform'].e)
        domain['colsize'] = abs(src.profile['transform'].a)
        domain['nrow'] = int(src.profile['height'])
        domain['ncol'] = int(src.profile['width'])
        domain['west'] = src.profile['transform'].c
        
        domain['east'] = src.profile['transform'].c + (src.profile['width']-1) * abs(src.profile['transform'].a)
        domain['north'] = src.profile['transform'].f
        domain['south'] = src.profile['transform'].f - (src.profile['height']-1) * abs(src.profile['transform'].e)
    return mbas
  
def load_aquifer_coeff(self, var, nlay, mult = 1, mask = None, splt = ";"):
    # can be used only after loading a mask
    
    varList = parseArray(cbinding(var), splt = splt)  # ['0.08', 'path...']
    emptyList = [None] * len(varList)
    lyrs = nlay
    if nlay < len(varList):
        raise ValueError(f'{var} has {len(varList)} parameters, but the model {nlay} layers')
    if nlay == len(varList):
        lyrs = 0
    for i in range(len(varList)):
        #print(i)
        val = varList[i].strip()
        
        #print(val)
        if is_float(val):
           # create map from float
           emptyList[i] = map_from_param(self = self, var = float(val), nlay_ = lyrs, mult = mult)
          # print(emptyList[i].shape)
        else:
           emptyList[i] = load_lyr_map(var = val, mask = mask, mult = mult, parsePath = False)
           #print(emptyList[i].shape)
    if nlay == 1:
        return np.array([np.array(emptyList)[0]])
        
    return np.array(emptyList)
    
def load_lyr_map(var, mask = None, parsePath = True, mult = 1):
    rasterio = importlib.import_module("rasterio", package=None)
    pth = var
    if parsePath:
        pth = cbinding(pth)
    with rasterio.open(pth, 'r') as src:
        lyrmap = src.read(1).astype(np.float32)
        if not mask is None:
            lyrmap[mask == False] = np.nan
        lyrmap = lyrmap * mult
    return lyrmap

def map_from_param(self, var, nlay_, mult = 1):
    if not is_float(var):
        var = float(cbinding(var))
    prm = var * mult
    #print(var)
    #print(mult)
    #print(prm)
    if nlay_ == 0:
        return np.full((self.domain['nrow'], self.domain['ncol']), prm)
    return np.full((nlay_, self.domain['nrow'], self.domain['ncol']), prm)
    
class groundwater_modflow:

    """
    **Functions**
    """
    
    def __init__(self, model):
        self.var = model.var
        self.model = model
    
    def get_corrected_modflow_cell_area(self):
        return np.bincount(
            self.indices['ModFlow_index'],
            weights=np.invert(self.var.mask.astype(bool)).ravel()[self.indices['CWatM_index']] * self.indices['area'],
            minlength=self.modflow.basin.size
        ).reshape(self.modflow.basin.shape)

    def get_corrected_cwatm_cell_area(self):
        return (self.var.cellArea_uncompressed.ravel() - np.bincount(
            self.indices['CWatM_index'],
            weights=np.invert(self.modflow.basin).ravel()[self.indices['ModFlow_index']] * self.indices['area'],
            minlength=self.var.mask.size
        )).reshape(self.var.mask.shape)
        
    def CWATM2modflow(self, variable, correct_boundary=False):
        # Converting flow [L/T] from 2D CWatM map to 2D ModFlow map
        if correct_boundary:
            self.var.modflow_cell_area = self.corrected_modflow_cell_area.ravel()
        else:
            self.var.modflow_cell_area = self.domain['rowsize'] * self.domain['colsize']  # in m2

        array = (np.bincount(
            self.indices['ModFlow_index'],
            variable.ravel()[self.indices['CWatM_index']] * self.indices['area'],
            minlength=self.domain['nrow'] * self.domain['ncol']
        ) / self.var.modflow_cell_area).reshape((self.domain['nrow'], self.domain['ncol'])).astype(variable.dtype)
        return array

    def modflow2CWATM(self, variable, correct_boundary=False):
        # Converting flow [L/T] from 2D ModFLow map to 2D CWatM map
        basin_map = self.modflow.basin.copy()
        basin_map = basin_map[0]
        variable_copy = variable.copy()
        variable_copy[basin_map== False] = 0
        assert not (np.isnan(variable_copy).any())
        assert basin_map.dtype == bool
        if correct_boundary:
            cwatm_cell_area = self.corrected_cwatm_cell_area.ravel()
        else:
            cwatm_cell_area = self.var.cellArea.ravel()  # in m2

        array = (np.bincount(
            self.indices['CWatM_index'],
            weights=variable_copy.ravel()[self.indices['ModFlow_index']] * self.indices['area'],
            minlength=maskinfo['shape'][0] * maskinfo['shape'][1]) / decompress(cwatm_cell_area, nanvalue=0)).reshape(
            maskinfo['shape']).astype(variable_copy.dtype)

        array[maskinfo['mask'] == 1] = np.nan

        return array

    def modflow2CWATMbis(self, variable, correct_boundary=False):
        """
        Converting the 2D ModFLow capillary rise map into the fraction of area where capillary rise occurs
        in the 2D CWatM maps. return a fraction for each CWatM cell, input is the ModFlow capillary rise map
        the returned array is used to apply leakage from water bodies to the ModFlow layer
        """
        variable_copy = variable.copy()
        variable_copy[self.modflow.basin == False] = 0
        # Each ModFlow cell is distinguished between producing or not producing capillary rise
        variable_copy[variable_copy > 0] = 1
        assert not (np.isnan(variable_copy).any())
        assert self.modflow.basin.dtype == bool
        if correct_boundary:
            cwatm_cell_area = self.corrected_cwatm_cell_area.ravel()
        else:
            cwatm_cell_area = self.var.cellArea.ravel()  # in m2
        array = (np.bincount(
            self.indices['CWatM_index'],
            weights=variable_copy.ravel()[self.indices['ModFlow_index']] * self.indices['area'],
            minlength=maskinfo['shape'][0] * maskinfo['shape'][1]
        ) / decompress(cwatm_cell_area, nanvalue=0)).reshape(maskinfo['shape']).astype(variable_copy.dtype)
        array[maskinfo['mask'] == 1] = np.nan
        return array
    
    def calcSaturatedCellFraction(self, lyr, head, omega = 10**-6):
    
        ''' 
        Applied from 'Documentaton for the MODFLOW 6 Groundwater Flow Model | Ch. 55 of Section A, Groundwater, Book6, Modeling Techniques'
        https://pubs.usgs.gov/tm/06/a55/tm6a55.pdf
           
        '''
        
        Aomega = 1 / (1-omega)
        # calculate cell saturated thickness - pp.62-63; eq.  4-4
        dv_n = np.maximum(np.minimum(head[lyr], self.layer_boundaries[lyr]) - self.layer_boundaries[lyr + 1], 0.)
        
        # calculate cell saturated fraction - pp.62-63; eq.  4-6
        sfn = divideArrays(dv_n, self.layer_boundaries[lyr] - self.layer_boundaries[lyr + 1])
        
        # use quadratic smooth function - pp.62-63; eq.  4-5
        
        sfn_quad = np.where(sfn < 0,  0., 
        np.where(sfn < omega, (Aomega / (2 * omega)) * sfn ** 2 , 
        np.where(sfn < 1 - omega, Aomega * sfn + 0.5 * (1 -Aomega) , 
        np.where(sfn < 1, 1 - ((Aomega/(2*omega)) * (1 -sfn) ** 2), 1.))))
        
        # recalculate cell saturated thickness - pp. 64 eq. 4-8

        #v_n_new = sfn_quad * (self.layer_boundaries[lyr] - self.layer_boundaries[lyr + 1])
        return(sfn_quad)
        
    def initial(self):
        

        # check if we we are using steady state option. Not yet implemented in this version, the user should provide an
        # estimate of the initial water table ("head" variable in the initial part)

        # ModFlow6 version is only daily currently
        self.var.modflowsteady = False  # returnBool('modflow_steadystate')
        self.var.modflow_timestep = 1  # int(loadmap('modflow_timestep'))

        self.var.Ndays_steady = 0  # int(loadmap('Ndays_steady'))

        if self.var.modflow:
            
            
            verboseGW = False
            if 'verbose_GW' in binding:
                verboseGW = returnBool('verbose_GW')

            # Define the size of the ModFlow time step - means that ModFlow is called each "modflow_timestep"
            self.var.modflow_timestep = int(loadmap('modflow_timestep'))
            if verboseGW:
                print('ModFlow is activated')
                print('ModFlow timestep is : ', self.var.modflow_timestep, ' days\n')

            #modflow_directory = cbinding('PathGroundwaterModflow')
            modflow_directory_output = cbinding('PathGroundwaterModflowOutput')
            directory_mf6dll = cbinding('path_mf6dll')
            if not (os.path.isdir(directory_mf6dll)):
                msg = "Error 222: Path to Modflow6 files does not exists "
                raise CWATMDirError(directory_mf6dll, msg, sname='path_mf6dll')
            
            nlay = 1
            if 'nlay' in binding:
                nlay = int(loadmap('nlay'))
                
            modflow_basin = load_modflow_basin(self, var = 'modflow_basin')
            
           
            nr = modflow_basin.shape[0]
            nc = modflow_basin.shape[1]
            self.modflow_basin = np.full((nlay, nr, nc), modflow_basin)
          
            print("NUMBER OF ACTIVE CELLS IN MODEL: " + str(np.nansum(self.modflow_basin)))
            
            # load topography
            topography = load_lyr_map(var = 'topo_modflow', mask = self.modflow_basin[0])  
            
            # load thickness
            self.thickness = load_aquifer_coeff(self, var = 'thickness', nlay = nlay)
            
            # set minimum values for missing values
            self.thickness[np.isnan(self.thickness)] = 50
            self.thickness[self.thickness == 0] = 50
            self.thickness[self.thickness < 0] = 50
            
            # currently not available
            self.confinedAquifer_flags = self.modflow_basin.copy()
            self.confinedAquifer_flags[self.confinedAquifer_flags > 0] = 1.
            
            # load channel ratio
            self.var.channel_ratio = load_lyr_map(var = 'chanRatio')
            self.var.waterdemandFixed = False
            # if a correcting add river_percent_factor is in settingsfile add this to self.var.channel_ratio,
            # but 0<= self.var.channel_ratio <= 1 correcting the channellratio from the value indicated in the
            # settings file
            factor_channelratio = 0.
            if "river_percent_factor" in binding:
                factor_channelratio = float(binding['river_percent_factor'])
                self.var.channel_ratio = np.where(self.var.channel_ratio > 0,
                                            np.where(self.var.channel_ratio + factor_channelratio > 1, 1,
                                                    np.where(self.var.channel_ratio + factor_channelratio < 0, 0,
                                                            self.var.channel_ratio + factor_channelratio)), 0)
                                                            
            # Coef to multiply transmissivity and storage coefficient (because ModFlow convergence is better if
            # aquifer's thicknes is big and permeability is small)
            self.coefficient = 1
            
            # load permeability
            self.permeability = load_aquifer_coeff(self, var = 'permeability', nlay = nlay,  mult =  24 * 3600 / self.coefficient)
            
            # set minimum values for missing values
            self.permeability[np.isnan(self.permeability)] = 1e-6 * 86400 / self.coefficient
            self.permeability[self.permeability == 0] = 1e-6 * 86400 / self.coefficient
            self.permeability[self.permeability < 0] = 1e-6 * 86400 / self.coefficient
       
  
            if 'permeability_vertical' in binding:
                self.permeability_v = load_aquifer_coeff(self, var = 'permeability_vertical', nlay = nlay,  mult =  24 * 3600 / self.coefficient)
                
                # set minimum values for missing values
                self.permeability_v[np.isnan(self.permeability_v)] = 1e-6 * 86400 / self.coefficient
                self.permeability_v[self.permeability_v == 0] = 1e-6 * 86400 / self.coefficient
                self.permeability_v[self.permeability_v < 0] = 1e-6 * 86400 / self.coefficient
            else:
                self.permeability_v = self.permeability.copy()
            
            
            # load porosity map/parameter  
            self.porosity = load_aquifer_coeff(self, var = 'poro', nlay = nlay)
            
            # set minimum values for missing values
            self.porosity[np.isnan(self.porosity)] = 0.02
            self.porosity[self.porosity == 0] = 0.02
            self.porosity[self.porosity < 0] = 0.02
            
            # load specific_yield
            self.s_yield = self.porosity.copy()
            if 'specific_yield' in binding:
                self.s_yield = load_aquifer_coeff(self = self, var = 'specific_yield', nlay = nlay, mask = self.modflow_basin[0])
            
            # load specific_storage
            self.s_stor = self.modflow_basin.copy()
            self.s_stor[self.s_stor > 0] = 0. 

            if 'specific_storage' in binding:
                self.s_stor = load_aquifer_coeff(self = self, var = 'specific_storage', nlay = nlay, mask = self.modflow_basin[0])
            # uploading arrays allowing to transform 2D arrays from ModFlow to CWatM and conversely
            modflow_x = np.load(os.path.join(cbinding('cwatm_modflow_indices'), 'modflow_x.npy'))
            modflow_y = np.load(os.path.join(cbinding('cwatm_modflow_indices'), 'modflow_y.npy'))
            cwatm_x = np.load(os.path.join(cbinding('cwatm_modflow_indices'), 'cwatm_x.npy'))
            cwatm_y = np.load(os.path.join(cbinding('cwatm_modflow_indices'), 'cwatm_y.npy'))

            self.indices = {
                'area': np.load(os.path.join(cbinding('cwatm_modflow_indices'), 'area.npy')),
                'ModFlow_index': np.array(modflow_y * self.domain['ncol'] + modflow_x),
                'CWatM_index': np.array(cwatm_y * maskinfo['shape'][1] + cwatm_x)
            }

            indices_cell_area = np.bincount(self.indices['CWatM_index'], weights=self.indices['area'],
                                            minlength=maskinfo['mapC'][0])
            area_correction = (decompress(self.var.cellArea, nanvalue=0) / indices_cell_area)[self.indices['CWatM_index']]
            self.indices['area'] = self.indices['area'] * (area_correction + (1-area_correction) * 0.5)

            ## load reduced drainange parameter
            self.reducedDrainage = 0.
            if 'reduceDrainange' in binding:
                self.reducedDrainage = load_aquifer_coeff(self, var = 'reduceDrainange', nlay = nlay, mask = self.modflow_basin[0])
                #print(self.reducedDrainage.shape)
                self.reducedDrainage[np.isnan(self.reducedDrainage)] = 0.
                
            # Converting the CWatM soil thickness into ModFlow map, then soil thickness will be removed from
            # topography if there is a lake or a reservoir soil depth should be replace (instead of 0) by the
            # averaged soil depth (if not the topography under the lake is above neighboring cells)
            soildepth_as_GWtop = False
            if 'use_soildepth_as_GWtop' in binding:
                soildepth_as_GWtop = returnBool('use_soildepth_as_GWtop')
            correct_depth_underlakes = False
            if 'correct_soildepth_underlakes' in binding:
                correct_depth_underlakes = returnBool('correct_soildepth_underlakes')

            if soildepth_as_GWtop:  # topographic minus soil depth map is used as groundwater upper boundary
                # in some regions or models soil depth is around zeros under lakes, so it should be similar to neighboring cells
                if correct_depth_underlakes:
                    if verboseGW:
                        print('Topography minus soil depth is the upper limit of groundwater. Correcting depth under lakes.')
                    
                    waterBodyID_temp = loadmap('waterBodyID').astype(np.int64)
                    soil_depth_temp = np.where(waterBodyID_temp != 0,
                                               np.nanmedian(self.var.soildepth12) - loadmap('depth_underlakes'),
                                               self.var.soildepth12)
                    soil_depth_temp = np.where(self.var.soildepth12 < 0.4, np.nanmedian(self.var.soildepth12),
                                               self.var.soildepth12)  # some cells around lake have small soil depths
                    soildepth_modflow = self.CWATM2modflow(decompress(soil_depth_temp)) + 0.05
                    soildepth_modflow[np.isnan(soildepth_modflow)] = 0
                else:
                    if verboseGW:
                        print('Topography minus soil depth is used as upper limit of groundwater. No correction of depth under lakes')
                    soildepth_modflow = self.CWATM2modflow(decompress(self.var.soildepth12)) + 0.05
                    soildepth_modflow[np.isnan(soildepth_modflow)] = 0
            else:  # topographic map is used as groundwater upper boundary
                if correct_depth_underlakes:  # we make a manual correction
                    if verboseGW:
                        print('Topography is used as upper limit of groundwater. Correcting depth under lakes. It can make ModFlow difficulties to converge')
                  
                    waterBodyID_temp = loadmap('waterBodyID').astype(np.int64)
                    
                    soil_depth_temp = np.where(waterBodyID_temp != 0, loadmap('depth_underlakes'), 0)
                   
                    soildepth_modflow = self.CWATM2modflow(decompress(soil_depth_temp))
                    soildepth_modflow[np.isnan(soildepth_modflow)] = 0
                else:
                    if verboseGW:
                        print('=> Topography is used as upper limit of groundwater. No correction of depth under lakes')
                    soildepth_modflow = np.zeros((self.domain['nrow'], self.domain['ncol']), dtype=np.float32)
           
            # defining the top of the ModFlow layer
            self.layer_boundaries = np.empty((nlay + 1, self.domain['nrow'], self.domain['ncol']), dtype=np.float32)
            self.layer_boundaries[0] = topography - soildepth_modflow
            
            self.layer_boundaries[0] = np.where(np.isnan(self.layer_boundaries[0]), 0, self.layer_boundaries[0])
            # defining the bottom of the ModFlow layer
            for lyr in range(nlay):            
                self.layer_boundaries[lyr + 1] = self.layer_boundaries[lyr] - self.thickness[lyr]

            if not correct_depth_underlakes:  # we make a manual correction
                waterBodyID_temp = loadmap('waterBodyID').astype(np.int64)
            lake_modf = np.where(waterBodyID_temp != 0, 1, 0)
            lake_modf = self.CWATM2modflow(decompress(lake_modf))
         
            soildepth_modflow[np.isnan(lake_modf)] = 0
            self.layer_boundaries[0] = np.where(lake_modf <= 0, np.where(self.var.channel_ratio > 0,
                                                                         self.layer_boundaries[0] - 1,
                                                                         self.layer_boundaries[0]),
                                                self.layer_boundaries[0])
  
            self.layer_boundaries[0] = self.layer_boundaries[0] - self.reducedDrainage[0]
            # saving soil thickness at modflow resolution to compute water table depth in postprocessing
            self.var.modflowtotalSoilThickness = soildepth_modflow

            # defining the initial water table map (it can be a hydraulic head map previously simulated)
            self.var.load_init_water_table = False
            if 'load_init_water_table' in binding:
                self.var.load_init_water_table = returnBool('load_init_water_table')

            if self.var.load_init_water_table:
                if verboseGW:
                    print('=> Initial water table depth is uploaded from ', cbinding('init_water_table'))
                watertable = cbinding('init_water_table')
                if watertable.split(".")[-1] == "npy":
                    head = np.load(cbinding('init_water_table'))
                    head =  np.where(np.isnan(head), 0., head)             
                else:
                    # MODIFIED DOR FRIDMAN
                    varList = parseArray(watertable, splt = ';')  # ['0.08', 'path...']
                    
                    emptyList = [None] * len(varList)
                    
                    # error raise if len(varList) < nlay or  >
                    for i in range(nlay):
                        emptyList[i] = self.CWATM2modflow(self, loadmap(varList[i]), fname = True)
                    head = np.array(emptyList)
            else:
                start_watertabledepth = load_aquifer_coeff(self, var = 'initial_water_table_depth',  nlay = nlay)
                
                start_watertabledepthIntFlag = False
                if is_float(parseArray(cbinding('initial_water_table_depth'), splt = ';')[0]):
                    start_watertabledepthIntFlag = True
                
                if start_watertabledepthIntFlag:
                    if verboseGW:
                        print('=> Water table depth is - ', start_watertabledepth, ' m at the begining')
                else:
                    if verboseGW:
                        print('=> The average water table depth is - ', np.nanmean(start_watertabledepth), ' m at the begining')
                # MODIFIED DOR FRIDMAN
                head = self.layer_boundaries[0:nlay] - start_watertabledepth   
                
            # Defining potential leakage under rivers and lakes or reservoirs
            self.var.leakageriver_factor = 0
            if 'leakageriver_permea' in binding:
                self.var.leakageriver_factor = loadmap('leakageriver_permea')  # in m/day
                if verboseGW:
                    print('=> Potential groundwater inflow from rivers is ', self.var.leakageriver_factor, ' m/day')
            
            self.var.leakagelake_factor = 0
            if 'leakagelake_permea' in binding:
                self.var.leakagelake_factor = loadmap('leakagelake_permea')  # in m/day
                if verboseGW:
                    print('=> Groundwater inflow from lakes/reservoirs is ', self.var.leakagelake_factor, ' m/day')
            
            ## Recharge mask
            self.var.rch_index = []
            self.rch_mask = np.copy(self.modflow_basin[0])

            index_modflowcell = 0            
            for ir in range(self.domain['nrow']):
                for ic in range(self.domain['ncol']):
                    if self.modflow_basin[0][ir][ic] == 1:
                        self.rch_mask[ir][ic] = True
                        self.var.rch_index.append(index_modflowcell)
                    else:
                        self.rch_mask[ir][ic] = False
                    index_modflowcell += 1  
            ## END BUILDING RECHARGE MASK
            
            # test if ModFlow pumping is used as defined in settings file
            Groundwater_pumping = False
            if 'Groundwater_pumping' in binding:
                self.var.GW_pumping = returnBool('Groundwater_pumping')
            
            if verboseGW:
                print('=> Groundwater pumping should be deactivated if includeWaterDemand is False')

            self.var.use_complex_solver_for_modflow = False
            if 'use_complex_solver_for_modflow' in option:
                if checkOption('use_complex_solver_for_modflow'):
                    self.var.use_complex_solver_for_modflow = True

            self.var.use_super_complex_solver_for_modflow = False
            if 'use_super_complex_solver_for_modflow' in option:
                if checkOption('use_super_complex_solver_for_modflow'):
                    self.var.use_super_complex_solver_for_modflow = True

            self.var.availableGWStorageFraction = 0.15
            if self.var.GW_pumping:
                self.var.correctPumpingDiscrepancy = False
                if 'correctPumpingDiscrepancy' in binding:
                    self.var.correctPumpingDiscrepancy = returnBool('correctPumpingDiscrepancy')
                   
                if verboseGW:
                    print('THE PUMPING MAP SHOULD BE DEFINED (In transient.py ALSO LINE 420) BEFORE TO RUN THE MODEL AND BE THE SAME FOR ALL THE SIMULATION')
                # creating a mask to set up pumping wells, TO DO MANUALLY HERE OR TO IMPORT AS A MAP, because running the model with zero pumping rates every cells is consuming
                self.wells_mask = np.copy(self.modflow_basin)
                wells_mask_from_file = np.copy(self.modflow_basin)
                # for Burgenland we set up only one well for each CWATM cell (1*1km), thus only one well every ten ModFlow cells
                # TODO NumOfModflowWellsInCWatMCell = 1 #Must be even or 1
                self.var.wells_index = []
                index_modflowcell = 0
                for layer in range(nlay):
                    index_modflowcell = 0
                    for ir in range(self.domain['nrow']):
                        for ic in range(self.domain['ncol']):
                            if self.modflow_basin[layer][ir][ic] == 1 & wells_mask_from_file[layer][ir][ic] == 1: #and int((ir+5.0)/10.0) - (ir+5.0)/10.0 == 0 and int((ic+5.0)/10.0) - (ic+5.0)/10.0 == 0:
                                #if ir != 0 and ic != 0 and ir != self.domain['nrow']-1 and ic != self.domain['ncol']-1:
                                self.wells_mask[layer][ir][ic] = True
                                self.var.wells_index.append(index_modflowcell)
                            else:
                                self.wells_mask[layer][ir][ic] = False
                            index_modflowcell += 1
                if 'water_table_limit_for_pumping' in binding:
                    # if available storage is too low, no pumping in this cell
                    self.var.availableGWStorageFraction = np.maximum(
                        np.minimum(0.15, loadmap('water_table_limit_for_pumping')),
                        0)  # if 85% of the ModFlow cell is empty, we prevent pumping in this cell
                if verboseGW:
                    print('=> Pumping in the ModFlow layer is prevented if water table is under',
                          int(100 * (self.var.availableGWStorageFraction)), '% of the layer capacity')
                # initializing the ModFlow6 model
                self.modflow = ModFlowSimulation(
                    'transient',
                    modflow_directory_output,
                    directory_mf6dll,
                    ndays=globals.dateVar['intEnd'],
                    timestep=self.var.modflow_timestep,
                    specific_storage=self.s_stor,
                    specific_yield=self.s_yield,
                    nlay=nlay,
                    nrow=self.domain['nrow'],
                    ncol=self.domain['ncol'],
                    rowsize=self.domain['rowsize'],
                    colsize=self.domain['colsize'],
                    top=self.layer_boundaries[0],
                    bottom=self.layer_boundaries[1:],
                    basin=self.modflow_basin,
                    confined_only=self.confinedAquifer_flags,                               
                    head=head,
                    topography=self.layer_boundaries[0],
                    permeability=self.permeability,
                    permeability_vertical=self.permeability_v,
                    load_from_disk=returnBool('load_modflow_from_disk'),
                    setpumpings=True,
                    pumpingloc=self.wells_mask,
                    verbose=verboseGW,
                    complex_solver=self.var.use_complex_solver_for_modflow,
                    super_complex_solver=self.var.use_super_complex_solver_for_modflow)



            else:  # no pumping
                # initializing the ModFlow6 model
                self.modflow = ModFlowSimulation(
                    'transient',
                    modflow_directory_output,
                    directory_mf6dll,
                    ndays=globals.dateVar['intEnd'],
                    timestep=self.var.modflow_timestep,
                    specific_storage=self.s_stor,
                    specific_yield=self.s_yield,
                    nlay=nlay,
                    nrow=self.domain['nrow'],
                    ncol=self.domain['ncol'],
                    rowsize=self.domain['rowsize'],
                    colsize=self.domain['colsize'],
                    top=self.layer_boundaries[0],
                    bottom=self.layer_boundaries[1:],
                    basin=self.modflow_basin,
                    confined_only=self.confinedAquifer_flags,                                                              
                    head=head,
                    topography=self.layer_boundaries[0],
                    permeability=self.permeability,
                    permeability_vertical=self.permeability_v,
                    load_from_disk=returnBool('load_modflow_from_disk'),
                    setpumpings=False,
                    pumpingloc=None,
                    verbose=verboseGW,
                    complex_solver=self.var.use_complex_solver_for_modflow,
                    super_complex_solver=self.var.use_super_complex_solver_for_modflow)

            # initializing arrays
            self.var.capillar = globals.inZero.copy()
            self.var.baseflow = globals.inZero.copy()
            self.var.depth = globals.inZero.copy()

            self.modflow_watertable = np.copy(head)  # water table will be also saved at modflow resolution

            # sumed up groundwater recharge for the number of days
            self.var.sumed_sum_gwRecharge = globals.inZero
            self.var.modflow_compteur = 0  # Usefull ?

            self.var.modflow_watertable = np.copy(head)  # water table will be also saved at modflow resolution

            # initial water table map is converting into CWatM map
            self.var.head = compressArray(self.modflow2CWATM(head[0]))
            self.var.head = np.array([self.var.head] * nlay)
            for lyr in range(nlay)[1:]:
                self.var.head[lyr,:] = compressArray(self.modflow2CWATM(head[lyr]))

            self.var.writeerror = False
            if 'writeModflowError' in binding:
                self.var.writeerror = returnBool('writeModflowError')
            if self.var.writeerror:
                # This one is to check model's water balance between ModFlow and CwatM exchanges ModFlow discrepancy
                # for each time step can be extracted from the listing file (.lst file) at the end of the simulation
                # as well as the actual pumping rate applied in ModFlow (ModFlow automatically reduces the pumping
                # rate once the ModFlow cell is almost saturated)
                if verboseGW:
                    print('=> ModFlow-CwatM water balance is checked\nModFlow discrepancy for each time step can be extracted from the listing file (.lst file) at the end of the simulation,\nas well as the actual pumping rate applied in ModFlow (ModFlow automatically reduces the pumping rate once the ModFlow cell is almost saturated)')
                self.var.modflowdiscrepancy = np.zeros(
                    dateVar['intEnd'])  # Will contain percentage discrepancy error for ModFLow simulation
                self.var.modflow_cell_area = self.domain['rowsize'] * self.domain['colsize']  # in m2
            else:
                if verboseGW:
                    print(
                    '=> ModFlow-CwatM water balance is not checked\nModFlow discrepancy for each time step can be extracted from the listing file (.lst file) at the end of the simulation,\nas well as the actual pumping rate applied in ModFlow (ModFlow automatically reduces the pumping rate once the ModFlow cell is almost saturated)')

            # then, we got the initial groundwater storage map at ModFlow resolution (in meter)
            self.groundwater_storage_n_layer = head.copy()
            for lyr in range(nlay):   
                
                '''
                 Calculate storage as the flow from storage if head was to drop to zero. 
                 Following 'Documentaton for the MODFLOW 6 Groundwater Flow Model | Ch. 5 of Section A, Groundwater, Book6, Modeling Techniques'
                 https://pubs.usgs.gov/tm/06/a55/tm6a55.pdf
                 
                 Q_sto =  Q_ss + Q_sy
                 
                 Q_ss = SS * A * (TOP - BOT) * (SF * ht) : SFt+1 * ht+1 = 0
                 Q_sy = SY * A * (TOP - BOT) * (SF)  : SFt+1 = 0    
                 
                 Q_sto = [A * (TOP -BOT) * SF] * [SS * ht + Sy]
                 
                 Whereas:
                 Q - flow of water from storage in m^3
                 A - grid cell area
                 SS\SY - specific storage/specific yield
                 TOP/BOT - top/bottom of the aquifer in meters
                 SF - Saturation fraction as calculated by: self.calcSaturatedCellFraction(lyr = lyr, head = head)
                 ht - head
                 
                 Here we calculate storage in meters so we do not account for the A (grid cell area). So:
                 Q_sto = [(TOP -BOT) * SF] * [SS * ht + Sy]                 
                 
                '''
                
                satFrac = self.calcSaturatedCellFraction(lyr = lyr, head = head)
                self.groundwater_storage_n_layer[lyr] = (self.layer_boundaries[lyr] - self.layer_boundaries[lyr + 1]) * satFrac * (self.s_stor[lyr] * head[lyr] + self.s_yield[lyr] * (self.confinedAquifer_flags[lyr] > 0))
            
            self.var.groundwater_storage_total = compressArray(self.modflow2CWATM(np.nansum(self.groundwater_storage_n_layer, axis = 0)))  
            #self.var.groundwater_storage_top_layer = (self.layer_boundaries[lyr] - self.layer_boundaries[lyr + 1]) *\
            #    satFrac * (self.s_stor[lyr] * head[lyr] + self.s_yield[lyr] * (self.confinedAquifer_flags[lyr] > 0))             
             
            # actual pumping output - Zero if no pumping
            self.var.modfPumpingM_actual = globals.inZero.copy()
           
            # calculate groundwater storage available for CWATM
            
            if self.var.GW_pumping:
                self.gwavailable_n_lyrs = head.copy()
                for lyr in range(nlay):
                    satFrac = self.calcSaturatedCellFraction(lyr = lyr, head = head)
                    satFrac_min = (self.var.availableGWStorageFraction * self.modflow.basin)[lyr]
                    head_min = self.layer_boundaries[lyr + 1] + satFrac_min * (self.layer_boundaries[lyr] - self.layer_boundaries[lyr + 1])
                    self.gwavailable_n_lyrs[lyr] =  (self.layer_boundaries[lyr] - self.layer_boundaries[lyr + 1]) * ((self.s_stor[lyr] * np.maximum(head[lyr] * satFrac - head_min * satFrac_min ,0)) + (self.s_yield[lyr] * np.maximum(satFrac - satFrac_min, 0)) * (self.confinedAquifer_flags[lyr] > 0))
                    self.gwavailable_n_lyrs[lyr] = np.where(self.gwavailable_n_lyrs[lyr]  < 0, 0., self.gwavailable_n_lyrs[lyr])
                self.var.groundwater_storage_available = compressArray(self.modflow2CWATM(np.nansum(self.gwavailable_n_lyrs  * self.wells_mask, axis = 0)))  # used in water demand module then
                self.groundwater_storage_available = np.nansum(self.gwavailable_n_lyrs * self.wells_mask, axis = 0)
                self.var.gwstorage_full = self.var.groundwater_storage_available.copy()
            else:
                self.var.groundwater_storage_available = self.var.groundwater_storage_total.copy()
                self.groundwater_storage_available = np.nansum(self.groundwater_storage_n_layer, axis = 0).copy()
                self.var.gwstorage_full = self.var.groundwater_storage_available.copy()
            # self.modflowGroupByCWATM(self.groundwater_storage_available) --> SEE IF CAN BE FIXED
            self.gwAvail_weights = np.minimum(divideArrays(self.groundwater_storage_available, self.CWATM2modflow(decompress(self.var.groundwater_storage_available))), 1.0)
            
                                                                        

            '''
            self.var.groundwater_storage_top_layer = (np.where(head > self.layer_boundaries[0],
                                                               self.layer_boundaries[0], head) - self.layer_boundaries[
                                                          1]) * self.porosity[0]
            
            
            # converting the groundwater storage from ModFlow to CWatM map (in meter)
            self.var.groundwater_storage_available = compressArray(
                self.modflow2CWATM(self.var.groundwater_storage_top_layer))  # used in water demand module then
            self.var.gwstorage_full = compressArray(self.modflow2CWATM(
                (self.layer_boundaries[0] - self.layer_boundaries[1]) * self.porosity[
                    0]))  # use in water damnd to limit pumping
            '''

            # permeability need to be translated into CWatM map to compute leakage from surface water bodies
            self.var.permeability = compressArray(self.modflow2CWATM(self.permeability_v[0])) * self.coefficient

        else:
            ii = 1
            # print('=> ModFlow coupling is not used')

    def dynamic(self):

        # Sumed recharge is re-initialized here for water budget computing purpose
        self.var.modfPumpingM_actual = globals.inZero.copy()  # compressArray(self.modflow2CWATM(self.permeability[0]*0))

        if self.var.modflow_timestep == 1 or ((dateVar['curr'] - int(dateVar[
                                                                         'curr'] / self.var.modflow_timestep) * self.var.modflow_timestep) == 1):  # if it is the first step of the week,months...
            # ? Can the above line be replaced with: (dateVar['curr']  % self.var.modflow_timestep) == 0:
            # setting sumed up recharge again to 7 (or 14 or 30...), will be sumed up for the following 7 days (or 14 or 30...)
            self.var.sumed_sum_gwRecharge = 0
            
        # Adding recharge of the day to the weekly (or bi-weekly, or monthly...) sum
        self.var.sumed_sum_gwRecharge = self.var.sumed_sum_gwRecharge + self.var.sum_gwRecharge
       
        
        # Every modflow timestep (e.g. 7,14,30... days)
        # print('dateVarcurr', dateVar['curr'])
        if dateVar['curr'] == 1 or (dateVar['curr'] % self.var.modflow_timestep) == 0:

            self.var.modflow_compteur += 1

            # converting the CWatM recharge into ModFlow recharge (in meter)
            # we avoid recharge on saturated ModFlow cells, thus CWatM recharge is concentrated  on unsaturated cells
            corrected_recharge = np.where(self.var.capriseindex == 1, 0,
                                          self.var.sumed_sum_gwRecharge / (1 - self.var.capriseindex))
            groundwater_recharge_modflow = self.CWATM2modflow(decompress(corrected_recharge, nanvalue=0),
                                                              correct_boundary=False)
            groundwater_recharge_modflow = np.where(self.var.modflow_watertable - self.layer_boundaries[0] >= 0,
                                                    0, groundwater_recharge_modflow)
            
            nlay_dyn = self.modflow.basin.shape[0]

            # recharge only for the top layer            
            zero_recharge = np.where(self.modflow.basin == 0, 0, 0)
            groundwater_recharge_modflow[1:, :, :] = zero_recharge[1:, :, :] 
            dryModflowCells = (self.modflow.decompress(self.modflow.head.astype(np.float32))  - self.layer_boundaries[1:, :, :]) < 0

            # Correct treatment of multiple layers 
            groundwater_recharge_modflow = np.where(np.nansum(dryModflowCells, axis = 0) == nlay_dyn, 0, groundwater_recharge_modflow)
            # give the information to ModFlow
            self.modflow.set_recharge(groundwater_recharge_modflow)
            
            
            ## ACTUAL RECHARGE
            
            actual_recharge = self.modflow.recharge / ( self.domain['rowsize'] * self.domain['colsize'] )
            actual_recharge_modflow_array = np.zeros((nlay_dyn, self.domain['nrow'], self.domain['ncol']), dtype=np.float32)
                                                                                            

            # modified for layers
            active_cells = np.nansum(self.modflow.basin) / nlay_dyn
            lyr = 0
            # apply separately for each layer 
            wghts = actual_recharge[int(lyr * active_cells):int((lyr + 1) * active_cells)] # active_cells - 1 -> active_cells
            rch_index = self.var.rch_index[int(lyr * active_cells):int((lyr + 1) * active_cells)] # active_cells - 1 -> active_cells
                                                                           

            actual_recharge_modflow_array[lyr, :, :] = np.bincount(rch_index, weights=wghts,
                                        minlength=int(self.modflow.nrow * self.modflow.ncol)).reshape((self.modflow.nrow, self.modflow.ncol))
    
            self.var.sum_gwRecharge_actualM = compressArray(self.modflow2CWATM(actual_recharge_modflow_array[0]))
            
            sum_gwRechargeOld = self.var.sum_prefFlow + self.var.sum_perc3toGW + self.var.pitLatrinToGW + self.var.riverbedExchangeM + self.var.leakageIntoGw + self.var.lakebedExchangeM 
            self.var.prefFlow_GW = divideValues(self.var.sum_prefFlow, sum_gwRechargeOld) * self.var.sum_gwRecharge_actualM
            self.var.perc3toGW_GW = divideValues(self.var.sum_perc3toGW, sum_gwRechargeOld) * self.var.sum_gwRecharge_actualM
            self.var.pitLartinToRunoff = self.var.pitLatrinToGW - divideValues(self.var.pitLatrinToGW, sum_gwRechargeOld) * self.var.sum_gwRecharge_actualM
            self.var.pitLatrinToGW = divideValues(self.var.pitLatrinToGW, sum_gwRechargeOld) * self.var.sum_gwRecharge_actualM
            self.var.riverbedExchangeM = divideValues(self.var.riverbedExchangeM, sum_gwRechargeOld) * self.var.sum_gwRecharge_actualM
            self.var.leakageIntoGw = divideValues(self.var.leakageIntoGw, sum_gwRechargeOld) * self.var.sum_gwRecharge_actualM
            self.var.lakebedExchangeM = divideValues(self.var.lakebedExchangeM, sum_gwRechargeOld) * self.var.sum_gwRecharge_actualM
            

            # update interflow an add rejected_recharge
            rejected_recharge = np.maximum(self.var.sum_gwRecharge - self.var.sum_gwRecharge_actualM, 0.)
            # share of landcover specific rejected runoff
            lc_rechargeShare = divideArrays(self.var.gwRecharge, self.var.sum_gwRecharge) * self.var.fracVegCover
            # update landcover specific interflow
            self.var.interflow += divideArrays(lc_rechargeShare * rejected_recharge, self.var.fracVegCover)
            #update collective gridcell interflow
            self.var.sum_interflow += np.nansum(lc_rechargeShare * rejected_recharge, axis = 0)
            
            # calculate rejected recharge that was not returned to interflow (i.e., from bed exchange - should be sent to runoff)
            self.var.bedExchangeLeakToRunoff = rejected_recharge - np.nansum(lc_rechargeShare * rejected_recharge, axis = 0)
            self.var.sum_directRunoff += self.var.bedExchangeLeakToRunoff 

            self.var.sum_runoff = self.var.sum_directRunoff + self.var.sum_interflow + self.var.leakageIntoRunoff
            
            # update cwatm sum_geRecharge
            self.var.sum_gwRecharge -= rejected_recharge
            
             # update sum_soilFlows
            self.var.sum_prefFlow = self.var.prefFlow_GW.copy()
            self.var.sum_perc3toGW = self.var.perc3toGW_GW.copy()

            
            
            #print(np.nansum(self.var.sum_gwRecharge) / np.nansum(self.var.prefFlow_GW  + self.var.perc3toGW_GW +  self.var.pitLatrinToGW + self.var.riverbedExchangeM + self.var.leakageIntoGw + self.var.lakebedExchangeM))
            
            ## INSTALLING WELLS
            if self.var.GW_pumping:
                ## Groundwater demand from CWatM installs wells in each Modflow cell

                # Groundwater pumping demand from the CWatM waterdemand module, will be decompressed to 2D array
                # CWatM 2D groundwater pumping array is converted into Modflow 2D array
                # Pumping is given to ModFlow in m3 and < 0
                # print('mean modflow pumping [m]: ', np.nanmean(self.var.modfPumpingM))
                groundwater_abstraction = - self.CWATM2modflow(decompress(self.var.Pumping_daily)) * domain['rowsize'] * \
                                          domain['colsize']
                #groundwater_abstraction = - self.CWATM2modflow(decompress(self.var.modfPumpingM)) * domain['rowsize'] * \
                #                          domain['colsize']
                
                #wellsMask = (np.nansum(self.wells_mask, axis = 0) > 0)
                wellsMask = 1
                groundwater_abstraction2 = np.array([groundwater_abstraction  * wellsMask] * nlay_dyn) 
            
                #self.gwAvail_weights = np.minimum(divideArrays(self.groundwater_storage_available, self.CWATM2modflow(decompress(self.var.groundwater_storage_available))), 1.0)
               
               
                #groundwater_abstraction2 = groundwater_abstraction2 * self.gwAvail_weights 

                if self.var.correctPumpingDiscrepancy:
                    '''
                    Correct pumping discrepancy between MODFLOW and CWatM. Due to conversion errors, and inaccurate storage estiamtes - it often occurs that CWatM allocated 
                    more groundwater for consumption relative to the volume that is actually being pumped. So np.nansum(self.var.modfPumpingM) > np.nansum(self.var.modfPumpingM_actual)
                
                    The correction suggested below, inflate requested pumping in all grid cells by a constant, calculated as CWatM pumping request divided by Modflow adjusted pumping request. 
                    The resulting pumping reduces the discrepancy between the two, but causes a spatial mismtach between demand and pumping locations. Therefore, this practice is discourage in 
                    coarse resolution appliactions.
                
                    User can set it on by setting correctPumpingDiscrepancy = True in the Settings file.
                    '''
                    #print(np.nansum(groundwater_abstraction))
                    #print(np.nansum(groundwater_abstraction2))
                    ### Huge differencve between gw_requested and groundwater_abstraction2 - when does it happen? why does it happen? is there prior restriction on gw_request from modflow to cwatm?
                    #print(np.nansum(groundwater_abstraction)/ np.nansum(groundwater_abstraction2))
                    groundwater_abstraction2 =  groundwater_abstraction2 * np.nansum(groundwater_abstraction)/ np.nansum(groundwater_abstraction2)
                self.var.modfPumpingM = globals.inZero.copy()
                #calculate available water per layer and split abstraction between layers. 
                StorageInAbstractionCells = self.gwavailable_n_lyrs * wellsMask  
            
                proportional_groundwater_storage = StorageInAbstractionCells / np.nansum(StorageInAbstractionCells, axis = 0)
                proportional_groundwater_storage = np.where(np.isnan(proportional_groundwater_storage), 1, proportional_groundwater_storage)
                groundwater_abstraction2 = groundwater_abstraction2 * proportional_groundwater_storage
                # * 100 because applied only in one ModFlow cell in Burgenland
                # give the information to ModFlow
                self.modflow.set_groundwater_abstraction(groundwater_abstraction2)

            # running ModFlow
            self.modflow.step()
            # self.modflow.finalize()
            # extracting the new simulated hydraulic head map
            head = self.modflow.decompress(self.modflow.head.astype(np.float32))
            # self.var.head = compressArray(self.modflow2CWATM(head))
            # MODIF LUCA
            if self.var.writeerror:
                # copying the previous groundwater storage at ModFlow resolution (in meter)
                groundwater_storage_n_layer0 = np.nansum(self.groundwater_storage_n_layer, axis = 0)  # for water balance computation; groundwater_storage_n_layer[0] -> groundwater_storage_top_layer
        
            # MODIFIED DOR FRIDMAN
            # computing the new groundwater storage map at ModFlow resolution (in meter)
            self.groundwater_storage_n_layer = head.copy()
            for lyr in range(nlay_dyn): 
                satFrac = self.calcSaturatedCellFraction(lyr = lyr, head = head)
                self.groundwater_storage_n_layer[lyr] = (self.layer_boundaries[lyr] - self.layer_boundaries[lyr + 1]) * satFrac * (self.s_stor[lyr] * head[lyr] + self.s_yield[lyr] * (self.confinedAquifer_flags[lyr] > 0))
            self.var.groundwater_storage_total = compressArray(self.modflow2CWATM(np.nansum(self.groundwater_storage_n_layer, axis = 0)))  

            if self.var.GW_pumping:
                self.gwavailable_n_lyrs = head.copy()
                for lyr in range(nlay_dyn):
                    satFrac = self.calcSaturatedCellFraction(lyr = lyr, head = head)
                    satFrac_min = self.var.availableGWStorageFraction
                    head_min = self.layer_boundaries[lyr + 1] + satFrac_min * (self.layer_boundaries[lyr] - self.layer_boundaries[lyr + 1])
                    self.gwavailable_n_lyrs[lyr] =  (self.layer_boundaries[lyr] - self.layer_boundaries[lyr + 1]) * ((self.s_stor[lyr] * np.maximum(head[lyr] * satFrac - head_min * satFrac_min ,0)) + (self.s_yield[lyr] * np.maximum(satFrac - satFrac_min, 0)) * (self.confinedAquifer_flags[lyr] > 0))
                    self.gwavailable_n_lyrs[lyr] = np.where(self.gwavailable_n_lyrs[lyr] < 0, 0. , self.gwavailable_n_lyrs[lyr])
                
                
                wellsMask = (np.nansum(self.wells_mask, axis = 0) > 0)
            
                # calculate groundwater storage available for MODFLOW & gwAvailable allocation weights
                self.groundwater_storage_available = np.nansum(self.gwavailable_n_lyrs  * wellsMask, axis = 0)
                self.var.groundwater_storage_available = compressArray(self.modflow2CWATM(np.nansum(self.gwavailable_n_lyrs  * wellsMask, axis = 0)))  # used in water demand module then
            else:
                self.var.groundwater_storage_available = self.var.groundwater_storage_total.copy()
                self.groundwater_storage_available = np.nansum(self.groundwater_storage_n_layer, axis = 0).copy()
            
            
        
            self.var.gwstorage_full = self.var.groundwater_storage_available.copy()
            assert self.permeability.ndim == 3
           

            # computing the groundwater outflow by re-computing water outflowing the aquifer through the DRAIN ModFlow package
            
    
            groundwater_outflow = np.where(head - self.layer_boundaries[0] >= 0,
                                           (head - self.layer_boundaries[0]) * self.coefficient * self.permeability[0],
                                           0)  # in m/day per ModFlow cell
            groundwater_outflow2 = np.where(head - self.layer_boundaries[0] >= 0, 1.0,
                                            0.0)  # For the next step, it will prevent recharge where ModFlow cells are saturated, even if there is no capillary rise (where h==topo)
             
            # capillary rise and baseflow from groundwater are allocated in function of the river percentage of each ModFlow cell
            capillar = groundwater_outflow * (1 - self.var.channel_ratio)  # We are still in ModFlow coordinate
            baseflow = groundwater_outflow * self.var.channel_ratio  # We are still in ModFlow coordinate
            
            if self.var.GW_pumping:

                    # extracting actual ModFlow pumping
                    actual_pumping = self.modflow.actualwell_rate.astype(np.float32)
                    # the size of actual pumping corresponds to the number of wells is masked array 'wellsloc'
                    actual_pumping_modflow_array = np.bincount(self.var.wells_index, weights=actual_pumping,
                                                               minlength=int(
                                                                   self.modflow.nrow * self.modflow.ncol)).reshape(
                        (self.modflow.nrow, self.modflow.ncol))
                    
                    self.var.modfPumpingM_actual = - compressArray(
                        self.modflow2CWATM(actual_pumping_modflow_array) / (domain['rowsize'] * domain['colsize']))
             
            if self.var.writeerror:
                # Check the water balance in the aquifer: recharge = capillary + baseflow + storage change (- pumping)

                if self.var.GW_pumping:

                    mid_gwflow = np.nansum(((groundwater_recharge_modflow + (
                                capillar + baseflow) * self.var.modflow_timestep) * self.var.modflow_cell_area - actual_pumping_modflow_array) / 2)
                    # in m3/timestep, for water balance computation
                    Budget_ModFlow_error = np.round(100 * (np.nansum(
                        (groundwater_recharge_modflow - (capillar + baseflow) * self.var.modflow_timestep -
                         (
                                     self.var.groundwater_storage_top_layer - groundwater_storage_top_layer0)) * self.var.modflow_cell_area + actual_pumping_modflow_array) / mid_gwflow),
                                                    2)
                    sumModFlowout = np.nansum((capillar + baseflow) * self.var.modflow_timestep * self.var.modflow_cell_area - groundwater_abstraction)  # m3/timestep, for water balance computation

                else:

                    mid_gwflow = np.nansum((groundwater_recharge_modflow + (
                                capillar + baseflow) * self.var.modflow_timestep) * self.var.modflow_cell_area / 2)  # in m3/timestep, for water balance computation
                    Budget_ModFlow_error = np.round(100 * (np.nansum(
                        (groundwater_recharge_modflow - (capillar + baseflow) * self.var.modflow_timestep -
                         (
                                     self.var.groundwater_storage_top_layer - groundwater_storage_top_layer0)) * self.var.modflow_cell_area) / mid_gwflow),
                                                    2)

                # to check CWatM-ModFlow projection is working well
                    sumModFlowout = np.nansum((capillar + baseflow) * self.var.modflow_timestep * self.var.modflow_cell_area)  # m3/timestep, for water balance computation

                #  saving modflow discrepancy, it will be written a text file at the end of the simulation
                self.var.modflowdiscrepancy = Budget_ModFlow_error
                # print('ModFlow discrepancy : ', Budget_ModFlow_error, ' % (if pumping, it considers pumping demand is satisfied)')

                # to check CWatM-ModFlow projection is working well
                sumModFlowin = np.nansum(groundwater_recharge_modflow * self.var.modflow_cell_area)  # m3, for water balance computation

            # converting flows from ModFlow to CWatM domain
            
            self.var.capillar = compressArray(self.modflow2CWATM(capillar[0,:,:]))
            baseflow = baseflow.astype('float64')  # required for runoff concentration (lib2.runoffConc)
            self.var.baseflow = compressArray(self.modflow2CWATM(baseflow[0,:,:]))                  
            
            # computing saturated fraction of each CWatM cells (where water table >= soil bottom)
            self.var.capriseindex = compressArray(self.modflow2CWATMbis(
                groundwater_outflow2))  # initialized in landcoverType module, self.var.capriseindex is the fraction of saturated ModFlow cells in each CWatM cell

            # head = self.modflow.decompress(self.modflow.head.astype(np.float32))
            # updating water table maps both at CWatM and ModFlow resolution
            for lyr in range(nlay_dyn):
                setattr(self.var, 'head' + str(lyr), compressArray(self.modflow2CWATM(head[lyr])))
                setattr(self.var, 'gwdepth' + str(lyr), compressArray(self.modflow2CWATM(self.layer_boundaries[lyr] - head[lyr])))

          
            if 'gw_depth_observations' in binding:
                self.var.gwdepth_difference_sim_obs = self.var.gwdepth - self.var.gwdepth_observations
            if 'gw_depth_sim_obs' in binding:
                self.var.gwdepth_adjusted = np.maximum(self.var.gwdepth - self.var.gwdepth_adjuster, 0)
                head_adjusted = self.layer_boundaries[0] - self.CWATM2modflow(decompress(self.var.gwdepth_adjusted))

                self.var.modflow_head_adjusted = np.copy(head_adjusted)
            
            self.var.modflow_watertable = np.copy(head)

            # Re-initializing the weekly (orbi-weekly, or monthly...) sum of groundwater pumping
            self.var.modfPumpingM = globals.inZero.copy()

            # to check CWatM-ModFlow projection is working well
            if self.var.writeerror:
                # computing the total recharge rate provided by CWatM for this time step (in m3/timestep)
                total_recharge_m3_cwatm = (self.var.sumed_sum_gwRecharge * self.var.cellArea).sum()  # m3/timestep, for water balance computation
           
                print('ModFlow-CWatM input conversion error : ', np.round(100 * (total_recharge_m3_cwatm - sumModFlowin)/sumModFlowin, 2), ' %')
                if self.var.GW_pumping:
                    print('ModFlow-CWatM output conversion error : ', np.round(100 * (np.nansum((self.var.capillar+self.var.baseflow) * self.var.modflow_timestep *
                                                                                         self.var.cellArea) - np.nansum(groundwater_abstraction) - sumModFlowout) / sumModFlowout, 2), ' %')
                    # Check crossed models:
                    print('ModFlow discrepancy crossed: ', np.round(100 * (total_recharge_m3_cwatm -
                                                                            np.nansum(((capillar + baseflow) * self.var.modflow_timestep +
                                                                                      (self.var.groundwater_storage_top_layer - groundwater_storage_top_layer0)) *
                                                                                     self.var.modflow_cell_area - groundwater_abstraction)) / mid_gwflow, 2), ' %')
            
                else:
                    print('ModFlow-CWatM output conversion error : ', np.round(100 * (np.nansum((self.var.capillar + self.var.baseflow) * self.var.modflow_timestep *
                                    self.var.cellArea) - sumModFlowout) / sumModFlowout, 2), ' %')
                    print('ModFlow discrepancy crossed: ', np.round(100 * (total_recharge_m3_cwatm -
                                                                    np.nansum(((capillar + baseflow) * self.var.modflow_timestep + (
                                                                                                 self.var.groundwater_storage_top_layer -
                                                                                                 groundwater_storage_top_layer0)) * self.var.modflow_cell_area)) / mid_gwflow, 2), ' %')

            # Writing ModFlow discrepancies at the end of simulation:
            if dateVar['currDate'] == dateVar['dateEnd']:
                if self.var.writeerror:
                    discrep_filename = cbinding('PathOut') + '/' + 'ModFlow_DiscrepancyError.txt'
                    discrep_file = open(discrep_filename, "w")
                    sum_modflow_errors = 0
                    threeshold_modflow_error = 0.01  # in percentage

                    for tt in range(dateVar['intEnd']):
                        if self.var.modflowdiscrepancy > threeshold_modflow_error:  # if error is higer than threeshold in %, we print it.
                            discrep_file.write(
                                "ModFlow stress period " + str(tt + 1) + " : percentage error in ModFlow is ")
                            discrep_file.write(str(self.var.modflowdiscrepancy))
                            discrep_file.write("\n")
                            sum_modflow_errors += 1
                    if sum_modflow_errors == 0:
                        discrep_file.write(
                            "ModFlow error was always below " + str(
                                threeshold_modflow_error) + ' % during the simulation')
                    discrep_file.close()
