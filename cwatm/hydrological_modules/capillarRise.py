# -------------------------------------------------------------------------
# Name:        Capillar Rise module
# Purpose: Capillary rise calculation module for groundwater-surface water interaction.
# Determines cell fractions influenced by capillary rise based on groundwater depth.
# Connects groundwater processes with soil moisture through capillary action.
#
# Author:      PB
# Created:     20/07/2016
# CWatM is licensed under GNU GENERAL PUBLIC LICENSE Version 3.
# -------------------------------------------------------------------------

from cwatm.management_modules.data_handling import *


class capillarRise(object):
    """
    Capillary rise module for groundwater-surface interaction.

    This class calculates the cell fraction influenced by capillary rise based on
    groundwater depth and relative elevation within each grid cell. It determines
    areas where groundwater can reach the surface through capillary action.

    Attributes
    ----------
    var : object
        Model variables container
    model : object
        CWatM model instance
    =====================================  ======================================================================  =====
    Variable [self.var]                    Description                                                             Unit 
    =====================================  ======================================================================  =====
    dzRel0100                              relative elevation above flood plains (max elevation above plain)       m    
    dzRel0090                              relative elevation above flood plains (90% elevation above plain)       m    
    dzRel0080                              relative elevation above flood plains (80% elevation above plain)       m    
    dzRel0070                              relative elevation above flood plains (70% elevation above plain)       m    
    dzRel0060                              relative elevation above flood plains (60% elevation above plain)       m    
    dzRel0050                              relative elevation above flood plains (median elevation above plain)    m    
    dzRel0040                              relative elevation above flood plains (40% elevation above plain)       m    
    dzRel0030                              relative elevation above flood plains (30% elevation above plain)       m    
    dzRel0020                              relative elevation above flood plains (20% elevation above plain)       m    
    dzRel0010                              relative elevation above flood plains (10% elevation above plain)       m    
    dzRel0005                              relative elevation above flood plains (5% elevation above plain)        m    
    dzRel0001                              relative elevation above flood plains (1% elevation above plain)        m    
    capRiseFrac                            fraction of a grid cell where capillar rise may happen                  m    
    modflow                                Flag: True if modflow_coupling = True in settings file                  --   
    storGroundwater                        Groundwater storage (non-fossil). This is primarily used when not usin  m    
    specificYield                          groundwater reservoir parameters (if ModFlow is not used) used to comp  m    
    maxGWCapRise                           influence of capillary rise above groundwater level                     m    
    =====================================  ======================================================================  =====

    **Functions**

    """

    def __init__(self, model):
        """
        Initialize the capillary rise module.

        Parameters
        ----------
        model : object
            CWatM model instance containing variables and configuration
        """
        self.var = model.var
        self.model = model

    def dynamic(self):
        """
        Calculate capillary rise dynamics for the current time step.

        Computes the cell fraction influenced by capillary rise based on the
        approximate height of groundwater and relative elevation within each grid cell.
        This determines areas where groundwater contributes to surface processes
        through capillary action.
        """

        if checkOption('CapillarRise') and not (self.var.modflow):

            # approximate height of groundwater table and corresponding reach of cell under influence of capillary rise
            dzGroundwater = self.var.storGroundwater / self.var.specificYield + self.var.maxGWCapRise
            #       10  9   8   7   6   5    4   3   2  1    0
            vvv = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01]
            i = 10
            for vv in vvv:
                if i > 0:
                    h = (vv - (self.var.dzRel[i] - dzGroundwater) * 0.1 / 
                         np.maximum(1e-3, self.var.dzRel[i] - self.var.dzRel[i-1]))
                else:
                    h = (vv - (self.var.dzRel[i] - dzGroundwater) * 0.1 / 
                         np.maximum(1e-3, self.var.dzRel[i]))
                CRFRAC = np.where(dzGroundwater < self.var.dzRel[i], h, CRFRAC)
                i -= 1
                
            self.var.capRiseFrac = np.maximum(0.0, np.minimum(1.0, CRFRAC))
        else:
            self.var.capRiseFrac = 0.
