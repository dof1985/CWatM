# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Name:        Water quality - erosion-sediment module (EROSED)
# Purpose:     simulate total suspended solids (TSS) in rivers
#
# Author:      TT, FS, PB, MS, DF
#
# Created:     20/01/2022
# Copyright:   (c) TT, FS, PB, MS, DF 2022
# -------------------------------------------------------------------------


from cwatm.management_modules.data_handling import *
from cwatm.management_modules.globals import *

class waterquality_erosed(object):
    """
        WATER QUALITY - EROSION-SEDIMENT MODULE TREATMENT

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
        """
                INITIAL PART OF THE EROSED MODULE
        """
        print('erosed init')

    def dynamic(self):
        """
        Dynamic part of EROSED module
        """
        print("erosed dyn")