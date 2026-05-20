import os
import sys
import numpy as np
import importlib
from cwatm.management_modules.data_handling import *

# =====================================================================
# MODFLOW PURE STEADY-STATE SUBCLASS (NO PUMPING)
# =====================================================================
class ModFlowSteadySimulation:
    """
    A dedicated standalone constructor class that forces Flopy 
    to compile and solve a pure, un-stressed Steady-State system.
    """
    def __init__(
        self, name, folder, path_mf6dll, specific_storage, specific_yield,
        nlay, nrow, ncol, rowsize, colsize, top, bottom, basin, confined_only, 
        head, topography, permeability, permeability_vertical, ndays=1.0, verbose=False
    ):
        flopy = importlib.import_module("flopy", package=None)

        self.name = name.upper()
        self.folder = folder
        self.dir_mf6dll = path_mf6dll
        self.nrow = nrow
        self.ncol = ncol
        self.rowsize = rowsize
        self.colsize = colsize
        self.basin = basin
        self.verbose = verbose
        self.working_directory = os.path.join(folder, 'wd_steady')
        
        if not os.path.exists(self.working_directory):
            os.makedirs(self.working_directory)

        # Build clear absolute path execution string to avoid lookup errors
        exe_path = os.path.join(path_mf6dll, 'mf6')
        if sys.platform.startswith("win") and not exe_path.endswith(".exe"):
            exe_path += ".exe"

        if self.verbose:
            print("Compiling Standalone Pre-Development Steady-State MODFLOW Grid Structure...")
            print(f"--> Target MODFLOW-6 Binary Verified At: {exe_path}")

        sim = flopy.mf6.MFSimulation(
            sim_name=self.name, version='mf6', exe_name=exe_path,
            sim_ws=self.working_directory, memory_print_option='all'
        )
        
        flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(float(ndays), 1, 1.0)])

        flopy.mf6.ModflowIms(
            sim, print_option=None, complexity='COMPLEX', linear_acceleration='BICGSTAB',
            under_relaxation='SIMPLE', under_relaxation_gamma=0.1, backtracking_number=5,
            backtracking_tolerance=10 ** 5, backtracking_reduction_factor=0.3, backtracking_residual_limit=150
        )

        gwf = flopy.mf6.ModflowGwf(sim, modelname=self.name, newtonoptions='newton under_relaxation', print_input=False, print_flows=False)

        flopy.mf6.ModflowGwfdis(
            gwf, nlay=nlay, nrow=self.nrow, ncol=self.ncol, delr=self.rowsize, delc=self.colsize,
            top=top, botm=bottom, idomain=self.basin, nogrb=True
        )

        flopy.mf6.ModflowGwfic(gwf, strt=head)
        flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=confined_only, k=permeability)
        flopy.mf6.ModflowGwfoc(gwf, head_filerecord=f'{self.name}.hds', saverecord=[('HEAD', 'ALL')])
 
        flopy.mf6.ModflowGwfsto(
            gwf, save_flows=False, iconvert=confined_only, ss=specific_storage, sy=specific_yield,
            steady_state=True, transient=False
        )

        self.recharge_cells = np.zeros((self.basin.sum(), 4), dtype=np.int32)
        recharge_locations = np.where(self.basin == True)
        self.recharge_cells[:, 0] = recharge_locations[0]
        self.recharge_cells[:, 1] = recharge_locations[1]
        self.recharge_cells[:, 2] = recharge_locations[2]
        
        self.rch_package = flopy.mf6.ModflowGwfrch(
            gwf, fixed_cell=False, print_input=False, print_flows=False, save_flows=False,
            boundnames=None, maxbound=self.basin.sum(), stress_period_data=self.recharge_cells.tolist()
        )

        drainage_cells = np.zeros((self.basin[0].sum(), 5))
        drn_locations = np.where(self.basin[0] == True)
        drainage_cells[:, 1] = drn_locations[0]
        drainage_cells[:, 2] = drn_locations[1]
        drainage_cells[:, 3] = topography[drn_locations] 
        drainage_cells[:, 4] = permeability[0][self.basin[0] == True] * self.rowsize * self.colsize
        
        drainage_list = [[int(i), int(j), int(k), l, m] for i, j, k, l, m in drainage_cells.tolist()]
        flopy.mf6.ModflowGwfdrn(
            gwf, maxbound=self.basin[0].sum(), stress_period_data=drainage_list,
            print_input=False, print_flows=False, save_flows=False
        )
        
        self.sim = sim

    def set_steady_recharge(self, recharge_forcing_array):
        """Maps computed 3D long term recharge directly onto Flopy stress period fields"""
        recharge_list = self.recharge_cells.copy().astype(object)
        for idx in range(len(recharge_list)):
            lyr, row, col = recharge_list[idx, 0], recharge_list[idx, 1], recharge_list[idx, 2]
            recharge_list[idx, 3] = recharge_forcing_array[lyr, row, col] * (self.rowsize * self.colsize)
            
        self.rch_package.stress_period_data.set_data(recharge_list.tolist(), key=0)

    def execute_and_solve(self):
        self.sim.write_simulation()
        self.sim.run_simulation()


# =====================================================================
# CWATM WRAPPER PIPELINE ROUTINE
# =====================================================================
def run_standalone_steady_state(
    transient_instance, folder_out, path_mf6dll, nlay, nrow, ncol,
    rowsize, colsize, layer_top, layer_bottom, basin_mask, confined_flags,
    initial_head, ss_storage, sy_yield, verbose_gw=False
):
    """
    Bypasses transient execution loops. Consumes prepared variables directly from 
    the transient module setup parameters to resolve spatial groundwater equilibrium matrices.
    """
    print("\n" + "="*60)
    print("  CWATM-MODFLOW-6 INTERCEPT: RUNNING INITIAL STEADY-STATE")
    print("="*60 + "\n")

    # Access calibration window tracking footprint length if present
    ndays_steady = 1.0
    if 'Ndays_steady' in binding:
        ndays_steady = float(loadmap('Ndays_steady'))

    # Initialize pure steady state execution wrapper using passed values
    model_runner = ModFlowSteadySimulation(
        name='steady', folder=folder_out, path_mf6dll=path_mf6dll,
        specific_storage=ss_storage, specific_yield=sy_yield, nlay=nlay, nrow=nrow, ncol=ncol,
        rowsize=rowsize, colsize=colsize, top=layer_top, bottom=layer_bottom, basin=basin_mask,
        confined_only=confined_flags, head=initial_head, topography=layer_top, 
        permeability=transient_instance.permeability, permeability_vertical=transient_instance.permeability_v,
        ndays=ndays_steady, verbose=verbose_gw
    )

    # Core Exogenous Recharge Forcing Loader
    if 'longterm_recharge' in binding:
        print("--> Extracting user-defined long-term recharge raster...")
        cwatm_rch_map = loadmap('longterm_recharge') 
        modflow_rch_2d = transient_instance.CWATM2modflow(decompress(cwatm_rch_map, nanvalue=0), correct_boundary=False)
    else:
        print("--> WARNING: 'longterm_recharge' missing in bindings. Initializing with flat floor value (0.001 m/day).")
        modflow_rch_2d = np.full((nrow, ncol), 0.001, dtype=np.float32)

    longterm_recharge_3d = np.zeros((nlay, nrow, ncol), dtype=np.float32)
    longterm_recharge_3d[0, :, :] = np.where(basin_mask[0] == 1, modflow_rch_2d, 0.0)
    
    # Pack array map specifications into model memory space pointers
    model_runner.set_steady_recharge(longterm_recharge_3d)

    print("--> Invoking mathematical matrix updates via MODFLOW-6...")
    model_runner.execute_and_solve()

    # Extract calibration configurations out to runtime context files
    flopy = importlib.import_module("flopy", package=None)
    head_file_obj = flopy.utils.HeadFile(os.path.join(model_runner.working_directory, 'STEADY.hds'))
    stable_heads = head_file_obj.get_data(totim=float(ndays_steady))
    
    export_save_path = os.path.join(folder_out, "steady_state_equilibrium_heads.npy")
    np.save(export_save_path, stable_heads)

    print("\n" + "="*60)
    print("  ✔ SUCCESS: EQUILIBRIUM SOLVED FROM LONG-TERM DRIVERS")
    print(f"  Target File Saved: {export_save_path}")
    print("  Forcing system process shutdown. Preventing CWatM dynamic loops.")
    print("="*60 + "\n")
    
    sys.exit(0)