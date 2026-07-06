import os
import numpy as np
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
import rasterio

def process_plot_and_export(hds_path, elevation_tif_path, output_dir="output_rasters"):
    # 1. Verify files exist
    if not os.path.exists(hds_path):
        print(f"Error: The head file '{hds_path}' does not exist.")
        return
    if not os.path.exists(elevation_tif_path):
        print(f"Error: The elevation file '{elevation_tif_path}' does not exist.")
        return
        
    os.makedirs(output_dir, exist_ok=True)

    # 2. Read Spatial Profile and data from the Elevation GeoTIFF
    with rasterio.open(elevation_tif_path) as src:
        elevation_array = src.read(1)  # Read band 1
        meta_profile = src.profile.copy()
        elev_nodata = src.nodata

    # 3. Load the binary head file
    try:
        headobj = bf.HeadFile(hds_path, precision='double')
        times = headobj.get_times()
    except Exception:
        headobj = bf.HeadFile(hds_path, precision='single')
        times = headobj.get_times()
    
    latest_time = times[-1] 
    head_array = headobj.get_data(totim=latest_time)
    num_layers = head_array.shape[0]

    # Verify grid dimensions match
    if elevation_array.shape != head_array[0].shape:
        print(f"Warning: Grid dimension mismatch! Elevation is {elevation_array.shape}, Head layer is {head_array[0].shape}")

    # Update metadata profile for our output files
    nodata_value = -9999.0
    meta_profile.update(dtype=rasterio.float32, count=1, nodata=nodata_value)

    # 4. Loop through each layer to compute, export, and plot
    for layer_idx in range(num_layers):
        layer_head = head_array[layer_idx, :, :].astype(np.float32)
        
        # --- MASKING ---
        # Mask out the inactive (1e30) and dry (-1e30 or -888) cells for plotting
        masked_head = np.ma.masked_greater(layer_head, 1e29)
        masked_head = np.ma.masked_less(masked_head, -1e29)
        masked_head = np.ma.masked_equal(masked_head, -888.0) 

        if masked_head.count() == 0:
            print(f"Warning: Layer {layer_idx + 1} contains no active data cells. Skipping.")
            continue

        # Create clean array for GeoTIFF export (assigning standard GIS nodata value)
        export_head = np.where(masked_head.mask, nodata_value, layer_head).astype(np.float32)

        # Calculate SWL (Depth to Water): Elevation - Head
        # Mask out locations where either head or elevation is invalid
        swl_array = np.where(
            (~masked_head.mask) & (elevation_array != elev_nodata),
            elevation_array - layer_head,
            nodata_value
        ).astype(np.float32)
        
        masked_swl = np.ma.masked_equal(swl_array, nodata_value)

        # --- EXPORT TO GEOTIFF ---
        head_out_path = os.path.join(output_dir, f"Layer_{layer_idx + 1}_Head.tif")
        with rasterio.open(head_out_path, 'w', **meta_profile) as dst:
            dst.write(export_head, 1)
            
        swl_out_path = os.path.join(output_dir, f"Layer_{layer_idx + 1}_SWL.tif")
        with rasterio.open(swl_out_path, 'w', **meta_profile) as dst:
            dst.write(swl_array, 1)
            
        print(f"Layer {layer_idx + 1}: Rasters exported successfully.")

        # --- PLOTTING ---
        # Plot 1: Hydraulic Head
        plt.figure(figsize=(8, 6))
        im1 = plt.imshow(masked_head, cmap='viridis')
        plt.colorbar(im1, label='Hydraulic Head (m)')
        plt.title(f"Layer {layer_idx + 1} - Active Hydraulic Heads")
        plt.xlabel("Column Index")
        plt.ylabel("Row Index")
        plt.show()

        # Plot 2: Static Water Level (SWL)
        plt.figure(figsize=(8, 6))
        # 'plasma' or 'YlGnBu_r' work nicely for depth maps
        im2 = plt.imshow(masked_swl, cmap='plasma') 
        plt.colorbar(im2, label='Depth to Water / SWL (m)')
        plt.title(f"Layer {layer_idx + 1} - Static Water Level (SWL)")
        plt.xlabel("Column Index")
        plt.ylabel("Row Index")
        plt.show()

    print(f"\nAll operations complete. Rasters saved to: '{os.path.abspath(output_dir)}'")

# --- Example Usage ---
hds_file = "P:/watmodel/CWATM/model/cwatm_modflow_layers/Toolkit/Calibration_modflow/runs_calibration/Run_165/165/wd_steady/STEADY.hds" 
elevation_file = "P:/watmodel/CWATM/modelruns/uppernile_wmz_G4DR/inputs/modf_out/5000.0m/elevation_modflow.tif"

process_plot_and_export(hds_file, elevation_file)