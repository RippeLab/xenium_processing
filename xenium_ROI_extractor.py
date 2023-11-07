import pandas as pd
import json
import os
import re
import tiffile
import argparse
import h5py

def get_arguments():	
    """
    Parses and checks command line arguments, and provides an help text.
    Assumes 3 and returns 3 positional command line arguments:
    """
    parser = argparse.ArgumentParser(description = "Split transcripts and images according to the given ROIs")
    parser.add_argument("xenium_folder", help = "path to the xenium folder")
    parser.add_argument("roi_folder", help = "path to the folder with the ROI images")   
    parser.add_argument("output_folder", help = "path to the output folder")
    args = parser.parse_args()
    return args.xenium_folder, args.roi_folder, args.output_folder

if __name__ == "__main__":
    xenium_folder, roi_folder, ouptut_folder = get_arguments()
	
    # Expects standard xenium filenames: https://www.10xgenomics.com/support/in-situ-gene-expression/documentation/steps/onboard-analysis/at-a-glance-xenium-output-files
    xenium_file_path = os.path.join(xenium_folder, "experiment.xenium")
    transcript_file_path = os.path.join(xenium_folder, "transcripts.parquet")
    cells_file_path = os.path.join(xenium_folder, "cells.parquet")
    cell_boundaries_file_path = os.path.join(xenium_folder, "cell_boundaries.parquet")
    mip_tiff_file_path = os.path.join(xenium_folder, "morphology_mip.ome.tif")
    print(f"Xenium file: {xenium_file_path}")
    print(f"Transcript file: {transcript_file_path}")
    print(f"Cell file: {cells_file_path}")
    print(f"Cell Boundary file: {cell_boundaries_file_path}")
    print(f"Maximum Projection Image File: {mip_tiff_file_path}")
    
    # Reads the pixel size in um
    with open(xenium_file_path) as metadata_file:
        d = json.load(metadata_file)
        pixel_size = d["pixel_size"]

    # Reads the transcripts and add locations in pixels
    print(f"Reading Transcript file: {transcript_file_path}")
    transcripts = pd.read_parquet(transcript_file_path)
    transcripts["x_location_px"] = transcripts["x_location"] / pixel_size
    transcripts["y_location_px"] = transcripts["y_location"] / pixel_size
    transcripts["z_location_px"] = transcripts["z_location"] / pixel_size
    
    # Reads the cells and add locations in pixels
    print(f"Reading Cell file: {cells_file_path}")
    cells = pd.read_parquet(cells_file_path)
    cells["x_centroid_px"] = cells["x_centroid"] / pixel_size
    cells["y_centroid_px"] = cells["y_centroid"] / pixel_size
    
    # Reads the cells and add locations in pixels
    print(f"Reading Cell Boundary file: {cell_boundaries_file_path}")
    cell_boundaries = pd.read_parquet(cell_boundaries_file_path)
    cell_boundaries["vertex_x_px"] = cell_boundaries["vertex_x"] / pixel_size
    cell_boundaries["vertex_y_px"] = cell_boundaries["vertex_y"] / pixel_size
    
    # Reads ROIS, expects _Xnnn_Ynnnn_Wnnn_Hnnnn_ with n = 0-9
    # makes a lists of lists [x, y, w, h]
    print(f"Looking for ROIs in: {roi_folder}")
    ROIs = os.listdir(roi_folder)
    ROIs = filter(lambda filename: re.match(".*_X[0-9]+_Y[0-9]+_W[0-9]+_H[0-9]+.*", filename), ROIs)
    ROIs = [ re.match(".*_X([0-9]+)_Y([0-9]+)_W([0-9]+)_H([0-9]+).*", filename).groups() for filename in ROIs]
    ROIs = [list(int(n) for n in i) for i in ROIs]
    print(f"Found ROIs: {ROIs}")
    
    # Read image 
    print(f"Reading Maximum Projection Image File: {mip_tiff_file_path}")
    image = tiffile.imread(mip_tiff_file_path)
    
    os.makedirs(ouptut_folder, exist_ok = True)
    
    for r in ROIs:
        roi_string = "X_" + str(r[0]) +  "_Y_" + str(r[0]) +  "_W_" + str(r[2]) +  "_H_"+ str(r[2])
        print(f"Processing: {roi_string}")
        
        sub_transcripts = transcripts.loc[(transcripts['x_location_px'] >= r[0]) & (transcripts['x_location_px'] <= (r[0] + r[2])) & (transcripts['y_location_px'] >= r[1]) & (transcripts['y_location_px'] <= (r[1] + r[3]))]
        sub_image = image[r[1] : r[1] + r[3], r[0] : r[0] + r[2]]
        
        # We exlude all cells not fully contained in the ROI!
        cell_boundaries["out"] = (cell_boundaries['vertex_x_px'] >= r[0]) & (cell_boundaries['vertex_x_px'] <= (r[0] + r[2])) & (cell_boundaries['vertex_y_px'] >= r[1]) & (cell_boundaries['vertex_y_px'] <= (r[1] + r[3]))
        filter_df = cell_boundaries[["cell_id", "out"]].groupby("cell_id").filter(lambda x: all(x["out"])).drop_duplicates()
            
        sub_cell_boundaries = pd.merge(cell_boundaries, filter_df["cell_id"], on = "cell_id", how = "inner")
        sub_cell_boundaries = sub_cell_boundaries.drop('out', axis=1)
        sub_cells = pd.merge(cells, filter_df["cell_id"], on = "cell_id", how = "inner")
            
        
        print(f"Writing {roi_string}")
        subfolder= os.path.join(ouptut_folder, roi_string)
        os.makedirs(subfolder, exist_ok = True)
        
        tiffile.imwrite(os.path.join(subfolder, roi_string + "-morphology_mip.ome.tif"),data = sub_image, ome = True)
        sub_transcripts.to_parquet(os.path.join(subfolder, roi_string + "-transcritps.parquet"))
        sub_transcripts.to_csv(os.path.join(subfolder, roi_string + "-transcritps.csv"))
        sub_cells.to_parquet(os.path.join(subfolder, roi_string + "-cells.parquet"))
        sub_cells.to_csv(os.path.join(subfolder, roi_string + "-cells.csv"))
        sub_cell_boundaries.to_parquet(os.path.join(subfolder, roi_string + "-cell_boundaries.parquet"))
        sub_cell_boundaries.to_csv(os.path.join(subfolder, roi_string + "-cell_boundaries.csv"))
        
