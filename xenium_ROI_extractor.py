import pandas as pd
import numpy as np
import json
import os
import re
import tifffile
import argparse
import h5py
import skimage
import scipy

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
    features_file_path = os.path.join(xenium_folder, "cell_feature_matrix.h5")
    mip_tiff_file_path = os.path.join(xenium_folder, "morphology_mip.ome.tif")
    print(f"Xenium file: {xenium_file_path}")
    print(f"Transcript file: {transcript_file_path}")
    print(f"Cell file: {cells_file_path}")
    print(f"Cell Boundary file: {cell_boundaries_file_path}")
    print(f"Feature Matrix file: {features_file_path}")
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
    
    # Read features
    # Features are equivalent to: transcripts.loc[(transcripts["qv"] > 20) & (transcripts["cell_id"] != "UNASSIGNED")]
    print(f"Processing features: {features_file_path}")
    with h5py.File(features_file_path, "r") as f:
        feature_attrs = dict(f.attrs.items())
        feature_h5 = {}
        feature_h5['barcodes'] = f["matrix"].get("barcodes")[:]
        feature_h5['data'] = f["matrix"].get("data")[:]
        feature_h5['indices'] = f["matrix"].get("indices")[:]
        feature_h5['indptr'] = f["matrix"].get("indptr")[:]
        feature_h5['shape'] = f["matrix"].get("shape")[:]
        feature_h5['features'] = {}
        feature_h5['features']["_all_tag_keys"] = f["matrix"].get("features").get("_all_tag_keys")[:]
        feature_h5['features']["feature_type"] = f["matrix"].get("features").get("feature_type")[:]
        feature_h5['features']["genome"] = f["matrix"].get("features").get("genome")[:]
        feature_h5['features']["id"] = f["matrix"].get("features").get("id")[:]
        feature_h5['features']["name"] = f["matrix"].get("features").get("name")[:]
    features = scipy.sparse.csr_matrix((feature_h5["data"], feature_h5["indices"], feature_h5["indptr"]), dtype=int)
    features = pd.DataFrame(features.toarray(), columns = feature_h5['features']["name"])
    features.columns = features.columns.astype(str)
    features["cell_id"] = cells["cell_id"]
    features = features[list(features.columns)[-1:] + list(features.columns)[:-1]]

    # Read image 
    print(f"Reading Maximum Projection Image File: {mip_tiff_file_path}")
    image = tifffile.imread(mip_tiff_file_path)
    
    os.makedirs(ouptut_folder, exist_ok = True)
    
    for r in ROIs:
        roi_string = "X_" + str(r[0]) +  "_Y_" + str(r[1]) +  "_W_" + str(r[2]) +  "_H_"+ str(r[3])
        print(f"Processing: {roi_string}")
        
        sub_transcripts = transcripts.loc[(transcripts['x_location_px'] >= r[0]) & (transcripts['x_location_px'] <= (r[0] + r[2])) & (transcripts['y_location_px'] >= r[1]) & (transcripts['y_location_px'] <= (r[1] + r[3]))]
        sub_transcripts.loc[:, "x_location_px_ROI"] = sub_transcripts["x_location_px"] - r[0]
        sub_transcripts.loc[:, "y_location_px_ROI"] = sub_transcripts["y_location_px"] - r[1]
        
        sub_image = image[r[1] : r[1] + r[3], r[0] : r[0] + r[2]]
        
        # We exlude all cells not fully contained in the ROI!
        cell_boundaries["out"] = (cell_boundaries['vertex_x_px'] >= r[0]) & (cell_boundaries['vertex_x_px'] <= (r[0] + r[2])) & (cell_boundaries['vertex_y_px'] >= r[1]) & (cell_boundaries['vertex_y_px'] <= (r[1] + r[3]))
        filter_df = cell_boundaries[["cell_id", "out"]].groupby("cell_id").filter(lambda x: all(x["out"])).drop_duplicates()
            
        sub_cell_boundaries = pd.merge(cell_boundaries, filter_df["cell_id"], on = "cell_id", how = "inner")
        sub_cell_boundaries = sub_cell_boundaries.drop('out', axis=1)
        sub_cell_boundaries["vertex_x_px_ROI"] = sub_cell_boundaries["vertex_x_px"] - r[0]
        sub_cell_boundaries["vertex_y_px_ROI"] = sub_cell_boundaries["vertex_y_px"] - r[1]
        
        sub_cells = pd.merge(cells, filter_df["cell_id"], on = "cell_id", how = "inner")
        sub_cells["x_centroid_px_ROI"] = sub_cells["x_centroid_px"] - r[0]
        sub_cells["y_centroid_px_ROI"] = sub_cells["y_centroid_px"] - r[1]
        
        # Make a cell mask for the ROI by getting the coordinates of every pixel inside each cell and then drawing them on ask of the size of the ROI
        # The cells are numbered in the same order as they appear in the cell_boundaries.parquet file
        cell_coords = sub_cell_boundaries.groupby("cell_id").apply(lambda x: skimage.draw.polygon(x["vertex_x_px_ROI"], x["vertex_y_px_ROI"], (r[2], r[3])))
        cell_mask = np.zeros((r[3], r[2]), dtype=np.int32) 
        for n, i in enumerate(cell_coords):
            cell_mask[i[1], i[0]] = n
        
        # Process features
        sub_features = pd.merge(features, filter_df["cell_id"], on = "cell_id", how = "inner")
        fmatrix = sub_features.drop("cell_id", axis=1)
        fmatrix = scipy.sparse.csr_array(fmatrix)
        
        # Make transcript table for reding with the Polylux viewer (Resolve Biosciences)
        polylux_tx = sub_transcripts[["x_location_px_ROI", "y_location_px_ROI", "z_location_px", "feature_name", "qv"]]
        polylux_tx["x_location_px_ROI"] = polylux_tx["x_location_px_ROI"].astype(np.int32)
        polylux_tx["y_location_px_ROI"] = polylux_tx["y_location_px_ROI"].astype(np.int32)
        polylux_tx["z_location_px"] = polylux_tx["z_location_px"].astype(np.int32)
    
        # Output
        print(f"Writing {roi_string}")
        subfolder= os.path.join(ouptut_folder, roi_string)
        os.makedirs(subfolder, exist_ok = True)
        
        # Make h5 file
        with h5py.File(os.path.join(subfolder, roi_string + "-cell_feature_matrix.h5"), "w") as f:
            for k,v in feature_attrs.items():
                f.attrs.create(k, v)
            m_to_w = f.create_group("matrix")
            m_to_w.create_dataset("barcodes", data=np.array(sub_features["cell_id"]))
            m_to_w.create_dataset("data", data=fmatrix.data)
            m_to_w.create_dataset("indices", data=fmatrix.indices)
            m_to_w.create_dataset("indptr", data=fmatrix.indptr)
            m_to_w.create_dataset("shape", data=fmatrix.shape)
            f_to_write = m_to_w.create_group("features")
            f_to_write.create_dataset("_all_tag_keys", data=feature_h5["features"]["_all_tag_keys"])
            f_to_write.create_dataset("feature_type", data=feature_h5["features"]["feature_type"])
            f_to_write.create_dataset("genome", data=feature_h5["features"]["genome"])
            f_to_write.create_dataset("id", data=feature_h5["features"]["id"]) 
            f_to_write.create_dataset("name", data=feature_h5["features"]["name"])
        
        # Write other output
        tifffile.imwrite(os.path.join(subfolder, roi_string + "-morphology_mip.ome.tif"), data = sub_image, ome = True, compression = "lzw")
        sub_transcripts.to_parquet(os.path.join(subfolder, roi_string + "-transcritps.parquet"))
        sub_transcripts.to_csv(os.path.join(subfolder, roi_string + "-transcritps.csv"))
        polylux_tx.to_csv(os.path.join(subfolder, roi_string + "-Polylux_transcritps.txt"), header = None, index = False, sep = "\t")
        sub_cells.to_parquet(os.path.join(subfolder, roi_string + "-cells.parquet"))
        sub_cells.to_csv(os.path.join(subfolder, roi_string + "-cells.csv"))
        sub_cell_boundaries.to_parquet(os.path.join(subfolder, roi_string + "-cell_boundaries.parquet"))
        sub_cell_boundaries.to_csv(os.path.join(subfolder, roi_string + "-cell_boundaries.csv"))
        sub_features.to_csv(os.path.join(subfolder, roi_string + "-feature_matrix.csv"))       
        tifffile.imwrite(os.path.join(subfolder, roi_string + "-xenium_cell_mask.ome.tif"), data = cell_mask, ome = True, compression = "lzw")
        
        