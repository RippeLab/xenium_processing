# xenium_processing

## Xenium data documentation:

Links:
- General Overview: https://www.10xgenomics.com/support/in-situ-gene-expression/documentation/steps/onboard-analysis/at-a-glance-xenium-output-files
- Output file specifications: https://www.10xgenomics.com/support/in-situ-gene-expression/documentation/steps/onboard-analysis/understanding-xenium-outputs
- Zarr outputs: https://www.10xgenomics.com/support/in-situ-gene-expression/documentation/steps/onboard-analysis/xenium-outputs-zarr
- Analyisis performed by the Xenium machine: https://www.10xgenomics.com/support/in-situ-gene-expression/documentation/steps/onboard-analysis/xenium-algorithms-overview
- Archiving Xenium Data: https://www.10xgenomics.com/support/in-situ-gene-expression/documentation/steps/onboard-analysis/archiving-xenium-data
- Analysis Guide: https://www.10xgenomics.com/resources/analysis-guides?query=&page=1&refinementList%5Btopics%5D=


## xenium_ROI-extractor.py
Extracts the cells, transcripts and maximum projection image from a Xenium dataset according to one or more ROIs.

### Input:
- `xenium_folder` = Path to the xenium folder, it expects the standard Xenium filenames.
- `roi_folder` = Path to the folder where to look for ROIs. It expects one or more files with names in this form: `_Xnnn_Ynnnn_Wnnn_Hnnnn_` with n = 0-9. The files are not used, only the names are needed to get the position and size of theROIs.
- `output_folder` = Path for output, the script will create it plus all the required subfolders.


### Processing:
1. Read the xenium metadata file to get the pixel size in µm.
2. Read the transcripts file, the cell file and the cell boundaries file, and add columns with the positions in pixels.
3. Look through the ROI folder for ROIs to extract.
4. Read the maximum projection image file.
5. For each ROI:
    - 5.1. Subsets the transcripts, cells and cell boundaries. Only cells fully included in the ROI are retained. 
    - 5.2. Adds columns to the transcripts, cells, and cell boundaries table with the positions in pixels relative to the top left corner of the ROI.
    - 5.3. Generates a semantic segmentation mask using the cell boundaries.
    - 5.4. Wites out the outputs. 

### Output:

For each ROI:
- Transcripts
    - Same columns as the original, plus:
        - `x_location_px` -> Pixel X coordinate.
        - `y_location_px` -> Pixel Y coordinate.
        - `x_location_px_ROI` -> Pixel X coordinate relative to the ROI.
        - `y_location_px_ROI` -> Pixel Y coordinate relative to the ROI.
    - Files:
        - parquet: `ROI-transcritps.parquet`
        - csv:  `ROI-transcritps.parquet`
- Maximum Projection Image: `ROI-morphology_mip.ome.tif`
- Cells
    - Same columns as the original, plus:
        - `x_centroid_px` -> Pixel X coordinate.
        - `y_centroid_px` -> Pixel X coordinate.
        - `x_centroid_px_ROI` -> Pixel X coordinate relative to the ROI.
        - `y_centroid_px_ROI` -> Pixel Y coordinate relative to the ROI.
    - Files
        - parquet: `ROI-cells.parquet`
        - csv: `ROI-cells.csv`
- Cell Boundaries
    - Same columns as the original, plus:
        - `vertex_x_px` -> Pixel X coordinate.
        - `vertex_y_px` -> Pixel Y coordinate.
        - `vertex_x_px_ROI` -> Pixel X coordinate relative to the ROI.
        - `vertex_y_px_ROI` -> Pixel Y coordinate relative to the ROI.
    - Files
        - parquet: `ROI-cell_boundaries.parquet`
        - csv:  `ROI-cell_boundaries.csv`
- Cell Mask: `ROI-xenium_cell_mask.ome.tif`

