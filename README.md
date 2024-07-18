# Matlab Toolbox for the Quantification of Cerebrovascular Tortuosity (TQTV)

(Legacy code - This was a graduate internship project from 2015)

The purpose of this project was to develop a tool to quantify the tortuosity ("bendiness") of cerebral blood vessels towards better understanding of age-related stroke. Various metrics exist in the literature, however no unified testing tool was built to validated and compare across research studies. 

### Highlights
- Open-source cerebrovascular tortuosity quantification tool developed.
- Validated with simulated blood vessels and additional parameter exploration toolbox.
- Tested in hypertensive and normal-tensive human subjects.

## Step-by-Step Guide

### Opening the toolbox
1. Open Matlab R2013a (untested on newer versions) and type “guide” in the command window
2. Set the TQTV toolbox folder as a path on Matlab, Environment Tab. Also set the additional functions as a path
3. Select “Browse” on the “Open Existing GUI” tab and within the folder “TQTV Toolbox”, open “Tortuosity GUI”
4. Press the Green arrow on the toolbar to start the GUI (or press ctrl +T) 

### Importing data
5. Make a selection of MRA DICOM images (Select first image, hold shift and press last image for multiple slice selection)
6. Set threshold (Maximum and Minimum) based on trial and error judgement (default threshold is usually sufficient)
7. Select full brain projection slice (i.e. where you can see all vessels on a single slice).

### ROI selection
8. Zoom into region of interest using Magnifier tool in toolbar.
9. Define the ROI by pressing the “Select ROI” button and then select polygon vertices, creating a closed polygon (Double-click when the end point is on top of the start point).
10. The vessel automatically is plotted, then adjust minimum and maximum z slice. Hint: Select “Rotate 3D” and right click on the “z-axis limit of ROI” plot and select either X-Z or Y-Z view for easier Z-axis limit selection. Here try to encompass as much of the vessel of interest as possible, extra vessels can be removed in the next step.

### Tortuosity calculations and centreline extraction
11. Once happy with ROI selection (Re-select ROI on first panel if not), press the “plot Centreline” button within the “Vessel parameter” box and the code will extract the centreline. This may take up to a minute given the size of ROI selection.
12. Observe the resulting centreline and decide whether cleaning tools are required (most probably) – This is the tricky part and involves lots of user-judgement unfortunately.
13. Firstly, defined a branch length minimum, setting it to maximum often produces a clean centreline.
14. You can also change the “Nearest neighbour” and “Minimum distance” parameters to optimise centreline, make take some “playing” to get used to these tools.
15. Finally, you can select the “Brush” from the toolbox and mark unwanted points and delete them with “Del” button.  Hint: again select the “Rotate 3D” button and view from all projections to unsure there are no unwanted branches remaining on the centreline.
16. Once finished cleaning, select, cleaning complete, this will calculate the tortuosity metrics for the updated centreline. Occasionally the re-stacking algorithm may miss a point and the final centreline with have a big green line going nowhere useful – change the branch length and re-clean the vessel – Hint: often brushing a few extra points from the bottom helps. If this happens, brush the centreline from many projection angles.
17. Read off the values 

Note 1: When vessel thins out (unknown signal loss) sometimes centrelines will have “bubbles” or clumps of points. Generally there is a seed branch, but it takes lots of brushing.

Note 2: Deviation from the above flow of the code may result in bugs, which can be resolved via restart. 

## Components

1. The vessel extraction and measurement tool

![Tool1](./Tortuosity%20Paper%20Images/Figure_2.jpg)

2. The simulated vessel metric validation tool:

![Tool](./Tortuosity%20Paper%20Images/Figure_4.jpg)

For an overview of the toolbox, you can view the [walthough mp4 video](https://youtu.be/z_a9z4SmE-A).

## Credits

The toolbox was developed in collaboration with CRIC at Bristol University under the supervision of Dr. Emma Hart.
