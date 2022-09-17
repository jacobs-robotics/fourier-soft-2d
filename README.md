# Fourier-SOFT-2D


Required packages:

* OpenCV
* Point Cloud Library
* OpenMP
* FFTW3

Start registration Voxel Based:
* `fourier_soft2D "scan1.jpg" "scan2.jpg" "outputFolderResults"`
* `fourier_soft2D inputFiles/voxelScan1.jpg inputFiles/voxelScan2.jpg tools/matlabVisualization/out`

Start registration PCL Based(worse results):
* `fourier_soft2D_pcl "scan1.pcd" "scan2.pcd" "outputFolderResults" "dimensionRegistration"`
* Example: `fourier_soft2D_pcl inputFiles/scan1.pcd  inputFiles/scan2.pcd tools/matlabVisualization/out 256`


Output can be printed with a small matlab script in `tools`