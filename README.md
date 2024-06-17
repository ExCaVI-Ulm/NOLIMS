NOLIMS is a open-source MATLAB simulation package for non-linear encoded MRI. An interface to magnetic field simulation software like CST is fairly easy via import of txt-files. A more detailed description of this repository is published in:
	-	

General
- Written in MATLAB
- Freely available under BSD License
- Validated simulation package for non-linear encoded MRI

- For visualization purposes arrayShow is used (Version 0.35): https://github.com/tsumpf/arrShow
- The following MATLAB Toolboxes ares used:
	- Signal Processing Toolbox V9.0
	- Image Processing Toolbox V11.5
	- Communications Toolbox V7.7
	- Parallel Computing Toolbox V7.6
	- MATLAB Parallel Server V7.6
	- Polypace Bug Finder V3.6

Installation
- Download / clone this repository
- Open NOLIMS.m in MATLAB
- Choose simulation settings and run the script

Example data
- In the folder FieldData, you will find:
	-- Example file demonstrating the file structure of CST exported field data, NOTE: In this file, the dominant magnetic field vector is in y-direction; in the 	  
           Simulation the dominant field vector should be in z-direction

	-- Folder ShimEncode100: Contains example data for (non-linear) encoding using up to second-order shim coils. Change to Branch Example_ShimEncode and run 	   
           NOLIMS.m. You should end up with an image as stored in the folder "ExampleResult_Encoding_Shimterms" (Slice 16 of last dimension in array: ":, :, 16").

Details on Code

Flowchart of most important steps:
Fig. 1 of paper

For simplicity, the following approximations were considered: 
- Only the z-component of the susceptibility-induced additional magnetic field is calculated assuming a constant (amplitude and direction) magnetic field across the phantom.
- For considering T_1-effects, a homogeneous RF-excitation and a magnetic field in z-direction is assumed. Furthermore, variations of T_1 due to local B_0 differences are considered to be negligible.
- The phase of the magnetization vector after RF-excitation and of the coil sensitivity maps is approximated by projection onto the x-y plane.
- Intravoxel Dephasing (IVD): Either choose matrixsize_signal > matrixsize_reco or choose RAM efficient option: Note, then only the magnetic field gradients pointing towards z are included. 

Citation
If you intend to use any part of the code in the repository, please cite
-

Acknowledgemts
NOLIMS was developed within the ExCaVI group of Ulm University and Ulm University Medical Center.
