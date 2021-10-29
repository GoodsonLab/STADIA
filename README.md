# STADIA - Statistical Tool for Automated Dynamic Instability Analysis

## Brief Introduction to STADIA
This repository is made available on behalf of the Goodson Lab the University of Notre Dame to serve the microtubule community interested in analyzing and quantifying microtubule length history data representative of dynamic instability behavior. The STADIA tool offers an automated process to help reduces human intervention in order to mitigate consistency and reproducibility between different labs and users. Additionally, classification stage in STADIA offers an opportunity to identify "stutter" behaviors, an intermediate dynamic instability phase shown to be strongly associated with catastrophe events.

## Getting Started
Users interested in using STADIA can download the corresponding files from this repository here, and follow the guidance in provided tutorials to gain a better understanding for what STADIA does and how to use it.


### Download STADIA files
To download the STADIA files, begin by clicking the green "Code" button near the top of this page to reveal a dropdown menu. You can use different options to download the STADIA files:
* Option 1: Download ZIP
A compressed ZIP file containing all of the content in this repository will be downloaded. Double-click on the corresponding file icon or follow the instructions for your machine to un-compress this file and make STADIA files accessible.

* Option 2: Clone the repository
Copy the repository URL and clone the repository on your local machine using the following terminal command line prompt:

	--- 	--- `git clone https://github.com/GoodsonLab/STADIA.git`


### Read the Tutorials
There are two tutorial files:
* STADIA-tutorial.pdf : tutorial in slides format that provide an overview of STADIA, the user input parameters, and the output content generated from running the code.
* STADIA_Diagnostic_Mode_Tutorial.pdf : tutorial to run STADIA in Diagnostic Mode, useful to help select the number of clusters to use during the classification stage. 


### Test STADIA with Example Data
You can test your version of STADIA by running it with the example length history data file provided inside the ZIP file along with the STADIA code. This example data files is named `length_13PF_10uM_MeanPF_3hr.txt`, which represents a detailed 13 protofilament MT behavior simulated over a 3 hour duration, and where MT length is measured by taking the mean protofilament length. Note that this data set differs from the 10 hour simulated data using the maximum protofilament length used in corresponding publications, which may generate output different that the reported results. However, the 3 hour data should still be ample enough to test the STADIA code, and to demonstrate similar results associated with identifying stutter phase segments.

To test STADIA with this example data, simply run the `Input_and_Run.m` file with the default parameters already provided. All required files should be where they need to be after downloading and un-compressing the ZIP file for the STADIA code.

### Contact Us
With questions, suggestions, or comments, please contact us via email:
* Shant Mahserejian - shant.mahserejian@pnnl.gov
* Holly Goodson - hgoodson@nd.edu

