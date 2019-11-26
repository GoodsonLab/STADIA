% This is the file that you run.

% This is a script to input parameters and run code to extract dynamic
% instability parameters from simulation data.

% This script calls the script Loop_Thru_Inputs_[date].m, which calls the
% script extract_DI_params_[date].m, which calls the function
% least_squares_regression_line.m. A total of six ".m"-files should be
% placed in the same folder/directory as the length history input files.

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT PARAMETERS

% Names of PARAMETERS are written in all uppercase letters.

% It is assumed that the input files have MT length in dimers.

% Are you running the program multiple times on the same input file?
SKIP_FILE_READ = 0; % 0-off; 1-on; When on, this option saves time by skipping the file read step
if(SKIP_FILE_READ == 0) % If off, clear variables and screen, and close figures
    clearvars;
    clear;
    close all; clc;
    SKIP_FILE_READ = 0;
end

% Name(s) of length history input file(s) to be analyzed. Each file name should be
% in single quotes, and the list of file names should be in {}. For
% example, FILE_NAME_INPUT = {'163.dat'} or FILE_NAME_INPUT = {'163.txt', '165.txt'};
FILE_NAME_INPUT = {'DetailedMTmodel_length_13PF_LongRun_10uM.dat'}; %{'BiPhase_Simulation_length1.csv'}; %

% What is the first row of the input file that contains data? (not column headers)
FIRST_DATA_ROW = 2;

% Which columns of the input file contain the microtubule lengths?
% For example, MT_LENGTH_COLUMN_INDICES = 2 or 
% MT_LENGTH_COLUMN_INDICES = 1:128;
% Or you can do some subset of the length columns, for example,
% MT_LENGTH_COLUMN_INDICES = [1:3, 5, 10:15];
%
% NOTE 1:
% If multiple columns of the input data file contain microtubule lengths,
% then they will be stitched together in the length history plots that are
% created. This is for visualization purposes only; only the first column
% of data will plotted corresponding to the actual time values in the input
% data, and the remaining columns will be shifted to further time steps. In
% any case, the time duration between input data points will remain
% unaffected.
% 
% NOTE 2: 
% If multiple input data files are being used, then all of the files must
% contain microtubule length values in the columns defined here.
MT_LENGTH_COLUMN_INDICES = 2;%:11;%[2:16];

% Which column contains the times?
% NOTE: only one column may contain time values.
TIME_COLUMN = 1;

% If the TIME_COLUMN of the input data file contains times in seconds, set
% TIME_CONVERSION_FACTOR = 1.
% If the TIME_COLUMN of the input data file contains times in minutes, set
% TIME_CONVERSION_FACTOR = 60, and so on. 
% If the TIME_COLUMN of the input data file contains frame/step numbers instead of time values, and
% there is a fixed time step between, then set TIME_CONVERSION_FACTOR equal to the
% frame rate, i.e. the length of each time step in seconds.
TIME_CONVERSION_FACTOR = 1;

% What is the delimeter character used in the input data files?
% The character used to seperate the entries in the length history data
% files must be consistent when using multiple input data files. 
% Example: INPUT_FILE_DELIMITER = ' '; or INPUT_FILE_DELIMITER = ', '; or INPUT_FILE_DELIMITER = ';';
INPUT_FILE_DELIMITER = ' ';

%%%%%%%%%%

% The following input parameters can be given as scalars or vectors. 
% For example, MIN_PEAK_HEIGHT_INPUT = 100 or 
% MIN_PEAK_HEIGHT_INPUT = [75, 100, 125];
% If an input parameter is a vector, the code will loop through each entry
% in the vector and run the script extract_DI_params_[date] in serial for
% each value of the input parameter in the vector.

% If any of the inputs are given as vectors, then the parameter
% INCLUDE_DATE_TIME should be set to 1 below, so that the date and time
% will be included in the output file names, to distinguish between the
% output files corresponding to each of the entries in the vector.

%%%%%
% Input Parameters for the Piecewise Linear Approximation to Length History
% Enforce a min time step separating those data points accepted into the
% approximation?
MIN_TIME_STEP_INPUT = 0.5; % Minimum time (in seconds) between segment vertices used in the approximation

% Point-wise Error Tolerance for linear segments approximating length history plot
%LINEAR_FIT_TOLERANCE_INPUT = 20;%8; %17; %12;  % Tolerance level: # of subunits margin of error for linear approximation
ERROR_TOLERANCE_LEVEL = 20;%8; %17; %12;  % Tolerance level: # of subunits margin of error for linear approximation
%%%%%

%%%%
% Input parameters for the DI Phase Classification using k-means
% The classification first assumes that Flat Stutters and Nucleation segments 
% have been removed using the following user defined thresholds
FLATSTUT_MAX_SEGMENTHEIGHT_THRESHOLD_INPUT = 1;%10; % Segment deltah abs value <= FLATSTUT_MAX_SEGMENTHEIGHT_THRESHOLD ---> Flat Stutters
FLATSTUT_MAX_SEGMENTSLOPE_THRESHOLD_INPUT = 0.5; %0.5;%1;% Segment slope abs value <= FLATSTUT_MAX_SEGMENTSLOPE_THRESHOLD ---> Flat Stutters
NUC_HEIGHT_THRESHOLD_INPUT = 75; % Segment endpoint height value <= NUC_Height_THRESHOLD ---> Nucleation

% DI Phase Classification has 3 diagnostic options, which allow the user to
% view the possible k-means classification options prior to commiting to
% them. When activated, seperate plots for the complete data set, positive
% slope segments only, and negative sloped segments only will be displayed.
% Additionally, the seperate k-means classification results will be shown
% using an automatically generated suggest k-value. Choosing the 
% appropriate k-values (number of clusters) s left to the user below. 
% Keeping the Diagnostic option activated will prevent the program from 
% continuing into the next stage of classification.

% Diagnostic Mode Options
KMEANS_DIAGNOSTIC_AllFLAG = 0; % 0-off; 1-on; Flag will generate plots to help develop k-means for the all segments (DEMONSTRATION PURPOSES ONLY)
KMEANS_DIAGNOSTIC_PosSlopeFLAG = 0; % 0-off; 1-on; Flag will generate plots to help develop k-means for the positive slope segments
KMEANS_DIAGNOSTIC_NegSlopeFLAG = 0; % 0-off; 1-on; Flag will generate plots to help develop k-means for the negative slope segments

% User-selected k-values to indication the number of cluster to be used
KMEANS_NumClust_PosSlope = 1; % # of clusters to identify for positive slope DI segments. Limit choices to 1, 2, or 3
KMEANS_NumClust_NegSlope = 1; % # of clusters to identify for negative slope DI segments. Limit choices to 1, 2, or 3

% KMEANS CLASSIFICATION OPTIONS ONLY IF k-value =2 CHOSEN ABOVE 
% Options for positive slope segments: 
% 'A' ==> use 2 growth phases (long and brief growth)
% 'B' ==> use 1 up-stutter phase, and 1 long growth phase
KMEANS_Pos2_Option = 'A'; 
% Options for negative slope segments: 
% 'A' ==> use 2 shortening phases (long and brief shortening)
% 'B' ==> use 1 down-stutter phase, and 1 long shortening phase
KMEANS_Neg2_Option = 'B'; 
%%%%%%%%%%


% Figure 1 shows the piecewise linear approximation to the original length history data
% Figure 2 shows the clustered results for the scaled and standardized pos & neg slope DI Segments
% Figure 3 shows the classification results within the DI phase variable space
% Figure 4 shows the plot in Fig1 with color labels of DI phases for each linear segment
% Figure 5 shows the average measurements for the different DI Phase Segments
% Figure 6 shows the total measurements for the different DI Phase Segments
% Figure 7 shows the resulting measurements for the possible changes in DI Phase

% Ths following parameters should be scalars, specifically 0 or 1, not vectors.

% Set PLOT_FIG_# = 1 to plot figure, set to 0 otherwise.
PLOT_FIG_1 = 0;
PLOT_FIG_2 = 1;
PLOT_FIG_3 = 1;
PLOT_FIG_4 = 1; 
PLOT_FIG_5 = 1;
PLOT_FIG_6 = 1; 
PLOT_FIG_7 = 1;  

% View Window for Figures 1 & 4 ( Time Values only. MT Height Window is
% automatically the data maximum)
FIG_WINDOW_START = 2000;
FIG_WINDOW_END = 3600;

% Set SAVE_FIG_# = 1 to save and close figure, set to 0 otherwise.
% SAVE_FIG_# can only be 1 if PLOT_FIG_# is also 1.
% If SAVE_MATLAB_WORKSPACE = 1 below, then all figures will be closed. So,
% in that case, the figures should be saved here if you want to see them.
SAVE_FIG_1 = 1;
SAVE_FIG_2 = 1;
SAVE_FIG_3 = 1;
SAVE_FIG_3_movie = 0; % A movie that spins the 3D plot for DI phase classification
SAVE_FIG_4 = 1; 
SAVE_FIG_5 = 1;
SAVE_FIG_6 = 1;
SAVE_FIG_7 = 1;
        
% Set SAVE_MATLAB_WORKSPACE = 1 to save the variables in the matlab
% workspace as a .mat file, 0 otherwise. 
% All figures will be closed before saving the workspace, so that they
% won't also be saved as part of the workspace.
% The time and length data from the input data files won't be saved as part
% of the .mat file because they take up a lot of space and they are in your
% input files anyway.
SAVE_MATLAB_WORKSPACE = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call script to loop through inputs and run code to extract DI parameters
% from simulation data for each set of inputs.

toc

Loop_Thru_Inputs_2019_07_16

toc

'done'


