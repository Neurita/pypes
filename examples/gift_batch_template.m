% Enter the values for the variables required for the ICA analysis.
% Variables are on the left and the values are on the right.
% Characters must be enterd in single quotes
%
% After entering the parameters, use icatb_batch_file_run(inputFile); 

% Input data directory name
input_directory_name = '{{ input_dir }}';

%% Enter directory to put results of analysis
outputDir = '{{ output_dir }}';

%% Enter location (full file path) of the image file to use as mask
% or use Default mask which is []
maskFile = '{{ mask_file }}';

%% Modality. Options are fMRI and EEG
modalityType = 'fMRI';

%% Type of analysis
% Options are 1, 2 and 3.
% 1 - Regular Group ICA
% 2 - Group ICA using icasso
% 3 - Group ICA using MST
which_analysis = 1;


%% Group ica type
% Options are spatial or temporal for fMRI modality. By default, spatial
% ica is run if not specified.
group_ica_type = 'spatial';


%% Parallel info
% enter mode serial or parallel. If parallel, enter number of
% sessions/workers to do job in parallel
parallel_info.mode = 'parallel';
parallel_info.num_workers = 2;


%% Group PCA performance settings. Best setting for each option will be selected based on variable MAX_AVAILABLE_RAM in icatb_defaults.m. 
% If you have selected option 3 (user specified settings) you need to manually set the PCA options. 
%
% Options are:
% 1 - Maximize Performance
% 2 - Less Memory Usage
% 3 - User Specified Settings
perfType = 2;


%% Design matrix selection
% Design matrix (SPM.mat) is used for sorting the components
% temporally (time courses) during display. Design matrix will not be used during the
% analysis stage except for SEMI-BLIND ICA.
% options are ('no', 'same_sub_same_sess', 'same_sub_diff_sess', 'diff_sub_diff_sess')
% 1. 'no' - means no design matrix.
% 2. 'same_sub_same_sess' - same design over subjects and sessions
% 3. 'same_sub_diff_sess' - same design matrix for subjects but different
% over sessions
% 4. 'diff_sub_diff_sess' - means one design matrix per subject.

keyword_designMatrix = 'no';

% specify location of design matrix here if you have selected 'same_sub_same_sess' or
% 'same_sub_diff_sess' option for keyword_designMatrix variable
%OnedesignMat = 'C:\MATLAB6p5p2\work\Example Subjects\Visuomotor_data\SPM.mat';


%% There are three ways to enter the subject data
% options are 1, 2, 3 or 4
dataSelectionMethod = 3;

%% Method 3 (Uses Regular expressions)

% Subject directory regular expression. This variable can have nested paths
% like Sub01_vis\Study1. To match this Sub\w+; Study\w+ regular expression can be used where semi-colon
% is used as a path separator. If there are no subject directories inside the input directory, leave it as empty like ''
subject_dir_regexp = '';

% Session directory regular expression. This variable cannot have nested
% paths. If there are no session directories inside subject directories, leave it as empty.
session_dir_regexp = '';

% Data file pattern. Use wild card for this and not regular expression.
data_file_pattern = '*.nii';

% File numbers to include. Leave it as empty if you want to include all of
% them.
file_numbers_to_include = [];

% SPM stats directory name relative to subject or session directories. Use this only when you specify
% 'diff_sub_diff_sess' as the value for keyword_designMatrix variable. GIFT
% will first search in the subject directories and later session
% directories to find SPM.mat files
spm_stats_dir = '';


%% Enter Name (Prefix) Of Output Files
prefix = '';

%% Group PCA Type. Used for analysis on multiple subjects and sessions.
% Options are 'subject specific' and 'grand mean'. 
%   a. Subject specific - Individual PCA is done on each data-set before group
%   PCA is done.
%   b. Grand Mean - PCA is done on the mean over all data-sets. Each data-set is
%   projected on to the eigen space of the mean before doing group PCA.
%
% NOTE: Grand mean implemented is from FSL Melodic. Make sure that there are
% equal no. of timepoints between data-sets.
%
group_pca_type = 'subject specific';

%% Back reconstruction type. Options are str and gica
backReconType = 'gica';

%% Data Pre-processing options
% 1 - Remove mean per time point
% 2 - Remove mean per voxel
% 3 - Intensity normalization
% 4 - Variance normalization
preproc_type = 3;


%% PCA Type. Also see options associated with the selected pca option. EM
% PCA options and SVD PCA are commented.
% Options are 1, 2, 3, 4 and 5.
% 1 - Standard 
% 2 - Expectation Maximization
% 3 - SVD
% 4 - MPOWIT
% 5 - STP
pcaType = 1;

%% PCA options (Standard)

% a. Options are yes or no
% 1a. yes - Datasets are stacked. This option uses lot of memory depending
% on datasets, voxels and components.
% 2a. no - A pair of datasets are loaded at a time. This option uses least
% amount of memory and can run very slower if you have very large datasets.
pca_opts.stack_data = 'no';

% b. Options are full or packed.
% 1b. full - Full storage of covariance matrix is stored in memory.
% 2b. packed - Lower triangular portion of covariance matrix is only stored in memory.
pca_opts.storage = 'packed';

% c. Options are double or single.
% 1c. double - Double precision is used
% 2c. single - Floating point precision is used.
pca_opts.precision = 'single';

% d. Type of eigen solver. Options are selective or all
% 1d. selective - Selective eigen solver is used. If there are convergence
% issues, use option all.
% 2d. all - All eigen values are computed. This might run very slow if you
% are using packed storage. Use this only when selective option doesn't
% converge.
pca_opts.eig_solver = 'selective';


% %% PCA Options (Expectation Maximization)
% % a. Options are yes or no
% % 1a. yes - Datasets are stacked. This option uses lot of memory depending
% % on datasets, voxels and components.
% % 2a. no - A pair of datasets are loaded at a time. This option uses least
% % amount of memory and can run very slower if you have very large datasets.
pca_opts.stack_data = 'no';

% 
% % b. Options are double or single.
% % 1b. double - Double precision is used
% % 2b. single - Floating point precision is used.
pca_opts.precision = 'single';

%
% % c. Stopping tolerance 
pca_opts.tolerance = 1e-4;

%
% % d. Maximum no. of iterations
pca_opts.max_iter = 1000;


% %% PCA Options (SVD)
% % a. Options are double or single.
% % 1a. double - Double precision is used
% % 2a. single - Floating point precision is used.
%pca_opts.precision = 'single';

% % b. Type of eigen solver. Options are selective or all
% % 1b. selective - svds function is used.
% % 2b. all - Economy size decomposition is used.
%pca_opts.solver = 'all';


%% Maximum reduction steps you can select is 2. Options are 1 and 2. For temporal ica, only one data reduction step is
% used.
% numReductionSteps = 2;

%% Batch Estimation. If 1 is specified then estimation of 
% the components takes place and the corresponding PC numbers are associated
% Options are 1 or 0
doEstimation = 1;


%% MDL Estimation options. This variable will be used only if doEstimation is set to 1.
% Options are 'mean', 'median' and 'max' for each reduction step. The length of cell is equal to
% the no. of data reductions used.
estimation_opts.PC1 = 'mean';
estimation_opts.PC2 = 'mean';

%% Number of pc to reduce each subject down to at each reduction step
% The number of independent components the will be extracted is the same as 
% the number of principal components after the final data reduction step.  
numOfPC1 = 20;
numOfPC2 = 16;


%% Scale the Results. Options are 0, 1, 2, 3 and 4
% 0 - Don't scale
% 1 - Scale to Percent signal change
% 2 - Scale to Z scores
% 3 - Normalize spatial maps using the maximum intensity value and multiply timecourses using the maximum intensity value
% 4 - Scale timecourses using the maximum intensity value and spatial maps using the standard deviation of timecourses
scaleType = 2;


%% 'Which ICA Algorithm Do You Want To Use';
% see icatb_icaAlgorithm for details or type icatb_icaAlgorithm at the
% command prompt.
% Note: Use only one subject and one session for Semi-blind ICA. Also specify atmost two reference function names

% 1 means infomax, 2 means fastICA, etc.
algoType = 1;

%% Specify atmost two reference function names if you select Semi-blind ICA algorithm.
% Reference function names can be acessed by loading SPM.mat in MATLAB and accessing 
% structure SPM.xX.name.
refFunNames = {'Sn(1) right*bf(1)', 'Sn(1) left*bf(1)'};


%% Specify spatial reference files for constrained ICA (spatial) or moo-icar
refFiles = {which('ref_default_mode.nii'), which('ref_left_visuomotor.nii'), which('ref_right_visuomotor.nii')};

%% ICA Options - Name by value pairs in a cell array. Options will vary depending on the algorithm.
%% See icatb_icaOptions for more details. Some options are shown below.
% Infomax -  {'posact', 'off', 'sphering', 'on', 'bias', 'on', 'extended', 0}
% FastICA - {'approach', 'symm', 'g', 'tanh', 'stabilization', 'on'}

icaOptions = {'posact', 'off', 'sphering', 'on', 'bias', 'on', 'extended', 0};
