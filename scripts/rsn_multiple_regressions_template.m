clear templateFiles; clear dispParameters;
clear all;

% this script is a modified version of icatb/icatb_scripts/icatb_example_spatial_sorting.m
% found in GroupICATv4.0a
% '/home/alexandre/data/thomas/ica_out/8mm/fmri_no-grptemplate_noWMcor_30ICs/'
output_dir = '{{ output_dir }}';

% 'Basal_Ganglia_21'
output_basename = '{{ output_basename }}';
template_index = {{ template_ic_index }};

% Specify ICA parameter file
% '/home/alexandre/data/thomas/ica_out/8mm/fmri_no-grptemplate_noWMcor_30ICs/_ica_parameter_info.mat';
param_file = '{{ ica_param_file }}';

% Specify the RSN template file and all indices
% '/home/alexandre/data/std_brains/resting_state/allen/baseline/ALL_HC_unthresholded_tmaps_resampled.nii';
template_file = '{{ rsn_atlas_file }}';

% build the list of RSN templates that is inside the NifTI file.
%template = load_nii(template_file);
%for rsn_idx = 1:length(template.img)
%    templateFiles(rsn_idx).name = [template_file, ',', int2str(rsn_idx)];
%end
templateFiles(1).name = [template_file, ',', int2str(template_index)];

% Options for selectedStr are:
% 1.'Same set of spatial templates for all data-sets'
% 2. 'Different set of spatial templates for sessions'
% 3. 'Different set of spatial templates for subjects and sessions'
selectedStr =  'Same set of spatial templates for all data-sets';

zscore_threshold = {{ zscore_threshold }};

%%%%%%%% End for specifying display parameters %%%%%%%%%%

icatb_defaults;

if ~exist('param_file', 'var')
    param_file = icatb_selectEntry('typeSelection', 'single', 'typeEntity', 'file', 'title', 'Select a valid ICA parameter file', ...
        'filter', ['*', PARAMETER_INFO_MAT_FILE, '*']);
end

load(param_file);

if ~exist('sesInfo', 'var')
    error(['file: ', param_file, ' is not a valid ica parameter file']);
end

% Get the output directory
[outputDir, fileName, extn] = fileparts(param_file);

if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

%%%%%%%%%% Get the required variables from sesInfo structure %%%%%%%%%%
% Number of subjects
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;

% Number of components
numComp = sesInfo.numComp;

dataType = sesInfo.dataType;

mask_ind = sesInfo.mask_ind;

% First scan
structFile = deblank(sesInfo.inputFiles(1).name(1, :));

structVol = icatb_get_vol_nifti(structFile);

DIM = structVol(1).dim(1:3);

% slices in mm
[parameters] = icatb_get_slice_def(structVol, 'axial');
slices_in_mm = parameters.slices; clear parameters;
%%%%%%%% End for getting the required vars from sesInfo %%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% Get component data %%%%%%%%%%%%%%%%%
% Get the ICA Output files
icaOutputFiles = sesInfo.icaOutputFiles;

[subjectICAFiles, meanICAFiles, ...
 tmapICAFiles, meanALL_ICAFile] = icatb_parseOutputFiles('icaOutputFiles', icaOutputFiles, 'numOfSub', numOfSub,  ...
                                                         'numOfSess', numOfSess, 'flagTimePoints', sesInfo.flagTimePoints);

% component data
numOfSess = 1; % this only works with 1 session per subject

compData = zeros(numComp, numOfSub, length(mask_ind));

for nSub = 1:numOfSub

    % Threshold, Include Image values and convert to z scores accordingly
    % to the template image
    disp(['Loading Subject ', num2str(nSub), ' Session ', num2str(numOfSess), '...']);
    compFiles = subjectICAFiles(nSub).ses(numOfSess).name;

    % component files
    compFiles = icatb_fullFile('directory', outputDir, 'files', compFiles);

    % Load the images
    % Apply Z-scores, threshold and convert to z criteria
    % load ICA images and time courses
    disp('Loading component data and applying display defaults ...');

    [icasig, HInfo, real_world_coords] = icatb_loadData(compFiles, 'real', [], [], ...
                                                        [1:numComp]);

    % Reshape icasig to components by voxels
    icasig = permute(icasig, [4 1 2 3]);

    % Structural volume
    HInfo.V = HInfo.V(1);

    img_dim = HInfo.DIM(1:3);

    % Reshape to 2d
    icasig = reshape(icasig, [numComp, prod(img_dim)]);

    % z-score spatial map
    icasig = icatb_applyDispParameters(icasig, 1, 2, zscore_threshold, img_dim, HInfo);

    % store in matrix
    compData(:, nSub, :) = icasig(:, mask_ind);

    clear icasig;
end

%%%%%%%%%%%%%%%%%%%%%%%%% End for getting component data %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Get template data %%%%%%%%%%%%%%%%%%%%%%%%%%
spatialTemplate = [];

tempV = icatb_get_vol_nifti(HInfo.V(1).fname);

% Handle flip
spatialTemplate = icatb_resizeData(tempV, templateFiles(1).name, 1);

% New dimensions for the spatial template
spatialTemplate = reshape(spatialTemplate, size(spatialTemplate, 1), prod(HInfo.DIM));

spatialTemplate = spatialTemplate(:, mask_ind);

clear tempV;
%%%%%%%%%%%%% End for getting template data %%%%%%%%%%%%%%%%%%%%%%%%%%


helpMsg = 'Calculating Regression ...';
disp(helpMsg)

%compDIMS = repmat(length(mask_ind), 1, numOfSub*numOfSess);

for nSub = 1:numOfSub
    % Multiple regression
    for nComp = 1:numComp
        disp(['Calculating Multiple regression for component ', num2str(nComp)]);

        compMap = compData(nComp, nSub, :);
        compMap = reshape(compMap, 1, size(compData, 3));
        % calculate regression
        %[comparison, regressCoeff, ModelIndices, otherIndices, linearRegress, removeTrend, ...
        %    icaTimecourse(:, nComp), sub_partial_corr, partialCorrSlopes] = ...
        %    icatb_multipleRegression(modelTimecourse, icaTimecourse(:, nComp), num_Regress, num_DataSets, diffTimePoints);

        % [comparison, regressCoeff, ModelIndices, otherIndices, linearRegress, removeTrend, ...
        % icaTimecourse(:, nComp), sub_partial_corr, partialCorrSlopes] = icatb_multipleRegression(spatialTemplates, compData(nComp, :), ...
        %                                                                                          numRegressors, numOfSub*numOfSess, compDIMS, 0);
        [comparison, regressCoeff, ModelIndices] = ...
                                icatb_multipleRegression(spatialTemplate, compMap);

        % detrend data for beta weights
        % tmpData = detrend(compData(nComp, :), 0);
        % tmpModel = detrend(multi_regress_data(spatialTemplates, :), 0);
        % [comparison, beta_weights, ModelIndices] = icatb_multipleRegression(tmpModel, tmpData, ...
        %                                                                     numRegressors, numOfSub*numOfSess, compDIMS, 0);

        % [comparison, regressCoeff, ModelIndices, otherIndices, linearRegress, removeTrend, ...
        % icaTimecourse(:, nComp), sub_partial_corr, partialCorrSlopes] = icatb_multipleRegression(tmpModel, tmpData, ...
        %                                                                                          numRegressors, numOfSub*numOfSess, compDIMS, 0);

        % store values from comparison
        %subject_partial_corr{nComp} = sub_partial_corr;
        %subject_partial_slopes{nComp} = partialCorrSlopes;

        % Truncate the outputs of regression
        % betaWeights{nComp} = beta_weights(ModelIndices);
        regressionCoeff{nComp, nSub} = regressCoeff(ModelIndices);
        clear  regressCoeff;
    end
end

clear spatialTemplate;

regression_values = cell2mat(regressionCoeff);
csvwrite([output_dir, output_basename, '_regression_values.csv'], regression_values);

% calculate the beta weights
% temp_regression = regression_values';
% for img_idx = 1:size(compData, 1)
%     betaWeights{img_idx} = temp_regression(:, img_idx)
% end
%
% beta_weights = cell2mat(betaWeights);
% csvwrite([output_dir, output_basename, '_beta_weights.csv'], beta_weights);

exit;
%%%%%%%% End for calculating Multiple Regression %%%%%%%%%%%
