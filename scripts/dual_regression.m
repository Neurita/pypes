%compFiles = '/home/alexandre/data/thomas/ica_out/8mm/gift_pet_grptemplate_no-pvc_30ICs/_sub01_component_ica_s1_.nii'
%outputDir = '/home/alexandre/data/thomas/ica_out/8mm/gift_pet_grptemplate_no-pvc_30ICs'

% PET
% Look in "tum/clinical/Resting-state ICA blob orgy.ipynb" on how these blobs were created.
% PET_fMRI_intersection_w_PET_values_on_PET
compFiles = '/home/alexandre/data/thomas/ica_out/8mm/spatio-temporal_regression/PET_fMRI_intersection_w_pet_values_on_pet/aDMN_RSN.nii'
outputDir = '/home/alexandre/data/thomas/ica_out/8mm/spatio-temporal_regression/PET_fMRI_intersection_w_pet_values_on_pet/aDMN_RSN'

% PET subjects input data
load('/home/alexandre/data/thomas/ica_out/8mm/gift_pet_grptemplate_no-pvc_30ICs/Subject.mat')

icatb_spatial_temp_regress(compFiles, struct2cell(files), 'outputDir', outputDir, 'format', '.nii')

exit
