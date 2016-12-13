#
#
# def dual_regression(y, X):
#     """ Use dual regression approach to compute the spatial maps
#     and time courses.
#
#     Parameters
#     ----------
#     y: np.ndarray
#         Observations in columns (Voxels by timepoints)
#
#     X: np.ndarray
#         design matrix (Voxels by components)
#
#     Returns
#     -------
#     tc: np.ndarray
#         Time courses (Timepoints by components)
#
#     spatial_maps: np.ndarray
#         Spatial maps (Components by voxels)
#     """
#     X =
#
# def spatial_maps_pairwise_regression(imgs1, imgs2, mask_file, method_num=0):
#     """ TODO!!
#
#     Look in icat_dual_regress:
#
#     function [tc, spatial_maps] = icatb_dual_regress(y, X)
#     %% Use dual regression approach to compute the spatial maps and time
#     % courses
#     %
#     % Inputs:
#     % 1. y - Observations in columns (Voxels by time points)
#     % 2. X - design matrix (Voxels by components)
#     %
#     % Outputs:
#     % 1. tc - Time courses (Timepoints by components)
#     % 2. spatial_maps - Spatial maps (Components by voxels)
#     %
#
#     %% First step. Fit model matrix to the data to get time courses.
#     X = icatb_remove_mean(X);
#
#     tc = pinv(X)*icatb_remove_mean(y);
#     tc = tc';
#     clear X;
#
#     % Store mean of timecourses
#     mean_tc = mean(tc);
#
#     % Remove mean of timecourse
#     tc = icatb_remove_mean(tc);
#
#     %% Second step. Fit Time courses at each voxel to get the spatial maps.
#     try
#
#         spatial_maps = pinv(tc)*icatb_remove_mean(y');
#
#     catch
#         %% Use less memory usage way to do regression
#
#         % Initialise spatial maps
#         spatial_maps = zeros(size(tc, 2), size(y, 1));
#
#
#         % Loop over voxels
#         for nVoxel = 1:size(y, 1)
#
#             bold_signal = detrend(y(nVoxel, :), 0);
#
#             spatial_maps(:, nVoxel) = pinv(tc)*bold_signal(:);
#
#             clear bold_signal;
#
#         end
#         % End of loop over voxels
#
#     end
#
#     clear y;
#
#     %% Add mean back to the timecourses
#     tc = tc + (repmat(mean_tc, size(tc, 1), 1));
#
#
#
#     Multiple Regression values of each RSN to each IC map in `ic_imgs` masked by `mask_file`.
#
#     Parameters
#     ----------
#     imgs1: list of niimg-like or 4D niimg-like
#
#     imgs2: list of niimg-like or 4D niimg-like
#
#     mask_file: niimg-like
#
#     method_name: str
#         Valid values for method_name are:
#         From scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan'].
#         From scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming',
#                                       'jaccard', 'kulsinski', 'mahalanobis', 'matching', 'minkowski',
#                                       'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
#                                       'sqeuclidean', 'yule']
#                                       See the documentation for scipy.spatial.distance for details on these metrics.
#
#     Returns
#     -------
#     corrs: np.ndarray
#         A matrix of shape MxN, where M is len(rsn_imgs) and N is len(ic_imgs).
#         It contains the correlation values.
#     """
#     rsn_img = niimg.load_img(imgs1)
#     ic_img  = niimg.load_img(imgs2)
#
#     n_rsns = rsn_img.shape[-1]
#     n_ics  =  ic_img.shape[-1]
#     corrs = np.zeros((n_rsns, n_ics), dtype=float)
#
#     mask_trnsf = niimg.resample_to_img(mask_file, niimg.index_img(ic_img, 0),
#                                        interpolation='nearest',
#                                        copy=True)
#
#     for rsn_idx, rsn in enumerate(niimg.iter_img(rsn_img)):
#         rsn_transf = niimg.resample_to_img(rsn, niimg.index_img(ic_img, 0), copy=True)
#         rsn_masked = nimask.apply_mask(rsn_transf, mask_trnsf)
#
#         for ic_idx, ic in enumerate(niimg.iter_img(ic_img)):
#             ic_masked = nimask.apply_mask(ic, mask_trnsf)
#             # dist = pairwise_distances(rsn_masked.reshape(1, -1),
#             #                            ic_masked.reshape(1, -1),
#             #                           metric=distance)
#
#             # -----------------------------------------------------------------------------------------------------
#             # CODE TO MODIFY
#
#             # http://stackoverflow.com/questions/17679140/multiple-linear-regression-with-python
#             # http://scikit-learn.org/stable/modules/linear_model.html#ordinary-least-squares
#             import pandas as pd
#             import numpy as np
#
#             path = 'DB2.csv'
#             data = pd.read_csv(path, header=None, delimiter=";")
#
#             data.insert(0, 'Ones', 1)
#             cols = data.shape[1]
#
#             X = data.iloc[:,0:cols-1]
#             y = data.iloc[:,cols-1:cols]
#
#             IdentitySize = X.shape[1]
#             IdentityMatrix= np.zeros((IdentitySize, IdentitySize))
#             np.fill_diagonal(IdentityMatrix, 1)
#
#             #For least squares method you use Numpy's numpy.linalg.lstsq. Here is Pyhton code
#             lamb = 1
#             th = np.linalg.lstsq(X.T.dot(X) + lamb * IdentityMatrix, X.T.dot(y))[0]
#
#             # Also you can use np.linalg.solve tool of numpy:
#             lamb = 1
#             XtX_lamb = X.T.dot(X) + lamb * IdentityMatrix
#             XtY = X.T.dot(y)
#             x = np.linalg.solve(XtX_lamb, XtY);
#
#             # For normal equation method use:
#             lamb = 1
#             xTx = X.T.dot(X) + lamb * IdentityMatrix
#             XtX = np.linalg.inv(xTx)
#             XtX_xT = XtX.dot(X.T)
#             theta = XtX_xT.dot(y)
#
#             # OLS
#             # >>> from sklearn import linear_model
#             # >>> reg = linear_model.LinearRegression()
#             # >>> reg.fit ([[0, 0], [1, 1], [2, 2]], [0, 1, 2])
#             # LinearRegression(copy_X=True, fit_intercept=True, n_jobs=1, normalize=False)
#             # >>> reg.coef_
#             # array([ 0.5,  0.5])
#
#             # Ridge Regression
#             # >>> from sklearn import linear_model
#             # >>> reg = linear_model.Ridge (alpha = .5)
#             # >>> reg.fit ([[0, 0], [0, 0], [1, 1]], [0, .1, 1])
#             # Ridge(alpha=0.5, copy_X=True, fit_intercept=True, max_iter=None,
#             #       normalize=False, random_state=None, solver='auto', tol=0.001)
#             # >>> reg.coef_
#             # array([ 0.34545455,  0.34545455])
#             # >>> reg.intercept_
#             # 0.13636...
#
#             # CODE TO MODIFY
#             # -----------------------------------------------------------------------------------------------------
#
#             # since this is a scalar value
#             dist = dist[0][0]
#
#             # since this is a distance based on correlation, not a correlation value
#             corr = 1 - dist
#
#             # store it
#             corrs[rsn_idx, ic_idx] = corr
#
#     return corrs

