from neuro_pypes.interfaces.nilearn.canica import CanICAInterface
from neuro_pypes.interfaces.nilearn.connectivity import ConnectivityCorrelationInterface
from neuro_pypes.interfaces.nilearn.image import (
    concat_imgs,
    math_img,
    mean_img,
    smooth_img,
    copy_header,
    resample,
    resample_to_img,
    ni2file)
from neuro_pypes.interfaces.nilearn.plot import (
    plot_all_components,
    plot_ica_components,
    plot_multi_slices,
    plot_ortho_slices
)
from neuro_pypes.interfaces.nilearn.roi import spread_labels
