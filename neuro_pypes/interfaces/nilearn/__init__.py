from .image import (concat_imgs,
                    math_img,
                    mean_img,
                    smooth_img,
                    copy_header,
                    resample,
                    resample_to_img,
                    ni2file)

from .canica import CanICAInterface

from .connectivity import ConnectivityCorrelationInterface

from .plot import (plot_all_components,
                   plot_ica_components,
                   plot_multi_slices,
                   plot_ortho_slices)

from .roi import spread_labels