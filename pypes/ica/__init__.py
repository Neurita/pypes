
from .decompose import attach_canica, attach_concat_canica

from .plotting import ( plot_connectivity_matrix,
                        plot_ica_results,
                        ICAResultsPlotter,
                        CanICAResultsPlotter,
                        MIALABICAResultsPlotter,
                        GIFTICAResultsPlotter,
                        GIFTGroupICAResultsPlotter,
                        SBMICAResultsPlotter,)

from .utils import get_largest_blobs

from .loadings import (filter_ics,
                       add_groups_to_loadings_table,
                       build_raw_loadings_table)

from .rsn_compare import (RestingStateNetworks,
                          spatial_maps_pairwise_similarity,
                          spatial_maps_goodness_of_fit,)