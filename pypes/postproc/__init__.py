
from .decompose import attach_canica, attach_concat_canica

from .plotting import ( plot_connectivity_matrix,
                        plot_ica_results,
                        ICAResultsPlotter,
                        CanICAResultsPlotter,
                        MIALABICAResultsPlotter,
                        GIFTICAResultsPlotter,
                        SBMICAResultsPlotter,)

from .ica_loadings import ( get_largest_blobs,
                            filter_ics,
                            add_groups_to_loadings_table,
                            build_raw_loadings_table)

from .rsn_compare import (RestingStateNetworks,
                          spatial_maps_pairwise_similarity,
                          spatial_maps_goodness_of_fit,)
