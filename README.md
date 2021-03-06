# polarized_interneuron
Code and dataset to calculate interneuron degrees and second order motifs observed in the presubiculum.
Two different .mlx MATLAB livescripts (interneuron_degree.mlx, second_order_motif.mlx) allow calculation of interneuron degrees and second order motifs. PDFs of livescripts show code and plotted figures. The respective .m files show the same code in binary format.

interneuron_degree.mlx and second_order_motif.mlx import the cell table (tcell) and the connection table (tconnection).
second_order_motif.mlx allows the selection of interneuron subtype and motif configuration via dropdown menus. See the comments and text in the livescripts or PDF for further details.

overlap_connectivity.Rmd is a R notebook which imports IN_overlap_connect.xlsx to analyse overlap scores for FS/NFS pairs.

polarity_eccentricity.xlsx contains the anatomically calculated data used for comparison of polarity and eccentricity of PrS/mEC FS/NFS neurons.
