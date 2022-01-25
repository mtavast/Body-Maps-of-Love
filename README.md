# Body Maps of Love
This repository contains analysis code for the Body Maps of Love -project.

Due to ethical restrictions, we cannot share the individual level data. Thus, most of the results of the manuscript can't be reproduced directly using the files available in this repository. However, from the code, one can review how the results of the manuscript were produced. 

Ethical restrictions do not apply to group-level data, and they are made available (output/body/bspm_ttest.mat, output/mind/mind_mds.mat, and output/sim/sim_mds.mat). These can be used to reproduce, for example, the clustering analyses and the Figures 1, 3 and 4 of the manuscript. 

mind_mds.mat -file contains the means and medians of the dimension ratings (based on both the raw ratings and Z transformed ratings), sim_mds.mat -file contains the mean distance matrix from experiment 3, and the bspm_ttest.mat -file contains the group level body maps.

Most of the code used here was repurposed from a previous project (https://version.aalto.fi/gitlab/eglerean/sensations/-/tree/master/) following MIT license. Thus, if you reuse the code, please cite:

Nummenmaa, L., Hari, R., Hietanen, J. K., & Glerean, E. (2018). Maps of subjective feelings. Proceedings of the National Academy of Sciences, 115(37), 9198-9203. https://doi.org/10.1073/pnas.1807390115
