# EWStimeCorrelatedNoise
This repository contains the code used to generate the results in "Warning Signs for Non-Markovian Bifurcations: Color Blindness and Scaling Laws", which is available at https://arxiv.org/abs/2106.08374. Please cite the preprint when using the material in this repository.

The code is written in MATLAB R2019a/b. Depending on the chosen noise (BM, fBM, color or Rosenblatt) by e.g. setting "noise = 'fBM'", sample paths are generated and the reduced Stommel Cessi SDE driven by the specified noise is solved numerically. Sample paths are drawn along the corresponding critical manifold. Furthermore, the scaling law is calculated and shown in a log-log-plot.

Please note that, due to memory allocation issues, the settings for the Rosenblatt process are slightly different. Uncomment and comment corresponding parameter choices as indicated in Section "Climate Tipping Numerics".
