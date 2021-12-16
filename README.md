# EWStimeCorrelatedNoise
This repository contains the code used to generate the results in "Warning Signs for Non-Markovian Bifurcations: Color Blindness and Scaling Laws", which is available at https://arxiv.org/abs/2106.08374. Please cite the preprint when using the material in this repository.

The code is written in MATLAB R2019a/b. Depending on the chosen noise (BM, fBM, color or Rosenblatt) by e.g. setting "noise = 'fBM'", sample paths are generated and the reduced Stommel Cessi SDE driven by the specified noise is solved numerically. Sample paths move along the corresponding critical manifold. Furthermore, the scaling law is calculated and shown in a log-log-plot.

Please note that, due to memory allocation issues, the settings for the Rosenblatt process are slightly different. Uncomment and comment corresponding parameter choices as indicated in Section "Climate Tipping Numerics". Furthermore, note that also for the fractional Brownian motion requested arrays may exceed maximum array size preference depending on you machine. In this case, you may switch to the parameter choices for the Rosenblatt process as these require less memory allocation.

Below you find a rough estimate of how long you could expect the code to run for each of the noise types. Note however that the code has not been tuned for performance and that the computing time heavily depends on the machine (rough estimates for a Desktop PC (64 Bit-operating system Windows 7, Intel® Core ™ i7-3770 CPU with 3.4 GHz, 16GB RAM))
<ul>
      <li> Brownian motion: 1 minute </li>
      <li> Colored noise process:  </li>
      <li> Fractional Brownian motion: </li>
       <li> Rosenblatt process: </li>
    </ul>
