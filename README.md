The initial commit of this code is the current version of code used in Lexy von Diezmann, Chloe Bristow, and Ofer Rog, "Diffusion within the synaptonemal complex can account for signal transduction along meiotic chromosomes,"
https://www.biorxiv.org/content/10.1101/2024.05.22.595404v1 with minor edits for clarity.

A repository for code used to simulate single-molecule diffusion in a confined space. This code works in one dimension.

brownianMotionInSC_generateTracks.m generates trajectories based on user-specified parameters.
This can used to simulate escape from an initial region, or capture at a small region, by setting a manual option (1 or 2).
This code can accept a set of empirical diffusion coefficients formatted as a .mat.

brownianMotionInSC_plotTracks.m selects individual trajectories for plotting and generates images of kymographs.

We have also included a compiled set of kymographs, and a short script that plots an overlay of derived quantities (exit and entrance rates into subregions).
