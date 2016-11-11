# Uplift
Matlab software for analysing particle displacement in TIF data from 
CET IIB project on particle uplift. 

EJR 2016 

Mainly rough work so far.

## tracking
This folder contains a set of GPL particle tracking software for Matlab, based
on Crocker + Grier's method, modified by DLB et al.
It is used by the scripts below

## trackuplift.m
This script applies particle finding, tracking, and analysis to a sample data file. 
It is hard-coded to analyse regions in the video file 'Exp1.MP4' and would need 
modifications to work on other data.

## trackuplift_v2.m
This script improves the robustness of particle finding + tracking, and 
generates overlays etc. in a cleaner way. It also allows (and requires, by default) the 
user to define a rectangular region of interest, so it can be applied to more data.
The particle tracking is adjusted to work with image data with not red lines on the 
glass. 

## timeresolved_uplift
Post-processes data in workspace from trackuplift_v2, and produces 
time-resolved visualisations of the analysis.

## timeresolved_distortion.m 
Post-processes data in workspace from trackuplift_v2, and produces 
time-resolved visualisations of the distortion of a square mesh.