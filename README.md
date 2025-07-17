# Kalman Filter Track Reconstruction

This project aims to simulate charged particle track reconstruction as it happens in ATLAS.

## Introduction

The Inner detectors at ATLAS and CMS aim to reconstruct charged particle trajectories using hits on pixel detectors and matching those hits using a combinatorial Kalman filter to reconstruct the track paramaters \cite{trackseeding}. A track is defined by 5 parameters, namely: The transverse and longitudinal distances of closest approach to a reference point(d0 and z0), the polar and azimuthal angle  of the momentum vector at the reference point, and the charge to momentum ratio (q/p)

## Implementation

For this project, we aim to write a python package that recreates a toy model of ATLAS track reconstruction. We will use numpy to create a 3d meshgrid of cylindrical coordinates to act as our detector
model, and simulate charged particle trajectories that register hits on them. Then we will use a Kalman filter algorithm that we will either write ourselves, or use a package, to reconstruct the tracks using the data from hits. We will use object oriented programming design in python to create classes for tracks, the detector, and for hits on the detector. We will use matplotlib to help visualize our tracks and see how the reconstructed tracks from the Kalman filter lines up with the original track.
###Inner Detector Simulation
The ATLAS Inner Detector reconstructs charged particle trajectories by detecting hits on several concentric tracking layers — pixels, silicon strips (SCT), and transition radiation trackers (TRT). Tracks are characterized by parameters such as impact parameters ($d_0$, $z_0$), angles ($\theta$, $\phi$), and momentum over charge ($q/p$).

This project aims to simulate a simplified version of the ATLAS Inner Detector and generate charged particle trajectories in a uniform magnetic field. Hits are registered where particle helices intersect detector layers. These hits form the basis for track reconstruction algorithms like the Kalman Filter.
####Implementation
Geometry: We construct a 3D cylindrical meshgrid representing 8 concentric radial detector layers, each spanning 1000 bins in the longitudinal ($z$) direction and 500 bins in azimuthal angle ($\phi$).
