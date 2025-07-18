# Kalman Filter Track Reconstruction

This project aims to simulate charged particle track reconstruction as it happens in ATLAS.

## Introduction

The Inner detectors at ATLAS and CMS aim to reconstruct charged particle trajectories using hits on pixel detectors and matching those hits using a combinatorial Kalman filter to reconstruct the track paramaters \cite{trackseeding}. A track is defined by 5 parameters, namely: The transverse and longitudinal distances of closest approach to a reference point(d0 and z0), the polar and azimuthal angle  of the momentum vector at the reference point, and the charge to momentum ratio (q/p)

## Implementation

For this project, we aim to write a python package that recreates a toy model of ATLAS track reconstruction. We will use numpy to create a 3d meshgrid of cylindrical coordinates to act as our detector
model, and simulate charged particle trajectories that register hits on them. Then we will use a Kalman filter algorithm that we will either write ourselves, or use a package, to reconstruct the tracks using the data from hits. We will use object oriented programming design in python to create classes for tracks, the detector, and for hits on the detector. We will use matplotlib to help visualize our tracks and see how the reconstructed tracks from the Kalman filter lines up with the original track.

# Inner Detector Simulation Module

### 📄 File: `src/Inner_Detector_simulation.py`

This module provides a detailed 3D simulation of the **ATLAS Inner Detector**, including the **Pixel Detector**, the **Semiconductor Tracker (SCT)**, and the **Transition Radiation Tracker (TRT)**.

It is composed of:

- `InnerDetector`: defines the geometry of the Inner Detector.
- `DetectorSimulation`: renders and saves the visualization using `matplotlib`.

---

## Features

- **Pixel Detector** layers at: `33.25 mm`, `50.5 mm`, `88.5 mm`, `122.5 mm`
- **SCT/TRT** layers: automatically generated from `200 mm` to `1050 mm`
- 3D rendering using `matplotlib` with labeled axes and legend
- Can be run as a standalone script or imported by other code (e.g. by Aissata)

---

## Requirements

Install dependencies with:

```bash
pip install matplotlib numpy


This project aims to simulate a simplified version of the ATLAS Inner Detector and generate charged particle trajectories in a uniform magnetic field. Hits are registered where particle helices intersect detector layers. These hits form the basis for track reconstruction algorithms like the Kalman Filter.

## Implementation of the Inner Detector Simulation

Geometry: We construct a 3D cylindrical meshgrid representing 8 concentric radial detector layers, each spanning 1000 bins in the longitudinal ($z$) direction and 500 bins in azimuthal angle ($\phi$).
##Generate and save a static 3D visualization
from src.Inner_Detector_simulation import DetectorSimulation

sim = DetectorSimulation()
sim.draw_static("Inner_Detector.png")

This will save the detector figure as Inner_Detector.png
##Import the detector in another script (Aissata's use case)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.Inner_Detector_simulation import InnerDetector

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

detector = InnerDetector()
detector.draw(ax)

ax.legend(handles=detector.get_legend(), loc='upper left')
plt.title("ATLAS Inner Detector Visualization")
plt.show()
This allows other team members to integrate the detector into larger figures or animations.
