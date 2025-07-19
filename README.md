# Kalman Filter Track Reconstruction

This project aims to simulate charged particle track reconstruction as it happens in ATLAS.

## Introduction

The Inner detectors at ATLAS and CMS aim to reconstruct charged particle trajectories using hits on pixel detectors and matching those hits using a combinatorial Kalman filter to reconstruct the track parameters \cite{trackseeding}. A track is defined by 5 parameters: the transverse and longitudinal distances of closest approach to a reference point (d0 and z0), the polar and azimuthal angle of the momentum vector at the reference point, and the charge to momentum ratio (q/p).

## Implementation

We aim to write a Python package that recreates a toy model of ATLAS track reconstruction. We simulate charged particle trajectories through a 3D cylindrical detector and reconstruct the tracks using a Kalman filter. We use object-oriented design to define `Track`, `Detector`, and `Hit` classes, and `matplotlib` for visualization.

## Steps to Run the Analysis

This repository is meant as a group project for SE4Sci Summer 2025. The key steps in this project include:

### 1. Generate Truth Track

Tracks are generated with high granularity using `TrackGenerator.py`. This file must be used as a module.

```python
from src import TrackGenerator

initial_guess = {"x": 2, "y": 1, "z": 0, "px": 30, "py": 40, "pz": 60, "charge": -1}
track = TrackGenerator.TrackGenerator(**initial_guess)
2. Check Detector Reachability
Check whether the track hits all pixel layers:
reachability = track.check_layer_reachability()
print(f"Reaches all layers: {reachability['all_reachable']}")
If not, generate a new track:


if not reachability["all_reachable"]:
    hits = track.generate_track_with_all_hits(max_attempts=100)
3. Extract Hit Coordinates

hits = track.find_layer_intersections()
These coordinates are used by the Kalman filter algorithm.

Visualization of Particle Trajectories
📄 File: particules_tracks.py

This script performs 3D visualization of a charged particle inside the Inner Detector. It includes:

Initialization of the particle

Simulation of its trajectory (5000 time steps)

Rendering of the cylindrical detector

Overlay of the track on the detector

Sample Code

from src import TrackGenerator
from src.Inner_Detector_simulation import DetectorSimulation

initial_guess = {"x": 2, "y": 1, "z": 0, "px": 30, "py": 40, "pz": 60, "charge": -1}
track = TrackGenerator.TrackGenerator(**initial_guess)

reachability = track.check_layer_reachability()
if not reachability["all_reachable"]:
    track.generate_track_with_all_hits()

x_vals, y_vals, z_vals, _, _ = track.evolve_track(time_steps=5000)

sim = DetectorSimulation()
sim.setup_detector_geometry()
sim.plot_detector()

sim.ax.plot(z_vals, x_vals, y_vals, color='red', linewidth=2, label='Particle Trajectory')

plt.savefig("char_partucle_tracks.png", dpi=300)
plt.show()
Output
The result is a 3D image of the particle’s path across ATLAS’s pixel layers:

```


## Future Improvements
Automate track generation for larger datasets

Improve robustness of track parameters

Integrate Kalman filter fully with hit extraction

