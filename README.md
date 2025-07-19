# Kalman Filter Track Reconstruction

This project aims to simulate charged particle track reconstruction as it happens in ATLAS.

## Introduction

The Inner detectors at ATLAS and CMS aim to reconstruct charged particle trajectories using hits on pixel detectors and matching those hits using a combinatorial Kalman filter to reconstruct the track paramaters \cite{trackseeding}. A track is defined by 5 parameters, namely: The transverse and longitudinal distances of closest approach to a reference point(d0 and z0), the polar and azimuthal angle  of the momentum vector at the reference point, and the charge to momentum ratio (q/p)

## Implementation

For this project, we aim to write a python package that recreates a toy model of ATLAS track reconstruction. We will use numpy to create a 3d meshgrid of cylindrical coordinates to act as our detector
model, and simulate charged particle trajectories that register hits on them. Then we will use a Kalman filter algorithm that we will either write ourselves, or use a package, to reconstruct the tracks using the data from hits. We will use object oriented programming design in python to create classes for tracks, the detector, and for hits on the detector. We will use matplotlib to help visualize our tracks and see how the reconstructed tracks from the Kalman filter lines up with the original track.


## Visualization of Particle Trajectories
📄 File:particules_tracks.py
As the final step of the project, we implemented a 3D visualization of a simulated charged particle traversing the layers of the ATLAS Inner Detector. This script integrates the full simulation pipeline: particle initialization, trajectory generation, and visualization on the detector model.

📦 Modules Used
```python
from src import TrackGenerator
Imports the class that handles particle initialization and track propagation under magnetic field influence.
```
``` python
from src.Inner_Detector_simulation import DetectorSimulation
#Imports the class that simulates and renders the ATLAS Inner Detector geometry in 3D.
```
These modules were developed earlier in the project and are reused here to build the final, integrated visualization.

Summary of Key Steps
1. Initialize the Particle
A charged particle is initialized with position and momentum components, plus electric charge:

```python
initial_guess = {"x": 2, "y": 1, "z": 0, "px": 30, "py": 40, "pz": 60, "charge": -1}
track = TrackGenerator.TrackGenerator(**initial_guess)
```
2. Check Detector Reachability
The simulation checks whether the particle reaches all detector layers. If not, the trajectory is adjusted:

```python
reachability = track.check_layer_reachability()
if not reachability["all_reachable"]:
    track.generate_track_with_all_hits()
```
3. Simulate the Trajectory
The particle's motion in 3D is computed over 5000 time steps:

```python
x_vals, y_vals, z_vals, _, _ = track.evolve_track(time_steps=5000)
```
4. Render the Detector
We generate and visualize the cylindrical geometry of the Inner Detector:

Pixel layers (colored in purple)

SCT and TRT layers (colored in blue)

5. Overlay the Particle Track
The computed trajectory is plotted as a red line on top of the detector geometry:

```python
sim.ax.plot(z_vals, x_vals, y_vals, color='red', linewidth=2, label='Particle Trajectory')
```
6. Save and Display the Figure
The complete visualization is saved as a high-resolution image and displayed:

```python
plt.savefig("char_partucle_tracks.png", dpi=300)
plt.show()
```
## Output
The result is a comprehensive 3D figure that shows the particle's full path as it moves through the layers of the ATLAS Inner Detector. This final visualization integrates both physics modeling and software design, and offers a strong validation tool for our simulation.
![OutputFiger](char_partucle_tracks.png)
