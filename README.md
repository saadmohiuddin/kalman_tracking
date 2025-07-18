# Kalman Filter Track Reconstruction

This project aims to simulate charged particle track reconstruction as it happens in ATLAS.

## Introduction

The Inner detectors at ATLAS and CMS aim to reconstruct charged particle trajectories using hits on pixel detectors and matching those hits using a combinatorial Kalman filter to reconstruct the track paramaters \cite{trackseeding}. A track is defined by 5 parameters, namely: The transverse and longitudinal distances of closest approach to a reference point(d0 and z0), the polar and azimuthal angle  of the momentum vector at the reference point, and the charge to momentum ratio (q/p)

## Implementation

For this project, we aim to write a python package that recreates a toy model of ATLAS track reconstruction. We will use numpy to create a 3d meshgrid of cylindrical coordinates to act as our detector
model, and simulate charged particle trajectories that register hits on them. Then we will use a Kalman filter algorithm that we will either write ourselves, or use a package, to reconstruct the tracks using the data from hits. We will use object oriented programming design in python to create classes for tracks, the detector, and for hits on the detector. We will use matplotlib to help visualize our tracks and see how the reconstructed tracks from the Kalman filter lines up with the original track.

## Steps to run the analysis
This repository is meant as a group project for SE4Sci Summer 2025 course. In this project, we try to reconstruct the track of charged particle as is done by the Inner Detector at ATLAS.
However, everything that we use here is built by group members on pure python. We have not used any of the already advanced libraries(ex Kalman filter algorithm) or ACTS(ATLAS Common Tracking Software). Thus, almost all content in this project is preliminary as we have tried to use what we learnt from the class with some tracking basics. The main steps that follow in this analysis are as below.

TLDR - Our goal is to reconstruct track of a charged particles based on the intersection of particle and pixel layers(hits).
## Generate Truth Track
First, we start with generating a 'Truth' track which has some granularity(higher than 4 layers of pixel). This part is handled by the file - TrackGenerator.py. This file can only be used when imported as module and not as a standalone script. This file can be interpreted as something that generates Monte Carlo tracks.
```python
from src import TrackGenerator
```
Currently, this works only on branch dev_suman, but we will eventually merge all four branches and use it.
After importing, you should give some initial parameters for the track.
```python
initial_guess = {"x": 2, "y": 1, "z": 0, "px": 30, "py": 40, "pz": 60, "charge": -1}
```

Now, create an instance of TrackGenerator class,
```python
track = TrackGenerator.TrackGenerator(**initial_guess)
```

This generates a track. But we need to make sure if this track hits all four layers of Pixel. This is needed because we need the intersection of track with layers(4 points) for Kalman filter which Saad is implementing. To check the track reachability,
```python
reachability = track.check_layer_reachability()
print(f"Reaches all Layer : {reachability['all_reachable']}")
```
This prints, either True or False. If True, we go ahead with this Track, if False we generate another track by varying its momentum or position coordinates.

### If the track does not reach all layers
```python
if not reachability["all_reachable"]:
    print("-" * 25, "Generating Track that reaches all layer", "-" * 25)
    # track with 100 points of X/Y/Z coordinates#
    hits = track.generate_track_with_all_hits(max_attempts=100)
    print(track)
    reachability = track.check_layer_reachability()
```

Now this track reaches all layers.

The next step is to get the (X,Y,Z) coordinates of the hit.

```python
# this track hits all layers
radius = (33.25, 50.5, 88.5, 122.5)  # radius of pixel layer
hits = track.find_layer_intersections()
hits_layer1 = hits[radius[0]]
hits_layer2 = hits[radius[1]]
hits_layer3 = hits[radius[2]]
hits_layer4 = hits[radius[3]]
print("-" * 50)
print(f"Reaches all Layer : {reachability['all_reachable']}")
print("-" * 50)
print("XYZ coordinate of hit in layer 1", hits_layer1)  # in mm
print("XYZ coordinate of hit in layer 2", hits_layer2)
print("XYZ coordinate of hit in layer 3", hits_layer3)
print("XYZ coordinate of hit in layer 4", hits_layer4)
```

The above command prints this in my terminal
```bash
XYZ coordinate of hit in layer 1 (np.float64(-22.741508113381812), np.float64(-24.256555761727547), np.float64(-7.894214353368604))
XYZ coordinate of hit in layer 2 (np.float64(-37.212922095389395), np.float64(-34.138661197197486), np.float64(-11.79834850637275))
XYZ coordinate of hit in layer 3 (np.float64(-73.68594825338614), np.float64(-49.01634744957389), np.float64(-20.59963779225072))
XYZ coordinate of hit in layer 4 (np.float64(-110.52897284511363), np.float64(-52.81619028780778), np.float64(-28.871911437826135))
```

Next, Saad uses this points for his Kalman Filter algorithm input. The momentum and charge are not needed as they can be reconstructed using Kalman filter from only position coordinates.

### Future Improvements
There are some areas of improvements in Track Generating, for example, automatically generate hundreds of tracks, more robust track parameter changing, etc. Given the short deadline, all of them can not be implemented.
