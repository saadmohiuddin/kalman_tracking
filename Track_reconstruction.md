## Track Reconstruction



This project aims to reconstruct a toy model of track reconstruction at High energy physics experiment. We are taking ATLAS as our inspiration for building this algorithm. We divided this work into 4 parts: Detector modelling, track generation and recording hit points on detector, track reconstruction using track fitting algorithms, and visualization. My part in the project was to write a track fitting algorithm. Track reconstruction at ATLAS utilizes an extended Kalman Filter to estimate track parameters in the presence of uncertainties on track hits. This algorithm relies on a seeded track.

### How to run:
The code is implemented in src/track_reconstruction. It takes as input the ```hits_layers``` coordinates and reconstructs a track trajectory from that. This is done in two parts.

1. Radius of curvature:
To calculate the radius of curvature, we use a least squares method to get the coordinates of the centre of the circle and its radius.

We fit a circle equation: (x-h)² + (y-k)² = r²

Expanding: x² - 2hx + h² + y² - 2ky + k² = r²

Rearranging: x² + y² = 2hx + 2ky + (r² - h² - k²)

So our system is:

A = [2x, 2y, 1] (coefficients)
b = x² + y² (target)
sol = [h, k, r² - h² - k²] (unknowns)
Where h, k are the circle center and the third term helps calculate the radius.

the function ```fit_circle_least_squares(points)``` takes as input the hit coordinates and outputs the radius of the fitted circle.

2. Track reconstruction usign 4th order Runge-Kutta method.
We use first hit point as the initial position and the normalized displacement between the first and second hit points mulitplied by the charge q, the magnetic field B, and the radius of curvature calculated in the first part. This follows from the relationshiop between the transverse momentum and the radius of curvature in a particle trajectory in a uniform B field $|p_T| =q*B*R$. The Lorentz force is calculated as the cross product $\frac{dp}{dt} = q\vec{v} \times \vec{B}$ The trajectory is updated backwards in time so it can build towards the beginning of the trajectory





The ```recreate_track(hits_final, B=2.0, q=q)``` function implements both the least squares radius calculation and the trajectory reconstruction so the code can be run as follows: 
```python
from src import track_reconstruction as tr


hits_final = np.array([hits_layer1, hits_layer2, hits_layer3, hits_layer4])
print("Hits:", hits_final)
# Simulate trajectory
x = hits_final[0] # initial position
p = tr.estimate_initial_momentum(hits_final, B=2.0, charge=q)  # initial momentum (GeV/c)
hits_final = np.array([hits_layer1, hits_layer2, hits_layer3, hits_layer4])
trajectory = tr.recreate_track(hits_final, B=2.0, q=q)
print(np.array(trajectory)[:, 0], np.array(trajectory)[:, 1],np.array(trajectory)[:, 2])
x = np.array(trajectory)[:, 0]
y = np.array(trajectory)[:, 1]
z = np.array(trajectory)[:, 2]
```