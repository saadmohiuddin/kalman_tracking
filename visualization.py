from src import TrackGenerator
from src.Inner_Detector_simulation import DetectorSimulation
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np

# 1. Define initial particle parameters
initial_guess = {"x": 2, "y": 1, "z": 0, "px": 30, "py": 40, "pz": 60, "charge": -1}
track = TrackGenerator.TrackGenerator(**initial_guess)

# 2. Check if the particle reaches all layers
reachability = track.check_layer_reachability()
if not reachability["all_reachable"]:
    track.generate_track_with_all_hits()

# 3. Get the full trajectory (X, Y, Z points)
x_vals, y_vals, z_vals, _, _ = track.evolve_track(time_steps=345)

# 4. Initialize the detector simulation
sim = DetectorSimulation()

# 5. Plot only the 4 pixel layers as transparent cylinders
z = np.linspace(-sim.layer.length / 2, sim.layer.length / 2, sim.layer.LONGITUDINAL_STEPS)
theta = np.linspace(0, 2 * np.pi, sim.layer.ANGULAR_STEPS)
Z, Theta = np.meshgrid(z, theta)

for i in range(4):  # Only the 4 pixel layers
    r = sim.layer.radii[i]
    X = r * np.sin(Theta)
    Y = r * np.cos(Theta)
    sim.ax.plot_surface(Z, Y, X, color='mediumorchid', alpha=0.2, edgecolor='none')

# 6. Set up the 3D axes
sim._setup_axes()

# 7. Plot the particle trajectory
sim.ax.plot(z_vals, x_vals, y_vals, color='red', linewidth=2, label='Particle Trajectory')

# 8. Add legend
legend_patches = [
    Patch(facecolor='mediumorchid', edgecolor='mediumorchid', label='Pixel Detector'),
    Patch(color='red', label='Particle Trajectory')
]
sim.ax.legend(handles=legend_patches, loc='upper left', bbox_to_anchor=(0.05, 0.95))

# 9. Display and save the figure
plt.tight_layout()
plt.savefig("trajectory.png", dpi=300)
plt.show()
