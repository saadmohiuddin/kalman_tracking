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

# 3. Get the full particle trajectory
x_vals, y_vals, z_vals, _, _ = track.evolve_track(time_steps=5000)

# 4. Initialize detector simulation
sim = DetectorSimulation()

# 5. Plot all layers of the inner detector (Pixel, SCT, TRT)
z = np.linspace(-sim.layer.length / 2, sim.layer.length / 2, sim.layer.LONGITUDINAL_STEPS)
theta = np.linspace(0, 2 * np.pi, sim.layer.ANGULAR_STEPS)
Z, Theta = np.meshgrid(z, theta)

# Define colors for each sub-detector group
colors = {
    "Pixel": "mediumorchid",
    "SCT": "royalblue",
    "TRT": "royalblue"
}

# Draw all cylindrical layers
for i, r in enumerate(sim.layer.radii):
    if i < 4:
        color = colors["Pixel"]
    elif i < 8:
        color = colors["SCT"]
    else:
        color = colors["TRT"]
    
    X = r * np.sin(Theta)
    Y = r * np.cos(Theta)
    sim.ax.plot_surface(Z, Y, X, color=color, alpha=0.2, edgecolor='none')

# 6. Set up 3D axes
sim._setup_axes()

# 7. Plot the particle trajectory
sim.ax.plot(z_vals, x_vals, y_vals, color='red', linewidth=2, label='Particle Trajectory')

# 8. Legend
legend_patches = [
    Patch(facecolor=colors["Pixel"], edgecolor=colors["Pixel"], label='Pixel Detector'),
    Patch(facecolor=colors["SCT"], edgecolor=colors["SCT"], label='SCT Detector'),
    Patch(facecolor=colors["TRT"], edgecolor=colors["TRT"], label='TRT Detector'),
    Patch(color='red', label='Particle Trajectory')
]
sim.ax.legend(handles=legend_patches, loc='upper left', bbox_to_anchor=(0.05, 0.95))

# 9. Show and save the plot
plt.tight_layout()
plt.savefig("char_particle_tracks.png", dpi=300)
plt.show()
