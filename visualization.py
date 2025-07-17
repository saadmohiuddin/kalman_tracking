import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch
from src import TrackGenerator

class InnerDetector:
    def __init__(self, length_mm=6200, diameter_mm=2100):
        self.length = length_mm
        self.diameter = diameter_mm
        self.radius = self.diameter / 2
        self.pixel_layers = (33.25, 50.5, 88.5, 122.5)  # mm
        # On ne conserve que les pixel layers, pas les SCT + TRT
        self.radii = list(self.pixel_layers)
        self.subsystems = [
            {'name': 'Pixel Detector', 'layers': range(0, 4), 'color': 'mediumorchid', 'alpha': 0.2}
        ]
        self.LONGITUDINAL_STEPS = 1000
        self.ANGULAR_STEPS = 500
        self._prepare_mesh()

    def _prepare_mesh(self):
        z = np.linspace(-self.length / 2, self.length / 2, self.LONGITUDINAL_STEPS)
        theta = np.linspace(0, 2 * np.pi, self.ANGULAR_STEPS)
        self.Z, self.Theta = np.meshgrid(z, theta)

    def draw(self, ax):
        for subsystem in self.subsystems:
            for i in subsystem['layers']:
                r = self.radii[i]
                X = r * np.sin(self.Theta)
                Y = r * np.cos(self.Theta)
                ax.plot_surface(self.Z, Y, X, color=subsystem['color'], alpha=subsystem['alpha'], edgecolor='none')

    def get_legend(self):
        return [
            Patch(facecolor='mediumorchid', edgecolor='mediumorchid', label='Pixel Detector')
        ]

class DetectorSimulation:
    def __init__(self):
        self.layer = InnerDetector()
        self.fig = plt.figure(figsize=(14, 12))
        self.ax = self.fig.add_subplot(111, projection='3d')

    def _setup_axes(self):
        axis_length = self.layer.radius * 1.2
        self.ax.scatter(0, 0, 0, color='black', s=5, label='Origin (0,0,0)')
        self.ax.quiver(0, 0, 0, axis_length, 0, 0, color='black', linewidth=1, arrow_length_ratio=0.05)
        self.ax.quiver(0, 0, 0, 0, axis_length, 0, color='black', linewidth=1, arrow_length_ratio=0.05)
        self.ax.quiver(0, 0, 0, 0, 0, axis_length, color='black', linewidth=1, arrow_length_ratio=0.05)
        self.ax.text(axis_length * 1.1, 0, 0, 'z', color='black', fontsize=14)
        self.ax.text(0, axis_length * 1.05, 0, 'x', color='black', fontsize=14)
        self.ax.text(0, 0, axis_length * 1.05, 'y', color='black', fontsize=14)
        self.ax.set_xlabel('Z (mm)', labelpad=10)
        self.ax.set_ylabel('X (mm)', labelpad=10)
        self.ax.set_zlabel('Y (mm)', labelpad=3)
        self.ax.set_box_aspect([
            self.layer.length,
            self.layer.radius * 2,
            self.layer.radius * 2
        ])
        self.ax.view_init(elev=10, azim=45)
        self.ax.set_title('ATLAS Pixel Detector with Track', fontsize=16)

    def draw_static_with_track(self, hits_dict, filename="Pixel_Detector_with_Track.png"):
        self.layer.draw(self.ax)
        self._setup_axes()

        # Hits
        hit_coords = list(hits_dict.values())
        xs, ys, zs = zip(*hit_coords)

        # Tracé amélioré de la trajectoire (ligne rouge)
        self.ax.plot(xs, ys, zs, color='crimson', linewidth=3, label='Particle Track')

        # Mise à jour de la légende
        legend_elements = self.layer.get_legend()
        legend_elements.append(Patch(facecolor='crimson', edgecolor='crimson', label='Particle Track'))

        self.ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0.05, 0.95))

        # Amélioration visuelle
        self.ax.set_facecolor('white')
        self.ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
        self.ax.view_init(elev=20, azim=30)

        plt.tight_layout()
        plt.savefig(filename, dpi=300)
        print(f"Figure saved as: {filename}")


if __name__ == "__main__":
    print("-" * 25, "Start", "-" * 25)
    initial_guess = {"x": 0, "y": 0, "z": 0, "px": 30, "py": 200, "pz": 0, "charge": -1}
    track = TrackGenerator.TrackGenerator(**initial_guess)

    reachability = track.check_layer_reachability()
    if not reachability["all_reachable"]:
        print("Track does not reach all layers. Generating a valid track...")
        track.generate_track_with_all_hits(max_attempts=100)
        reachability = track.check_layer_reachability()

    print(f"Reaches all Layers: {reachability['all_reachable']}")

    # Obtenir les hits
    hits = track.find_layer_intersections()

    # Dessiner
    sim = DetectorSimulation()
    sim.draw_static_with_track(hits, filename="Pixel_Detector_with_Track2A.png")
    print("-" * 25, "End", "-" * 25)
