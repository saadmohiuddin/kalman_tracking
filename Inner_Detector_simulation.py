import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch

class InnerDetectorVisualizer:
    def __init__(self):
        # Dimensions
        self.length = 6.2  # meters
        self.diameter = 2.1  # meters
        self.radius = self.diameter / 2
        
        # Resolution
        self.num_layers = 8
        self.z_steps = 1000
        self.theta_steps = 500

        # Grid
        self.z = np.linspace(-self.length / 2, self.length / 2, self.z_steps)
        self.theta = np.linspace(0, 2 * np.pi, self.theta_steps)
        self.Z, self.Theta = np.meshgrid(self.z, self.theta)

        # Subsystems
        self.subsystems = [
            {'name': 'Pixel Detector', 'layers': range(0, 4), 'color': 'mediumorchid', 'alpha': 0.1},
            {'name': 'SCT + TRT', 'layers': range(4, 8), 'color': 'dimgray', 'alpha': 0.1}
        ]


    def create_detector_layers(self, ax):
        radii = np.linspace(0, self.radius, self.num_layers)
        for subsystem in self.subsystems:
            for i in subsystem['layers']:
                r = radii[i]
                X = r * np.cos(self.Theta)
                Y = r * np.sin(self.Theta)
                ax.plot_surface(self.Z, Y, X,
                                color=subsystem['color'],
                                alpha=subsystem['alpha'],
                                edgecolor='none')


    def draw_axes(self, ax):
        origin_length = self.radius * 1.2

        # Origin
        ax.text(-0.07, -0.07, -0.12, '0', color='black', fontsize=12)
        ax.scatter(0, 0, 0, color='black', s=5)

        # Arrows
        ax.quiver(0, 0, 0, origin_length, 0, 0, color='black', linewidth=1, arrow_length_ratio=0.1)  
        ax.quiver(0, 0, 0, 0, origin_length, 0, color='black', linewidth=1, arrow_length_ratio=0.1)  
        ax.quiver(0, 0, 0, 0, 0, origin_length, color='black', linewidth=1, arrow_length_ratio=0.1) 

        # Labels
        ax.text(origin_length * 1.1, 0, 0, 'z', color='black', fontsize=14)
        ax.text(0, origin_length * 1.05, 0, 'x', color='black', fontsize=14)
        ax.text(0, 0, origin_length * 1.05, 'y', color='black', fontsize=14)
        ax.set_xlabel('Z (m)', labelpad=10)
        ax.set_ylabel('X (m)', labelpad=10)
        ax.set_zlabel('Y (m)', labelpad=3)


    def plot(self):
        fig = plt.figure(figsize=(14, 12))
        ax = fig.add_subplot(111, projection='3d')

        # Plot components
        self.create_detector_layers(ax)
        self.draw_axes(ax)

        # Aesthetics
        ax.set_box_aspect([self.length, self.radius * 2, self.radius * 2])
        ax.view_init(elev=10, azim=45)
        ax.set_title('ATLAS Inner Detector Structure', fontsize=16)

        # Legend
        legend_elements = [
            Patch(facecolor='mediumorchid', edgecolor='mediumorchid', label='Pixel Detector'),
            Patch(facecolor='dimgray', edgecolor='dimgray', label='Semiconductor Tracker (SCT) + Transition Radiation Tracker (TRT)')
        ]
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0.05, 0.95))

        plt.tight_layout()
       # plt.show() saad: plt.show doesnt work in terminal in a py file
       plt.savefig("detector.png")

# Utilization
if __name__ == "__main__":
    visualizer = InnerDetectorVisualizer()
    visualizer.plot()

