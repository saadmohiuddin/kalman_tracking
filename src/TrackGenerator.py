"""
This file produces the 'MC truth' track of a charged particle.
Other group members use it to get the 4 intersection layers
for Kalman Filtering.
"""

import numpy as np


def perigee_representation(x, y, z, px, py, pz, q, B=2):
    """
    Takes in cartesian coordinates and gives Perigee coordinates.
    These coordinates is used inside ATLAS
    """

    import numpy as np

    momentum_trans = np.hypot(px, py)
    momentum_mag = np.linalg.norm(x=[px, py, pz])

    theta_polar = np.arccos(pz / momentum_mag)
    azimuth_angle = np.arctan2(py, px)
    signed_curvature = q / momentum_trans
    radius_curvature = momentum_trans / (q * B)

    x_c = x + py / (q * B)
    y_c = y - px / (q * B)

    D = np.hypot(x_c, y_c)
    d0 = D - (np.absolute(radius_curvature))
    z0 = z - (pz / momentum_trans) * d0

    return d0, z0, azimuth_angle, theta_polar, signed_curvature


class TrackGenerator:
    def __init__(
        self,
        x: float,
        y: float,
        z: float,
        px: float,
        py: float,
        pz: float,
        charge: int,
    ):
        self.magfield = 2  # Tesla - ATLAS inner magnet
        self.x = x
        self.y = y
        self.z = z
        self.px = px
        self.py = py
        self.pz = pz
        self.charge = charge
        self.pixel_layers = (33.25, 50.5, 88.5, 122.5)  # millimeter

        self._pT = np.hypot(self.px, self.py)
        self._momentum_mag = np.linalg.norm([self.px, self.py, self.pz])
        self._radius = self._pT / (abs(self.charge) * self.magfield)
        self._x_center = self.x + self.py / (self.charge * self.magfield)
        self._y_center = self.y - self.px / (self.charge * self.magfield)
        self._angular_frequency = (self.charge * self.magfield) / self._pT

    def evolve_track(self, time_steps=100):
        """
        Evolve the particle's track in the magnetic field.

        Parameters
        ----------
            time_steps: int
            Number of time steps to simulate

        Returns
        ----------
            tuple
            Arrays of (x, y, z, phi, time) coordinates along the track
        """
        if time_steps <= 0:
            raise ValueError("time_steps must be positive")

        phi_0 = np.arctan2(self.y - self._y_center, self.x - self._x_center)

        times = np.arange(0, time_steps, 1)
        phi_values = phi_0 + self._angular_frequency * times

        x_values = self._x_center + self._radius * np.cos(phi_values)
        y_values = self._y_center + self._radius * np.sin(phi_values)
        z_values = self.z + (self.pz / self._pT) * (
            self._radius * self._angular_frequency * times
        )

        return x_values, y_values, z_values, phi_values, times

    def find_layer_intersections(self, verbose=False, max_time_steps=500):
        """
        Find intersections of the track with detector layers.

        Parameters
        ----------
            verbose: bool
            Print debugging information
        max_time_steps: int
            Maximum number of time steps to use for track evolution

        Returns
        ----------
            hits: dict
            Dictionary mapping layer radius to (x, y, z) intersection coordinates
        """
        # Start with default time steps, but increase if needed
        time_steps = 100
        hits = {}

        while time_steps <= max_time_steps:
            x_values, y_values, z_values, _, time_values = self.evolve_track(
                time_steps=time_steps
            )
            r_track = np.sqrt(x_values**2 + y_values**2)

            hits = {}
            missing_layers = []

            for layer_radius in self.pixel_layers:
                crossing_idx = self._find_crossing_index(r_track, layer_radius)

                if crossing_idx is None:
                    missing_layers.append(layer_radius)
                    if verbose:
                        print(
                            f"No crossing found for layer {layer_radius} mm with {time_steps} time steps"
                        )
                    continue

                if verbose:
                    print(
                        f"Found crossing at index {crossing_idx} for layer {layer_radius}"
                    )

                hit_coords = self._interpolate_intersection(
                    x_values,
                    y_values,
                    z_values,
                    r_track,
                    crossing_idx,
                    layer_radius,
                )
                hits[layer_radius] = hit_coords

            # If we found all layers or reached max time steps, break
            if len(missing_layers) == 0 or time_steps >= max_time_steps:
                break

            # Increase time steps and try again
            time_steps = min(time_steps * 2, max_time_steps)
            if verbose:
                print(
                    f"Increasing time steps to {time_steps} to find missing layers: {missing_layers}"
                )

        return hits

    def _find_crossing_index(self, r_track, layer_radius):
        """
        Find the index of time(or other X,Y,Z,Phi coordinate) where the track crosses a given layer radius
        Parameters
        ----------
            r_track: np.ndarray
                Radius of track
            layer_radius: float
                Radius of the layer(single number)

        Returns
        ----------
            crossing index: int | None
            Index of time/X/Y/Z coordinate array where particle hits pixel layers

        """
        delta_r = r_track - layer_radius

        for i in range(len(delta_r) - 1):
            if delta_r[i] * delta_r[i + 1] < 0:
                return i
        return None

    def _interpolate_intersection(self, x_vals, y_vals, z_vals, r_vals, idx, target_r):
        """
        Interpolate to find exact intersection coordinates

        Parameters
        ----------
            x_vals: np.ndarray
                    X coordinates of track
            y_vals: np.ndarray
                    Y coordinates of track
            z_vals: np.ndarray
                    Z coordinates of track
            r_vals: np.ndarray
                    Radius of track. TODO - Can be calculated in this function in future improvements.
            idx: int
                Position where the track intersects layer
            target_r:

        Returns
        ----------
            intersection coordinates: tuple[np.ndarray]
                                      tuple of X,Y,Z coordinates of intersection

        """
        r1, r2 = r_vals[idx], r_vals[idx + 1]

        if abs(r1 - r2) > 1e-10:
            alpha = (target_r - r1) / (r2 - r1)
            x_hit = x_vals[idx] + alpha * (x_vals[idx + 1] - x_vals[idx])
            y_hit = y_vals[idx] + alpha * (y_vals[idx + 1] - y_vals[idx])
            z_hit = z_vals[idx] + alpha * (z_vals[idx + 1] - z_vals[idx])
            return (x_hit, y_hit, z_hit)
        else:
            # Use the closer point[using radius]
            if abs(r1 - target_r) <= abs(r2 - target_r):
                return (x_vals[idx], y_vals[idx], z_vals[idx])
            else:
                return (x_vals[idx + 1], y_vals[idx + 1], z_vals[idx + 1])

    def check_layer_reachability(self):
        """
        Check if the track can geometrically reach all detector layers.

        Returns
        ----------
            dict: Information about which layers the track can reach
        """
        center_distance = np.hypot(self._x_center, self._y_center)
        min_radius = abs(center_distance - self._radius)
        max_radius = center_distance + self._radius

        reachable_layers = []
        for layer in self.pixel_layers:
            if min_radius <= layer <= max_radius:
                reachable_layers.append(layer)

        return {
            "reachable_layers": reachable_layers,
            "all_reachable": len(reachable_layers) == len(self.pixel_layers),
            "track_range": (min_radius, max_radius),
            "center": (self._x_center, self._y_center),
            "radius": self._radius,
        }

    def get_perigee_parameters(self):
        """
        Get the track parameters in perigee representation.

        Returns
        ----------
            tuple: (d0, z0, phi, theta, q/pT) perigee parameters
        Reference
        ----------
            https://atlassoftwaredocs.web.cern.ch/internal-links/tracking-tutorial/idoverview/
            https://indico.cern.ch/event/1271601/contributions/5349697/attachments/2625445/4540216/ActsDevMeeting.pdf

        """
        return perigee_representation(
            self.x,
            self.y,
            self.z,
            self.px,
            self.py,
            self.pz,
            self.charge,
            self.magfield,
        )

    def update_parameters(self, **kwargs):
        """
        Update track parameters and recalculate cached values.

        Parameters
        ----------
            **kwargs
            Any combination of x, y, z, px, py, pz, charge, magfield
        """
        # Update provided parameters
        for param, value in kwargs.items():
            if hasattr(self, param):
                setattr(self, param, value)
            else:
                raise ValueError(f"Invalid parameter: {param}")

        # Recalculate cached values
        if any(k in kwargs for k in ["px", "py", "charge", "magfield"]):
            self._pT = np.hypot(self.px, self.py)
            if self._pT == 0:
                raise ValueError("Transverse momentum cannot be zero")
            self._momentum_mag = np.linalg.norm([self.px, self.py, self.pz])
            self._radius = self._pT / (abs(self.charge) * self.magfield)
            self._angular_frequency = (self.charge * self.magfield) / self._pT

        if any(k in kwargs for k in ["x", "y", "px", "py", "charge", "magfield"]):
            self._x_center = self.x + self.py / (self.charge * self.magfield)
            self._y_center = self.y - self.px / (self.charge * self.magfield)

    def generate_track_with_all_hits(self, max_attempts=100, verbose=False):
        """
        Generate a track that hits all detector layers, adjusting parameters if needed.

        Parameters
        ----------
            max_attempts: int
                        Maximum number of adjustment attempts
            verbose: bool
                    Print debugging information

        Returns
        ----------
            dict | None:
            Dictionary of layer intersections, or None if the track does not intersect
        """
        reachability = self.check_layer_reachability()

        if reachability["all_reachable"]:
            return self.find_layer_intersections(verbose=verbose)

        if verbose:
            print("Track doesn't reach all layers, attempting adjustments...")

        original_params = {
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "px": self.px,
            "py": self.py,
            "pz": self.pz,
        }

        min_radius_needed = max(self.pixel_layers)
        min_pt_needed = min_radius_needed * abs(self.charge) * self.magfield

        for attempt in range(max_attempts):
            try:
                self._adjust_track_parameters(attempt, min_pt_needed, verbose)

                new_reachability = self.check_layer_reachability()
                if new_reachability["all_reachable"]:
                    if verbose:
                        print(f"Success on attempt {attempt + 1}")
                    return self.find_layer_intersections(verbose=verbose)

            except Exception as e:
                if verbose:
                    print(f"Attempt {attempt + 1} failed: {e}")
                continue

        self.update_parameters(**original_params)

        if verbose:
            print(f"Failed to find valid track after {max_attempts} attempts")

        return None

    def _adjust_track_parameters(self, attempt, min_pt_needed, verbose):
        """Adjust track parameters for a given attempt."""
        aggression_factor = 1.0 + (attempt * 0.1)

        if attempt % 3 == 0:
            # momentum scale
            if self._pT < min_pt_needed:
                scale = (min_pt_needed / self._pT) * (1.1 * aggression_factor)
            else:
                scale = 1.1 * aggression_factor

            self.update_parameters(px=self.px * scale, py=self.py * scale)

        elif attempt % 3 == 1:
            # position shift
            position_factor = 0.1 * aggression_factor
            offset_x = -self._x_center * position_factor
            offset_y = -self._y_center * position_factor

            self.update_parameters(x=self.x + offset_x, y=self.y + offset_y)

        else:
            # position and momentum shift
            if self._pT < min_pt_needed:
                scale = (min_pt_needed / self._pT) * (1.05 * aggression_factor)
            else:
                scale = 1.05 * aggression_factor

            position_factor = 0.05 * aggression_factor
            offset_x = -self._x_center * position_factor
            offset_y = -self._y_center * position_factor

            self.update_parameters(
                px=self.px * scale,
                py=self.py * scale,
                x=self.x + offset_x,
                y=self.y + offset_y,
            )

        if verbose:
            print(
                f"Attempt {attempt + 1}: pT={self._pT:.2f}, center=({self._x_center:.2f}, {self._y_center:.2f})"
            )

    def __repr__(self):
        """
        String representation of the track.
        """
        return (
            f"TrackGenerator(x={self.x:.2f}, y={self.y:.2f}, z={self.z:.2f}, "
            f"px={self.px:.2f}, py={self.py:.2f}, pz={self.pz:.2f}, "
            f"charge={self.charge}, pT={self._pT:.2f})"
        )


if __name__ == "__main__":
    print("This code can only be used when imported by another file")
    print("Import and use it")
    print("Quitting....")
    raise RuntimeError("This module is not intended to be run as a script.")
