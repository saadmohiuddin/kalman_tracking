import numpy as np
import matplotlib.pyplot as plt
import TrackGenerator


# Constants



#track = TrackGenerator.TrackGenerator(hit_coordinates, B, q=q, h=h, R_t=R_t, d0=d0, phi0=phi0, z0=z0, theta=theta)

def fit_circle_least_squares(points):
    """
    Fits a circle to 4 (or more) 2D points using least squares.
    Returns the radius.
    """
    x = np.array([p[0] for p in points])
    y = np.array([p[1] for p in points])
    
    # Formulating matrix for least squares
    A = np.c_[2*x, 2*y, np.ones(len(points))]
    b = x**2 + y**2
    
    # Solving for circle parameters: Ax + By + C = r^2
    sol = np.linalg.lstsq(A, b, rcond=None)[0]
    center_x, center_y = sol[0], sol[1]
    radius = np.sqrt(sol[2] + center_x**2 + center_y**2)
    
    return radius

def estimate_initial_momentum(hits, B=2.0, charge=1):
    r1, r2, r3, r4 = hits[:4]
    direction = r2 - r1
    direction /= np.linalg.norm(direction)

    # Crude estimate of radius (use proper circle fit for accuracy)
    R = fit_circle_least_squares([r1[:2], r2[:2], r3[:2],r4[:2]])

    pt =  charge * B * R  # GeV
    momentum = pt * direction
    return momentum

#estimate_initial_momentum(hit_coordinates, B=2.0, charge=-1)

def lorentz_force(q, p, B):
    # Velocity from momentum (non-relativistic approx)
    v = p #/ np.linalg.norm(p)  # Normalize momentum to get direction
    return q * np.array([v[1], -v[0], 0]) * B  # Cross product in 3D

def rk4_step(x, p, h=0.01, B=2.0,q=1):
    # RK4 integration for position and momentum
    def dp_dt(p): return lorentz_force(q,p, B)
    def dx_dt(p): return p   # direction of motion

    k1_x = h * dx_dt(p)
    k1_p = h * dp_dt(p)

    k2_x = h * dx_dt(p + 0.5 * k1_p)
    k2_p = h * dp_dt(p + 0.5 * k1_p)

    k3_x = h * dx_dt(p + 0.5 * k2_p)
    k3_p = h * dp_dt(p + 0.5 * k2_p)

    k4_x = h * dx_dt(p + k3_p)
    k4_p = h * dp_dt(p + k3_p)

    x_new = x + (1/6) * (k1_x + 2*k2_x + 2*k3_x + k4_x)
    p_new = p + (1/6) * (k1_p + 2*k2_p + 2*k3_p + k4_p)

    return x_new, p_new

def recreate_track(hits, B=2.0, q=1,h=-0.001,time_steps=5000):
    # Estimate initial momentum

    p = estimate_initial_momentum(hits, B, q)
    x = np.array(hits[0][:3])  # Start from the first hit position
    trajectory = [x]
    
    for i in range(time_steps):
        x, p = rk4_step(x, p,h=h,q=q)
        trajectory.append(x)
    return trajectory

# Example usage
def get_impact_parameter(hits):
    """
    Calculate the impact parameter from the first two hits.
    """
    r1, r2 = hits[0][:2], hits[1][:2]
    return np.linalg.norm(r1 - r2) / 2  # Simple average distance

if __name__ == "__main__":
    # Initial conditions
    initial_guess = {"x": 0.0, "y": 0.0, "z": 0.0, "px": 0.2, "py": 0.1, "pz": 0.3, "charge": q}
    track = TrackGenerator.TrackGenerator(**initial_guess)

    print("-" * 25, "Original Track", "-" * 25)
    print(track)
    print("-" * 50)

    # Generate track with RK4
    hit_coordinates = track.generate_track_with_all_hits(max_attempts=100)
    print("Hit Coordinates:", hit_coordinates)

    # Simulate trajectory
    x = np.array([0.0, 0.0, 0.0])  # initial position
    p = np.array([0.2, 0.1, 0.3])  # initial momentum (GeV/c)

    trajectory = recreate_track(hit_coordinates, B=2.0, charge=q)

    

    print(np.array(trajectory)[:, 0], np.array(trajectory)[:, 1],np.array(trajectory)[:, 2])
    x = np.array(trajectory)[:, 0]
    y = np.array(trajectory)[:, 1]
    z = np.array(trajectory)[:, 2]
    print(x)
    print(y)
    print(z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z)
    ax.set_xlabel('X Position (m)')
    ax.set_ylabel('Y Position (m)')
    ax.set_zlabel('Z Position (m)')
    ax.set_title('Particle Trajectory in Magnetic Field')

    fig.savefig('particle_trajectory.png')
        