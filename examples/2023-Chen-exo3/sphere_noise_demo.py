import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

def generate_small_scale_noise(n_lat = 400, n_lon = 800, n_waves = 15, scale=100.0):
    """Generate small-scale noise pattern on a sphere.
    
    Args:
        n_waves: Number of waves
        scale: Number of waves per unit sphere
        seed: Random seed for reproducibility
    
    Returns:
        theta: Latitude angles in radians
        phi: Longitude angles in radians
        pattern: 2D array of noise values
    """
    # Create grid of points
    theta = np.linspace(0, np.pi, n_lat)
    phi = np.linspace(0, 2*np.pi, n_lon)
    theta_grid, phi_grid = np.meshgrid(theta, phi)
    
    # Initialize pattern
    pattern = np.zeros((n_lon, n_lat))
    
    # Generate random positions on sphere
    # Using fibonacci sphere method for even distribution
    golden_ratio = (1 + np.sqrt(5)) / 2
    i = np.arange(n_waves)
    
    # Add random perturbations to positions
    theta_centers = np.arccos(1 - 2*(i+0.5)/n_waves) + np.random.normal(0, 0.1, n_waves)
    phi_centers = 2*np.pi * i/golden_ratio + np.random.normal(0, 0.1, n_waves)
    
    # Random amplitudes and phase shifts for each wave
    amplitudes = np.random.uniform(0.5, 1.5, n_waves)
    phase_shifts = np.random.uniform(0, 2*np.pi, n_waves)
    
    for i in range(n_waves):
        angular_dist = np.arccos(
            np.sin(theta_centers[i])*np.cos(phi_centers[i])*np.sin(theta_grid)*np.cos(phi_grid) +
            np.sin(theta_centers[i])*np.sin(phi_centers[i])*np.sin(theta_grid)*np.sin(phi_grid) +
            np.cos(theta_centers[i])*np.cos(theta_grid)
        )
        pattern += amplitudes[i] * np.sin(scale * angular_dist + phase_shifts[i])
    
    pattern = pattern / n_waves
    
    return theta, phi, pattern

def init_centers(n_waves):    
    # Generate random positions on sphere
    # Using fibonacci sphere method for even distribution
    golden_ratio = (1 + np.sqrt(5)) / 2
    i = np.arange(n_waves)
    
    # Add random perturbations to positions
    theta_centers = np.arccos(1 - 2*(i+0.5)/n_waves) + np.random.normal(0, 0.1, n_waves)
    phi_centers = 2*np.pi * i/golden_ratio + np.random.normal(0, 0.1, n_waves)

    return theta_centers, phi_centers

def generate__noise_given_centers(centers_lat, centers_lon, scale, n_lat = 400, n_lon = 800): 
    theta = np.linspace(0, np.pi, n_lat)
    phi = np.linspace(0, 2*np.pi, n_lon)
    theta_grid, phi_grid = np.meshgrid(theta, phi)
    
    # Initialize pattern
    pattern = np.zeros((n_lon, n_lat))
    
    for i in range(n_waves):
        angular_dist = np.arccos(
            np.sin(centers_lat[i])*np.cos(centers_lon[i])*np.sin(theta_grid)*np.cos(phi_grid) +
            np.sin(centers_lat[i])*np.sin(centers_lon[i])*np.sin(theta_grid)*np.sin(phi_grid) +
            np.cos(centers_lat[i])*np.cos(theta_grid)
        )
        pattern += np.sin(scale * angular_dist)
    
    n_waves = len(centers_lat)
    pattern = pattern / n_waves
    
    return pattern

class SmallScaleNoiseField:
    def __init__(self, n_lat=400, n_lon=800, n_waves=15, scale=100.0, alpha=0.001):
        """
        Initialize the small scale noise field with EMA
        
        Args:
            n_lat: Number of latitude points
            n_lon: Number of longitude points 
            n_waves: Number of waves for noise generation
            scale: Scale of the noise pattern
            alpha: EMA decay factor (0-1), smaller means slower decay
        """
        self.n_lat = n_lat
        self.n_lon = n_lon
        self.n_waves = n_waves
        self.scale = scale
        self.alpha = alpha
        
        # Initialize first field
        self.current_field = np.zeros((self.n_lon, self.n_lat))
        
        # Store field history for animation
        self.field_history = [self.current_field.copy()]
        
    def step(self):
        # Generate new noise field
        _, _, new_field = generate_small_scale_noise(
            n_lat=self.n_lat,
            n_lon=self.n_lon,
            n_waves=self.n_waves,
            scale=self.scale
        )

        # Check if new field is the same as the previous field
        if np.allclose(new_field, self.current_field):
            print("New field is the same as the previous field")
            return
        
        # Update EMA
        self.current_field = (1 - self.alpha) * self.current_field + self.alpha * new_field
        
        # Store for animation
        self.field_history.append(self.current_field.copy())
        
    def animate(self, filename='ema_noise_evolution.gif', steps=100, fps=10):
        """
        Create animation from field history
        
        Args:
            filename: Output gif filename
            fps: Frames per second
        """
        fig, ax = plt.subplots(figsize=(10, 5))

        for i in tqdm(range(steps)):
            self.step()
        
        def update(frame):
            ax.clear()
            im = ax.imshow(self.field_history[frame].T, 
                         cmap='RdBu',
                         aspect='equal',
                         extent=[0, 2*np.pi, 0, np.pi])
            ax.set_title(f'Frame {frame}')
            return [im]
            
        anim = FuncAnimation(fig, update, 
                           frames=len(self.field_history),
                           interval=1000/fps,
                           blit=True)
        
        anim.save(filename, writer='pillow', fps=fps)
        plt.close()


class DynamicForcingField:
    def __init__(self, L, Li, N, M, noise_std, total_steps):
        self.L = L
        self.Li = Li
        self.N = N
        self.M = M
        self.noise_std = noise_std
        self.total_steps = total_steps
        
        # Create spatial grid
        self.dx = L / N
        x = np.linspace(0, L, N)
        y = np.linspace(0, L, N)
        self.X, self.Y = np.meshgrid(x, y)
        
        # Initialize forcing parameters
        self.A = np.random.uniform(-0.5, 0.5, M)  # Strength
        self.B = np.random.uniform(0, 2*np.pi, M)  # Phase
        self.C = np.random.uniform(0, 2*np.pi, M)  # Direction
        
        # Store field history for animation
        self.field_history = []

    def compute_field(self, time_step):
        """Compute the forcing field at current time step"""
        field = np.zeros_like(self.X)
        for i in range(self.M):
            field += self.A[i] * np.cos(
                (2*np.pi/self.Li) * 
                (self.X*np.cos(self.C[i]) + self.Y*np.sin(self.C[i])) + 
                self.B[i] * time_step * 0.001
            )
        return field

    def evolve_parameters(self):
        """Evolve the forcing parameters using normal noise"""
        self.A += np.random.normal(0, self.noise_std * 3, self.M)
        self.B += np.random.normal(0, self.noise_std * 2 * np.pi, self.M)
        self.C += np.random.normal(0, self.noise_std * 2 * np.pi, self.M)

    def simulate(self):
        """Run the simulation and store field history"""
        for step in range(self.total_steps):
            field = self.compute_field(step)
            self.field_history.append(field)
            self.evolve_parameters()

    def create_animation(self, filename='forcing_field.gif'):
        """Create and save animation"""
        fig, ax = plt.subplots(figsize=(8, 8))
        
        def update(frame):
            ax.clear()
            im = ax.imshow(self.field_history[frame], 
                         extent=[0, self.L, 0, self.L],
                         cmap='RdBu_r',
                         aspect='equal')
            ax.set_title(f'Time step: {frame}')
            
        ani = FuncAnimation(fig, update, frames=self.total_steps, 
                          interval=100, blit=False)
        
        ani.save(filename, writer='pillow')
        plt.close()

# Example usage
if __name__ == "__main__":

    # EMA noise field
    n_lat = 400
    n_lon = 800
    n_waves = 15
    scale = 100.0
    alpha = 0.1
    noise_field = SmallScaleNoiseField(n_lat, n_lon, n_waves, scale, alpha)
    noise_field.animate()
