import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

class DynamicForcingField:
    def __init__(self, L, Li, N, M, noise_std, total_steps):
        """
        Initialize the forcing field parameters
        L: domain size
        Li: forcing scale
        N: number of grid points in each direction
        M: number of forcing components
        noise_std: standard deviation for parameter evolution
        total_steps: number of time steps to simulate
        """
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
    # Set parameters
    L = 1.0E8          # Domain size
    Li = 4.0E6         # Forcing scale
    N = 1024           # Grid resolution
    M = 100             # Number of forcing components
    noise_std = 0.001  # Standard deviation for parameter evolution
    total_steps = 80  # Total number of time steps

    # Create and run simulation
    simulator = DynamicForcingField(L, Li, N, M, noise_std, total_steps)
    simulator.simulate()
    simulator.create_animation()