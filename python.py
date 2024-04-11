import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters
L = 1.0  # Length of the casting (m)
Nx = 100  # Number of grid points in the casting
dx = L / Nx  # Grid spacing
Nt = 100  # Number of time steps
dt = 0.1  # Time step (s)

# Material properties
density = 7000  # Density of the material (kg/m^3)
specific_heat = 500  # Specific heat capacity of the material (J/kg*K)
thermal_conductivity = 50  # Thermal conductivity of the material (W/m*K)

# Initial temperature distribution (initial condition)
T_initial = 300  # Initial temperature (K)
T = np.ones((Nt+1, Nx + 1)) * T_initial  # Initial temperature array

# Boundary conditions
T_left = 1000  # Temperature at the left boundary (K)
T_right = 300  # Temperature at the right boundary (K)

# Function to update temperature using finite difference method
def update_temperature(T):
    for t in range(1, Nt+1):
        for i in range(1, Nx):
            T[t, i] = T[t-1, i] + dt * thermal_conductivity / (density * specific_heat * dx**2) * (T[t-1, i-1] - 2*T[t-1, i] + T[t-1, i+1])

        # Apply boundary conditions
        T[t, 0] = T_left
        T[t, -1] = T_right

    return T

# Perform FEA simulation
T = update_temperature(T)

# Create meshgrid for 3D plotting
X, Y = np.meshgrid(np.linspace(0, L, Nx + 1), np.arange(Nt + 1) * dt)
Z = T

# Plot 3D surface
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap='coolwarm', edgecolor='none')
ax.set_title('Temperature Distribution in Casting Over Time')
ax.set_xlabel('Position (m)')
ax.set_ylabel('Time (s)')
ax.set_zlabel('Temperature (K)')
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
plt.show()
