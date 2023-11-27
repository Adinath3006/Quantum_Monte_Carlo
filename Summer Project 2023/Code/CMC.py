import numpy as np

def initialize_spins(N):
    return np.random.choice([-1, 1], size=N)

def calculate_energy(spin_path, J):
    N = spin_path.shape[0]
    energy = 0
    for i in range(N):
        energy += -J * spin_path[i] * spin_path[(i + 1) % N]
    return energy

def metropolis_hastings_update(spin_path, J, beta):
    N = spin_path.shape[0]
    for i in range(N):
        flipped_spin_path = spin_path.copy()
        flipped_spin_path[i] *= -1

        delta_E = 2 * J * spin_path[i] * (spin_path[(i - 1) % N] + spin_path[(i + 1) % N])
        transition_prob = np.exp(-beta * delta_E)
        
        if np.random.rand() < transition_prob:
            spin_path = flipped_spin_path

    return spin_path

def simulate_ising_model(N, J, beta, num_steps):
    spin_path = initialize_spins(N)
    energies = []
    for step in range(num_steps):
        spin_path = metropolis_hastings_update(spin_path, J, beta)
        energy = calculate_energy(spin_path, J)
        energies.append(energy)
    return spin_path, energies

# Example usage
N = 100  # Number of lattice sites
J = 1    # Coupling constant
beta = 1 # Inverse temperature
num_steps = 10000

spin_path, energies = simulate_ising_model(N, J, beta, num_steps)

print("Final spin path:", spin_path)
print("Energies:", energies)
