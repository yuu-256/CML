import numpy as np
import matplotlib.pyplot as plt

# Load output file
fnm = 'output/simnx20ny10nsteps5000.txt'
fnm = 'output/simnx80ny40nsteps5000.txt'
# fnm = 'output/simnx80ny40nsteps1000.txt'
# fnm = 'output/simnx160ny80nsteps5000.txt'
# fnm = 'output/simnx40ny20nsteps5000.txt'
# fnm = 'output/simnx80ny40nsteps10000.txt'
# fnm = 'output/simnx40ny20nsteps5000.txt'
# fnm = 'output/simnx80ny40nsteps200.txt'
# fnm = 'output/simnx80ny40nsteps500.txt'
with open(fnm, 'r') as f:
    lines = f.readlines()

# Get nx and ny
for line in lines:
    if 'nx' in line and 'ny' in line:
        parts = line.split()
        nx = int(parts[2])
        ny = int(parts[5])
        break

# Horizontal velocity
start_index = None
for i, line in enumerate(lines):
    if 'Horizontal velocity(u)' in line:
        start_index = i + 1
        hv = np.array([list(map(float, line.split())) for line in lines[start_index:start_index+ny]])
    if 'Vertical velocity(v)' in line:
        start_index = i + 1
        vv = np.array([list(map(float, line.split())) for line in lines[start_index:start_index+ny]])
    if 'Internal energy' in line:
        start_index = i + 1
        ie = np.array([list(map(float, line.split())) for line in lines[start_index:start_index+ny]])
    if 'Liquid water' in line:
        start_index = i + 1
        lw = np.array([list(map(float, line.split())) for line in lines[start_index:start_index+ny]])
    if 'Water vapor' in line:
        start_index = i + 1
        wv = np.array([list(map(float, line.split())) for line in lines[start_index:start_index+ny]])

# Get data and transform to 2D array
plt.imshow(lw, origin='lower', extent=[0, nx, 0, ny], cmap='Greys_r')
plt.title('Liquid water mass')
plt.colorbar()
plt.clim(0, 0.0001)
plt.savefig('image/ie.png')
plt.show()
