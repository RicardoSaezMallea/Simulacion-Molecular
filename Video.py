import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, PillowWriter
from tqdm import tqdm

# Leer el fichero
filename = 'C:/Users/ricar/Desktop/Pruebas VSC/Prubea/Moleculas_Trayectoria_N64_eta1.0_h0.0100_k1000_b1_it1000.txt'
data = np.loadtxt(filename)

# Número de monómeros
N = 64  # Cambia esto al número real de monómeros

# Reshape de los datos para tener la forma (frames, monomeros, 3)
data = data.reshape(-1, N, 3)

# Crear la figura y el eje 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Ajustar los límites de los ejes basándonos en los datos
x_min, x_max = np.min(data[:, :, 0]), np.max(data[:, :, 0])
y_min, y_max = np.min(data[:, :, 1]), np.max(data[:, :, 1])
z_min, z_max = np.min(data[:, :, 2]), np.max(data[:, :, 2])

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_zlim(z_min, z_max)

# Inicializar la posición de los monómeros
scat = ax.scatter(data[0, :, 0], data[0, :, 1], data[0, :, 2])

# Añadir las líneas de los muelles
lines = [ax.plot([], [], [], 'r-')[0] for _ in range(N-1)]

# Función de actualización para la animación
def update(frame):
    # Actualizar posiciones de los monómeros
    scat._offsets3d = (data[frame, :, 0], data[frame, :, 1], data[frame, :, 2])
    
    # Actualizar los muelles
    for i in range(N-1):
        lines[i].set_data(data[frame, i:i+2, 0], data[frame, i:i+2, 1])
        lines[i].set_3d_properties(data[frame, i:i+2, 2])
    return scat, *lines

# Crear la animación con barra de progreso
frames = tqdm(range(data.shape[0]), desc="Generating animation")
ani = FuncAnimation(fig, update, frames=frames, interval=50, blit=False)

# Guardar la animación
ani.save('polimero_tryout.gif', writer=PillowWriter(fps=20))

# Mostrar la animación
plt.show()
