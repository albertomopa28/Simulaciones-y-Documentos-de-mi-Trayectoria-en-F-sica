import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Nombre del archivo con datos
archivo = "cohete_data.dat"

# Leer datos del archivo
def leer_datos(archivo):
    datos = []
    with open(archivo, "r") as f:
        frame = []
        for linea in f:
            if linea.strip():
                frame.append([float(x) for x in linea.split(",")])
            else:
                datos.append(np.array(frame))
                frame = []
        if frame:
            datos.append(np.array(frame))
    return datos

# Procesar datos
datos = leer_datos(archivo)

# Crear figura y ejes
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_xlim(-200000, 200000)
ax.set_ylim(-200000, 200000)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Movimiento del cohete, Tierra, Marte y Sol")

# Dibujar el Sol (posición fija en el origen)
sol, = ax.plot(0, 0, 'yo', markersize=10, label="Sol")

# Inicializar trayectorias
cohete_line, = ax.plot([], [], 'b-', label="Cohete")
tierra_line, = ax.plot([], [], 'orange', label="Tierra")
marte_line, = ax.plot([], [], 'g-', label="Marte")

# Inicializar puntos
cohete_point, = ax.plot([], [], 'bo')
tierra_point, = ax.plot([], [], 'o', color="orange")
marte_point, = ax.plot([], [], 'go')

# Leyenda
ax.legend()

# Función de inicialización
def init():
    cohete_line.set_data([], [])
    tierra_line.set_data([], [])
    marte_line.set_data([], [])
    cohete_point.set_data([], [])
    tierra_point.set_data([], [])
    marte_point.set_data([], [])
    return cohete_line, tierra_line, marte_line, cohete_point, tierra_point, marte_point

# Función de actualización
def actualizar(frame):
    cohete = datos[frame][:1, :]
    tierra = datos[frame][1:2, :]
    marte = datos[frame][2:3, :]

    cohete_line.set_data(cohete[:, 0], cohete[:, 1])
    tierra_line.set_data(tierra[:, 0], tierra[:, 1])
    marte_line.set_data(marte[:, 0], marte[:, 1])

    cohete_point.set_data(cohete[:, 0], cohete[:, 1])
    tierra_point.set_data(tierra[:, 0], tierra[:, 1])
    marte_point.set_data(marte[:, 0], marte[:, 1])

    return cohete_line, tierra_line, marte_line, cohete_point, tierra_point, marte_point

# Crear animación
anim = FuncAnimation(fig, actualizar, frames=len(datos), init_func=init, blit=True)

# Guardar animación como video
anim.save("planetas_cohete.mp4", writer="ffmpeg", fps=30)

# Mostrar gráfica inicial
plt.show()
