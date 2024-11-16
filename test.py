import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Filtrar la parte superior del perfil (mitad inicial de perfil_x y perfil_y)
mitad = len(perfil_x) // 2
perfil_x_superior = perfil_x[:mitad]
perfil_y_superior = perfil_y[:mitad]

# Parámetros de la simulación
U = 1.0  # velocidad de corriente libre
c = 1.0  # longitud de la cuerda
epsilon = 1e-5  # valor pequeño para evitar divisiones por cero
fuente_x, fuente_y = 0.08, 0.0

# Función para calcular la velocidad en los puntos de estancamiento
def velocidad_punto_estancamiento(params):
    Lambda, n = params
    n = int(n)  # n debe ser un entero
    sumideros_x = np.linspace(0.1, c, n)
    sumideros_y = np.zeros(n)
    intensidad_sumidero = -Lambda / n

    # Velocidad en el punto de estancamiento delantero (0, 0)
    velocidad_1 = U + (Lambda / (2 * np.pi)) * (
        (0 - fuente_x) / ((0 - fuente_x) ** 2 + (0 - fuente_y) ** 2 + epsilon)
    ) + np.sum(
        (intensidad_sumidero / (2 * np.pi)) * (
            (0 - sumideros_x) / ((0 - sumideros_x) ** 2 + (0 - sumideros_y) ** 2 + epsilon)
        )
    )

    # Velocidad en el punto de estancamiento trasero (c, 0)
    velocidad_2 = U + (Lambda / (2 * np.pi)) * (
        (c - fuente_x) / ((c - fuente_x) ** 2 + (0 - fuente_y) ** 2 + epsilon)
    ) + np.sum(
        (intensidad_sumidero / (2 * np.pi)) * (
            (c - sumideros_x) / ((c - sumideros_x) ** 2 + (0 - sumideros_y) ** 2 + epsilon)
        )
    )

    # Queremos minimizar la suma de los cuadrados de las velocidades en los puntos de estancamiento
    return velocidad_1**2 + velocidad_2**2

# Valores iniciales para Lambda y n
valores_iniciales = [0.4, 10]

# Realizar la optimización
resultado = minimize(velocidad_punto_estancamiento, valores_iniciales, bounds=[(0.1, 2.0), (1, 20)])

# Obtener los valores óptimos de Lambda y n
Lambda_opt, n_opt = resultado.x
n_opt = int(n_opt)
sumideros_x = np.linspace(0.1, c, n_opt)
sumideros_y = np.zeros(n_opt)
intensidad_sumidero = -Lambda_opt / n_opt

# Cálculo de la distribución de velocidad sobre la parte superior del perfil
velocidad_x_superior = U + (Lambda_opt / (2 * np.pi)) * (
    (perfil_x_superior - fuente_x) / ((perfil_x_superior - fuente_x) ** 2 + (perfil_y_superior - fuente_y) ** 2 + epsilon)
) + np.sum(
    (intensidad_sumidero / (2 * np.pi)) * (
        (perfil_x_superior[:, np.newaxis] - sumideros_x) / ((perfil_x_superior[:, np.newaxis] - sumideros_x) ** 2 + (perfil_y_superior[:, np.newaxis] - sumideros_y) ** 2 + epsilon)
    ),
    axis=1
)

velocidad_y_superior = (Lambda_opt / (2 * np.pi)) * (
    (perfil_y_superior - fuente_y) / ((perfil_x_superior - fuente_x) ** 2 + (perfil_y_superior - fuente_y) ** 2 + epsilon)
) + np.sum(
    (intensidad_sumidero / (2 * np.pi)) * (
        (perfil_y_superior[:, np.newaxis] - sumideros_y) / ((perfil_x_superior[:, np.newaxis] - sumideros_x) ** 2 + (perfil_y_superior[:, np.newaxis] - sumideros_y) ** 2 + epsilon)
    ),
    axis=1
)

velocidad_total_superior = np.sqrt(velocidad_x_superior**2 + velocidad_y_superior**2)

# Cálculo del coeficiente de presión Cp sobre la parte superior del perfil
Cp_superior = 1 - (velocidad_total_superior / U) ** 2

# Cálculo de la velocidad en los puntos de estancamiento con los valores óptimos
velocidad_punto_estancamiento_1 = U + (Lambda_opt / (2 * np.pi)) * (
    (0 - fuente_x) / ((0 - fuente_x) ** 2 + (0 - fuente_y) ** 2 + epsilon)
) + np.sum(
    (intensidad_sumidero / (2 * np.pi)) * (
        (0 - sumideros_x) / ((0 - sumideros_x) ** 2 + (0 - sumideros_y) ** 2 + epsilon)
    )
)
velocidad_punto_estancamiento_2 = U + (Lambda_opt / (2 * np.pi)) * (
    (c - fuente_x) / ((c - fuente_x) ** 2 + (0 - fuente_y) ** 2 + epsilon)
) + np.sum(
    (intensidad_sumidero / (2 * np.pi)) * (
        (c - sumideros_x) / ((c - sumideros_x) ** 2 + (0 - sumideros_y) ** 2 + epsilon)
    )
)

print(f"Valores óptimos: Lambda = {Lambda_opt:.4f}, n = {n_opt}")
print(f"Velocidad en el punto de estancamiento delantero (0, 0): {velocidad_punto_estancamiento_1:.4f}")
print(f"Velocidad en el punto de estancamiento trasero (c, 0): {velocidad_punto_estancamiento_2:.4f}")

# Graficar distribución de velocidad sobre la parte superior del perfil
plt.figure()
plt.plot(perfil_x_superior, velocidad_total_superior, "b-", linewidth=2)
plt.xlabel("Posición x sobre el perfil")
plt.ylabel("Velocidad sobre la parte superior del perfil")
plt.title("Distribución de Velocidad sobre la parte superior del perfil NACA 0024")

# Graficar Cp sobre la parte superior del perfil
plt.figure()
plt.plot(perfil_x_superior, Cp_superior, "r-", linewidth=2)
plt.xlabel("Posición x sobre el perfil")
plt.ylabel("Coeficiente de presión Cp")
plt.title("Distribución del Coeficiente de Presión Cp sobre la parte superior del perfil NACA 0024")
plt.gca().invert_yaxis()  # Invertir el eje y para visualizar Cp correctamente

plt.show()
