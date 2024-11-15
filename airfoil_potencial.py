import numpy as np
import matplotlib.pyplot as plt

# Datos del perfil NACA 0024
perfil_x = np.array(
    [
        1.0000,
        0.9500,
        0.9000,
        0.8000,
        0.7000,
        0.6000,
        0.5000,
        0.4000,
        0.3000,
        0.2500,
        0.2000,
        0.1500,
        0.1000,
        0.0750,
        0.0500,
        0.0250,
        0.0125,
        0.0000,
        0.0125,
        0.0250,
        0.0500,
        0.0750,
        0.1000,
        0.1500,
        0.2000,
        0.2500,
        0.3000,
        0.4000,
        0.5000,
        0.6000,
        0.7000,
        0.8000,
        0.9000,
        0.9500,
        1.0000,
    ]
)
perfil_y = np.array(
    [
        0.00252,
        0.01613,
        0.02896,
        0.05247,
        0.07328,
        0.09127,
        0.10588,
        0.11607,
        0.12004,
        0.11883,
        0.11475,
        0.10691,
        0.09365,
        0.08400,
        0.07109,
        0.05229,
        0.03788,
        0.00000,
        -0.03788,
        -0.05229,
        -0.07109,
        -0.08400,
        -0.09365,
        -0.10691,
        -0.11475,
        -0.11883,
        -0.12004,
        -0.11607,
        -0.10588,
        -0.09127,
        -0.07328,
        -0.05247,
        -0.02896,
        -0.01613,
        -0.00252,
    ]
)

# Parámetros de la simulación
U = 1.0  # velocidad de corriente libre
Lambda = 0.4  # intensidad de la fuente
n = 10  # número de sumideros
c = 1.0  # longitud de la cuerda

# Posiciones de la fuente y los sumideros
fuente_x, fuente_y = 0.05, 0.0
sumideros_x = np.linspace(0.1, c, n)
sumideros_y = np.zeros(n)
intensidad_sumidero = -Lambda / n

# Generar una malla de puntos para calcular las funciones de corriente
x = np.linspace(-1.0, 2.0, 200)
y = np.linspace(-1.0, 1.0, 200)
X, Y = np.meshgrid(x, y)

# Función de corriente para la corriente uniforme
psi_uniforme = U * Y

# Función de corriente para la fuente
psi_fuente = (Lambda / (2 * np.pi)) * np.arctan2(Y - fuente_y, X - fuente_x)

# Función de corriente para los sumideros
psi_sumideros = np.zeros_like(X)
for i in range(n):
    psi_sumideros += (intensidad_sumidero / (2 * np.pi)) * np.arctan2(
        Y - sumideros_y[i], X - sumideros_x[i]
    )

# Superposición de las funciones de corriente
psi_total = psi_uniforme + psi_fuente + psi_sumideros

# Graficar el perfil NACA, líneas de corriente y puntos de estancamiento
plt.figure(figsize=(10, 6))
plt.contour(
    X,
    Y,
    psi_total,
    levels=np.linspace(-2, 2, 50),
    colors="blue",
    linestyles="solid",
)
plt.plot(
    perfil_x, perfil_y, "k-", linewidth=2, label="Perfil NACA 0024"
)  # Perfil NACA

# Marcar puntos de estancamiento (aproximados)
plt.plot(0, 0, "ro", label="Punto de Estancamiento (delantero)")
plt.plot(c, 0, "ro", label="Punto de Estancamiento (trasero)")

# Identificar zonas de presión estática
plt.fill_between(
    perfil_x,
    perfil_y,
    color="lightgray",
    alpha=0.3,
    label="Zona de Presión Estática",
)

# Título y etiquetas
plt.title(
    "Líneas de corriente alrededor del perfil NACA 0024 con zonas de presión estática"
)
plt.xlabel("x")
plt.ylabel("y")
plt.legend()

# Cálculo de la distribución de velocidad sobre el perfil
velocidad_x = U - (Lambda / (2 * np.pi)) / perfil_x
velocidad_y = np.zeros_like(
    perfil_x
)  # simplificación para flujo potencial sin ángulo
velocidad_total = np.sqrt(velocidad_x * 2 + velocidad_y * 2)

# Cálculo del coeficiente de presión Cp
Cp = 1 - (velocidad_total / U) ** 2

# Graficar distribución de velocidad sobre el perfil
plt.figure()
plt.plot(perfil_x, velocidad_total, "b-", linewidth=2)
plt.xlabel("Posición x sobre el perfil")
plt.ylabel("Velocidad sobre el perfil")
plt.title("Distribución de Velocidad sobre el perfil NACA 0024")

# Graficar Cp sobre el perfil
plt.figure()
plt.plot(perfil_x, Cp, "r-", linewidth=2)
plt.xlabel("Posición x sobre el perfil")
plt.ylabel("Coeficiente de presión Cp")
plt.title(
    "Distribución del Coeficiente de Presión Cp sobre el perfil NACA 0024"
)
plt.gca().invert_yaxis()  # Invertir el eje y para visualizar Cp correctamente

plt.show()
