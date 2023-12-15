#!/usr/bin/env python3

#Pendulo elastico
#version 3

import numpy as np
from runge_kutta import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sys import exit as cerrar

print("El siguiente programa simula un pendulo que cuelga de un resorte dadas unas condiciones iniciales dadas.")
#propiedades del pendulo
largo_natural: float = float(input("Ingrese el largo natural del resorte en metros: "))
masa: float = float(input("Ingrese la masa del pendulo en kilogramos: "))
constante_elastica: float = float(input("Ingrese la constante elastica del resorte: "))
g: float = 9.8

for i in (largo_natural, masa):
    if i <= 0: 
        print(">:(")
        cerrar()

largo_corregido: float = largo_natural + masa * g / constante_elastica
print(f"El largo del pendulo en la posición de equilibrio es {largo_corregido}")


def pendulo_elastico(t, mov):

    """Este es el sistema de ecuaciones diferenciales asociado a las ecuaciones
    de movimiento del pendulo elastico
    Aquí que t es el tiempo y mov es un arreglo que contiene el angulo, velocidad angular,
    velocidad radial del pendulo y la deformacion del resorte
    En general se escribira mov = [theta, omega, deformacion, Vradial] """

    theta = mov[0]
    omega = mov[1]
    deformacion = mov[2]
    Vradial = mov[3]

    var_1 = omega
    var_2 = -2 * Vradial * omega / (largo_corregido + deformacion) - g * np.sin(theta) / (largo_corregido + deformacion)
    var_3 = Vradial
    var_4 = (largo_corregido + deformacion) * np.power(omega, 2) - constante_elastica * deformacion / masa + g *np.cos(theta)

    return np.array([var_1, var_2, var_3, var_4])

#condiciones iniciales
print("Ingrese las condiciones iniciales para el pendulo elastico")
theta_0: float = float(input("Ingrese el angulo inicial en radianes: "))
omega_0: float = float(input("Ingrese la velocidad angular inicial en radianes: "))
deformacion_inicial: float = float(input("Ingrese la deformación inicial del resorte en metros: "))
Vradial_0: float = float(input("Ingrese la velocidad radial inicial en metros por segundo: "))

#Aquí se obtienen los valores para simular el pendulo con el metodo de Runge Kutta
puntos_t, puntos_y = RK6(pendulo_elastico, 0, [theta_0, omega_0, deformacion_inicial, Vradial_0], 200)

#Configurar la animación
fig, ax = plt.subplots()
ax.set_xlim(-3 * largo_corregido, 3 * largo_corregido)
ax.set_ylim(-3 * largo_corregido, 3 * largo_corregido)
linea, = ax.plot([], [], lw = 2) 
bola = ax.scatter([], [])

# Listas para almacenar las posiciones anteriores del pendulo
estela_x = []
estela_y = []

#Funcion que actualiza los frames
def update(frame):
    x = [0, (largo_corregido + puntos_y[frame, 2] ) * np.sin(puntos_y[frame, 0])]
    y = [0, - (largo_corregido + puntos_y[frame, 2]) * np.cos(puntos_y[frame, 0])]
    linea.set_data(x, y)

    estela_x.append(x[1]) #Esto acumula las posiciones por las que pasa el pendulo
    estela_y.append(y[1])

    bola_x = (largo_corregido + puntos_y[frame, 2] ) * np.sin(puntos_y[frame, 0])
    bola_y = - (largo_corregido + puntos_y[frame, 2]) * np.cos(puntos_y[frame, 0])
    bola.set_offsets([[bola_x, bola_y]])

    estela, = ax.plot(estela_x, estela_y, color='blue', alpha=0.5)
    
    #aquí cambio los colores
    linea.set_color('blue')
    bola.set_color('blue')

    return linea, bola, estela

animacion = FuncAnimation(fig, update, frames = len(puntos_t), interval = 50, blit=True)


plt.show()





