#!/usr/bin/env python3

#Pendulo doble
#Versión 4

import numpy as np
from runge_kutta import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sys import exit as cerrar

#Propiedades del pendulo
print("El siguiente programa simula un pendulo doble dadas unas condiciones iniciales dadas.")
print("Ingrese las propiedades del pendulo")
print("Intente no exagerar")

masa_1: float = float(input("Ingrese el valor de la masa 1 en kilogramos: "))
masa_2: float = float(input("Ingrese el valor de la masa 2 en kilogramos: "))
largo_1: float = float(input("Ingrese el largo del pendulo 1 en metros: "))
largo_2: float = float(input("Ingrese el largo del pendulo 2 en metros: "))
g: float = 9.8

for i in (masa_1, masa_2, largo_1, largo_2):
    if i <= 0: 
        print(">:(")
        cerrar()


def pendulo_doble(t, mov):

    """Este es el sistema de ecuaciones diferenciales de primer orden asociado a
    las ecuaciones de movimiento del pendulo doble
    Aquí t es el tiempo y mov es un arreglo que contiene los angulos y velocidades
    angulares de ambos pendulos
    En general, debería ser, mov = [theta1, omega1, theta2, omega2]"""

    theta_1 = mov[0]
    omega_1 = mov[1]
    theta_2 = mov[2]
    omega_2 = mov[3]

    #Esto es solo para simplificar notación
    masa_f = (masa_1 + masa_2) / masa_2
    delta = theta_1 - theta_2
    div = masa_f - np.power(np.cos(delta), 2)
    
    #Esta es la parte importante
    var_1 = omega_1
    var_2 = (g * (np.sin(theta_2) * np.cos(delta) - masa_f * np.sin(theta_1) ) - np.sin(delta) * (largo_1 * np.power(omega_1, 2) * np.cos(delta) + largo_2 * np.power(omega_2, 2) ) ) / (largo_1 * div)
    var_3 = omega_2
    var_4 = (g * masa_f * (np.sin(theta_1) * np.cos(delta) - np.sin(theta_2) ) + np.sin(delta) * (largo_1 * masa_f * np.power(omega_1, 2) + largo_2 * np.power(omega_2, 2) * np.cos(delta) ) ) / (largo_2 * div)

    return np.array([var_1, var_2, var_3, var_4 ])


#condiciones iniciales
print("\n...\n")
print("Ingrese las condiciones iniciales para el pendulo doble")

theta_1_0: float = float(input("Ingrese el angulo inicial del primer pendulo en radianes: "))
theta_2_0: float = float(input("Ingrese el angulo inicial del segundo pendulo en radianes: "))
omega_1_0: float = float(input("Ingrese la velocidad angular inicial del primer pendulo en radianes: "))
omega_2_0: float = float(input("Ingrese la velocidad angular inicial del segundo pendulo en radianes: "))

#Aquí se obtienen los valores para simular el pendulo con el metodo de Runge Kutta
puntos_t, puntos_y = RK(pendulo_doble, 0, [theta_1_0, omega_1_0, theta_2_0, omega_2_0], 200)

#Configurar la animación
fig, ax = plt.subplots()
ax.set_xlim(-1.2 * (largo_1 + largo_2) , 1.2 * (largo_1 + largo_2))
ax.set_ylim(-1.2 * (largo_1 + largo_2) , 1.2 * (largo_1 + largo_2))
linea_1, = ax.plot([], [], lw = 2) 
bola_1 = ax.scatter([], [])
linea_2, = ax.plot([], [], lw = 2)
bola_2 = ax.scatter([], [])

# Listas para almacenar las posiciones anteriores del pendulo
estela_1_x = []
estela_1_y = []
estela_2_x = []
estela_2_y = []

#Función que actualiza los frames
def update(frame):
    #pendulo 1
    x_1 = [0, largo_1 * np.sin(puntos_y[frame, 0])]
    y_1 = [0, - largo_1 * np.cos(puntos_y[frame, 0])]
    linea_1.set_data(x_1, y_1)

    estela_1_x.append(x_1[1])
    estela_1_y.append(y_1[1])

    bola_1_x = largo_1 * np.sin(puntos_y[frame, 0]) 
    bola_1_y = -largo_1 * np.cos(puntos_y[frame, 0])
    bola_1.set_offsets([[bola_1_x, bola_1_y]])

    #pendulo 2
    x_2 = [largo_1 * np.sin(puntos_y[frame, 0]), largo_1 * np.sin(puntos_y[frame, 0]) + largo_2 * np.sin(puntos_y[frame, 2])]
    y_2 = [- largo_1 * np.cos(puntos_y[frame, 0]), - largo_1 * np.cos(puntos_y[frame, 0]) - largo_2 * np.cos(puntos_y[frame, 2])]
    linea_2.set_data(x_2, y_2)

    estela_2_x.append(x_2[1])
    estela_2_y.append(y_2[1])

    bola_2_x = largo_1 * np.sin(puntos_y[frame, 0]) + largo_2 * np.sin(puntos_y[frame, 2])
    bola_2_y = - largo_1 * np.cos(puntos_y[frame, 0]) - largo_2 * np.cos(puntos_y[frame, 2])
    bola_2.set_offsets([[bola_2_x, bola_2_y]])

    estela_1, = ax.plot(estela_1_x, estela_1_y, color = 'cyan', alpha= 0.5)
    estela_2, = ax.plot(estela_2_x, estela_2_y, color = 'blue', alpha = 0.5)
    
    #aquí cambio los colores
    linea_1.set_color('blue')
    bola_1.set_color('blue')
    linea_2.set_color('blue')
    bola_2.set_color('blue')

    return linea_1, bola_1, linea_2, bola_2, estela_1, estela_2


animacion = FuncAnimation(fig, update, frames = len(puntos_t), interval = 50, blit=True)


plt.show()
















