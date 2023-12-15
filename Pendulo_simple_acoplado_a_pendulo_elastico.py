#!/usr/bin/env python3

#Pendulo simple acoplado a pendulo elastico
#Versión 1

import numpy as np
from runge_kutta import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sys import exit as cerrar

#Propiedades del pendulo
print("El siguiente programa simula un pendulo simple acoplado a un pendulo elastico dadas unas condiciones iniciales dadas.")
print("Ingrese las propiedades del pendulo")
print("Intente no exagerar")


m_1: float = float(input("Ingrese el valor de la masa 1 en kilogramos: "))
m_2: float = float(input("Ingrese el valor de la masa 2 en kilogramos: "))
l_0: float = float(input("Ingrese el largo del pendulo 1 en metros: "))
l_2: float = float(input("Ingrese el largo del pendulo 2 en metros: "))
k: float = float(input("inserte el valor de la constante elastica: ")) #constante elastica
g: float = 9.8

for i in (m_1, m_2, l_0, l_2, k):
    if i <= 0: 
        print(">:(")
        cerrar()

l_1 = l_0 + ((m_1 + m_2)/k)*g

print(f"El largo corregido del resorte es: {l_1}")

def pendulo_doble(t, mov):

    """Este es el sistema de ecuaciones diferenciales de primer orden asociado a
    las ecuaciones de movimientos de movimiento del sistema formado por un pendulo simple acoplado a un pendulo elastico
    Aquí t es el tiempo y mov es un arreglo que contiene los angulos y velocidades
    angulares de ambos pendulos
    En general, debería ser, mov = [theta1, omega1, theta2, omega2]"""

    theta_1 = mov[0]
    omega_1 = mov[1]
    theta_2 = mov[2]
    omega_2 = mov[3]
    deformacion = mov[4]
    V_radial = mov[5]
    #Esto es solo para simplificar notación
    
    #Esta es la parte importante
    var_1 = omega_1
    var_2 = (-2*V_radial*m_1**2*omega_1 + 2*V_radial*m_1*m_2*omega_1*np.sin(theta_1 - theta_2)**2 + 2*V_radial*m_1*m_2*omega_1*np.cos(theta_1 - theta_2)**2 - 4*V_radial*m_1*m_2*omega_1 + 2*V_radial*m_2**2*omega_1*np.sin(theta_1 - theta_2)**2 + 2*V_radial*m_2**2*omega_1*np.cos(theta_1 - theta_2)**2 - 2*V_radial*m_2**2*omega_1 - g*m_1**2*np.sin(theta_1) + g*m_1*m_2*(np.sin(theta_1 - 2*theta_2) + np.sin(3*theta_1 - 2*theta_2))/4 + g*m_1*m_2*np.sin(theta_1)*np.sin(theta_1 - theta_2)**2 - 2*g*m_1*m_2*np.sin(theta_1) + g*m_1*m_2*np.sin(theta_2)*np.cos(theta_1 - theta_2) + g*m_2**2*(np.sin(theta_1 - 2*theta_2) + np.sin(3*theta_1 - 2*theta_2))/4 + g*m_2**2*np.sin(theta_1)*np.sin(theta_1 - theta_2)**2 - g*m_2**2*np.sin(theta_1) + g*m_2**2*np.sin(theta_2)*np.cos(theta_1 - theta_2) - k*m_2*deformacion*np.sin(2*theta_1 - 2*theta_2)/2 - l_2*m_1*m_2*omega_2**2*np.sin(theta_1 - theta_2) + l_2*m_2**2*omega_2**2*np.sin(theta_1 - theta_2)**3 + l_2*m_2**2*omega_2**2*np.sin(theta_1 - theta_2)*np.cos(theta_1 - theta_2)**2 - l_2*m_2**2*omega_2**2*np.sin(theta_1 - theta_2))/(m_1*(l_1*m_1 + l_1*m_2 + m_1*deformacion + m_2*deformacion))
    var_3 = omega_2
    var_4 = ((k * deformacion)/(m_1 * l_2)) * np.sin(theta_1-theta_2) 
    var_5 = V_radial
    var_6 = (g*m_1**2*np.cos(theta_1) - g*m_1*m_2*(np.cos(theta_1 - 2*theta_2) - np.cos(3*theta_1 - 2*theta_2))/4 + g*m_1*m_2*np.sin(theta_2)*np.sin(theta_1 - theta_2) - g*m_1*m_2*np.cos(theta_1)*np.cos(theta_1 - theta_2)**2 + 2*g*m_1*m_2*np.cos(theta_1) - g*m_2**2*(np.cos(theta_1 - 2*theta_2) - np.cos(3*theta_1 - 2*theta_2))/4 + g*m_2**2*np.sin(theta_2)*np.sin(theta_1 - theta_2) - g*m_2**2*np.cos(theta_1)*np.cos(theta_1 - theta_2)**2 + g*m_2**2*np.cos(theta_1) - k*m_1*deformacion + k*m_2*deformacion*np.cos(theta_1 - theta_2)**2 - k*m_2*deformacion + l_1*m_1**2*omega_1**2 - l_1*m_1*m_2*omega_1**2*np.sin(theta_1 - theta_2)**2 - l_1*m_1*m_2*omega_1**2*np.cos(theta_1 - theta_2)**2 + 2*l_1*m_1*m_2*omega_1**2 - l_1*m_2**2*omega_1**2*np.sin(theta_1 - theta_2)**2 - l_1*m_2**2*omega_1**2*np.cos(theta_1 - theta_2)**2 + l_1*m_2**2*omega_1**2 + l_2*m_1*m_2*omega_2**2*np.cos(theta_1 - theta_2) - l_2*m_2**2*omega_2**2*np.sin(theta_1 - theta_2)**2*np.cos(theta_1 - theta_2) - l_2*m_2**2*omega_2**2*np.cos(theta_1 - theta_2)**3 + l_2*m_2**2*omega_2**2*np.cos(theta_1 - theta_2) + m_1**2*omega_1**2*deformacion - m_1*m_2*omega_1**2*deformacion*np.sin(theta_1 - theta_2)**2 - m_1*m_2*omega_1**2*deformacion*np.cos(theta_1 - theta_2)**2 + 2*m_1*m_2*omega_1**2*deformacion - m_2**2*omega_1**2*deformacion*np.sin(theta_1 - theta_2)**2 - m_2**2*omega_1**2*deformacion*np.cos(theta_1 - theta_2)**2 + m_2**2*omega_1**2*deformacion)/(m_1*(m_1 + m_2))                         



    return np.array([var_1, var_2, var_3, var_4, var_5, var_6])


#condiciones iniciales
print("\n...\n")
print("Ingrese las condiciones iniciales para el pendulo supremo")

theta_1_0: float = float(input("Ingrese el angulo inicial del primer pendulo en radianes: "))
theta_2_0: float = float(input("Ingrese el angulo inicial del segundo pendulo en radianes: "))
omega_1_0: float = float(input("Ingrese la velocidad angular inicial del primer pendulo en radianes: "))
omega_2_0: float = float(input("Ingrese la velocidad angular inicial del segundo pendulo en radianes: "))
deformacion_0: float = float(input("Ingrese la deformacion inicial: ")) 
v_radial_0: float = float(input("Ingrese la velocidad radial inicial: ")) 


#Aquí se obtienen los valores para simular el pendulo con el metodo de Runge Kutta
puntos_t, puntos_y = RK6(pendulo_doble, 0, [theta_1_0, omega_1_0, theta_2_0, omega_2_0, deformacion_0, v_radial_0], 200)

#Configurar la animación
fig, ax = plt.subplots()
ax.set_xlim(-1.2 * (l_1 + l_2) , 1.2 * (l_1 + l_2))
ax.set_ylim(-1.2 * (l_1 + l_2) , 1.2 * (l_1 + l_2))
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
    x_1 = [0, (l_1 + puntos_y[frame, 4]) * np.sin(puntos_y[frame, 0])]
    y_1 = [0, - (l_1 + puntos_y[frame, 4]) * np.cos(puntos_y[frame, 0])]
    linea_1.set_data(x_1, y_1)

    estela_1_x.append(x_1[1])
    estela_1_y.append(y_1[1])

    bola_1_x = (l_1 + puntos_y[frame, 4])* np.sin(puntos_y[frame, 0]) 
    bola_1_y = -(l_1 + puntos_y[frame, 4]) * np.cos(puntos_y[frame, 0])
    bola_1.set_offsets([[bola_1_x, bola_1_y]])

    #pendulo 2
    x_2 = [(l_1 + puntos_y[frame, 4]) * np.sin(puntos_y[frame, 0]), (l_1 + puntos_y[frame, 4]) * np.sin(puntos_y[frame, 0]) + l_2 * np.sin(puntos_y[frame, 2])]
    y_2 = [- (l_1 + puntos_y[frame, 4]) * np.cos(puntos_y[frame, 0]), - (l_1 + puntos_y[frame, 4]) * np.cos(puntos_y[frame, 0]) - l_2 * np.cos(puntos_y[frame, 2])]
    linea_2.set_data(x_2, y_2)

    estela_2_x.append(x_2[1])
    estela_2_y.append(y_2[1])

    bola_2_x = (l_1 + puntos_y[frame, 4]) * np.sin(puntos_y[frame, 0]) + l_2 * np.sin(puntos_y[frame, 2])
    bola_2_y = - (l_1 + puntos_y[frame, 4]) * np.cos(puntos_y[frame, 0]) - l_2 * np.cos(puntos_y[frame, 2])
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
