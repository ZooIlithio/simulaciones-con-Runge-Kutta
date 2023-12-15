#!/usr/bin/env python3

#pendulo simple

import numpy as np
from runge_kutta import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#propiedades del pendulo
largo = float(input("Ingrese el largo del pendulo en metros: "))
g = 9.8
masa = 0.1 #Esto no afecta en nada xd, pero por costumbre

def pendulo_simple(t, mov):

    """Este es el sistema de ecuaciones diferenciales de primer orden asociado a la ecuacion del pendulo simple
    aqui mov debe ser un arreglo que contenga el angulo y la velocidad angular"""

    theta = mov[0]
    omega = mov[1]

    
    var1 = omega
    var2 = - g / largo * np.sin(theta)
    return np.array([var1, var2])


#condiciones iniciales
theta0 = float(input("Ingrese el angulo inicial en radianes: "))
omega0 = float(input("Ingrese la velocidad angular inicial en radianes: "))

#Aqui se obtienen valores para simular el pendulo utilizando el metodo de Runge Kutta
puntos_t, puntos_y = RK(pendulo_simple, 0, [theta0, omega0], 200)



# Configurar la animación
fig, ax = plt.subplots()   #Esto es para crear la figura y los ejes
ax.set_xlim(-1.5*largo, 1.5*largo)  #Aqui se ajusta el tamaño de los ejes para que dependa del largo del pendulo 
ax.set_ylim(-1.5*largo, 1.5*largo)  
linea, = ax.plot([], [], lw=2) #Esto crea la linea que une el centro de las coordenadas con lo que sería la posición de la esfera, lw es para el grosor
bola = ax.scatter([], []) #Esto crea un punto de dispersión que simula ser la esfera en el pendulo, inicialmente los puntos estan vacios

# Listas para almacenar las posiciones anteriores del pendulo
estela_x = []
estela_y = []

#Aquí se define una función que actualiza los frames
def update(frame):
    x = [0, largo * np.sin(puntos_y[frame, 0])]  # Aquí se convierte a coordenadas cartesianas.
    y = [0, - largo * np.cos(puntos_y[frame, 0])]
    linea.set_data(x,y)   #Esto actualiza la posición de la linea

    estela_x.append(x[1]) #Esto acumula las posiciones por las que pasa el pendulo
    estela_y.append(y[1])

    bola_x = largo * np.sin(puntos_y[frame, 0]) 
    bola_y = -largo * np.cos(puntos_y[frame, 0])
    bola.set_offsets([[bola_x, bola_y]])  #Esto actualiza la posición del punto de dispersión

    estela, = ax.plot(estela_x, estela_y, color='blue', alpha=0.5)
    
    #aquí cambio los colores
    linea.set_color('blue')
    bola.set_color('blue')

    return linea, bola, estela


#interval es 0.05 para simular a 20fps
#blit = True es para realizar una tecnica de blitting, solo actualiza las partes que han cambiado.
animacion = FuncAnimation(fig, update, frames = len(puntos_t), interval = 0.05, blit=True)


plt.show()
