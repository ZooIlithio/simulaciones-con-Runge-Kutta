#!/usr/bin/env python3

#runge kutta orden 4

import numpy as np

def RK(funcion, t0, y0, tf, h = 0.05):

    """Metodo de Runge kutta de orden 4 para sistemas de ecuaciones diferenciales
    t0 y tf representan el tiempo inicial y final, y0 es el vector con condiciones iniciales
    h es el paso del tiempo, por defecto es 0.05 segundos para simular a 20 fps"""

    nvariables = len(y0) #para generalizar el número de variables en la edo.
    intervalos = int((tf - t0) / h)   #Número de intervalos
    puntost = np.linspace(t0, tf, intervalos + 1)   # Este es un arreglo que contiene los valores de tiempo a evaluar
    puntosy = np.zeros((intervalos + 1, nvariables))  #Este es un arreglo con elementos de dimensión igual a la dimensión de y0 

    y = np.array(y0) #Este es un arreglo con las condiciones iniciales para iterar en Runge Kutta  

    for i in range(intervalos + 1):
        puntosy[i] = y   #Esto va rellenando el arreglo con valores de  y
        k1 = h * np.array(funcion(puntost[i], y))
        k2 = h * np.array(funcion(puntost[i] + h/2, y + k1 * h / 2))
        k3 = h * np.array(funcion(puntost[i] + h/2, y + k2 * h / 2))
        k4 = h * np.array(funcion(puntost[i] + h, y + k3 * h))
        y = y.astype(float)   #Esto en sí no es necesario, pero si se quita a veces hay errores.
        y += (k1 + 2 * k2 + 2 *k3 + k4) * h / 6

        #print(puntosy[i])
    return puntost, puntosy

def RK6(funcion, t0, y0, tf, h = 0.05):

    """Metodo de Runge kutta de orden 6 para sistemas de ecuaciones diferenciales
    t0 y tf representan el tiempo inicial y final, y0 es el vector con condiciones iniciales
    h es el paso del tiempo, por defecto es 0.05 segundos para simular a 20 fps"""

    nvariables = len(y0) #para generalizar el número de variables en la edo.
    intervalos = int((tf - t0) / h)   #Número de intervalos
    puntost = np.linspace(t0, tf, intervalos + 1)   # Este es un arreglo que contiene los valores de tiempo a evaluar
    puntosy = np.zeros((intervalos + 1, nvariables))  #Este es un arreglo con elementos de dimensión igual a la dimensión de y0 

    y = np.array(y0) #Este es un arreglo con las condiciones iniciales para iterar en Runge Kutta  

    for i in range(intervalos + 1):
        puntosy[i] = y   #Esto va rellenando el arreglo con valores de  y
        k1 = h * np.array(funcion(puntost[i], y))
        k2 = h * np.array(funcion(puntost[i] + h, y + k1 * h ))
        k3 = h * np.array(funcion(puntost[i] + h/2, y + (3 * k1 + k2) * h / 8  ) )
        k4 = h * np.array(funcion(puntost[i] + 2 * h / 3, y + (8 * k1 + 2 * k2 + 8 * k3) * h / 27  ) )
        k5 = h * np.array(funcion(puntost[i] + (7 - np.sqrt(21) ) * h / 14 , y + (3 * ( 3 * np.sqrt(21) - 7) * k1 - 8 * (7 - np.sqrt(21)) * k2 + 48 * (7 - np.sqrt(21)) * k3 - 3 * (21 - np.sqrt(21)) * k3 ) * h / 392    ) )
        k6 = h * np.array(funcion(puntost[i] + (7 - np.sqrt(21)) * h / 14 , y + (-5 * (231 + 51 * np.sqrt(21)) * k1 - 40 * ( 7 + np.sqrt(21)) * k2 - 320 * np.sqrt(21) * k3 + 3 * (21 + 121 * np.sqrt(21)) * k4 + 392 * (6 + np.sqrt(21)) * k5  ) * h / 1960  ))
        k7 = h * np.array(funcion(puntost[i] + h , y + (15 * (22 + 7 * np.sqrt(21)) * k1 + 120 * k2 + 40 * ( 7 * np.sqrt(21) - 5 ) * k3 - 63 * (3 * np.sqrt(21) - 2 ) * k4 - 14 * (49 + 9 * np.sqrt(21)) * k5 + 70 * (7 - np.sqrt(21)) * k6   ) * h / 180   ) )

        y = y.astype(float)   
        y += h * (9 * k1 + 64 * k3 + 49 * k5 + 49 * k6 + 9 * k7 )   / 180

        #print(puntosy[i])
    return puntost, puntosy


