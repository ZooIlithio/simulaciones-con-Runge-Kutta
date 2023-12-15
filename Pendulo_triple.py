#!/usr/bin/env python3

#Pendulo triple
#Versión 3

import numpy as np
from runge_kutta import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sys import exit as cerrar


#Propiedades del pendulo
print("El siguiente programa simula un pendulo triple dadas unas condiciones iniciales dadas.")
print("Ingrese las propiedades del pendulo")

m_1: float = float(input("Ingrese el valor de la masa 1 en kilogramos: "))
m_2: float = float(input("Ingrese el valor de la masa 2 en kilogramos: "))
m_3: float = float(input("Ingrese el valor de la masa 3 en kilogramos: "))
l_1: float = float(input("Ingrese el largo del pendulo 1 en metros: "))
l_2: float = float(input("Ingrese el largo del pendulo 2 en metros: "))
l_3: float = float(input("Ingrese el largo del pendulo 3 en metros: "))
g: float = 9.8

for i in (m_1, m_2, m_3, l_1, l_2, l_3):
    if i <= 0: 
        print(">:(")
        cerrar()

def pendulo_triple(t, mov):

    """Este es el sistema de ecuaciones diferenciales de primer orden asociado a
    las ecuaciones de movimiento del pendulo triple
    Aquí t es el tiempo y mov es un arreglo que contiene los angulos y velocidades
    angulares de ambos pendulos
    En general, debería ser, mov = [theta_1, omega_1, theta_2, omega_2, theta_3, omega_3]"""

    theta_1 = mov[0]
    omega_1 = mov[1]
    theta_2 = mov[2]
    omega_2 = mov[3]
    theta_3 = mov[4]
    omega_3 = mov[5]

    #Esto es solo para simplificar notación
    masa_t = m_1 + m_2 + m_3
    div = (1 - (m_3 / (m_2 + m_3) * np.power(np.cos(theta_2 - theta_3) ,2)) )

    
    #Esta es la parte importante
    var_1 = omega_1
    var_2 = (g*m_1*m_2*np.sin(theta_1) - g*m_1*m_3*np.sin(theta_1)*np.cos(theta_2 - theta_3)**2 + g*m_1*m_3*np.sin(theta_1) + g*m_2**2*np.sin(theta_1) - g*m_2**2*np.sin(theta_2)*np.cos(theta_1 - theta_2) - g*m_2*m_3*np.sin(theta_1)*np.cos(theta_2 - theta_3)**2 + 2*g*m_2*m_3*np.sin(theta_1) - 2*g*m_2*m_3*np.sin(theta_2)*np.cos(theta_1 - theta_2) + g*m_2*m_3*np.sin(theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) + g*m_2*m_3*np.sin(theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_2 - theta_3) - g*m_2*m_3*np.sin(theta_3)*np.cos(theta_1 - theta_3) - g*m_3**2*np.sin(theta_1)*np.cos(theta_2 - theta_3)**2 + g*m_3**2*np.sin(theta_1) - g*m_3**2*np.sin(theta_2)*np.cos(theta_1 - theta_2) + g*m_3**2*np.sin(theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) + g*m_3**2*np.sin(theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_2 - theta_3) - g*m_3**2*np.sin(theta_3)*np.cos(theta_1 - theta_3) + l_1*m_2**2*omega_1**2*np.sin(2*theta_1 - 2*theta_2)/2 - l_1*m_2*m_3*omega_1**2*np.sin(theta_1 - theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) - l_1*m_2*m_3*omega_1**2*np.sin(theta_1 - theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_2 - theta_3) + l_1*m_2*m_3*omega_1**2*np.sin(2*theta_1 - 2*theta_2) + l_1*m_2*m_3*omega_1**2*np.sin(2*theta_1 - 2*theta_3)/2 - l_1*m_3**2*omega_1**2*np.sin(theta_1 - theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) - l_1*m_3**2*omega_1**2*np.sin(theta_1 - theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_2 - theta_3) + l_1*m_3**2*omega_1**2*np.sin(2*theta_1 - 2*theta_2)/2 + l_1*m_3**2*omega_1**2*np.sin(2*theta_1 - 2*theta_3)/2 + l_2*m_1*m_2*omega_2**2*np.sin(theta_1 - theta_2) - l_2*m_1*m_3*omega_2**2*np.sin(theta_1 - theta_2)*np.cos(theta_2 - theta_3)**2 + l_2*m_1*m_3*omega_2**2*np.sin(theta_1 - theta_2) + l_2*m_2**2*omega_2**2*np.sin(theta_1 - theta_2) - l_2*m_2*m_3*omega_2**2*(-np.sin(theta_1 - 3*theta_2 + 2*theta_3) + np.sin(theta_1 + theta_2 - 2*theta_3))/4 - l_2*m_2*m_3*omega_2**2*np.sin(theta_1 - theta_2)*np.cos(theta_2 - theta_3)**2 + l_2*m_2*m_3*omega_2**2*np.sin(theta_1 - theta_2) + l_2*m_2*m_3*omega_2**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_3) - l_2*m_3**2*omega_2**2*(-np.sin(theta_1 - 3*theta_2 + 2*theta_3) + np.sin(theta_1 + theta_2 - 2*theta_3))/4 + l_2*m_3**2*omega_2**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_3) + l_3*m_2*m_3*omega_3**2*np.sin(theta_1 - theta_3) - l_3*m_2*m_3*omega_3**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_2) + l_3*m_3**2*omega_3**2*(-np.sin(theta_1 - 2*theta_2 + theta_3) + np.sin(theta_1 + 2*theta_2 - 3*theta_3))/4 - l_3*m_3**2*omega_3**2*np.sin(theta_1 - theta_3)*np.cos(theta_2 - theta_3)**2 + l_3*m_3**2*omega_3**2*np.sin(theta_1 - theta_3) - l_3*m_3**2*omega_3**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_2))/(l_1*(-m_1*m_2 + m_1*m_3*np.cos(theta_2 - theta_3)**2 - m_1*m_3 + m_2**2*np.cos(theta_1 - theta_2)**2 - m_2**2 + 2*m_2*m_3*np.cos(theta_1 - theta_2)**2 - 2*m_2*m_3*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) + m_2*m_3*np.cos(theta_1 - theta_3)**2 + m_2*m_3*np.cos(theta_2 - theta_3)**2 - 2*m_2*m_3 + m_3**2*np.cos(theta_1 - theta_2)**2 - 2*m_3**2*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) + m_3**2*np.cos(theta_1 - theta_3)**2 + m_3**2*np.cos(theta_2 - theta_3)**2 - m_3**2))
    var_3 = omega_2
    var_4 = (-g*m_1*m_2*np.sin(theta_1)*np.cos(theta_1 - theta_2) + g*m_1*m_2*np.sin(theta_2) - g*m_1*m_3*np.sin(theta_1)*np.cos(theta_1 - theta_2) + g*m_1*m_3*np.sin(theta_1)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) + g*m_1*m_3*np.sin(theta_2) - g*m_1*m_3*np.sin(theta_3)*np.cos(theta_2 - theta_3) - g*m_2**2*np.sin(theta_1)*np.cos(theta_1 - theta_2) + g*m_2**2*np.sin(theta_2) - 2*g*m_2*m_3*np.sin(theta_1)*np.cos(theta_1 - theta_2) + g*m_2*m_3*np.sin(theta_1)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) - g*m_2*m_3*np.sin(theta_2)*np.cos(theta_1 - theta_3)**2 + 2*g*m_2*m_3*np.sin(theta_2) + g*m_2*m_3*np.sin(theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3) - g*m_2*m_3*np.sin(theta_3)*np.cos(theta_2 - theta_3) - g*m_3**2*np.sin(theta_1)*np.cos(theta_1 - theta_2) + g*m_3**2*np.sin(theta_1)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) - g*m_3**2*np.sin(theta_2)*np.cos(theta_1 - theta_3)**2 + g*m_3**2*np.sin(theta_2) + g*m_3**2*np.sin(theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3) - g*m_3**2*np.sin(theta_3)*np.cos(theta_2 - theta_3) - l_1*m_1*m_2*omega_1**2*np.sin(theta_1 - theta_2) - l_1*m_1*m_3*omega_1**2*np.sin(theta_1 - theta_2) + l_1*m_1*m_3*omega_1**2*np.sin(theta_1 - theta_3)*np.cos(theta_2 - theta_3) - l_1*m_2**2*omega_1**2*np.sin(theta_1 - theta_2) - l_1*m_2*m_3*omega_1**2*(-np.sin(-3*theta_1 + theta_2 + 2*theta_3) + np.sin(theta_1 + theta_2 - 2*theta_3))/4 + l_1*m_2*m_3*omega_1**2*np.sin(theta_1 - theta_2)*np.cos(theta_1 - theta_3)**2 - 2*l_1*m_2*m_3*omega_1**2*np.sin(theta_1 - theta_2) + l_1*m_2*m_3*omega_1**2*np.sin(theta_1 - theta_3)*np.cos(theta_2 - theta_3) - l_1*m_3**2*omega_1**2*(-np.sin(-3*theta_1 + theta_2 + 2*theta_3) + np.sin(theta_1 + theta_2 - 2*theta_3))/4 + l_1*m_3**2*omega_1**2*np.sin(theta_1 - theta_2)*np.cos(theta_1 - theta_3)**2 - l_1*m_3**2*omega_1**2*np.sin(theta_1 - theta_2) + l_1*m_3**2*omega_1**2*np.sin(theta_1 - theta_3)*np.cos(theta_2 - theta_3) - l_2*m_1*m_2*omega_2**2*np.sin(2*theta_1 - 2*theta_2)/2 + l_2*m_1*m_3*omega_2**2*np.sin(theta_1 - theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) - l_2*m_1*m_3*omega_2**2*np.sin(2*theta_1 - 2*theta_2)/2 + l_2*m_1*m_3*omega_2**2*np.sin(2*theta_2 - 2*theta_3)/2 - l_2*m_2**2*omega_2**2*np.sin(2*theta_1 - 2*theta_2)/2 + l_2*m_2*m_3*omega_2**2*np.sin(theta_1 - theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) - l_2*m_2*m_3*omega_2**2*np.sin(2*theta_1 - 2*theta_2)/2 - l_2*m_2*m_3*omega_2**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3) + l_2*m_2*m_3*omega_2**2*np.sin(2*theta_2 - 2*theta_3)/2 - l_2*m_3**2*omega_2**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3) + l_2*m_3**2*omega_2**2*np.sin(2*theta_2 - 2*theta_3)/2 + l_3*m_1*m_3*omega_3**2*np.sin(theta_2 - theta_3) - l_3*m_2*m_3*omega_3**2*np.sin(theta_1 - theta_3)*np.cos(theta_1 - theta_2) + l_3*m_2*m_3*omega_3**2*np.sin(theta_2 - theta_3) + l_3*m_3**2*omega_3**2*(-np.sin(-2*theta_1 + theta_2 + theta_3) + np.sin(2*theta_1 + theta_2 - 3*theta_3))/4 - l_3*m_3**2*omega_3**2*np.sin(theta_1 - theta_3)*np.cos(theta_1 - theta_2) - l_3*m_3**2*omega_3**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_3)**2 + l_3*m_3**2*omega_3**2*np.sin(theta_2 - theta_3))/(l_2*(-m_1*m_2 + m_1*m_3*np.cos(theta_2 - theta_3)**2 - m_1*m_3 + m_2**2*np.cos(theta_1 - theta_2)**2 - m_2**2 + 2*m_2*m_3*np.cos(theta_1 - theta_2)**2 - 2*m_2*m_3*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) + m_2*m_3*np.cos(theta_1 - theta_3)**2 + m_2*m_3*np.cos(theta_2 - theta_3)**2 - 2*m_2*m_3 + m_3**2*np.cos(theta_1 - theta_2)**2 - 2*m_3**2*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) + m_3**2*np.cos(theta_1 - theta_3)**2 + m_3**2*np.cos(theta_2 - theta_3)**2 - m_3**2)) 
    var_5 = omega_3
    var_6 = (g*m_1*m_2*np.sin(theta_1)*np.cos(theta_1 - theta_2)*np.cos(theta_2 - theta_3) - g*m_1*m_2*np.sin(theta_1)*np.cos(theta_1 - theta_3) - g*m_1*m_2*np.sin(theta_2)*np.cos(theta_2 - theta_3) + g*m_1*m_2*np.sin(theta_3) + g*m_1*m_3*np.sin(theta_1)*np.cos(theta_1 - theta_2)*np.cos(theta_2 - theta_3) - g*m_1*m_3*np.sin(theta_1)*np.cos(theta_1 - theta_3) - g*m_1*m_3*np.sin(theta_2)*np.cos(theta_2 - theta_3) + g*m_1*m_3*np.sin(theta_3) + g*m_2**2*np.sin(theta_1)*np.cos(theta_1 - theta_2)*np.cos(theta_2 - theta_3) - g*m_2**2*np.sin(theta_1)*np.cos(theta_1 - theta_3) + g*m_2**2*np.sin(theta_2)*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3) - g*m_2**2*np.sin(theta_2)*np.cos(theta_2 - theta_3) - g*m_2**2*np.sin(theta_3)*np.cos(theta_1 - theta_2)**2 + g*m_2**2*np.sin(theta_3) + 2*g*m_2*m_3*np.sin(theta_1)*np.cos(theta_1 - theta_2)*np.cos(theta_2 - theta_3) - 2*g*m_2*m_3*np.sin(theta_1)*np.cos(theta_1 - theta_3) + 2*g*m_2*m_3*np.sin(theta_2)*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3) - 2*g*m_2*m_3*np.sin(theta_2)*np.cos(theta_2 - theta_3) - 2*g*m_2*m_3*np.sin(theta_3)*np.cos(theta_1 - theta_2)**2 + 2*g*m_2*m_3*np.sin(theta_3) + g*m_3**2*np.sin(theta_1)*np.cos(theta_1 - theta_2)*np.cos(theta_2 - theta_3) - g*m_3**2*np.sin(theta_1)*np.cos(theta_1 - theta_3) + g*m_3**2*np.sin(theta_2)*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3) - g*m_3**2*np.sin(theta_2)*np.cos(theta_2 - theta_3) - g*m_3**2*np.sin(theta_3)*np.cos(theta_1 - theta_2)**2 + g*m_3**2*np.sin(theta_3) + l_1*m_1*m_2*omega_1**2*np.sin(theta_1 - theta_2)*np.cos(theta_2 - theta_3) - l_1*m_1*m_2*omega_1**2*np.sin(theta_1 - theta_3) + l_1*m_1*m_3*omega_1**2*np.sin(theta_1 - theta_2)*np.cos(theta_2 - theta_3) - l_1*m_1*m_3*omega_1**2*np.sin(theta_1 - theta_3) - l_1*m_2**2*omega_1**2*(-np.sin(-3*theta_1 + 2*theta_2 + theta_3) + np.sin(theta_1 - 2*theta_2 + theta_3))/4 + l_1*m_2**2*omega_1**2*np.sin(theta_1 - theta_2)*np.cos(theta_2 - theta_3) + l_1*m_2**2*omega_1**2*np.sin(theta_1 - theta_3)*np.cos(theta_1 - theta_2)**2 - l_1*m_2**2*omega_1**2*np.sin(theta_1 - theta_3) - l_1*m_2*m_3*omega_1**2*(-np.sin(-3*theta_1 + 2*theta_2 + theta_3) + np.sin(theta_1 - 2*theta_2 + theta_3))/2 + 2*l_1*m_2*m_3*omega_1**2*np.sin(theta_1 - theta_2)*np.cos(theta_2 - theta_3) + 2*l_1*m_2*m_3*omega_1**2*np.sin(theta_1 - theta_3)*np.cos(theta_1 - theta_2)**2 - 2*l_1*m_2*m_3*omega_1**2*np.sin(theta_1 - theta_3) - l_1*m_3**2*omega_1**2*(-np.sin(-3*theta_1 + 2*theta_2 + theta_3) + np.sin(theta_1 - 2*theta_2 + theta_3))/4 + l_1*m_3**2*omega_1**2*np.sin(theta_1 - theta_2)*np.cos(theta_2 - theta_3) + l_1*m_3**2*omega_1**2*np.sin(theta_1 - theta_3)*np.cos(theta_1 - theta_2)**2 - l_1*m_3**2*omega_1**2*np.sin(theta_1 - theta_3) + l_2*m_1*m_2*omega_2**2*(-np.sin(-2*theta_1 + theta_2 + theta_3) + np.sin(2*theta_1 - 3*theta_2 + theta_3))/4 - l_2*m_1*m_2*omega_2**2*np.sin(theta_1 - theta_2)*np.cos(theta_1 - theta_3) - l_2*m_1*m_2*omega_2**2*np.sin(theta_2 - theta_3) + l_2*m_1*m_3*omega_2**2*(-np.sin(-2*theta_1 + theta_2 + theta_3) + np.sin(2*theta_1 - 3*theta_2 + theta_3))/4 - l_2*m_1*m_3*omega_2**2*np.sin(theta_1 - theta_2)*np.cos(theta_1 - theta_3) - l_2*m_1*m_3*omega_2**2*np.sin(theta_2 - theta_3) + l_2*m_2**2*omega_2**2*(-np.sin(-2*theta_1 + theta_2 + theta_3) + np.sin(2*theta_1 - 3*theta_2 + theta_3))/4 - l_2*m_2**2*omega_2**2*np.sin(theta_1 - theta_2)*np.cos(theta_1 - theta_3) + l_2*m_2**2*omega_2**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_2)**2 - l_2*m_2**2*omega_2**2*np.sin(theta_2 - theta_3) + l_2*m_2*m_3*omega_2**2*(-np.sin(-2*theta_1 + theta_2 + theta_3) + np.sin(2*theta_1 - 3*theta_2 + theta_3))/4 - l_2*m_2*m_3*omega_2**2*np.sin(theta_1 - theta_2)*np.cos(theta_1 - theta_3) + 2*l_2*m_2*m_3*omega_2**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_2)**2 - 2*l_2*m_2*m_3*omega_2**2*np.sin(theta_2 - theta_3) + l_2*m_3**2*omega_2**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_2)**2 - l_2*m_3**2*omega_2**2*np.sin(theta_2 - theta_3) - l_3*m_1*m_3*omega_3**2*np.sin(2*theta_2 - 2*theta_3)/2 + l_3*m_2*m_3*omega_3**2*np.sin(theta_1 - theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_2 - theta_3) - l_3*m_2*m_3*omega_3**2*np.sin(2*theta_1 - 2*theta_3)/2 + l_3*m_2*m_3*omega_3**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3) - l_3*m_2*m_3*omega_3**2*np.sin(2*theta_2 - 2*theta_3)/2 + l_3*m_3**2*omega_3**2*np.sin(theta_1 - theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_2 - theta_3) - l_3*m_3**2*omega_3**2*np.sin(2*theta_1 - 2*theta_3)/2 + l_3*m_3**2*omega_3**2*np.sin(theta_2 - theta_3)*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3) - l_3*m_3**2*omega_3**2*np.sin(2*theta_2 - 2*theta_3)/2)/(l_3*(-m_1*m_2 + m_1*m_3*np.cos(theta_2 - theta_3)**2 - m_1*m_3 + m_2**2*np.cos(theta_1 - theta_2)**2 - m_2**2 + 2*m_2*m_3*np.cos(theta_1 - theta_2)**2 - 2*m_2*m_3*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) + m_2*m_3*np.cos(theta_1 - theta_3)**2 + m_2*m_3*np.cos(theta_2 - theta_3)**2 - 2*m_2*m_3 + m_3**2*np.cos(theta_1 - theta_2)**2 - 2*m_3**2*np.cos(theta_1 - theta_2)*np.cos(theta_1 - theta_3)*np.cos(theta_2 - theta_3) + m_3**2*np.cos(theta_1 - theta_3)**2 + m_3**2*np.cos(theta_2 - theta_3)**2 - m_3**2))

    return np.array([var_1, var_2, var_3, var_4, var_5, var_6 ])

#condiciones iniciales
print("\n...\n")
print("Ingrese las condiciones iniciales para el pendulo triple")

theta_1_0: float = float(input("Ingrese el angulo inicial del primer pendulo en radianes: "))
theta_2_0: float = float(input("Ingrese el angulo inicial del segundo pendulo en radianes: "))
theta_3_0: float = float(input("Ingrese el angulo inicial del tercer pendulo en radianes: "))
omega_1_0: float = float(input("Ingrese la velocidad angular inicial del primer pendulo en radianes: "))
omega_2_0: float = float(input("Ingrese la velocidad angular inicial del segundo pendulo en radianes: "))
omega_3_0: float = float(input("Ingrese la velocidad angular inicial del tercer pendulo en radianes: "))

#Aquí se obtienen los valores para simular el pendulo con el metodo de Runge Kutta
puntos_t, puntos_y = RK6(pendulo_triple, 0, [theta_1_0, omega_1_0, theta_2_0, omega_2_0, theta_3_0, omega_3_0], 200)

#Configurar la animación
fig, ax = plt.subplots()
ax.set_xlim(-1.2 * (l_1 + l_2 + l_3) , 1.2 * (l_1 + l_2 + l_3))
ax.set_ylim(-1.2 * (l_1 + l_2 + l_3) , 1.2 * (l_1 + l_2 + l_3))
linea_1, = ax.plot([], [], lw = 2) 
bola_1 = ax.scatter([], [])
linea_2, = ax.plot([], [], lw = 2)
bola_2 = ax.scatter([], [])
linea_3, = ax.plot([], [], lw = 2)
bola_3 = ax.scatter([], [])

# Listas para almacenar las posiciones anteriores del pendulo
estela_1_x = []
estela_1_y = []
estela_2_x = []
estela_2_y = []
estela_3_x = []
estela_3_y = []

#Función que actualiza los frames
def update(frame):
    #pendulo 1
    x_1 = [0, l_1 * np.sin(puntos_y[frame, 0])]
    y_1 = [0, - l_1 * np.cos(puntos_y[frame, 0])]
    linea_1.set_data(x_1, y_1)

    estela_1_x.append(x_1[1])
    estela_1_y.append(y_1[1])

    bola_1_x = l_1 * np.sin(puntos_y[frame, 0]) 
    bola_1_y = -l_1 * np.cos(puntos_y[frame, 0])
    bola_1.set_offsets([[bola_1_x, bola_1_y]])

    #pendulo 2
    x_2 = [l_1 * np.sin(puntos_y[frame, 0]), l_1 * np.sin(puntos_y[frame, 0]) + l_2 * np.sin(puntos_y[frame, 2])]
    y_2 = [- l_1 * np.cos(puntos_y[frame, 0]), - l_1 * np.cos(puntos_y[frame, 0]) - l_2 * np.cos(puntos_y[frame, 2])]
    linea_2.set_data(x_2, y_2)

    estela_2_x.append(x_2[1])
    estela_2_y.append(y_2[1])

    bola_2_x = l_1 * np.sin(puntos_y[frame, 0]) + l_2 * np.sin(puntos_y[frame, 2])
    bola_2_y = - l_1 * np.cos(puntos_y[frame, 0]) - l_2 * np.cos(puntos_y[frame, 2])
    bola_2.set_offsets([[bola_2_x, bola_2_y]])

    #pendulo 3
    x_3 = [l_1 * np.sin(puntos_y[frame, 0]) + l_2 * np.sin(puntos_y[frame, 2]), l_1 * np.sin(puntos_y[frame, 0]) + l_2 * np.sin(puntos_y[frame, 2]) + l_3 * np.sin(puntos_y[frame, 4])]
    y_3 =  [- l_1 * np.cos(puntos_y[frame, 0]) - l_2 * np.cos(puntos_y[frame, 2]), - l_1 * np.cos(puntos_y[frame, 0]) - l_2 * np.cos(puntos_y[frame, 2]) - l_3 * np.cos(puntos_y[frame, 4])]
    linea_3.set_data(x_3, y_3)

    estela_3_x.append(x_3[1])
    estela_3_y.append(y_3[1])

    bola_3_x = l_1 * np.sin(puntos_y[frame, 0]) + l_2 * np.sin(puntos_y[frame, 2]) + l_3 * np.sin(puntos_y[frame, 4])
    bola_3_y = - l_1 * np.cos(puntos_y[frame, 0]) - l_2 * np.cos(puntos_y[frame, 2]) - l_3 * np.cos(puntos_y[frame, 4])
    bola_3.set_offsets([[bola_3_x, bola_3_y]])


    estela_1, = ax.plot(estela_1_x, estela_1_y, color = 'cyan', alpha= 0.5)
    estela_2, = ax.plot(estela_2_x, estela_2_y, color = 'blue', alpha = 0.5)
    estela_3, = ax.plot(estela_3_x, estela_3_y, color = 'gray', alpha = 0.5)
    
    #aquí cambio los colores
    linea_1.set_color('blue')
    bola_1.set_color('blue')
    linea_2.set_color('blue')
    bola_2.set_color('blue')
    linea_3.set_color('blue')
    bola_3.set_color('blue')

    return linea_1, bola_1, linea_2, bola_2, linea_3, bola_3, estela_1, estela_2, estela_3


animacion = FuncAnimation(fig, update, frames = len(puntos_t), interval = 50, blit=True)


plt.show()




