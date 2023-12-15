#!/usr/bin/env python3

#3 pendulos triples
#Versión 1

import numpy as np
from runge_kutta import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sys import exit as cerrar


#Propiedades del pendulo
print("El siguiente programa simula 3 pendulos triples dadas unas condiciones iniciales dadas variando el tercer angulo inicial de todos en 0.000001 radianes")
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
puntos_1_t, puntos_1_y = RK6(pendulo_triple, 0, [theta_1_0, omega_1_0, theta_2_0, omega_2_0, theta_3_0, omega_3_0], 200)
puntos_2_t, puntos_2_y = RK6(pendulo_triple, 0, [theta_1_0, omega_1_0, theta_2_0, omega_2_0, theta_3_0 + 0.000001, omega_3_0], 200)
puntos_3_t, puntos_3_y = RK6(pendulo_triple, 0, [theta_1_0, omega_1_0, theta_2_0, omega_2_0, theta_3_0 + 0.000002, omega_3_0], 200)


#Configurar la animación
fig, ax = plt.subplots()
ax.set_xlim(-1.2 * (l_1 + l_2 + l_3) , 1.2 * (l_1 + l_2 + l_3))
ax.set_ylim(-1.2 * (l_1 + l_2 + l_3) , 1.2 * (l_1 + l_2 + l_3))
#1
linea_1_1, = ax.plot([], [], lw = 2) 
bola_1_1 = ax.scatter([], [])
linea_1_2, = ax.plot([], [], lw = 2)
bola_1_2 = ax.scatter([], [])
linea_1_3, = ax.plot([], [], lw = 2)
bola_1_3 = ax.scatter([], [])
#2
linea_2_1, = ax.plot([], [], lw = 2) 
bola_2_1 = ax.scatter([], [])
linea_2_2, = ax.plot([], [], lw = 2)
bola_2_2 = ax.scatter([], [])
linea_2_3, = ax.plot([], [], lw = 2)
bola_2_3 = ax.scatter([], [])
#3
linea_3_1, = ax.plot([], [], lw = 2) 
bola_3_1 = ax.scatter([], [])
linea_3_2, = ax.plot([], [], lw = 2)
bola_3_2 = ax.scatter([], [])
linea_3_3, = ax.plot([], [], lw = 2)
bola_3_3 = ax.scatter([], [])


# Listas para almacenar las posiciones anteriores del pendulo
#1
estela_1_1_x = []
estela_1_1_y = []
estela_1_2_x = []
estela_1_2_y = []
estela_1_3_x = []
estela_1_3_y = []
#2
estela_2_1_x = []
estela_2_1_y = []
estela_2_2_x = []
estela_2_2_y = []
estela_2_3_x = []
estela_2_3_y = []
#3
estela_3_1_x = []
estela_3_1_y = []
estela_3_2_x = []
estela_3_2_y = []
estela_3_3_x = []
estela_3_3_y = []


#Función que actualiza los frames
def update(frame):
    #1
    #pendulo 1
    x_1_1 = [0, l_1 * np.sin(puntos_1_y[frame, 0])]
    y_1_1 = [0, - l_1 * np.cos(puntos_1_y[frame, 0])]
    linea_1_1.set_data(x_1_1, y_1_1)

    estela_1_1_x.append(x_1_1[1])
    estela_1_1_y.append(y_1_1[1])

    bola_1_1_x = l_1 * np.sin(puntos_1_y[frame, 0]) 
    bola_1_1_y = -l_1 * np.cos(puntos_1_y[frame, 0])
    bola_1_1.set_offsets([[bola_1_1_x, bola_1_1_y]])

    #pendulo 2
    x_1_2 = [l_1 * np.sin(puntos_1_y[frame, 0]), l_1 * np.sin(puntos_1_y[frame, 0]) + l_2 * np.sin(puntos_1_y[frame, 2])]
    y_1_2 = [- l_1 * np.cos(puntos_1_y[frame, 0]), - l_1 * np.cos(puntos_1_y[frame, 0]) - l_2 * np.cos(puntos_1_y[frame, 2])]
    linea_1_2.set_data(x_1_2, y_1_2)

    estela_1_2_x.append(x_1_2[1])
    estela_1_2_y.append(y_1_2[1])

    bola_1_2_x = l_1 * np.sin(puntos_1_y[frame, 0]) + l_2 * np.sin(puntos_1_y[frame, 2])
    bola_1_2_y = - l_1 * np.cos(puntos_1_y[frame, 0]) - l_2 * np.cos(puntos_1_y[frame, 2])
    bola_1_2.set_offsets([[bola_1_2_x, bola_1_2_y]])

    #pendulo 3
    x_1_3 = [l_1 * np.sin(puntos_1_y[frame, 0]) + l_2 * np.sin(puntos_1_y[frame, 2]), l_1 * np.sin(puntos_1_y[frame, 0]) + l_2 * np.sin(puntos_1_y[frame, 2]) + l_3 * np.sin(puntos_1_y[frame, 4])]
    y_1_3 =  [- l_1 * np.cos(puntos_1_y[frame, 0]) - l_2 * np.cos(puntos_1_y[frame, 2]), - l_1 * np.cos(puntos_1_y[frame, 0]) - l_2 * np.cos(puntos_1_y[frame, 2]) - l_3 * np.cos(puntos_1_y[frame, 4])]
    linea_1_3.set_data(x_1_3, y_1_3)

    estela_1_3_x.append(x_1_3[1])
    estela_1_3_y.append(y_1_3[1])

    bola_1_3_x = l_1 * np.sin(puntos_1_y[frame, 0]) + l_2 * np.sin(puntos_1_y[frame, 2]) + l_3 * np.sin(puntos_1_y[frame, 4])
    bola_1_3_y = - l_1 * np.cos(puntos_1_y[frame, 0]) - l_2 * np.cos(puntos_1_y[frame, 2]) - l_3 * np.cos(puntos_1_y[frame, 4])
    bola_1_3.set_offsets([[bola_1_3_x, bola_1_3_y]])

    estela_1_1, = ax.plot(estela_1_1_x, estela_1_1_y, color = 'cyan', alpha= 0.5)
    estela_1_2, = ax.plot(estela_1_2_x, estela_1_2_y, color = 'blue', alpha = 0.5)
    estela_1_3, = ax.plot(estela_1_3_x, estela_1_3_y, color = 'gray', alpha = 0.5)

    if len(estela_1_1_x) >= 100:
        estela_1_1_x.remove(estela_1_1_x[0])
        estela_1_2_x.remove(estela_1_2_x[0])
        estela_1_3_x.remove(estela_1_3_x[0])
        estela_1_1_y.remove(estela_1_1_y[0])
        estela_1_2_y.remove(estela_1_2_y[0])
        estela_1_3_y.remove(estela_1_3_y[0])




    #2
    #pendulo 1
    x_2_1 = [0, l_1 * np.sin(puntos_2_y[frame, 0])]
    y_2_1 = [0, - l_1 * np.cos(puntos_2_y[frame, 0])]
    linea_2_1.set_data(x_2_1, y_2_1)

    estela_2_1_x.append(x_2_1[1])
    estela_2_1_y.append(y_2_1[1])

    bola_2_1_x = l_1 * np.sin(puntos_2_y[frame, 0]) 
    bola_2_1_y = -l_1 * np.cos(puntos_2_y[frame, 0])
    bola_2_1.set_offsets([[bola_2_1_x, bola_2_1_y]])

    #pendulo 2
    x_2_2 = [l_1 * np.sin(puntos_2_y[frame, 0]), l_1 * np.sin(puntos_2_y[frame, 0]) + l_2 * np.sin(puntos_2_y[frame, 2])]
    y_2_2 = [- l_1 * np.cos(puntos_2_y[frame, 0]), - l_1 * np.cos(puntos_2_y[frame, 0]) - l_2 * np.cos(puntos_2_y[frame, 2])]
    linea_2_2.set_data(x_2_2, y_2_2)

    estela_2_2_x.append(x_2_2[1])
    estela_2_2_y.append(y_2_2[1])

    bola_2_2_x = l_1 * np.sin(puntos_2_y[frame, 0]) + l_2 * np.sin(puntos_2_y[frame, 2])
    bola_2_2_y = - l_1 * np.cos(puntos_2_y[frame, 0]) - l_2 * np.cos(puntos_2_y[frame, 2])
    bola_2_2.set_offsets([[bola_2_2_x, bola_2_2_y]])

    #pendulo 3
    x_2_3 = [l_1 * np.sin(puntos_2_y[frame, 0]) + l_2 * np.sin(puntos_2_y[frame, 2]), l_1 * np.sin(puntos_2_y[frame, 0]) + l_2 * np.sin(puntos_2_y[frame, 2]) + l_3 * np.sin(puntos_2_y[frame, 4])]
    y_2_3 =  [- l_1 * np.cos(puntos_2_y[frame, 0]) - l_2 * np.cos(puntos_2_y[frame, 2]), - l_1 * np.cos(puntos_2_y[frame, 0]) - l_2 * np.cos(puntos_2_y[frame, 2]) - l_3 * np.cos(puntos_2_y[frame, 4])]
    linea_2_3.set_data(x_2_3, y_2_3)

    estela_2_3_x.append(x_2_3[1])
    estela_2_3_y.append(y_2_3[1])

    bola_2_3_x = l_1 * np.sin(puntos_2_y[frame, 0]) + l_2 * np.sin(puntos_2_y[frame, 2]) + l_3 * np.sin(puntos_2_y[frame, 4])
    bola_2_3_y = - l_1 * np.cos(puntos_2_y[frame, 0]) - l_2 * np.cos(puntos_2_y[frame, 2]) - l_3 * np.cos(puntos_2_y[frame, 4])
    bola_2_3.set_offsets([[bola_2_3_x, bola_2_3_y]])

    estela_2_1, = ax.plot(estela_2_1_x, estela_2_1_y, color = 'darksalmon', alpha= 0.5)
    estela_2_2, = ax.plot(estela_2_2_x, estela_2_2_y, color = 'rosybrown', alpha = 0.5)
    estela_2_3, = ax.plot(estela_2_3_x, estela_2_3_y, color = 'firebrick', alpha = 0.5)

    if len(estela_2_1_x) >= 100:
        estela_2_1_x.remove(estela_2_1_x[0])
        estela_2_2_x.remove(estela_2_2_x[0])
        estela_2_3_x.remove(estela_2_3_x[0])
        estela_2_1_y.remove(estela_2_1_y[0])
        estela_2_2_y.remove(estela_2_2_y[0])
        estela_2_3_y.remove(estela_2_3_y[0])


    #3

    #pendulo 1
    x_3_1 = [0, l_1 * np.sin(puntos_3_y[frame, 0])]
    y_3_1 = [0, - l_1 * np.cos(puntos_3_y[frame, 0])]
    linea_3_1.set_data(x_3_1, y_3_1)

    estela_3_1_x.append(x_3_1[1])
    estela_3_1_y.append(y_3_1[1])

    bola_3_1_x = l_1 * np.sin(puntos_3_y[frame, 0]) 
    bola_3_1_y = -l_1 * np.cos(puntos_3_y[frame, 0])
    bola_3_1.set_offsets([[bola_3_1_x, bola_3_1_y]])

    #pendulo 2
    x_3_2 = [l_1 * np.sin(puntos_3_y[frame, 0]), l_1 * np.sin(puntos_3_y[frame, 0]) + l_2 * np.sin(puntos_3_y[frame, 2])]
    y_3_2 = [- l_1 * np.cos(puntos_3_y[frame, 0]), - l_1 * np.cos(puntos_3_y[frame, 0]) - l_2 * np.cos(puntos_3_y[frame, 2])]
    linea_3_2.set_data(x_3_2, y_3_2)

    estela_3_2_x.append(x_3_2[1])
    estela_3_2_y.append(y_3_2[1])

    bola_3_2_x = l_1 * np.sin(puntos_3_y[frame, 0]) + l_2 * np.sin(puntos_3_y[frame, 2])
    bola_3_2_y = - l_1 * np.cos(puntos_3_y[frame, 0]) - l_2 * np.cos(puntos_3_y[frame, 2])
    bola_3_2.set_offsets([[bola_3_2_x, bola_3_2_y]])

    #pendulo 3
    x_3_3 = [l_1 * np.sin(puntos_3_y[frame, 0]) + l_2 * np.sin(puntos_3_y[frame, 2]), l_1 * np.sin(puntos_3_y[frame, 0]) + l_2 * np.sin(puntos_3_y[frame, 2]) + l_3 * np.sin(puntos_3_y[frame, 4])]
    y_3_3 =  [- l_1 * np.cos(puntos_3_y[frame, 0]) - l_2 * np.cos(puntos_3_y[frame, 2]), - l_1 * np.cos(puntos_3_y[frame, 0]) - l_2 * np.cos(puntos_3_y[frame, 2]) - l_3 * np.cos(puntos_3_y[frame, 4])]
    linea_3_3.set_data(x_3_3, y_3_3)

    estela_3_3_x.append(x_3_3[1])
    estela_3_3_y.append(y_3_3[1])

    bola_3_3_x = l_1 * np.sin(puntos_3_y[frame, 0]) + l_2 * np.sin(puntos_3_y[frame, 2]) + l_3 * np.sin(puntos_3_y[frame, 4])
    bola_3_3_y = - l_1 * np.cos(puntos_3_y[frame, 0]) - l_2 * np.cos(puntos_3_y[frame, 2]) - l_3 * np.cos(puntos_3_y[frame, 4])
    bola_3_3.set_offsets([[bola_3_3_x, bola_3_3_y]])

    estela_3_1, = ax.plot(estela_3_1_x, estela_3_1_y, color = 'mediumturquoise', alpha= 0.5)
    estela_3_2, = ax.plot(estela_3_2_x, estela_3_2_y, color = 'mediumaquamarine', alpha = 0.5)
    estela_3_3, = ax.plot(estela_3_3_x, estela_3_3_y, color = 'mediumseagreen', alpha = 0.5)

    if len(estela_3_1_x) >= 100:
        estela_3_1_x.remove(estela_3_1_x[0])
        estela_3_2_x.remove(estela_3_2_x[0])
        estela_3_3_x.remove(estela_3_3_x[0])
        estela_3_1_y.remove(estela_3_1_y[0])
        estela_3_2_y.remove(estela_3_2_y[0])
        estela_3_3_y.remove(estela_3_3_y[0])



    #aquí cambio los colores
    #1
    linea_1_1.set_color('blue')
    bola_1_1.set_color('blue')
    linea_1_2.set_color('blue')
    bola_1_2.set_color('blue')
    linea_1_3.set_color('blue')
    bola_1_3.set_color('blue')
    #2
    linea_2_1.set_color('firebrick')
    bola_2_1.set_color('firebrick')
    linea_2_2.set_color('firebrick')
    bola_2_2.set_color('firebrick')
    linea_2_3.set_color('firebrick')
    bola_2_3.set_color('firebrick')
    #3
    linea_3_1.set_color('mediumseagreen')
    bola_3_1.set_color('mediumseagreen')
    linea_3_2.set_color('mediumseagreen')
    bola_3_2.set_color('mediumseagreen')
    linea_3_3.set_color('mediumseagreen')
    bola_3_3.set_color('mediumseagreen')




    return linea_1_1, bola_1_1, linea_1_2, bola_1_2, linea_1_3, bola_1_3, estela_1_1, estela_1_2, estela_1_3, linea_2_1, bola_2_1, linea_2_2, bola_2_2, linea_2_3, bola_2_3, estela_2_1, estela_2_2, estela_2_3, linea_3_1, bola_3_1, linea_3_2, bola_3_2, linea_3_3, bola_3_3, estela_3_1, estela_3_2, estela_3_3


animacion = FuncAnimation(fig, update, frames = len(puntos_1_t), interval = 50, blit=True, repeat = False)

plt.show()

#animacion.save('aaa.mp4', fps=20, metadata=dict(artist='aaaaaa'), bitrate=1200)

#print("listo")
