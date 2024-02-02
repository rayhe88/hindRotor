#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
    Author  : Raymundo Hernandez-Esparza
    Date    : March, 2023
    Project : hRotor

    This little program was created during my postdoc stay at Argonne. It will
    be part of the hRotor project, in the future I hope to write inside the main
    code.

    The code cleans the energies.dat file that was proporcionated by the group of
    FQOT-BUAP, sets to zero degrees the value of the minimum energy, and deletes
    if there are repetead values.

    for example:
       0.0 degrees -> 360.0 degrees,
    -180.0 degrees -> 180.0 degrees.

    If there are two values the program takes the mean value, and delete the last
    value in energy and theta.

    All the new information is dumped into a file with the extension 'ndatx'.
'''
import sys
import os
import numpy as np
import math
import matplotlib.pyplot as plt

def printing(theta, energy):
    ''' Printing function, only to check if the cleaning was correct.
    '''
    xtics = np.linspace(0,360,13)
    plt.figure()
    theta = np.append(theta,360)
    energy = np.append(energy,energy[0])
    plt.plot(theta, energy, 'bo', theta, energy, 'k')
    plt.xlabel("theta / degree")
    plt.ylabel("Energy / Hartree")
    plt.title("Rotational profile")
    plt.xlim(0,360)
    plt.xticks(xtics)
    plt.grid(True)
    plt.show()

def getMin(array):
    ''' Function that gets the minimun value of an array and also the index
    '''
    amin = 1.E9;

    for idx, x in enumerate(array):
        if x < amin:
            amin = x;
            imin = idx

    return imin, amin


def reorderData(theta, energy):
    ''' Function that check if there are repetead values and set to zero degrees
        the energy of the minimum.
    '''
    if len(theta) != len(energy):
        print("Error")
        raise(Error)

    ndat = len(theta)

    theta = theta - theta[0];

    # First part: Delete the repetead values
    if math.isclose(theta[-1], 360):
        print("\n\033[31;49;1m Repeated values!! It is necessary to remove it.\033[0m")
        val = (energy[0] + energy[-1])/2.0
        energy[0] = val;
        energy = energy[:-1]
        theta = theta[:-1]
        ndat = len(theta)
        print(" new ndat",ndat)

    imin, val = getMin(energy)

    # Get the index for the minimum energy
    print(" index for min: ",imin)
    energy2 = np.zeros(ndat)

    # Second part: Set to zero degrees the energy.
    for i in range(ndat):
        new_i = i + imin
        if new_i >= ndat:
            new_i -= ndat
        energy2[i] = energy[new_i]

    # Printing before and after set zero
    printing(theta,energy)
    printing(theta,energy2)

    return theta, energy2


def writeData(fname, theta, energy):
    ''' Function that dump the values of theta, and energy into a new file
        with extension 'ndatx'.
    '''
    subnames = fname.split('.')

    name = subnames[0] + '.ndatx'
    print("\033[32;49;1m The file {} will be created!\033[0m\n".format(name))
    with open(name, 'w') as fout:
        for i in range(len(theta)):
            print("{0:7.2f}   {1:14.9f}".format(theta[i],energy[i]),file=fout)

def cleanData(fname):
    ''' Main function here we read the original file '*-energies.dat',
        the values are store into two list, after convert to a numpy
        array. It calls the reorderData file with the numpy arrays for
        energy and theta, and finally calls the writeData function to dump
        into a new file.
    '''
    finp = open(fname, 'r')
    lines = finp.readlines()
    finp.close()

    theta = []
    energy = []

    for line in lines:
        line = line.rstrip()
        if line.count('#') == 0 :
            sublines = line.split()
            theta.append(float(sublines[0]))
            energy.append(float(sublines[1]))

    mat_theta = np.array(theta)
    mat_energy = np.array(energy)

    th, en = reorderData(mat_theta, mat_energy)

    writeData(fname, th, en)


if __name__ == "__main__":

    cleanData(sys.argv[1])

