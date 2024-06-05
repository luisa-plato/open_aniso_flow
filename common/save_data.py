#!/usr/bin/env python 


from fenics import *
import matplotlib.pyplot as plt
import numpy as np
import datetime
import os
import json

def create_dic(name):
    e = datetime.datetime.now()
    time_stamp = e.strftime("%d_%m_%Y_%H_%M_%S")
    path = './data/' + name + "_" + time_stamp

    if not os.path.exists(path):
        os.makedirs(path)
    
    return path
    
def create_vtk(path, function_names):

    if not os.path.exists(path + '/vtk'):
        os.makedirs(path + '/vtk')

    n = len(function_names)
    files = []
    for i in range(n):
        files.append(File(path + '/vtk/' + function_names[i] + '.pvd'))
    return files

def write_to_vtk(vtk_files, functions, t):
    n = len(vtk_files)
    for i in range(n):
        vtk_files[i] << (functions[i], t)

def save_plot(path, name, f, x_axis = []):
    if not os.path.exists(path + '/plots'):
        os.makedirs(path + '/plots')
    if isinstance(f, list):
        if x_axis != []:
            plt.plot(x_axis, f)
        else:
            plt.plot(f)
    else:
        plot(f)
    plt.savefig(path + '/plots/' + name + '.png')
    plt.close()

def write_param_file(path, parameters):
	with open(path + '/parameters.txt', 'w') as file:
		json.dump(parameters, file)

def read_param_file(path):
    with open(path + '/parameters.txt', 'r') as f:
        parameters = json.load(f)
    return parameters