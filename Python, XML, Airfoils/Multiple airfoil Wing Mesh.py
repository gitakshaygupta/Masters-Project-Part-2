# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 10:55:18 2016

@author: jmascolo
"""

import xmltodict
import math
from subprocess import Popen
import numpy as np

#Tolerance for nonzero line length
delta = 0

user_input = raw_input(
    'Enter the name of the XML file containing the wing definition, please\n')

with open(user_input+'.xml') as fd:
    doc = xmltodict.parse(fd.read())

#Create lists to store sectional variables
xle = []
yle = []
zle = []
theta = []
c = []
#tc = []
#camc = []
airfoil_files = []

#Read the number of panel columns and rows, the relative thickness of the reference airfoil and its coordinates file
cols = int(doc['WingData']['cols'])
rows = int(doc['WingData']['rows'])
#tc_ref = float(doc['WingData']['tcref'])
#camc_ref = float(doc['WingData']['camcref'])
# airfoil_file = str(doc['WingData']['airfoil']['@file'])

#Read the geometry variables from the XML file and store them into lists
for item in doc['WingData']['section']:
    xle.append(float(item['xle']))
    yle.append(float(item['yle']))
    zle.append(float(item['zle']))
    theta.append(float(item['theta']))
    c.append(float(item['c']))
    #tc.append(float(item['tc']))
    #camc.append(float(item['camc']))
    airfoil_files.append(str(item['airfoil']['@file']))

#Number of sections
sec_num = len(doc['WingData']['section'])

#Write gmsh file

#Write coordinates and splines of the airfoil sections
with open(user_input+'.geo', 'w') as f:
    for i in range(sec_num):

        #Coordinates of the reference airfoil's upper surface
        upper_coord = []

        #Coordinates of the reference airfoil's lower surface
        lower_coord = []

        with open(airfoil_files[i]) as airfoil_file:
            lines = airfoil_file.readlines()
            lines = [line.split() for line in lines]

            coord = False
            for line in lines:
                if len(line) == 4:
                    if line[0] == 'XU':
                        coord = True
                if coord == True and len(line) == 4 and line[0] != 'XU':
                    upper_coord.append((line[0], line[1]))
                    lower_coord.append((line[2], line[3]))

        #Number of points of one side of the airfoil
        foil_num = len(upper_coord)

        #Compute camber and half-thickness of the reference airfoil
        camber_ref = []
        half_thick_ref = []

        for k in range(foil_num):
            camber_ref.append((float(upper_coord[k][1]) + float(lower_coord[k][1]))/2)
            half_thick_ref.append(
                (float(upper_coord[k][1]) - float(lower_coord[k][1]))/2)
                
        #Write upper surface coordinates
        for j in range(foil_num):
            if (j == 0 or j == foil_num - 1):
                y_u = camber_ref[j] + \
                    half_thick_ref[j]
                f.write('Point('+str(j+1+2*i*len(upper_coord))+') = {'+str(
                    c[i]*float(upper_coord[j][0]))+', 0, '+str(c[i]*y_u+delta)+'};\n')
            else:
                y_u = camber_ref[j] + \
                    half_thick_ref[j]
                f.write('Point('+str(j+1+2*i*len(upper_coord)) +
                        ') = {'+str(c[i]*float(upper_coord[j][0]))+', 0, '+str(c[i]*y_u)+'};\n')

        #Write lower surface coordinates
        offset = len(upper_coord)
        for j in range(foil_num):
            y_l = camber_ref[j] - \
                half_thick_ref[j]
            f.write('Point('+str(j+1+offset+2*i*len(lower_coord)) +
                    ') = {'+str(c[i]*float(lower_coord[j][0]))+', 0, '+str(c[i]*y_l)+'};\n')
                    
#        #Write upper surface coordinates
#        for j in range(foil_num):
#            if (j == 0 or j == foil_num - 1):
#                y_u = camc[i]/camc_ref*camber_ref[j] + \
#                    tc[i]/tc_ref*half_thick_ref[j]
#                f.write('Point('+str(j+1+2*i*len(upper_coord))+') = {'+str(
#                    c[i]*float(upper_coord[j][0]))+', 0, '+str(c[i]*y_u+delta)+'};\n')
#            else:
#                y_u = camc[i]/camc_ref*camber_ref[j] + \
#                    tc[i]/tc_ref*half_thick_ref[j]
#                f.write('Point('+str(j+1+2*i*len(upper_coord)) +
#                        ') = {'+str(c[i]*float(upper_coord[j][0]))+', 0, '+str(c[i]*y_u)+'};\n')
#
#        #Write lower surface coordinates
#        offset = len(upper_coord)
#        for j in range(foil_num):
#            y_l = camc[i]/camc_ref*camber_ref[j] - \
#                tc[i]/tc_ref*half_thick_ref[j]
#            f.write('Point('+str(j+1+offset+2*i*len(lower_coord)) +
#                    ') = {'+str(c[i]*float(lower_coord[j][0]))+', 0, '+str(c[i]*y_l)+'};\n')

        #Create upper surface spline
        f.write('Spline('+str(2*i+1) +
                ') = {'+str(1+2*i*len(upper_coord))+':'+str((2*i+1)*len(upper_coord))+'};\n')

        #Create lower surface spline
        f.write('Spline('+str(2*i+2)+') = {'+str(len(upper_coord)+1+2*i*len(
            upper_coord))+':'+str(2*(i+1)*len(upper_coord))+'};\n\n')

#Translation of the airfoil sections
    for i in range(sec_num):
        f.write('Translate {'+str(xle[i])+', '+str(yle[i])+', '+str(
            zle[i])+'} {Duplicata{Line{'+str(2*i+1)+', '+str(2*i+2)+'};}}\n\n')

#Rotation of the airfoil sections
    for i in range(sec_num):
        f.write('Rotate {{0, 1, 0}, {'+str(xle[i])+', '+str(yle[i])+', '+str(zle[i])+'}, '+str(
            math.radians(theta[i]))+'} {Line{'+str(2*i+1+2*sec_num)+', '+str(2*i+2+2*sec_num)+'};}\n\n')

#Create trailing edge lines (upper surface)
    for i in range(1, sec_num):
        f.write('Line('+str(4*sec_num+i)+') = {'+str(2*foil_num*sec_num+(i-1)*2*(
            foil_num+1)+foil_num)+', '+str(2*foil_num*sec_num+i*2*(foil_num+1)+foil_num)+'};\n')
    f.write('\n')

#Create leading edge lines (upper surface)
    for i in range(1, sec_num):
        f.write('Line('+str(5*sec_num+i-1)+') = {'+str(2*foil_num*sec_num+(
            2*(i-1)*(foil_num+1)+1))+', '+str(2*foil_num*sec_num+2*i*(foil_num+1)+1)+'};\n')
    f.write('\n')

#Create trailing edge lines (lower surface)
    for i in range(1, sec_num):
        f.write('Line('+str(4*sec_num+2*(sec_num-1)+i)+') = {'+str(2*foil_num*sec_num+i*2*(
            foil_num+1))+', '+str(2*foil_num*sec_num+(i+1)*2*(foil_num+1))+'};\n')
    f.write('\n')

#Create leading edge lines (lower surface)
    for i in range(1, sec_num):
        f.write('Line('+str(4*sec_num+3*(sec_num-1)+i)+') = {'+str(2*foil_num*sec_num+(i-1)*2*(
            foil_num+1)+foil_num+3)+', '+str(2*foil_num*sec_num+i*2*(foil_num+1)+foil_num+3)+'};\n')
    f.write('\n')

#Create trailing and leading edge wingtip lines
    f.write('Line('+str(4*sec_num+3*(sec_num-1)+i+1)+') = {'+str(2*foil_num*sec_num+i*2*(
        foil_num+1)+foil_num)+', '+str(2*foil_num*sec_num+(i+1)*2*(foil_num+1))+'};\n')
    f.write('Line('+str(4*sec_num+3*(sec_num-1)+i+2)+') = {'+str(2*foil_num*sec_num+i*2*(
        foil_num+1)+foil_num+3)+', '+str(2*foil_num*sec_num+2*i*(foil_num+1)+1)+'};\n')

#Transform trailing edges into transfinite lines (upper surface)
    for i in range(1, sec_num):
        f.write('Transfinite Line {'+str(4*sec_num+i) +
                '} = ('+str(cols+1)+') Using Progression 1;\n')

#Transform leading edges into transfinite lines (upper surface)
    for i in range(1, sec_num):
        f.write('Transfinite Line {'+str(5*sec_num+i-1) +
                '} = ('+str(cols+1)+') Using Progression 1;\n')
    f.write('\n')

#Transform trailing edges into transfinite lines (upper surface)
    for i in range(1, sec_num):
        f.write('Transfinite Line {'+str(4*sec_num+2*(sec_num-1)+i) +
                '} = ('+str(cols+1)+') Using Progression 1;\n')

#Transform leading edges into transfinite lines (upper surface)
    for i in range(1, sec_num):
        f.write('Transfinite Line {'+str(4*sec_num+3*(sec_num-1)+i) +
                '} = ('+str(cols+1)+') Using Progression 1;\n')
    f.write('\n')

#Transform wingtip lines into transfinite lines
    f.write('Transfinite Line {'+str(4*sec_num+3 *
                                     (sec_num-1)+i+1)+'} = (3) Using Progression 1;\n')
    f.write('Transfinite Line {'+str(4*sec_num+3 *
                                     (sec_num-1)+i+2)+'} = (3) Using Progression 1;\n')

#Transform airfoil splines into transfinite lines
    for i in range(2*sec_num):
        f.write(
            'Transfinite Line {'+str(2*sec_num+i+1)+'} = ('+str(rows+1)+');\n')
    f.write('\n')

#Create upper surface line loops
    for i in range(1, sec_num):
        f.write('Line Loop('+str(i)+') = {'+str(4*sec_num+i)+', '+str(-(
            2*sec_num+2*i+1))+', '+str(-(4*sec_num+i+sec_num-1))+', '+str(2*sec_num+2*i-1)+'};\n')
    f.write('\n')

#Create lower surface line loops
    for i in range(1, sec_num):
        f.write('Line Loop('+str(i+sec_num-1)+') = {'+str(4*sec_num+i+3*(sec_num-1))+', '+str(
            2*sec_num+2*(i+1))+', '+str(-(4*sec_num+i+2*(sec_num-1)))+', '+str(-(2*sec_num+2*i))+'};\n')
    f.write('\n')

#Create wingtip surface line loop
    f.write('Line Loop('+str(i+sec_num)+') = {'+str((4*sec_num+3*(sec_num-1)+i+1))+', '+str(-(
        2*sec_num+2*(i+1)))+', '+str((4*sec_num+3*(sec_num-1)+i+2))+', '+str((2*sec_num+2*i+1))+'};\n')
    f.write('\n')

#Create surfaces
    for i in range(1, 2*sec_num):
        f.write('Ruled Surface('+str(i)+') = {'+str(i)+'};\n')
        f.write('Transfinite Surface {'+str(i)+'};\n')
        f.write('Recombine Surface {'+str(i)+'};\n\n')

#Create 2D mesh
    f.write('Mesh 2;\n\n')

#Create physical surface from all surfaces
    f.write('Physical Surface(1) = {1:'+str(2*(sec_num-1)+1)+'};')

#Create the .msh file with GMSH
p = Popen('gmsh.exe '+user_input+'.geo -2 -o '+user_input+'.msh')
p.wait()

#Create the .wgs file from the newly created .msh file
input_file = user_input + '.msh'

nodes = {}

elm_list = []

elm_dict = {}

with open(input_file) as f:
    lines = f.readlines()
    nodes_beg = lines.index('$Nodes\n')
    nodes_end = lines.index('$EndNodes\n')
    elms_beg = lines.index('$Elements\n')
    elms_end = lines.index('$EndElements\n')
    lines = [i.split() for i in lines]

for i in range(len(lines)):
    if i > (nodes_beg+1) and i < nodes_end:
        nodes[int(lines[i][0])] = [float(lines[i][1]),
                                   float(lines[i][2]), float(lines[i][3])]

    if i > (elms_beg+1) and i < elms_end:
        elm_list.append(int(lines[i][0]))
        elm_dict[int(lines[i][0])] = [int(lines[i][5]), int(
            lines[i][6]), int(lines[i][7]), int(lines[i][8])]

#Redefine number of columns of the whole wing
cols = (sec_num - 1)*cols

#Table containing the elements of the upper network
elm_net_upper = np.zeros((rows, cols), dtype=np.dtype(int))

#Table containing the elements of the lower network
elm_net_lower = np.zeros((rows, cols), dtype=np.dtype(int))

#Table containing the elements of the wingtip surface
elm_net_tip = np.zeros((rows, 2), dtype=np.dtype(int))

for i in range(rows):
    for j in range(cols):
        elm_net_upper[i][j] = elm_list[j*rows+i]
        elm_net_lower[i][j] = elm_list[j*rows+i+rows*cols]

    for j in range(2):
        elm_net_tip[i][j] = elm_list[j*rows+i+2*rows*cols]

#Table containing the nodes of the upper network
nodes_net_upper = np.zeros((rows+1, cols+1))

#Table containing the nodes of the lower network
nodes_net_lower = np.zeros((rows+1, cols+1))

#Table containing the nodes of wing tip
nodes_net_tip = np.zeros((rows+1, 3))

for i in range(rows):
    for j in range(cols):
        elm_upper = elm_net_upper[i][j]
        nodes_net_upper[i][j] = elm_dict[elm_upper][0]
        nodes_net_upper[i][j+1] = elm_dict[elm_upper][1]
        nodes_net_upper[i+1][j+1] = elm_dict[elm_upper][2]
        nodes_net_upper[i+1][j] = elm_dict[elm_upper][3]

        elm_lower = elm_net_lower[i][j]
        nodes_net_lower[i][j] = elm_dict[elm_lower][0]
        nodes_net_lower[i][j+1] = elm_dict[elm_lower][1]
        nodes_net_lower[i+1][j+1] = elm_dict[elm_lower][2]
        nodes_net_lower[i+1][j] = elm_dict[elm_lower][3]

    for j in range(2):
        elm_tip = elm_net_tip[i][j]
        nodes_net_tip[i][j] = elm_dict[elm_tip][0]
        nodes_net_tip[i][j+1] = elm_dict[elm_tip][1]
        nodes_net_tip[i+1][j+1] = elm_dict[elm_tip][2]
        nodes_net_tip[i+1][j] = elm_dict[elm_tip][3]

#Remove gap between upper and lower surfaces at the trailing and leading edges
for j in range(cols+1):
    for i in range(3):
        nodes[nodes_net_upper[0][j]][i] = nodes[nodes_net_lower[rows][j]][i]
        nodes[nodes_net_upper[rows][j]][i] = nodes[nodes_net_lower[0][j]][i]

for j in range(2):
    for i in range(3):
        nodes[nodes_net_tip[0][j]][i] = nodes[nodes_net_tip[0][-1]][i]
        nodes[nodes_net_tip[-1][j]][i] = nodes[nodes_net_tip[-1][-1]][i]

#Write the wgs geometry file
with open(user_input+'.wgs', 'w') as f:
    f.write("'Wing'\n")
    f.write("'upper_surface'\n")
    f.write('1  '+str(cols+1)+' '+str(rows+1)+'  0  0 0 0  0 0 0  1 1 1  1\n')

    for j in range(cols+1):
        for i in range(rows+1):
            f.write(str(nodes[nodes_net_upper[i][j]][0])+'  '+str(
                nodes[nodes_net_upper[i][j]][1])+'  '+str(nodes[nodes_net_upper[i][j]][2])+'\n')

    f.write("'lower_surface'\n")
    f.write('2  '+str(cols+1)+' '+str(rows+1)+'  0  0 0 0  0 0 0  1 1 1  1\n')

    for j in range(cols+1):
        for i in range(rows+1):
            f.write(str(nodes[nodes_net_lower[i][j]][0])+'  '+str(
                nodes[nodes_net_lower[i][j]][1])+'  '+str(nodes[nodes_net_lower[i][j]][2])+'\n')

    f.write("'wingtip_surface'\n")
    f.write('3  '+str(3)+' '+str(rows+1)+'  0  0 0 0  0 0 0  1 1 1  1\n')

    for j in range(3):
        for i in range(rows+1):
            f.write(str(nodes[nodes_net_tip[i][j]][0])+'  '+str(
                nodes[nodes_net_tip[i][j]][1])+'  '+str(nodes[nodes_net_tip[i][j]][2])+'\n')
