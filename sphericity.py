# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 11:28:11 2021

@author: Jouni Kuusisto

* insert filename in if __name__ == '__main__': section at the bottom

*numpy-stl library displays error message of a not closed mesh when calling get_mass_properties() even 
 when a closed mesh (from FreeCAD v0.18 as .SLT file) is loaded. Volume and cog agree with FreeCAD and manual calculation.
"""

import numpy as np
import math
from stl import mesh
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

#Loading mesh and add vectors to the plot
def loadMesh(file: str):
    return mesh.Mesh.from_file(file)

def meshPlot(case_mesh: mesh.Mesh):
    figure = plt.figure()
    axes = mplot3d.Axes3D(figure)
    axes.add_collection3d(mplot3d.art3d.Line3DCollection(case_mesh.vectors))
    #Auto scale to mesh size
    scale = case_mesh.points.flatten()
    axes.auto_scale_xyz(scale, scale, scale)
    #Show the plot
    plt.show()


'''
Following function from: https://newbedev.com/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to-x-and-y
ax.set_box_aspect([1,1,1]) is also needed before this function is called
This sets axis equal in matplotlibs 3d plots
'''
def set_axes_equal(ax): 
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

#Strictly own code from here on

#Calculating center point of triangles. case_mesh.vectors returns a numpy array, where lines represent each triangle and columns coordinates
def triangleCenters(case_mesh: mesh.Mesh):
    tri_center_coordinates = np.zeros((len(case_mesh), 3))
    for count, tri in enumerate(case_mesh.vectors):
        tri_center_coordinates[count,:] = np.mean(tri, axis = 0)
    return tri_center_coordinates

#Calculating triangle areas by cross product of vectors
def triangleAreas(case_mesh: mesh.Mesh):
    areasArray = np.zeros((len(case_mesh), 1))
    for count, tri in enumerate(case_mesh.vectors):
        AB = tri[1,:] - tri[0,:]
        AC = tri[2,:] - tri[0,:]
        areasArray[count] = 0.5 * np.linalg.norm(np.cross(AB, AC))
    return areasArray
              
#Calculate triangle centers by weighing them by their area
def weightedCenters(case_mesh: mesh.Mesh):
    centers = triangleCenters(case_mesh)
    areas = triangleAreas(case_mesh)
    weightedCenters = centers.copy()
    i = 0
    for center, area in zip(centers, areas):
        weightedCenters[i,:] = area * center
        i += 1
    return weightedCenters

#Calculate mass center, this is done also with stl package method - just to be sure results agree (they did up to decimal level). This is done in accordance to paper by Bisbal
def massCenter(case_mesh: mesh.Mesh):
    return np.sum(weightedCenters(case_mesh), axis = 0) / np.sum(triangleAreas(case_mesh))


#Calculate radius of best fitted sphere. First calculate sum of all radii from center of mass to all triangles weighted by respective triangle areas. Then division by all surface area and number of triangles results AR.
def averageWeightedRadius(case_mesh: mesh.Mesh):
    LA_surface_area = np.sum(triangleAreas(case_mesh))
    LA_mass_center = massCenter(case_mesh)
    return np.sum(triangleAreas(case_mesh) * (np.linalg.norm((triangleCenters(case_mesh) - LA_mass_center), axis = 1))) / (LA_surface_area * len(case_mesh))

#Calculate standard deviation of the AR (Bisbal)
def sdAr(case_mesh: mesh.Mesh): 
    AR = averageWeightedRadius(case_mesh)
    LA_surface_area = np.sum(triangleAreas(case_mesh))
    LA_mass_center = LA_mass_center = massCenter(case_mesh)
    return ((np.sum(triangleAreas(case_mesh) * ( np.linalg.norm((triangleCenters(case_mesh) - LA_mass_center), axis = 1) - AR)**2)) / (LA_surface_area * len(case_mesh)))**0.5

#CVS stands for "coefficient of variation of the sphere" (Bisbal)
def cvs(case_mesh: mesh.Mesh): #CVS stands for "coefficient of variation of the sphere" (Bisbal)
    return sdAr(case_mesh)/averageWeightedRadius(case_mesh)

#Calculate LA sphericity as proposed by Bisbal
def laspBisbal(case_mesh: mesh.Mesh, percentage = False): #Method proposed by Bisbal
    lasp = 1-cvs(case_mesh)
    if percentage:
        return lasp*100
    else:
        return lasp
    
#Calculate own 3D sphericity
def threeDsKuusisto(case_mesh: mesh.Mesh, percentage = False): #Our own method
    volume, _, _ = case_mesh.get_mass_properties()
    threeDs = (math.pi**(1/3) * ((6*volume)**(2/3))) / np.sum(triangleAreas(case_mesh))
    if percentage:
        return threeDs*100
    else:
        return threeDs
    

def plotTriangleCentersAndCenterOfMass(case_mesh: mesh.Mesh, title = "", save = False, save_name=""):
    #LA_surface_area = np.sum(triangleAreas(case_mesh))
    LA_mass_center = massCenter(case_mesh)  
    tri_center_coordinates = triangleCenters(case_mesh)
    AR = averageWeightedRadius(case_mesh)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.view_init(90,0)
    #ax.scatter(tri_center_coordinates[:,0], tri_center_coordinates[:,1], tri_center_coordinates[:,2], c='r', marker="x", s=2)
    ax.scatter(LA_mass_center[0], LA_mass_center[1], LA_mass_center[2], c='r')
    
    ax.add_collection3d(mplot3d.art3d.Line3DCollection(case_mesh.vectors, linewidths=0.5))
        
    # draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    x *= AR
    y *= AR
    z *= AR
    x += LA_mass_center[0]
    y += LA_mass_center[1]
    z += LA_mass_center[2]
    ax.plot_wireframe(x, y, z, color="r", alpha=0.5, linewidths = 1)
    
    plt.title(title)
    ax.set_box_aspect([1,1,1])
    ax.grid(False)
    set_axes_equal(ax)
    
    if save:
        fig.savefig(save_name,dpi=300)    
    plt.show()
    
def printResults(case_mesh: mesh.Mesh):
    bisbal = laspBisbal(case_mesh)
    kuusisto = threeDsKuusisto(case_mesh)
    volume, _, _ = case_mesh.get_mass_properties()
    print("Volume:           " + str(volume/1000))
    print("LASP (Bisbal):    " + str(bisbal))
    print("3DS (Kuusisto):   " + str(kuusisto))
    return volume, bisbal, kuusisto


if __name__ == '__main__':
    
    case_mesh = loadMesh("dotu_part2.stl") # insert filename here
    meshPlot(case_mesh)
    plotTriangleCentersAndCenterOfMass(case_mesh)
    
    #Printing mesh characteristics
    volume, cog, inertia = case_mesh.get_mass_properties()
    print("Mesh Volume                             = {0}".format(volume))
    print("Position of the center of gravity (COG) = {0}".format(cog))
    
    printResults(case_mesh)