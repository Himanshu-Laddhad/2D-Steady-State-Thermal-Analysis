import matplotlib as mt
import matplotlib.pyplot as plt
import numpy as np

#Base parameter required for the simulation
l = float(input("length: "))
w = float(input("width: "))
rx = float(input("successive ratio in x direction: "))
ry = float(input("successive ratio in y direction: "))
mode = input("Enter N for specifying number of cells and L for specifying length of first cell: ")

#For specified number of cells
if mode == "N":

    #Number of cells 
    cx = int(input("cells in x direc:"))
    cy = int(input("cells in y direc:"))

    #calculating the length of first cell using GP
    if rx == 1:
        fclx = l/cx
    else:
        fclx = l * abs((1 - rx)/(1 - rx**cx))
    
    if ry == 1:
        fcly = w/cy
    else:
        fcly = w * abs((1 - ry)/(1 - ry**cy))
    
    #Determining x coordinates and y coordiantes of the TempM including boundary
    x_coo = np.zeros(1)
    t = 0
    for i in range(cx):
        t += fclx*(rx**i)
        x_coo = np.append(x_coo, t)
    x_cells = cx

    y_coo = np.zeros(1)
    t = 0
    for i in range(cy):
        t += fcly*(ry**i)
        y_coo = np.append(y_coo, t)
    y_cells = cy

    #plotting the TempM using x and y coordinates
    for i in x_coo:
        plt.vlines(i, ymin = 0, ymax = w)
    for j in y_coo:
        plt.hlines(j, xmin = 0, xmax = l)
    plt.show()

#For given length  of the first cell    
elif mode == "L":
    #Length of first cells
    fclx = float(input("length of first element in x direc: "))
    fcly = float(input("length of first element in y direc: "))

    #Determing x and y coordinates of the cell corners
    x_coo = np.zeros(1)
    t = fclx
    x_cells = 0
    while t <= l:
        x_coo = np.append(x_coo, t)
        t += fclx*(rx**x_cells)
        x_cells += 1
    x_cells -= 1

    y_coo = np.zeros(1)
    t = fcly
    y_cells = 0
    while t <= w:
        y_coo = np.append(y_coo, t)
        t += fcly*(ry**y_cells)
        y_cells += 1
    y_cells -= 1

#plotting the TempM 
    for i in x_coo:
        plt.vlines(i, ymin = 0, ymax = w)
    for j in y_coo:
        plt.hlines(j, xmin = 0, xmax = l)
    #top and rightmost boundary
    plt.hlines(w, xmin = 0, xmax = l) 
    plt.vlines(l, ymin = 0, ymax = w)
    plt.show()
else:
    quit("Invalid Input")

#Matrix of the mesh
Mesh = np.empty(shape = (y_cells + 1, x_cells + 1))

#boundary conditions
lb_t = float(input("left boundary temperature:"))
rb_t = float(input("right boundary temperature:"))
tb_t = float(input("top boundary temperature:"))
bb_t = float(input("bottom boundary temperature:"))
#real boundary 
rows = [i for i in range(0,y_cells)]
columns = [i for i in range(0,x_cells)]

#Initialising a cell center temperaure matrix
TempM = np.empty(shape = (y_cells, x_cells))

#Initialising the value of cell centers in the matrix
for i in rows:
    for j in columns:
        if (i - j) > 0:
            if (i - j) < (w - l):
                TempM[i,j] = lb_t
            elif (i - j) > (w - l):
                TempM[i,j] = tb_t
            else:
                TempM[i,j] = (lb_t + tb_t)/2
        elif (i - j) < 0:
            if (i - j) < (w - l):
                TempM[i,j] = bb_t
            elif (i - j) > (w - l):
                TempM[i,j] = rb_t
            else:
                TempM[i,j] = (bb_t + rb_t)/2
        else:
            if (i - j) < (w - l):
                TempM[i,j] = (bb_t + lb_t)/2
            elif (i - j) > (w - l):
                TempM[i,j] = (rb_t + tb_t)/2
            else:
                TempM[i,j] = (bb_t + rb_t + tb_t + lb_t)/4

wt = np.zeros((y_cells, x_cells, 5))

for i in rows:
    for j in columns:
        if (i >= 1 and i <= y_cells - 2) and (j >= 1 and j <= x_cells - 2):
#matrix ecluding boundary cells
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
#boundary cells
        elif j == x_cells - 1 and i in range(1, y_cells - 2):
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 1] - x_coo[j])
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif i == 0 and j in range(1, x_cells - 2):
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i])
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif i == y_cells-1 and j in range(1, x_cells - 2):
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 1] - y_coo[i])
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif j == 0 and i in range(1, y_cells - 1):
            dxw = (x_coo[j + 1] - x_coo[j])
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
#corners
        elif (i,j) == (0,0):
            dxw = (x_coo[j + 1] - x_coo[j])
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i])
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif (i,j) == (y_cells - 1, x_cells - 1):
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 1] - x_coo[j])
            dyn = (y_coo[i + 1] - y_coo[i])
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif (i,j) == (0, x_cells-1):
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 1] - x_coo[j])
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i])
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif (i,j) == (y_cells-1, 0):
            dxw = (x_coo[j + 1] - x_coo[j])
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 1] - y_coo[i])
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        #wt are saved in the format North-South-East-West-Summation
        wt[i,j,0] = a_n
        wt[i,j,1] = a_s
        wt[i,j,2] = a_e
        wt[i,j,3] = a_w
        wt[i,j,4] = a_p

iterations = int(input("iterations: "))
for k in range(iterations):
    for i in rows:
        for j in columns:
            if (i >= 1 and i <= y_cells - 2) and (j >= 1 and j <= x_cells - 2):
    #matrix ecluding boundary cells
                TempM[i,j] = (wt[i,j,0]*TempM[i+1,j] + wt[i,j,1]*TempM[i-1,j] + wt[i,j,2]*TempM[i,j+1] + wt[i,j,3]*TempM[i,j-1])*(-1)/wt[i,j,4]
    #boundary cells
            elif j == x_cells and i in range(1, y_cells - 1):#east
                TempM[i,j] = (wt[i,j,0]*TempM[i+1,j] + wt[i,j,1]*TempM[i-1,j] + wt[i,j,2]*rb_t*2 + wt[i,j,3]*TempM[i,j-1])*(-1)/(wt[i,j,4] - wt[i,j,2])
            elif i == 0 and j in range(1, x_cells - 1):#south
                TempM[i,j] = (wt[i,j,0]*TempM[i+1,j] + wt[i,j,1]*bb_t*2 + wt[i,j,2]*TempM[i,j+1] + wt[i,j,3]*TempM[i,j-1])*(-1)/(wt[i,j,4] - wt[i,j,1])
            elif i == y_cells and j in range(1, x_cells - 1):#north
                TempM[i,j] = (wt[i,j,0]*tb_t*2 + wt[i,j,1]*TempM[i-1,j] + wt[i,j,2]*TempM[i,j+1] + wt[i,j,3]*TempM[i,j-1])*(-1)/(wt[i,j,4] - wt[i,j,0])
            elif j == 0 and i in range(1, y_cells - 1):#west
                TempM[i,j] = (wt[i,j,0]*TempM[i+1,j] + wt[i,j,1]*TempM[i-1,j] + wt[i,j,2]*TempM[i,j+1] + wt[i,j,3]*lb_t*2)*(-1)/(wt[i,j,4] - wt[i,j,3])
    #corners
            elif (i,j) == (0,0):
                TempM[i,j] = (wt[i,j,0]*TempM[i+1,j] + wt[i,j,1]*bb_t*2 + wt[i,j,2]*TempM[i,j+1] + wt[i,j,3]*lb_t*2)*(-1)/(wt[i,j,4] - wt[i,j,3] - wt[i,j,1])            
            elif (i,j) == (y_cells, x_cells):
                TempM[i,j] = (wt[i,j,0]*tb_t*2 + wt[i,j,1]*TempM[i-1,j] + wt[i,j,2]*rb_t*2 + wt[i,j,3]*TempM[i,j-1])*(-1)/(wt[i,j,4] - wt[i,j,0] - wt[i,j,2])
            elif (i,j) == (0, x_cells):
                TempM[i,j] = (wt[i,j,0]*TempM[i+1,j] + wt[i,j,1]*bb_t*2 + wt[i,j,2]*rb_t*2 + wt[i,j,3]*TempM[i,j-1])*(-1)/(wt[i,j,4] - wt[i,j,2] - wt[i,j,1])
            elif (i,j) == (y_cells, 0):
                TempM[i,j] = (wt[i,j,0]*tb_t*2 + wt[i,j,1]*TempM[i-1,j] + wt[i,j,2]*TempM[i,j+1] + wt[i,j,3]*lb_t*2)*(-1)/(wt[i,j,4] - wt[i,j,0] -  wt[i,j,3])

print(TempM)




            









