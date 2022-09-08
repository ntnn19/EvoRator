from bs4 import BeautifulSoup as bs
import os
import re
import sys
import pandas as pd
#from chart_studio.plotly import plot, iplot
import plotly.graph_objs as go
import matplotlib.pyplot as plt

fp_xyz=open(sys.argv[1],'r')
fp_edgelist=open(sys.argv[2],'r')
df = pd.read_csv(sys.argv[3])
plt.set_cmap('autumn')


Xn,Yn,Zn=[],[],[]
coords=[]
for i in fp_xyz.readlines():
        Xn.append(float(i.split()[1]))
        Yn.append(float(i.split()[2]))
        Zn.append(float(i.split()[3]))
        coords.append([float(i.split()[1]),float(i.split()[2]),float(i.split()[3])])
fp_xyz.close()
edgelist=[]
for i in fp_edgelist.readlines():
        edgelist.append([int(i.split()[0]), int(i.split()[1])])


Xe=[]
Ye=[]
Ze=[]
for i in range(len(edgelist)):
        tempx,tempy,tempz=[Xn[edgelist[i][0]],Xn[edgelist[i][1]],None],[Yn[edgelist[i][0]],Yn[edgelist[i][1]],None],[Zn[edgelist[i][0]],Zn[edgelist[i][1]],None]
        Xe.append(tempx)
        Ye.append(tempy)
        Ze.append(tempz)

minimum = 0.0
maximum = 9.0
n_colours = 9
# Colours are calculated by dividing the RGB colours by 255
colours = [[16,200,209],[140,255,255],[215,255,255],[234,255,255],[255,255,255],
       [252,237,244],[250,201,222],[240,125,171],[160,37,96]]
# colours = [
#     [0.039215686, 0.490196078, 0.509803922],
#     [0.294117647, 0.68627451, 0.745098039],
#     [0.647058824, 0.862745098, 0.901960784],
#     [0.843137255, 0.941176471, 0.941176471],
#     [1, 1, 1],
#     [0.980392157, 0.921568627, 0.960784314],
#     [0.980392157, 0.784313725, 0.862745098],
#     [0.941176471, 0.490196078, 0.666666667],
#     [0.62745098, 0.156862745, 0.37254902]]
bin_size = (maximum - minimum) / n_colours
bin_size = 1

# Loop through colour intervals
score_type = ''
if 'Score' in df.columns.tolist():
    score_type= 'Score'
else:
    score_type= 'dScore'
for i in range(n_colours):
    lower = minimum + (i + 1) * bin_size
    upper = lower + bin_size
    colour = colours[i]
    if (i < n_colours):
        df.loc[(lower<=df[score_type]) & (upper>df[score_type]),'color1']=colour[0]
        df.loc[(lower<=df[score_type]) & (upper>df[score_type]),'color2']=colour[1]
        df.loc[(lower<=df[score_type]) & (upper>df[score_type]),'color3']=colour[2]
    else:
        df.loc[(lower <= df[score_type]) and not (upper < df[score_type]), 'color1'] = colour[0]
        df.loc[(lower <= df[score_type]) and not (upper < df[score_type]), 'color2'] = colour[1]
        df.loc[(lower <= df[score_type]) and not (upper < df[score_type]), 'color3'] = colour[2]
        # sel_string += " & ! b > " + str(upper)

# # Xe = [[-43.864, -42.333, None], [-43.864, -39.831, None], [-43.864, -38.831, None], [-43.864, -38.217, None], [-42.333, -39.831, None], [-42.333, -38.831, None], [-42.333, -38.217, None], [-42.333, -36.06, None], [-39.831, -38.831, None], [-39.831, -38.217, None], [-39.831, -36.06, None], [-39.831, -33.929, None], [-38.831, -38.217, None], [-38.831, -36.06, None], [-38.831, -33.929, None], [-38.831, -33.811, None], [-38.217, -36.06, None], [-38.217, -33.929, None], [-38.217, -33.811, None], [-38.217, -32.47, None], [-36.06, -33.929, None], [-36.06, -33.811, None], [-36.06, -32.47, None], [-36.06, -30.089, None], [-33.929, -33.811, None], [-33.929, -32.47, None], [-33.929, -30.089, None], [-33.929, -28.705, None], [-33.811, -32.47, None], [-33.811, -30.089, None], [-33.811, -28.705, None], [-33.811, -28.544, None], [-32.47, -30.089, None], [-32.47, -28.705, None], [-32.47, -28.544, None], [-32.47, -26.534, None], [-30.089, -28.705, None], [-30.089, -28.544, None], [-30.089, -26.534, None], [-30.089, -24.035, None], [-28.705, -28.544, None], [-28.705, -26.534, None], [-28.705, -24.035, None], [-28.705, -23.458, None], [-28.544, -26.534, None], [-28.544, -24.035, None], [-28.544, -23.458, None], [-28.544, -23.043, None], [-26.534, -24.035, None], [-26.534, -23.458, None], [-26.534, -23.043, None], [-26.534, -20.567, None], [-24.035, -23.458, None], [-24.035, -23.043, None], [-24.035, -20.567, None], [-24.035, -18.822, None], [-23.458, -23.043, None], [-23.458, -20.567, None], [-23.458, -18.822, None], [-23.458, -18.633, None], [-23.043, -20.567, None], [-23.043, -18.822, None], [-23.043, -18.633, None], [-23.043, -17.285, None], [-20.567, -18.822, None], [-20.567, -18.633, None], [-20.567, -17.285, None], [-20.567, -14.432, None], [-18.822, -18.633, None], [-18.822, -17.285, None], [-18.822, -14.432, None], [-18.822, -13.388, None], [-18.633, -17.285, None], [-18.633, -14.432, None], [-18.633, -13.388, None], [-18.633, -12.971, None], [-17.285, -14.432, None], [-17.285, -13.388, None], [-17.285, -12.971, None], [-17.285, -10.99, None], [-14.432, -13.388, None], [-14.432, -12.971, None], [-14.432, -10.99, None], [-14.432, -8.766, None], [-13.388, -12.971, None], [-13.388, -10.99, None], [-13.388, -8.766, None], [-13.388, -8.199, None], [-12.971, -10.99, None], [-12.971, -8.766, None], [-12.971, -8.199, None], [-12.971, -7.423, None], [-10.99, -8.766, None], [-10.99, -8.199, None], [-10.99, -7.423, None], [-10.99, -4.811, None], [-8.766, -8.199, None], [-8.766, -7.423, None], [-8.766, -4.811, None], [-8.766, -3.251, None], [-8.199, -7.423, None], [-8.199, -4.811, None], [-8.199, -3.251, None], [-8.199, -3.082, None], [-7.423, -4.811, None], [-7.423, -3.251, None], [-7.423, -3.082, None], [-7.423, -1.389, None], [-4.811, -3.251, None], [-4.811, -3.082, None], [-4.811, -1.389, None], [-4.811, 1.052, None], [-3.251, -3.082, None], [-3.251, -1.389, None], [-3.251, 1.052, None], [-3.251, 1.956, None], [-3.082, -1.389, None], [-3.082, 1.052, None], [-3.082, 1.956, None], [-3.082, 2.563, None], [-1.389, 1.052, None], [-1.389, 1.956, None], [-1.389, 2.563, None], [-1.389, 4.973, None], [1.052, 1.956, None], [1.052, 2.563, None], [1.052, 4.973, None], [1.052, 6.543, None], [1.956, 2.563, None], [1.956, 4.973, None], [1.956, 6.543, None], [1.956, 6.827, None], [2.563, 4.973, None], [2.563, 6.543, None], [2.563, 6.827, None], [2.563, 8.279, None], [4.973, 6.543, None], [4.973, 6.827, None], [4.973, 8.279, None], [4.973, 10.874, None], [6.543, 6.827, None], [6.543, 8.279, None], [6.543, 10.874, None], [6.543, 11.845, None], [6.827, 8.279, None], [6.827, 10.874, None], [6.827, 11.845, None], [6.827, 12.401, None], [8.279, 10.874, None], [8.279, 11.845, None], [8.279, 12.401, None], [8.279, 14.236, None], [10.874, 11.845, None], [10.874, 12.401, None], [10.874, 14.236, None], [10.874, 16.29, None], [11.845, 12.401, None], [11.845, 14.236, None], [11.845, 16.29, None], [11.845, 17.723, None], [12.401, 14.236, None], [12.401, 16.29, None], [12.401, 17.723, None], [12.401, 18.946, None], [14.236, 16.29, None], [14.236, 17.723, None], [14.236, 18.946, None], [16.29, 17.723, None], [16.29, 18.946, None], [16.29, 21.085, None], [17.723, 18.946, None], [17.723, 21.085, None], [18.946, 21.085, None]]
# Xe = Xe
#
#
# # Ye = [[27.976, 31.377, None], [27.976, 30.223, None], [27.976, 27.602, None], [27.976, 30.401, None], [31.377, 30.223, None], [31.377, 27.602, None], [31.377, 30.401, None], [31.377, 32.552, None], [30.223, 27.602, None], [30.223, 30.401, None], [30.223, 32.552, None], [30.223, 29.48, None], [27.602, 30.401, None], [27.602, 32.552, None], [27.602, 29.48, None], [27.602, 29.168, None], [30.401, 32.552, None], [30.401, 29.48, None], [30.401, 29.168, None], [30.401, 32.672, None], [32.552, 29.48, None], [32.552, 29.168, None], [32.552, 32.672, None], [32.552, 32.229, None], [29.48, 29.168, None], [29.48, 32.672, None], [29.48, 32.229, None], [29.48, 29.052, None], [29.168, 32.672, None], [29.168, 32.229, None], [29.168, 29.052, None], [29.168, 30.906, None], [32.672, 32.229, None], [32.672, 29.052, None], [32.672, 30.906, None], [32.672, 33.733, None], [32.229, 29.052, None], [32.229, 30.906, None], [32.229, 33.733, None], [32.229, 31.42, None], [29.052, 30.906, None], [29.052, 33.733, None], [29.052, 31.42, None], [29.052, 29.396, None], [30.906, 33.733, None], [30.906, 31.42, None], [30.906, 29.396, None], [30.906, 32.624, None], [33.733, 31.42, None], [33.733, 29.396, None], [33.733, 32.624, None], [33.733, 33.751, None], [31.42, 29.396, None], [31.42, 32.624, None], [31.42, 33.751, None], [31.42, 30.373, None], [29.396, 32.624, None], [29.396, 33.751, None], [29.396, 30.373, None], [29.396, 30.66, None], [32.624, 33.751, None], [32.624, 30.373, None], [32.624, 30.66, None], [32.624, 34.2, None], [33.751, 30.373, None], [33.751, 30.66, None], [33.751, 34.2, None], [33.751, 33.803, None], [30.373, 30.66, None], [30.373, 34.2, None], [30.373, 33.803, None], [30.373, 30.772, None], [30.66, 34.2, None], [30.66, 33.803, None], [30.66, 30.772, None], [30.66, 33.19, None], [34.2, 33.803, None], [34.2, 30.772, None], [34.2, 33.19, None], [34.2, 35.718, None], [33.803, 30.772, None], [33.803, 33.19, None], [33.803, 35.718, None], [33.803, 32.852, None], [30.772, 33.19, None], [30.772, 35.718, None], [30.772, 32.852, None], [30.772, 31.546, None], [33.19, 35.718, None], [33.19, 32.852, None], [33.19, 31.546, None], [33.19, 34.975, None], [35.718, 32.852, None], [35.718, 31.546, None], [35.718, 34.975, None], [35.718, 35.883, None], [32.852, 31.546, None], [32.852, 34.975, None], [32.852, 35.883, None], [32.852, 32.463, None], [31.546, 34.975, None], [31.546, 35.883, None], [31.546, 32.463, None], [31.546, 33.031, None], [34.975, 35.883, None], [34.975, 32.463, None], [34.975, 33.031, None], [34.975, 36.51, None], [35.883, 32.463, None], [35.883, 33.031, None], [35.883, 36.51, None], [35.883, 35.187, None], [32.463, 33.031, None], [32.463, 36.51, None], [32.463, 35.187, None], [32.463, 32.159, None], [33.031, 36.51, None], [33.031, 35.187, None], [33.031, 32.159, None], [33.031, 34.632, None], [36.51, 35.187, None], [36.51, 32.159, None], [36.51, 34.632, None], [36.51, 36.593, None], [35.187, 32.159, None], [35.187, 34.632, None], [35.187, 36.593, None], [35.187, 33.428, None], [32.159, 34.632, None], [32.159, 36.593, None], [32.159, 33.428, None], [32.159, 31.995, None], [34.632, 36.593, None], [34.632, 33.428, None], [34.632, 31.995, None], [34.632, 35.337, None], [36.593, 33.428, None], [36.593, 31.995, None], [36.593, 35.337, None], [36.593, 35.163, None], [33.428, 31.995, None], [33.428, 35.337, None], [33.428, 35.163, None], [33.428, 31.657, None], [31.995, 35.337, None], [31.995, 35.163, None], [31.995, 31.657, None], [31.995, 32.878, None], [35.337, 35.163, None], [35.337, 31.657, None], [35.337, 32.878, None], [35.337, 36.124, None], [35.163, 31.657, None], [35.163, 32.878, None], [35.163, 36.124, None], [35.163, 33.648, None], [31.657, 32.878, None], [31.657, 36.124, None], [31.657, 33.648, None], [31.657, 31.251, None], [32.878, 36.124, None], [32.878, 33.648, None], [32.878, 31.251, None], [32.878, 34.657, None], [36.124, 33.648, None], [36.124, 31.251, None], [36.124, 34.657, None], [33.648, 31.251, None], [33.648, 34.657, None], [33.648, 35.55, None], [31.251, 34.657, None], [31.251, 35.55, None], [34.657, 35.55, None]]
# Ye = Ye
# # Ze = [[17.706, 18.594, None], [17.706, 21.098, None], [17.706, 18.654, None], [17.706, 16.108, None], [18.594, 21.098, None], [18.594, 18.654, None], [18.594, 16.108, None], [18.594, 18.336, None], [21.098, 18.654, None], [21.098, 16.108, None], [21.098, 18.336, None], [21.098, 18.926, None], [18.654, 16.108, None], [18.654, 18.336, None], [18.654, 18.926, None], [18.654, 15.111, None], [16.108, 18.336, None], [16.108, 18.926, None], [16.108, 15.111, None], [16.108, 14.576, None], [18.336, 18.926, None], [18.336, 15.111, None], [18.336, 14.576, None], [18.336, 17.624, None], [18.926, 15.111, None], [18.926, 14.576, None], [18.926, 17.624, None], [18.926, 15.98, None], [15.111, 14.576, None], [15.111, 17.624, None], [15.111, 15.98, None], [15.111, 12.625, None], [14.576, 17.624, None], [14.576, 15.98, None], [14.576, 12.625, None], [14.576, 14.205, None], [17.624, 15.98, None], [17.624, 12.625, None], [17.624, 14.205, None], [17.624, 15.956, None], [15.98, 12.625, None], [15.98, 14.205, None], [15.98, 15.956, None], [15.98, 12.75, None], [12.625, 14.205, None], [12.625, 15.956, None], [12.625, 12.75, None], [12.625, 10.871, None], [14.205, 15.956, None], [14.205, 12.75, None], [14.205, 10.871, None], [14.205, 13.532, None], [15.956, 12.75, None], [15.956, 10.871, None], [15.956, 13.532, None], [15.956, 13.341, None], [12.75, 10.871, None], [12.75, 13.532, None], [12.75, 13.341, None], [12.75, 9.511, None], [10.871, 13.532, None], [10.871, 13.341, None], [10.871, 9.511, None], [10.871, 9.81, None], [13.532, 13.341, None], [13.532, 9.511, None], [13.532, 9.81, None], [13.532, 12.38, None], [13.341, 9.511, None], [13.341, 9.81, None], [13.341, 12.38, None], [13.341, 10.226, None], [9.511, 9.81, None], [9.511, 12.38, None], [9.511, 10.226, None], [9.511, 7.305, None], [9.81, 12.38, None], [9.81, 10.226, None], [9.81, 7.305, None], [9.81, 9.284, None], [12.38, 10.226, None], [12.38, 7.305, None], [12.38, 9.284, None], [12.38, 10.506, None], [10.226, 7.305, None], [10.226, 9.284, None], [10.226, 10.506, None], [10.226, 7.014, None], [7.305, 9.284, None], [7.305, 10.506, None], [7.305, 7.014, None], [7.305, 5.681, None], [9.284, 10.506, None], [9.284, 7.014, None], [9.284, 5.681, None], [9.284, 8.342, None], [10.506, 7.014, None], [10.506, 5.681, None], [10.506, 8.342, None], [10.506, 7.835, None], [7.014, 5.681, None], [7.014, 8.342, None], [7.014, 7.835, None], [7.014, 4.074, None], [5.681, 8.342, None], [5.681, 7.835, None], [5.681, 4.074, None], [5.681, 4.524, None], [8.342, 7.835, None], [8.342, 4.074, None], [8.342, 4.524, None], [8.342, 7.084, None], [7.835, 4.074, None], [7.835, 4.524, None], [7.835, 7.084, None], [7.835, 4.92, None], [4.074, 4.524, None], [4.074, 7.084, None], [4.074, 4.92, None], [4.074, 2.039, None], [4.524, 7.084, None], [4.524, 4.92, None], [4.524, 2.039, None], [4.524, 4.06, None], [7.084, 4.92, None], [7.084, 2.039, None], [7.084, 4.06, None], [7.084, 5.299, None], [4.92, 2.039, None], [4.92, 4.06, None], [4.92, 5.299, None], [4.92, 1.749, None], [2.039, 4.06, None], [2.039, 5.299, None], [2.039, 1.749, None], [2.039, 0.526, None], [4.06, 5.299, None], [4.06, 1.749, None], [4.06, 0.526, None], [4.06, 3.371, None], [5.299, 1.749, None], [5.299, 0.526, None], [5.299, 3.371, None], [5.299, 2.191, None], [1.749, 0.526, None], [1.749, 3.371, None], [1.749, 2.191, None], [1.749, -1.324, None], [0.526, 3.371, None], [0.526, 2.191, None], [0.526, -1.324, None], [0.526, -0.591, None], [3.371, 2.191, None], [3.371, -1.324, None], [3.371, -0.591, None], [3.371, 1.492, None], [2.191, -1.324, None], [2.191, -0.591, None], [2.191, 1.492, None], [2.191, -1.176, None], [-1.324, -0.591, None], [-1.324, 1.492, None], [-1.324, -1.176, None], [-1.324, -2.298, None], [-0.591, 1.492, None], [-0.591, -1.176, None], [-0.591, -2.298, None], [1.492, -1.176, None], [1.492, -2.298, None], [1.492, 0.865, None], [-1.176, -2.298, None], [-1.176, 0.865, None], [-2.298, 0.865, None]]
# Ze = Ze
print(len(Xe),len(Ze),len(Ye))
# Xn = x
# Xn = [-43.864, -42.333, -39.831, -38.831, -38.217, -36.06, -33.929, -33.811, -32.47, -30.089, -28.705, -28.544, -26.534, -24.035, -23.458, -23.043, -20.567, -18.822, -18.633, -17.285, -14.432, -13.388, -12.971, -10.99, -8.766, -8.199, -7.423, -4.811, -3.251, -3.082, -1.389, 1.052, 1.956, 2.563, 4.973, 6.543, 6.827, 8.279, 10.874, 11.845, 12.401, 14.236, 16.29, 17.723, 18.946, 21.085]
# Yn = y
# Yn = [27.976, 31.377, 30.223, 27.602, 30.401, 32.552, 29.48, 29.168, 32.672, 32.229, 29.052, 30.906, 33.733, 31.42, 29.396, 32.624, 33.751, 30.373, 30.66, 34.2, 33.803, 30.772, 33.19, 35.718, 32.852, 31.546, 34.975, 35.883, 32.463, 33.031, 36.51, 35.187, 32.159, 34.632, 36.593, 33.428, 31.995, 35.337, 35.163, 31.657, 32.878, 36.124, 33.648, 31.251, 34.657, 35.55]
# Zn = z
# Zn = [17.706, 18.594, 21.098, 18.654, 16.108, 18.336, 18.926, 15.111, 14.576, 17.624, 15.98, 12.625, 14.205, 15.956, 12.75, 10.871, 13.532, 13.341, 9.511, 9.81, 12.38, 10.226, 7.305, 9.284, 10.506, 7.014, 5.681, 8.342, 7.835, 4.074, 4.524, 7.084, 4.92, 2.039, 4.06, 5.299, 1.749, 0.526, 3.371, 2.191, -1.324, -0.591, 1.492, -1.176, -2.298, 0.865]
cl=df[['color1','color2','color3']].values.ravel().tolist()
cl=[cl[i:i+3] for i in range(0,len(cl),3)]
print(cl)
# exit()
# cl=[f'({cl[i:i+3][0]},{cl[i:i+3][1]},{cl[i:i+3][2]})' for i in range(len(cl))]

pos_names = df['pdb_position'].str.split("_").str[0].replace(r"[A-Z]+", '', regex=True)
pos_names = df['pdb_aa'] + pos_names.tolist()
print(pos_names)
trace1=go.Scatter3d(x=[item for sublist in Xe for item in sublist],
               y=[item for sublist in Ye for item in sublist],
               z=[item for sublist in Ze for item in sublist],
               mode='lines',
               line=dict(color='black', width=3),
               hoverinfo='none'
               )

trace2=go.Scatter3d(x=Xn,
               y=Yn,
               z=Zn,
               mode='markers',
               name='actors',
               marker=dict(symbol='circle',
                             size=10,
                             color=cl,
                             line=dict(color='black', width=0.5)
                             ),
               text=pos_names,
               hoverinfo='text'
               )

trace3=go.Scatter3d(x=Xn,
               y=Yn,
               z=Zn,
                mode='lines',
                line=dict(color='black', width=10),
                hoverinfo='none'

               )

axis=dict(showbackground=False,
          showline=False,
          zeroline=False,
          showgrid=False,
          showticklabels=False,
          showspikes= False,
          title=''
          )

layout = go.Layout(
         title="",
         width=1000,
         height=1000,
         showlegend=False,
         scene=dict(
             xaxis=dict(axis),
             yaxis=dict(axis),
             zaxis=dict(axis),
        ),
     margin=dict(
        t=100
    ),
    hovermode='closest',
    annotations=[
           dict(
           showarrow=False,

            text="" , #"Data source: <a href='http://bost.ocks.org/mike/miserables/miserables.json'>[1] miserables.json</a>"
            xref='paper',
            yref='paper',
            x=0,
            y=0.1,
            xanchor='left',
            yanchor='bottom',
            font=dict(
            size=14
            )
            )
        ],    )

data=[trace1, trace2,trace3]
fig=go.Figure(data=data, layout=layout)
fig.update_layout(xaxis_showgrid=False, yaxis_showgrid=False)
fig.write_html(sys.argv[4]+"/network.html",include_plotlyjs="cdn")

# Open the HTML in which you want to make changes
html = open(sys.argv[4]+'/network.html')

# Parse HTML file in Beautiful Soup
soup = bs(html, 'html.parser')

old_text = soup.find("div", {"class": "plotly-graph-div"})
print(old_text)
# Replace the already stored text with
# the new text which you wish to assign
old_text['style'] =   "margin: auto; background: white;  max-width: 1000px;"
new_text = soup.find("div", {"class": "plotly-graph-div"})

print(new_text)
# margin:auto;background:white;
with open(sys.argv[4]+'/network.html', "wb") as f_output:
    f_output.write(soup.prettify("utf-8"))



