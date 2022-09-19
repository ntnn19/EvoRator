import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D
import base64
import re
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from mpl_toolkits.mplot3d.proj3d import proj_transform
from matplotlib.text import Annotation


class Annotation3D(Annotation):
    '''Annotate the point xyz with text s'''

    def __init__(self, s, xyz, *args, **kwargs):
        Annotation.__init__(self,s, xy=(0,0), *args, **kwargs)
        self._verts3d = xyz

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.xy=(xs,ys)
        Annotation.draw(self, renderer)

def annotate3D(ax, s, *args, **kwargs):
    '''add anotation text s to to Axes3d ax'''

    tag = Annotation3D(s, *args, **kwargs)
    ax.add_artist(tag)

'''Program to view the 3d plot of the network
command line arguments:
1) xyz file: file containing the x,y,z coordinates of the residues obtained from pdb file.
	column 1: residueid
	column 2: x coordinate
	column 3: y coordinate
	column 4: z coordinate
2)edgelist file
'''

plt.set_cmap('autumn')

fp_xyz=open(sys.argv[1],'r')
x,y,z=[],[],[]
coords=[]
for i in fp_xyz.readlines():
        x.append(float(i.split()[1]))
        y.append(float(i.split()[2]))
        z.append(float(i.split()[3]))
        coords.append([float(i.split()[1]),float(i.split()[2]),float(i.split()[3])])
fp_xyz.close()
fp_edgelist=open(sys.argv[2],'r')
edgelist=[]
for i in fp_edgelist.readlines():
        edgelist.append([int(i.split()[0]), int(i.split()[1])])

#fig=plt.figure(figsize=(20,20))
fig=plt.figure(figsize=(20,20))
ax=Axes3D(fig)
df = pd.read_csv(sys.argv[3])
minimum = 0.0
maximum = 9.0
n_colours = 9
# Colours are calculated by dividing the RGB colours by 255
# RGB = [[16,200,209],[140,255,255],[215,255,255],[234,255,255],[255,255,255],
#        [252,237,244],[250,201,222],[240,125,171],[160,37,96]]
colours = [
    [0.039215686, 0.490196078, 0.509803922],
    [0.294117647, 0.68627451, 0.745098039],
    [0.647058824, 0.862745098, 0.901960784],
    [0.843137255, 0.941176471, 0.941176471],
    [1, 1, 1],
    [0.980392157, 0.921568627, 0.960784314],
    [0.980392157, 0.784313725, 0.862745098],
    [0.941176471, 0.490196078, 0.666666667],
    [0.62745098, 0.156862745, 0.37254902]]
bin_size = (maximum - minimum) / n_colours
# bin_size = 1

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

ax.scatter(x,y,z,s=40,c=df[['color1','color2','color3']].values,alpha=1, depthshade=False)


for i in range(len(edgelist)):
        tempx,tempy,tempz=[x[edgelist[i][0]],x[edgelist[i][1]]],[y[edgelist[i][0]],y[edgelist[i][1]]],[z[edgelist[i][0]],z[edgelist[i][1]]]

        ax.plot(tempx,tempy,tempz,c='gray',label=edgelist[i][0],lw=0.5)
ax.plot(x,y,z,c='black',lw=2)


for i,xyz_  in enumerate(coords):
    pos_names = df['pdb_position'].str.split("_").str[0].replace(r"[A-Z]+", '', regex=True).tolist()[i]
    annotate3D(ax, s=str((df['pdb_aa'] + pos_names).tolist()[i]), xyz=xyz_, fontsize=6, xytext=(0, 0),
               textcoords='offset points', ha='center', va='center', fontweight="bold")


ax.grid(b=None)
ax.axis('off')


plt.tight_layout()
from io import BytesIO
temp = BytesIO()

#fig.savefig(temp, format="png",dpi=150)
fig.savefig(sys.argv[4]+"/network.png", format="png",dpi=96)
exit()
fig_encode_bs64 = base64.b64encode(temp.getvalue()).decode('utf-8')

html_string = """
<img src = 'data:image/png;base64,{}', width="1600" height = 1200"/>
""".format(fig_encode_bs64)
with open(sys.argv[4]+"/network.html", "w") as f:
    f.write("""<style>
    #content {
        position: relative;
    }
    #content img {
        position: absolute;
        top: 0px;
        left: 0px;
    }
</style>
""")
    f.write(html_string)


