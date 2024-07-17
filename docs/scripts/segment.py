import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Function to plot a portion of a sphere
def plot_sphere(ax,center,radius,phi_range,theta_range):
    phi=np.linspace(phi_range[0],phi_range[1],1000)
    theta=np.linspace(theta_range[0],theta_range[1],1000)
    phi,theta=np.meshgrid(phi,theta)

    x=center[0]+radius*np.sin(phi)*np.cos(theta)
    y=center[1]+radius*np.sin(phi)*np.sin(theta)
    z=center[2]+radius*np.cos(phi)

    ax.plot_surface(x,y,z,color='b',alpha=.5,rstride=1,cstride=1)
    ##ax.plot_wireframe(x+.01,y+.01,z+.01,color='w',lw=2,alpha=.5,rstride=100,cstride=100)

# Function to plot normal vector at a point on the sphere
def plot_normal_vector(ax,point,normal):
    ax.quiver(point[0],point[1],point[2],normal[0],normal[1],normal[2],color='g',length=.2,normalize=True)

def plot_radial_vector(ax,point,normal):
    ax.plot([point[0],normal[0]],[point[1],normal[1]],[point[2],normal[2]],color='k')

# Parameters
center=np.array([0,0,0])
radius=1
phi_range=[0,np.pi/2]  # portion of the sphere
theta_range=[0,np.pi/2]  # portion of the sphere
point1=np.array([radius/np.sqrt(2),radius/np.sqrt(4),radius/np.sqrt(4)])  # first point on the sphere
normal1=point1/np.linalg.norm(point1)  # normal at the first point
point2=np.array([radius/np.sqrt(4),radius/np.sqrt(4),radius/np.sqrt(2)])  # second point on the sphere
normal2=point2/np.linalg.norm(point2)  # normal at the second point

# Plotting
fig=plt.figure()
ax=fig.add_subplot(111, projection='3d')
plot_sphere(ax,center,radius,phi_range,theta_range)
plot_normal_vector(ax,point1,normal1)
plot_normal_vector(ax,point2,normal2)
plot_radial_vector(ax,center,point1)
plot_radial_vector(ax,center,point2)

# Setting the aspect ratio of the plot to be equal
ax.set_box_aspect([1,1,1])

# Remove the frame and the grid
ax.grid(False)
ax.axis('off')
ax.view_init(elev=-15,azim=-36,roll=0)

fig.set_size_inches(16,16)
fig.savefig('fig_segment.png',format='png',dpi=300,bbox_inches="tight")

# Show plot
plt.show()
