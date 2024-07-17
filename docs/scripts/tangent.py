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
    ax.plot_wireframe(x+.01,y+.01,z+.01,color='w',lw=2,alpha=1,rstride=100,cstride=100)

# Function to plot the tangent plane
def plot_tangent_plane(ax,point,normal,size):
    d=-point.dot(normal)
    xx,yy=np.meshgrid(np.linspace(point[0]-size,point[0]+size,10),
                      np.linspace(point[1]-size,point[1]+size,10))
    zz=(-normal[0]*xx-normal[1]*yy-d)/normal[2]
    ax.plot_surface(xx,yy,zz,color='r',alpha=1)

# Function to plot the tangent disk
def plot_tangent_disk(ax,point,normal,radius):
    # Create orthogonal vectors to the normal
    if normal[2]!=0:
        a=np.array([-normal[1],normal[0],0])
    else:
        a=np.array([0,-normal[2], normal[1]])
    a/=np.linalg.norm(a)
    b=np.cross(normal,a)
    
    theta=np.linspace(0,2*np.pi,100)
    r=np.linspace(0,radius,50)
    theta,r=np.meshgrid(theta,r)
    
    # Parametric equations for the disk
    X=point[0]+r*(a[0]*np.cos(theta)+b[0]*np.sin(theta))
    Y=point[1]+r*(a[1]*np.cos(theta)+b[1]*np.sin(theta))
    Z=point[2]+r*(a[2]*np.cos(theta)+b[2]*np.sin(theta))
    
    ax.plot_surface(X,Y,Z,color='r',alpha=1,rstride=100,cstride=100)

# Parameters
center=np.array([0,0,0])
radius=1
phi_range=[0,np.pi/2]  # portion of the sphere
theta_range=[0,np.pi/2]  # portion of the sphere
point=np.array([radius/np.sqrt(3),radius/np.sqrt(3),radius/np.sqrt(3)])  # point on the sphere
normal=point/np.linalg.norm(point)  # normal at the point
#size=.1  # size of the tangent plane
disk_radius=.25  # radius of the tangent disk

# Plotting
fig=plt.figure()
ax=fig.add_subplot(111, projection='3d')
plot_sphere(ax,center,radius,phi_range,theta_range)
#plot_tangent_plane(ax,point,normal,size)
plot_tangent_disk(ax, point,normal,disk_radius)

# Plot the point of tangency
ax.scatter(point[0]+.02,point[1]+.02,point[2]+.02,color='k',s=50)

# Setting the aspect ratio of the plot to be equal
ax.set_box_aspect([1,1,1])

# Remove the frame and the grid
ax.grid(False)
ax.axis('off')
ax.view_init(elev=30,azim=45,roll=0)

fig.set_size_inches(16,16)
fig.savefig('fig_tangent.png',format='png',dpi=300,bbox_inches="tight")

# Show plot
plt.show()
