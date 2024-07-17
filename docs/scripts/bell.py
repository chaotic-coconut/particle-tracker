import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters
r=1.5
x=np.linspace(-3.,3.,600)
y=np.linspace(-3.,3.,600)
X,Y=np.meshgrid(x,y)
Z=np.exp(-(X**2+Y**2)/r**2)

# Apply circular boundary to limit the bottom part of the Gaussian bell
mask=np.sqrt(X**2+Y**2)<=3.  # Change this value to adjust the circular boundary radius
Z[~mask]=np.nan  # Set values outside the circular boundary to NaN

# Create a 3D plot
fig=plt.figure(figsize=(10,8))
ax=fig.add_subplot(111, projection='3d')

# Plot the vertical axis at (x=0, y=0)
ax.plot([0,0],[0,0],[0,.99],color='red',linewidth=4)

ax.plot([0,r*np.cos(-3*np.pi/8)],[0,r*np.sin(-3*np.pi/8)],[np.exp(-1),np.exp(-1)],color='green',linewidth=4)

# Generate the circle of radius r at the level corresponding to the width r
theta=np.linspace(0,2*np.pi,100)
circle_x=.99*r*np.cos(theta)
circle_y=.99*r*np.sin(theta)
circle_z=np.exp(-1)*np.ones_like(theta)  # Z level for radius r

# Plot the circle
ax.plot(circle_x,circle_y,circle_z,color='blue',linewidth=4,alpha=.9)

# Plot the Gaussian bell
ax.plot_surface(X,Y,Z,cmap='cool',alpha=.3)
#ax.plot_wireframe(X,Y,Z,color='yellow',alpha=.3)

# Remove the frame and the grid
ax.grid(False)
ax.axis('off')
ax.view_init(elev=30,azim=270,roll=0)

# Labels and title
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')

fig.set_size_inches(16,16)
fig.savefig('fig_bell.png',format='png',dpi=300,bbox_inches="tight")

# Show the plot
plt.show()
