import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rc

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times'],'size':30})
rc('text',usetex=True)
rc('axes',linewidth=2)

# Generate a realistic-looking 2D field using a combination of sine and cosine functions
x=np.linspace(-3.,3.,100)
y=np.linspace(-3.,3.,100)
X,Y=np.meshgrid(x,y)
Z=np.sin(X)*np.sin(Y)

# Create the filled contour plot
fig=plt.figure(figsize=plt.figaspect(1))
ax=fig.add_subplot(111)

contour=ax.contourf(X,Y,Z,cmap='coolwarm',levels=100)
#plt.colorbar(contour)

# Generate a rectangular set of dots
rect_x=np.linspace(-2.1,1.9,5)
rect_y=np.linspace(-2.1,1.9,5)
rect_X,rect_Y=np.meshgrid(rect_x,rect_y)

# Plot the dots on top of the contour plot
ax.scatter(rect_X,rect_Y,color='black')

# Add a circle of radius 1.5
circle=patches.Circle((.4,.4),radius=1.5,edgecolor='green',facecolor='none',linewidth=4)
plt.gca().add_patch(circle)
ax.scatter([.4],[.4],color='green',s=120)

# Add labels and title
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
#plt.title('Filled Contour Plot with Rectangular Set of Dots')

fig.set_size_inches(16,16)
fig.savefig('fig_field.png',format='png',dpi=300,bbox_inches="tight")

# Show the plot
plt.show()

