import matplotlib as mpl
# mpl.use("Qt5Agg")
# print(mpl.get_backend())
# valid strings are ['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 
# 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX',
# 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']

import matplotlib.pyplot as plt
import math
import atom
import rbn

def visualize(particle):
    """
    Set the scatter plot to have axes -100 up to 100 for both x and y.
    Calculate the angle between nodes as 2pi/N. Calculate the x, y coord
    """
    print("Visualize")

    n = particle.rbn.n
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)

    # Create the theta values as fractions of 2pi.
    thetas = [i * 2 * math.pi/n for i in range(n)]
    # All r values are just 1.
    radius = 10
    radii = n * [radius]

    # Plot the points.
    plt.polar(thetas, radii, 'o', markersize=15)

    for thisnode in particle.rbn.nodes:
        for thatnode in thisnode.in_edges:
            ax.annotate('', # A textless annotation is just an arrow.
                xy=(thetas[thisnode.loc_idx], radius),  # TO this point
                xytext=(thetas[thatnode.loc_idx], radius),    # FROM this point
                arrowprops=dict(arrowstyle='->', shrinkA=15, shrinkB=15)
                )


    plt.axis('off')
    # For some reason I can't use plt.show() because matplotlib thinks I'm headless.
    fig.savefig('temp.png')

if __name__ == "__main__":
    print("viz.py invoked as script...")
    rbn.NodeSpace(100)
    a = atom.Atom(0, 12, 2)
    print(a.rbn)
    visualize(a)
    

