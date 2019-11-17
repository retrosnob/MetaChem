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

    # for thisnode in particle.rbn.nodes:
    #     for thatnode in thisnode.in_edges:
    #         pass
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
    # Plot lines between points where edges exist.
    # plt.plot([thetas[0],thetas[1]],[radii[0],radii[1]],'k-')
    # plt.arrow(thetas[0], radii[0], thetas[1]-thetas[0], 0, head_length=0.05, length_includes_head=True, head_width=0.05, zorder=5)

    for thisnode in particle.rbn.nodes:
        for thatnode in thisnode.in_edges:
            # plt.arrow(thetas[thatnode.loc_idx], radii[0], thetas[thisnode.loc_idx]-thetas[thatnode.loc_idx], 0, head_length=0.05, length_includes_head=True, head_width=0.05, zorder=5)
            ax.annotate('', # A textless annotation is just an arrow.
                xy=(thetas[thisnode.loc_idx], radius),  # to this point
                xytext=(thetas[thatnode.loc_idx], radius),    # from this point
                arrowprops=dict(arrowstyle='->', shrinkA=15, shrinkB=15)
                )


    plt.axis('off')

    fig.savefig('temp.png')

if __name__ == "__main__":
    print("viz.py invoked as script...")
    rbn.NodeSpace(100)
    a = atom.Atom(0, 12, 2)
    print(a.rbn)
    visualize(a)
    

