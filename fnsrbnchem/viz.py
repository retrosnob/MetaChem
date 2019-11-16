import matplotlib as mpl
# mpl.use("Qt5Agg")
# print(mpl.get_backend())
# valid strings are ['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 
# 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX',
# 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']

import matplotlib.pyplot as plt
import atom
import rbn
print(mpl.get_backend())
plt.switch_backend("pdf")

def visualize(particle):
    print("Visualize")

    # for thisnode in particle.rbn.nodes:
    #     for thatnode in thisnode.in_edges:
    #         pass
    # fig = plt.figure()
    print(mpl.get_backend())
    plt.plot([1, 2, 3, 4, 5])
    # fig.show()

if __name__ == "__main__":
    print("viz.py invoked as script...")
    rbn.NodeSpace(100)
    a = atom.Atom(0, 12, 2)
    visualize(a)
    

