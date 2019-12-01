import rbn
import atom
import reaction
import viz

rbn.NodeSpace(100)

count = 0
for i in range(1000):
    a = atom.Atom(0, 12, 2)
    b = atom.Atom(0, 12, 2)
    a_site, b_site = reaction.particlescanbond(a, b)
    if a_site is not None:
        print("Bond possible")
        break
# print(a)
# print(b)
# print(a_site)
# print(b_site)
# c = reaction.bond(a, b, a_site, b_site)
c = atom.Atom(particle1=a, particle2=b)
viz.visualize(c)

