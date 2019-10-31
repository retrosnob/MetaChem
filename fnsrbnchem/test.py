import rbn
import atom
import reaction

rbn.NodeSpace(100)

count = 0
for i in range(100):
    a = atom.Atom(0, 12, 2)
    b = atom.Atom(0, 12, 2)
    if reaction.trybond(a, b):
        print("Bond possible")
        count += 1
print(str(count) + " bonds possible.")
