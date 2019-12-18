import rbn
import particle
import reaction
import viz

rbn.NodeSpace(100)

count = 0
for i in range(1000):
    a = particle.Particle.new(12, 2)
    b = particle.Particle.new(12, 2)
    a_site, b_site = reaction.particlescanbond(a, b)
    if a_site is not None:
        print("Bond possible")
        break
# print(a)
# print(b)
# print(a_site)
# print(b_site)
# c = reaction.bond(a, b, a_site, b_site)
c = Particle.compose(particle1=a, particle2=b)
viz.visualize(c)

