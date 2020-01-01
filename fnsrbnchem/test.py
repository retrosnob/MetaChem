import rbn
import particle
import reaction
import viz
import util

rbn.NodeSpace(100)

count = 0
for i in range(1000):
    a = particle.Particle.new(12, 2, 0)
    b = particle.Particle.new(12, 2, 0)
    a_site, b_site = reaction.particlescanbond(a, b)
    if a_site is not None:
        print("Bond possible")
        # Pickle the particles here
        util.create_pickled_particle(the_particle = a, filename = 'a')
        util.create_pickled_particle(the_particle = b, filename = 'b')
        break
print(a)
print(b)
print(a_site)
for node in a_site.nodelist:
    print(node)
print(b_site)
# c = reaction.bond(a, b, a_site, b_site)
# c = particle.Particle.compose(particle1=a, particle2=b)
# viz.visualize(c)

