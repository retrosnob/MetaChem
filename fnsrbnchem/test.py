import rbn
import particle
import reaction
import viz
import util

rbn.NodeSpace(100)

def createbondingpair():
    for i in range(1000):
        a = particle.Particle.new(12, 2, 0)
        b = particle.Particle.new(12, 2, 0)
        a_site, b_site = reaction.particlescanbond(a, b)
        if a_site is not None:
            print("Bond possible: ")
            print(a_site)
            print(b_site)
            # Pickle the particles here
            util.create_pickled_particle(the_particle = a, filename='c')
            util.create_pickled_particle(the_particle = b, filename='d')
            break

def loadbondingpair(filename1, filename2):
    a = util.load_pickled_particle(filename = filename1)
    b = util.load_pickled_particle(filename = filename2)
    return (a, b)

def test1():
    a, b = loadbondingpair('c', 'd')
    a_site, b_site = reaction.particlescanbond(a, b)
    print('Bonding sites: ')
    print(a_site)
    print(b_site)
    c = particle.Particle.compose(particle1=a, particle2=b, int_site1=a_site, int_site2=b_site)
    print(c.rbn.summarystring())
    viz.visualize(c)

#createbondingpair()
test1()