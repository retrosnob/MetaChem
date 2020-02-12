import rbn
import particle
import reaction
import viz
import util

def savebondingpair(a, b):
    # Pickle the particles here
    if not isinstance(a, particle.Particle) or not isinstance(b, particle.Particle):
        raise TypeError("Both arguments must be particles.")
    else:
        util.create_pickled_particle(the_particle = a)
        util.create_pickled_particle(the_particle = b)

def getbondingpair():
    for i in range(1000):
        a = particle.Particle.new(12, 2, 0)
        b = particle.Particle.new(12, 2, 0)
        a_site, b_site = reaction.particlescanbond(a, b)
        if a_site is not None:
            print("Bond possible: ")
            print(a_site)
            print(b_site)
            break
    return (a, b)

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

def test2():
    # Create a bonding pair and bond them.
    a, b = getbondingpair()    
    a_site, b_site = reaction.particlescanbond(a, b)
    print("Particle a:")
    print(a)
    print("Particle b:")
    print(b)
    print('Bonding sites: ')
    print(a_site)
    print(b_site)
    c = particle.Particle.compose(particle1=a, particle2=b, int_site1=a_site, int_site2=b_site)
    print(c.rbn.summarystring())
    print(c)
    viz.visualize(c)


#createbondingpair()
test2()