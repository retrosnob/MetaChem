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
        a_site, b_site = reaction.getbondablesites(a, b)
        if a_site is not None:
            # print("Bond possible: ")
            # print(a_site)
            # print(b_site)
            break
    return (a, a_site, b, b_site)

def loadbondingpair(filename1, filename2):
    a = util.load_pickled_particle(filename = filename1)
    b = util.load_pickled_particle(filename = filename2)
    return (a, b)

def getpairparticle():
    a, a_site, b, b_site, = getbondingpair()    
    return particle.Particle.compose(particle1=a, particle2=b, int_site1=a_site, int_site2=b_site)

def test1():
    a, b = loadbondingpair('c', 'd')
    a_site, b_site = reaction.getbondablesites(a, b)
    print('Bonding sites: ')
    print(a_site)
    print(b_site)
    c = particle.Particle.compose(particle1=a, particle2=b, int_site1=a_site, int_site2=b_site)
    print(c.rbn.summarystring())
    viz.visualize(c)

def test2():
    # Create a bonding pair and bond them.
    a, a_site, b, b_site = getbondingpair()    
    print("Particle a:")
    print(a)
    print("Particle b:")
    print(b)
    print('Bonding sites: ')
    print(a_site)
    print(b_site)
    c = particle.Particle.compose(particle1=a, particle2=b, int_site1=a_site, int_site2=b_site)
    print(a_site)
    print(b_site)
    print(c.rbn.summarystring())
    print(c)
    viz.visualize(c)

def test3():
    # Build a particle with three atoms in it.
    pairs = []
    triples = []
    for i in range(100):
        pairs.append(getpairparticle())
    while not triples:
        a = particle.Particle.new(12, 2, 0)
        for pair in pairs:
            a_site, pair_site = reaction.getbondablesites(a, pair)
            if a_site is not None and pair_site is not None:
                triples.append(particle.Particle.compose(particle1=a, particle2=pair, int_site1=a_site, int_site2=pair_site))
                break
    print(triples[0])

def test4():
    # Bond two particles and then unbond them to check the algorithm works.
    a, sitea, b, siteb = getbondingpair()
    print('Before')
    print(f'Site A: {sitea}')
    print(f'Site B: {siteb}')
    c = particle.Particle.compose(particle1=a, particle2=b, int_site1=sitea, int_site2=siteb)
    print('After')
    print(f'Site A: {sitea}')
    print(f'Site B: {siteb}')
    pass


test4()
