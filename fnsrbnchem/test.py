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

def test5():
    # Bond two particles and then unbond them to check the algorithm works.
    # When I first did this test it occasionally gave a different attractor for an atom
    # after it had been unbonded from the attractor before it was bonded. I found that 
    # RBN.get_cycle() was zeroing the values of the nodes properly. Test 6 tries to check
    # that this bug no longer exists.
    a, sitea, b, siteb = getbondingpair()
    # print('Before')
    print(f'Site A: {sitea}')
    # for strnode in [f'{str(node)}' for node in sitea.nodes]:
    #     print(strnode)
    # _, attractor = rbn.RBN.get_cycle(sitea.parent_atom.rbn.nodes)
    # print(attractor)
    # print()
    print(f'Site B: {siteb}')
    # for strnode in [f'{str(node)}' for node in siteb.nodes]:
    #     print(strnode)
    # _, attractor = rbn.RBN.get_cycle(siteb.parent_atom.rbn.nodes)
    # print(attractor)
    # print(f'Sites can bond: {reaction.sitescanbond(sitea, siteb)}')
    typea1 = sitea.spike_type()
    vala1 = sitea.spike_value()
    typeb1 = siteb.spike_type()
    valb1 = siteb.spike_value()
    print(f'{typea1}, {vala1}, {typeb1}, {valb1}')
    # print()


    reaction.do_edge_swaps(sitea, siteb)
    c = particle.Composite([a] + [b])
    # print('After')
    print(f'Site A: {sitea}')
    # for strnode in [f'{str(node)}' for node in sitea.nodes]:
    #     print(strnode)
    # print()
    print(f'Site B: {siteb}')
    # for strnode in [f'{str(node)}' for node in siteb.nodes]:
    #     print(strnode)
    # print(f'Sites can bond: {reaction.sitescanbond(sitea, siteb)}')
    # print()

    reaction.do_edge_swaps(sitea, siteb)
    sitea.parent_atom.parent_composite = None
    siteb.parent_atom.parent_composite = None

    # print('And back again...')
    print(f'Site A: {sitea}') # ! This was the cause of the error but the zeroing of the nodes states fixes it.
    # for strnode in [f'{str(node)}' for node in sitea.nodes]:
    #     print(strnode)
    # _, attractor = rbn.RBN.get_cycle(sitea.parent_atom.rbn.nodes)
    # print(attractor)
    # print()
    print(f'Site B: {siteb}')
    # for strnode in [f'{str(node)}' for node in siteb.nodes]:
    #     print(strnode)
    # _, attractor = rbn.RBN.get_cycle(siteb.parent_atom.rbn.nodes)
    # print(attractor)
    # print(f'Sites can bond: {reaction.sitescanbond(sitea, siteb)}')
    typea2 = sitea.spike_type()
    vala2 = sitea.spike_value()
    typeb2 = siteb.spike_type()
    valb2 = siteb.spike_value()
    print(f'{typea2}, {vala2}, {typeb2}, {valb2}')

    # Check
    if typea1 != typea2 or vala1 != vala2 or typeb1 != typeb2 or valb1 != valb2:
        print("Error")

def test6():
    # This checks that interaction sites are the same after bonding and unbonding.
    for i in range(1000):
        a, sitea, b, siteb = getbondingpair()

        # Get spikes before
        typea1 = sitea.spike_type()
        vala1 = sitea.spike_value()
        typeb1 = siteb.spike_type()
        valb1 = siteb.spike_value()
        # print(f'{typea1}, {vala1}, {typeb1}, {valb1}')
        
        # Now bond
        reaction.do_edge_swaps(sitea, siteb)
        c = particle.Composite([a] + [b])

        # Now unbond
        reaction.do_edge_swaps(sitea, siteb)
        sitea.parent_atom.parent_composite = None
        siteb.parent_atom.parent_composite = None
        
        # Get spikes after
        typea2 = sitea.spike_type()
        vala2 = sitea.spike_value()
        typeb2 = siteb.spike_type()
        valb2 = siteb.spike_value()
        # print(f'{typea2}, {vala2}, {typeb2}, {valb2}')

        # Check
        if typea1 != typea2 or vala1 != vala2 or typeb1 != typeb2 or valb1 != valb2:
            print("Error")
            break

def test7():
    # Get two particles that meeting the bonding criteria but fail the stability criteria.
    for i in range(100):
        a, sitea, b, siteb = getbondingpair()
        assert reaction.sitescanbond(sitea, siteb)
        reaction.do_edge_swaps(sitea, siteb)
        c = particle.Composite([a] + [b])
        stable = reaction.sitescanbond(sitea, siteb)
        if not stable:
            print(f'{i} Failed stability')
            break

test7()
