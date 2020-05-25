import rbn
import particle
import reaction
import viz
import util

def savebondingpair(a, b):
    # Pickle the Atoms here
    if not isinstance(a, particle.Atom) or not isinstance(b, particle.Atom):
        raise TypeError("Both arguments must be Atoms.")
    else:
        util.create_pickled_Atom(the_Atom = a)
        util.create_pickled_Atom(the_Atom = b)

def getbondingpair():
    for i in range(1000):
        a = particle.Atom.new(12, 2)
        b = particle.Atom.new(12, 2)
        a_site, b_site = reaction.getbondablesites(a, b)
        if a_site is not None:
            # print("Bond possible: ")
            # print(a_site)
            # print(b_site)
            break
    return (a, a_site, b, b_site)

def loadbondingpair(filename1, filename2):
    a = util.load_pickled_Atom(filename = filename1)
    b = util.load_pickled_Atom(filename = filename2)
    return (a, b)

def getpairAtom():
    a, a_site, b, b_site, = getbondingpair()    
    return Atom.Composite(int_site1=a_site, int_site2=b_site)

def get_composites_of_size(size):
    """
    Returns a list of composites of size size. Starts with 1000 2-composites. 
    At the moment, if size = 10 it tends to return only a handful of 10-composites.

    """
    if size < 2:
        print ('Size has to be >= 2')
        return 

    composites = []
    for i in range(10000):
        a, sitea, b, siteb = getbondingpair()
        assert reaction.sitescanbond(sitea, siteb)
        # reaction.do_edge_swaps(sitea, siteb)
        # c = Atom.Composite([a] + [b])
        # sitea.bondedto, siteb.bondedto = siteb, sitea
        c = reaction.bond(sitea, siteb)
        stable = reaction.sitescanbond(sitea, siteb)
        if stable:
            composites.append(c)
    for i in range(size - 2):
        print(len(composites))
        compositesold = composites
        composites = []
        a = particle.Atom.new(12, 2)
        for composite in compositesold:
            if len(composite.atoms) != len(particle.Composite.traverse(composite.atoms[0])):
                print("error")
            for atom in composite.atoms:
                sitea, siteb = reaction.getbondablesites(a, atom)
                if sitea and siteb:
                    # reaction.do_edge_swaps(sitea, siteb)
                    # sitea.bondedto, siteb.bondedto = siteb, sitea
                    # c = Atom.Composite([a] + composite.atoms)
                    c = reaction.bond(sitea, siteb)
                    composites.append(c)
                    a = particle.Atom.new(12, 2)
                    break
    for c in composites:    
        print([atom.id for atom in c.atoms])
        for i in range(size):
            print([atom.id for atom in particle.Composite.traverse(c.atoms[i])])
        if len(c.atoms) != len(particle.Composite.traverse(c.atoms[0])) or len(c.atoms) != len(particle.Composite.traverse(c.atoms[1])):
            print("error")
    print('done')
    return composites

def test1():
    a, b = loadbondingpair('c', 'd')
    a_site, b_site = reaction.getbondablesites(a, b)
    print('Bonding sites: ')
    print(a_site)
    print(b_site)
    c = particle.Atom.compose(Atom1=a, Atom2=b, int_site1=a_site, int_site2=b_site)
    print(c.rbn.summarystring())
    viz.visualize(c)

def test2():
    # Create a bonding pair and bond them.
    a, a_site, b, b_site = getbondingpair()    
    print("Atom a:")
    print(a)
    print("Atom b:")
    print(b)
    print('Bonding sites: ')
    print(a_site)
    print(b_site)
    c = particle.Atom.compose(Atom1=a, Atom2=b, int_site1=a_site, int_site2=b_site)
    print(a_site)
    print(b_site)
    print(c.rbn.summarystring())
    print(c)
    viz.visualize(c)

def test3():
    # Build a Atom with three atoms in it.
    pairs = []
    triples = []
    for i in range(100):
        pairs.append(getpairAtom())
    while not triples:
        a = particle.Atom.new(12, 2)
        for pair in pairs:
            a_site, pair_site = reaction.getbondablesites(a, pair)
            if a_site is not None and pair_site is not None:
                triples.append(particle.Atom.compose(Atom1=a, Atom2=pair, int_site1=a_site, int_site2=pair_site))
                break
    print(triples[0])

def test4():
    # Bond two Atoms and then unbond them to check the algorithm works.
    a, sitea, b, siteb = getbondingpair()
    print('Before')
    print(f'Site A: {sitea}')
    print(f'Site B: {siteb}')
    c = particle.Atom.compose(Atom1=a, Atom2=b, int_site1=sitea, int_site2=siteb)
    print('After')
    print(f'Site A: {sitea}')
    print(f'Site B: {siteb}')
    pass

def test5():
    # Bond two Atoms and then unbond them to check the algorithm works.
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
    # Get two Atoms that meeting the bonding criteria but fail the stability criteria.
    for i in range(100):
        a, sitea, b, siteb = getbondingpair()
        assert reaction.sitescanbond(sitea, siteb)
        reaction.do_edge_swaps(sitea, siteb)
        c = particle.Composite([a] + [b])
        stable = reaction.sitescanbond(sitea, siteb)
        if not stable:
            print(f'{i} Failed stability')
            break

def test8():
    """
    Build a composite of more than 2 atoms.
    Repeat 100:
        Get a composite c that is stable
        Add it to a list composites
        Get an atom a
        for each composite in composites:
            for each atom in composite:
                try to bond a to atom (which uses the composite's attractor)
    """
    composites = []
    for i in range(100):
        a, sitea, b, siteb = getbondingpair()
        assert reaction.sitescanbond(sitea, siteb)
        reaction.do_edge_swaps(sitea, siteb)
        c = particle.Composite([a] + [b])
        sitea.bondedto, siteb.bondedto = siteb, sitea
        stable = reaction.sitescanbond(sitea, siteb)
        if not stable:
            print(f'{i} Failed stability')
        else:
            composites.append(c)
            a = particle.Atom.new(12, 2)
            for composite in composites:
                for atom in composite.atoms:
                    sitea, siteb = reaction.getbondablesites(a, atom)
                    if sitea and siteb:
                        reaction.do_edge_swaps(sitea, siteb)
                        sitea.bondedto, siteb.bondedto = siteb, sitea
                        c = particle.Composite([a] + composite.atoms)
                        print(c.atoms)
                        break
                    else:
                        print("No bonding possible between composite and atom.")

def test9():
    # ! Problems here. The number of atoms in self.atoms disagrees with the atoms returned
    # ! by the Atom traversal. The traversal is correct, so bonds are happening where 
    # ! they shouldn't be.
    # * The problem was that I kept going through the 2-composites over and over and the 
    # * the second and subsequent times through some of them were already bonded.
    # Build some larger molecules to test the traversal defined in Composite.
    composites = []
    for i in range(10000):
        a, sitea, b, siteb = getbondingpair()
        assert reaction.sitescanbond(sitea, siteb)
        reaction.do_edge_swaps(sitea, siteb)
        c = particle.Composite([a] + [b])
        sitea.bondedto, siteb.bondedto = siteb, sitea
        stable = reaction.sitescanbond(sitea, siteb)
        if stable:
            composites.append(c)
    print(len(composites))
    composites3 = []
    a = particle.Atom.new(12, 2)
    for composite in composites:
        if len(composite.atoms) != len(particle.Composite.traverse(composite.atoms[0])):
            print("error")
        for atom in composite.atoms:
            sitea, siteb = reaction.getbondablesites(a, atom)
            if sitea and siteb:
                reaction.do_edge_swaps(sitea, siteb)
                sitea.bondedto, siteb.bondedto = siteb, sitea
                c = particle.Composite([a] + composite.atoms)
                composites3.append(c)
                a = particle.Atom.new(12, 2)
                break
    # for composite in composites3:
    #     print([atom.id for atom in particle.Composite.traverse(composite.atoms[0])])
    #     print([atom.id for atom in particle.Composite.traverse(composite.atoms[1])])
    #     print([atom.id for atom in particle.Composite.traverse(composite.atoms[2])])
    #     print()
    print(len(composites3))
    composites4 = []
    a = particle.Atom.new(12, 2)
    for composite in composites3:
        if len(composite.atoms) != len(particle.Composite.traverse(composite.atoms[0])):
            print("error")
        for atom in composite.atoms:
            sitea, siteb = reaction.getbondablesites(a, atom)
            if sitea and siteb:
                reaction.do_edge_swaps(sitea, siteb)
                sitea.bondedto, siteb.bondedto = siteb, sitea
                c = particle.Composite([a] + composite.atoms)
                composites4.append(c)
                a = particle.Atom.new(12, 2)
                break
    print(len(composites4))

def test10():
    # Create 100 2-atom composites and check that their len(self.atoms) agrees with 
    # the number of atoms returned by the traversal.
    composites = []
    for i in range(1000):
        a, sitea, b, siteb = getbondingpair()
        assert reaction.sitescanbond(sitea, siteb)
        reaction.do_edge_swaps(sitea, siteb)
        c = particle.Composite([a] + [b])
        sitea.bondedto, siteb.bondedto = siteb, sitea
        stable = reaction.sitescanbond(sitea, siteb)
        print([atom.id for atom in c.atoms])
        print([atom.id for atom in particle.Composite.traverse(c.atoms[0])])
        print([atom.id for atom in particle.Composite.traverse(c.atoms[1])])
        if len(c.atoms) != len(particle.Composite.traverse(c.atoms[0])) or len(c.atoms) != len(particle.Composite.traverse(c.atoms[1])):
            print("error")
    print('done')

def test11():
    # Check decomposition of Atom.
    # Start with 8-composite.
    composite = get_composites_of_size(10)[0]
    print(composite.atoms)
    print(composite.is_valid())
    decomposed = particle.Composite.decompose(composite.atoms[3])
    print(decomposed)

def test12():
    a = particle.Atom.new(12, 2)

def test13():
    # Muck around with some 6, 1 RBNs
    cyc_lens = []
    count = 1000
    for i in range(count):
        atom = rbn.RBN.new(6, 4)
        # print(atom.summarystring())
        cyc = rbn.RBN.get_cycle(atom.nodes)[1]
        # print(cyc)
        cyc_lens.append(len(cyc))
    print(sum(cyc_lens)/count)
    for i in range(1, max(cyc_lens) + 1):
        print(f'{i} {cyc_lens.count(i)}')

def num_unique_atoms():
    lst = []
    for i in range(2**6):
        st = unique_orientations(i)
        if not st in lst:
            lst.append(st)
    return lst

def list_to_int(lst):
    return sum([2**i * val for i, val in enumerate(lst[::-1])])

def int_to_list(i):
    return list(map(int, f'{bin(i)[2:]:>06}'))

def unique_orientations(i):
    current = int_to_list(i)
    config_list = []

    # y1 at top: x1, x2, y1, y2, z1, z2
    next = [current[0], current[1], current[2], current[3], current[4], current[5]]
    config_list.append(list_to_int(next))
    # rotate 3 times in new vertical axis = 0 -> 4 -> 1 -> 5 -> 0
    for _ in range(3):
        next = [next[5], next[4], next[2], next[3], next[0], next[1]]
        config_list.append(list_to_int(next))

    # y2 at top: z1, z2, y2, y1, x1, x2
    next = [current[4], current[5], current[3], current[2], current[0], current[1]]
    config_list.append(list_to_int(next))
    # rotate 3 times in new vertical axis = 0 -> 4 -> 1 -> 5 -> 0
    for _ in range(3):
        next = [next[5], next[4], next[2], next[3], next[0], next[1]]
        config_list.append(list_to_int(next))

    # x1 at top: y2, y1, x1, x2, z1, z2
    next = [current[3], current[2], current[0], current[1], current[4], current[5]]
    config_list.append(list_to_int(next))
    # rotate 3 times in new vertical axis = 0 -> 4 -> 1 -> 5 -> 0
    for _ in range(3):
        next = [next[5], next[4], next[2], next[3], next[0], next[1]]
        config_list.append(list_to_int(next))

    # x2 at top: z1, z2, x2, c1, y2, y1
    next = [current[4], current[5], current[1], current[0], current[3], current[2]]
    config_list.append(list_to_int(next))
    # rotate 3 times in new vertical axis = 0 -> 4 -> 1 -> 5 -> 0
    for _ in range(3):
        next = [next[5], next[4], next[2], next[3], next[0], next[1]]
        config_list.append(list_to_int(next))

    # z1 at top: x1, x2, z1, z2, y2, y1
    next = [current[0], current[1], current[4], current[5], current[3], current[2]]
    config_list.append(list_to_int(next))
    # rotate 3 times in new vertical axis = 0 -> 4 -> 1 -> 5 -> 0
    for _ in range(3):
        next = [next[5], next[4], next[2], next[3], next[0], next[1]]
        config_list.append(list_to_int(next))

    # z2 at top: y2, y1, z2, z1, x1, x2
    next = [current[3], current[2], current[5], current[4], current[0], current[1]]
    config_list.append(list_to_int(next))
    # rotate 3 times in new vertical axis = 0 -> 4 -> 1 -> 5 -> 0
    for _ in range(3):
        next = [next[5], next[4], next[2], next[3], next[0], next[1]]
        config_list.append(list_to_int(next))

    return set(config_list)
    
test13()