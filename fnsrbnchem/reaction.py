import particle

def getbondablesites(p1, p2):
    """
    Checks if two Atoms can bond by checking all combinations
    of their available interaction sites. There is an ambiguity 
    here because the Atoms might be able to bond using more
    than one combination of sites but no explicit decision is made
    as to which combination is chosen. We just take the first we
    find.

    Parameters: The two candidate Atoms.

    Returns: The two interactions sites as a tuple or (None, None)
    if no bond is possible. Note that the sites are returned in the
     same order as the Atoms to which they belong.

    """
    for site1 in p1.interaction_sites:
        if not site1.bondedto:
            for site2 in p2.interaction_sites:
                if not site2.bondedto:
                    # Test the bonding criteria
                    if sitescanbond(site1, site2):
                        return site1, site2
    return (None, None)

def sitescanbond(site1, site2):
    """
    Checks if two interaction sites can bond.

    Parameters: The two candidate interaction sites.

    Returns: True if sites can bond, false otherwise.

    """
    # * Watson rules for now
    if len(site1) < 2 or len(site2) < 2:
        return False
    else:
        spike_sum = abs(site1.spike_value() + site2.spike_value())
        spike_type1 = site1.spike_type()
        spike_type2 = site2.spike_type()

        if ((spike_sum == 0 and spike_type1 == 1 and spike_type2 == 1)
        or (spike_sum <= 1 and (spike_type1 + spike_type2) > 2 and (spike_type1 + spike_type2) < 6)
        or (spike_sum <= 2 and spike_type1 == 3 or spike_type2 == 3)):
            return True
        else:
            return False

def bond(site1, site2):
    assert site1.bondedto is None and site2.bondedto is None, "Trying to bonded sites that are already bonded."
    do_edge_swaps(site1, site1)
    site1.bondedto, site2.bondedto = site2, site1
    return particle.Composite(particle.Composite.traverse(site1.parent_atom))
    # Need to set new parent composites here otherwise wrong attractor will be calculated.
    # Create composite from each atom by traversing.
    # Then return two composites.

def break_bond(site1, site2):
    assert site1.bondedto is site2 and site2.bondedto is site1, "Trying to break bond between unbonded Atoms."
    do_edge_swaps(site1, site1)
    site1.bondedto, site2.bondedto = None, None
    # Need to set new parent composites here otherwise wrong attractor will be calculated.
    # Create composite from each atom by traversing.
    # Then return two composites.
    return particle.Composite(particle.Composite.traverse(site1.parent_atom)), particle.Composite(particle.Composite.traverse(site2.parent_atom))

def is_stable(site1, site2):
    # For the time being....
    return sitescanbond(site1, site2)

def do_edge_swaps(int_site1, int_site2):
    """
    Parameters:
    Description:
    
    Swaps nodes to form a pair of bonded interaction sites.
    
    a <- b <- c <- d <- e
    w <- x <- y <- z
    becomes
    a <- x <- c <- z <- e
    
    w <- b <- y <- d
    Parameters:
    int_site_1: The interaction site from one Atom.
    
    int_site_2: The interaction site from the other Atom.
    """
    # THIS IS WHERE BONDING TAKES PLACE
    # Find the nodes at the odd indices running up to the length of the shorter of the two interaction sites involved in the bond, then switch the edges, e.g in _switch_edges(2, 1, 4, 3):

    # 1 <- 2   becomes    1 <- 4
    # 3 <- 4              3 <- 2
    indices = [i for i in range(1, min(len(int_site1.nodes), len(int_site2.nodes)))]
    for index in indices:
        # Make appropriate changes to the inward edges lists for each node.s
        _switch_edges(int_site1.nodes[index], int_site1.nodes[index-1], int_site2.nodes[index], int_site2.nodes[index-1])
    for index in indices:
        if (index) % 2 == 1: # ie 1, 3, 5, 7, ...
            # Swap every other node in the interaction sites themselves. Have to operate on index - 1 node so as not to 
            # affect the edge switching.
            int_site1.nodes[index], int_site2.nodes[index] = int_site2.nodes[index], int_site1.nodes[index] 
            # int_site1.nodes[index], int_site2.nodes[index] = int_site2.nodes[index], int_site1.nodes[index] 

def _switch_edges(from_node1, to_node1, from_node2, to_node2):
    """
    In an RBN in which there is an edge from from_node1 to to_node1 and from from_node2 to
    to_node2, this method creates an edge from from_node1 to to_node2 and from from_node2 
    to to_node1, removing the original edges.
    It is important to remember that this method changes in_edges and out_edges, not any interaction site. 
    Spike details should be recalculated after a called to switch_edges.
    Parameters:
    from_node1
    to_node1
    from_node2
    to_node2
    """
    
    # This method preserves k and so doesn't corrupt the RBN and can be called from outside the class.

    # Assert that the nodes are where they should be to start with.
    assert from_node1 in to_node1.in_edges, f'Node {from_node1} is not in in_edges {[str(node) for node in to_node1.in_edges]}'
    assert to_node1 in from_node1.out_edges, f'Node {to_node1} is not in out_edges {[str(node) for node in from_node1.out_edges]}'
    assert from_node2 in to_node2.in_edges, f'AAA Node {from_node2} is not in in_edges {[str(node) for node in to_node2.in_edges]}'
    assert to_node2 in from_node2.out_edges, f'BBB Node {to_node2} is not in out_edges {[str(node) for node in from_node2.out_edges]}'

    # Add nodes to the corresponding in_edges and out_edges lists
    to_node1.in_edges.insert(to_node1.in_edges.index(from_node1), from_node2)
    to_node2.in_edges.insert(to_node2.in_edges.index(from_node2), from_node1)
    from_node1.out_edges.insert(from_node1.out_edges.index(to_node1), to_node2)
    from_node2.out_edges.insert(from_node2.out_edges.index(to_node2), to_node1)

    # Remove nodes from the corresponding in_edges and out_edges lists
    to_node1.in_edges.remove(from_node1)
    to_node2.in_edges.remove(from_node2)
    from_node1.out_edges.remove(to_node1)
    from_node2.out_edges.remove(to_node2)

    # Assert that the nodes are where they should be to start with.
    assert(from_node1 in to_node2.in_edges)
    assert(to_node1 in from_node2.out_edges)
    assert(from_node2 in to_node1.in_edges)
    assert(to_node2 in from_node1.out_edges)      

if __name__ == '__main__':
    print('reaction.py invoked as main')
    pass