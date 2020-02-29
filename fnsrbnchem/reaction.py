def getbondablesites(p1, p2):
    """
    Checks if two particles can bond by checking all combinations
    of their available interaction sites. There is an ambiguity 
    here because the particles might be able to bond using more
    than one combination of sites but no explicit decision is made
    as to which combination is chosen. We just take the first we
    find.

    Parameters: The two candidate particles.

    Returns: The two interactions sites as a tuple or (None, None)
    if no bond is possible. Note that the sites are returned in the
     same order as the particles to which they belong.

    """
    for site1 in p1.interaction_sites:
        if site1.available:
            for site2 in p2.interaction_sites:
                if site2.available:
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
        diff = abs(len(site1) - len(site2))
        spike_sum = site1.spike_type() + site2.spike_type()
        # print("diff: " + str(diff))

        if ((diff == 0 and site1.spike_type == 1 and site2.spike_type == 1)
        or (diff <= 1 and spike_sum > 2 and spike_sum < 6)
        or (diff <= 2 and site1.spike_type == 3 or site2.spike_type == 3)):
            return True
        else:
            return False

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
    int_site_1: The interaction site from one particle.
    
    int_site_2: The interaction site from the other particle.
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
    assert from_node1 in to_node1.in_edges, f'Node {from_node1} is not in in_edges {to_node1.in_edges}'
    assert(to_node1 in from_node1.out_edges)
    assert(from_node2 in to_node2.in_edges)
    assert(to_node2 in from_node2.out_edges)

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