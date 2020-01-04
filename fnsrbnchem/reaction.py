def printatom(atom):
    print(atom)


def particlescanbond(p1, p2):
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
                    if _sitescanbond(site1, site2):
                        return site1, site2
    return (None, None)

def _sitescanbond(site1, site2):
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
        spike_sum = site1.spike_type + site2.spike_type
        print("diff: " + str(diff))

        if ((diff == 0 and site1.spike_type == 1 and site2.spike_type == 1)
        or (diff <= 1 and spike_sum > 2 and spike_sum < 6)
        or (diff <= 2 and site1.spike_type == 3 or site2.spike_type == 3)):
            return True
        else:
            return False
