def printatom(atom):
    print(atom)


def trybond(p1, p2):
    for site1 in p1.interaction_sites:
        if site1.available:
            for site2 in p2.interaction_sites:
                if site2.available:
                    # Test the bonding criteria
                    if canbond(site1, site2):
                        return True
    return False

def canbond(site1, site2):
    # * Watson rules for now
    diff = abs(len(site1.nodelist) + len(site2.nodelist))
    spike_sum = site1.spike_type + site2.spike_type

    if ((diff == 0 and site1.spike_type == 1 and site2.spike_type == 1)
    or (diff <= 1 and spike_sum > 2 and spike_sum < 6)
    or (diff <= 2 and site1.spike_type == 3 or site2.spike_type == 3)):
        return True
    else:
        return False
    