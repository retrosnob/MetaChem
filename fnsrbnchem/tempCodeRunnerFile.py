
    for strnode in [f'{str(node)}' for node in sitea.nodes]:
        print(strnode)
    _, attractor = rbn.RBN.get_cycle(sitea.parent_atom.rbn.nodes)