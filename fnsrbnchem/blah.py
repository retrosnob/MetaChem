        for index in [i for i in range(1, min(len(nodes1), len(nodes2)))]:
            obj._switch_edges(nodes1[index], nodes1[index-1], nodes2[index], nodes2[index-1])