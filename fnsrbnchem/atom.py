import rbn
from collections import defaultdict
from operator import add

class Atom:

    """ Represents an Atom in the Frozen-Node Spiky RBN Chemistry.

    To do:

    Should bonding be a method in the Atom or in the Chemistry itself?

    """

    def __init__(self, id, n, k):

        """ The Atom class constructor.

        Parameters:

            id = An id number for the atom.
            n = The number of nodes in the Atom. This is passed directly to the RBN constructor.
            k = The indegree of each node. Passed to the RBN constructor.

        Returns: None

        """
        self.id = id
        self.rbn = rbn.RBN(n, k, 0, rbn.NodeSpace(100))
        self.ILs = self.calculate_interaction_lists()

    def __str__(self):
        s = 'Atom ' + str(self.id) + '\n'
        s += str(self.interaction_lists) + '\n'
        s += str(self.interaction_list_spikes) + '\n'
        return s

    def calculate_spikes(self):
        """
        Each interaction list is given a "spike" value, calculated over the length of the attractor cycle
        This is done by adding 1 for each true state and -1 for each false state for each node over the length of 
        the attractor cycle.

        Returns:

        interaction_list_spikes (list of integers): A list of the spike values for each interaction list.
        # for node_list in self.rbn.attractor_cycle:
        #     print(node_list)
        zipped_cycle = list(zip(*self.rbn.attractor_cycle))
        spikes = [sum(tpl) for tpl in zipped_cycle]
        # TO DO: Make each 0 contribute -1
        # This can be calculated since the number of zeroes in the spike is equal to
        #  attractor cycle length - current node spike value
        # so if we multiply this value by -1 and add it to the current node spike value then we're done
        node_spikes = list(map(lambda x: (len(self.rbn.attractor_cycle) - x) * -1 + x, spikes))
        # self.spikes now contains the spike value for each node
        # In order to get the spikes for each IL, we have to sum up the spikes for each node on that IL
        total = 0
        interaction_list_spikes = []
        for interaction_list in self.interaction_lists:
            for node in interaction_list:
                total += node_spikes[node]
            interaction_list_spikes.append(total)
            total = 0
        # self.interaction_list_spikes = [reduce(add, self.node_spikes[node]) for interaction_list in self.interaction_lists for node in interaction_list]
        return interaction_list_spikes   
                """
         

    def calculate_interaction_lists(self):
        """
        Creates and calculates the interaction list for an Atom.

        Returns:
          
        ILs (List of lists of RBNNode objects): A list of the interaction lists for this Atom.
        """
        ILs = []
        # Sort nodes in ascending order of influence.
        sortednodes = sorted([node for node in self.rbn.nodes], key = lambda node: len(node.out_edges))
        while sortednodes:
            node = sortednodes.pop(0)
            IL = []
            IL.append(node)
            nextnodes = [nextnode for nextnode in node.in_edges if nextnode in sortednodes]
            while nextnodes:
                nextnode = nextnodes[0]
                sortednodes.remove(nextnode)
                IL.append(nextnode)
                node = nextnode
                nextnodes = [nextnode for nextnode in node.in_edges if nextnode in sortednodes]
            ILs.append(IL)
        
        return ILs

        # Create interaction lists
        # Sort nodes ascending on number of outgoing edges
        # while N is not empty do
        #   Remove first node n from N
        #   Create new Interaction List ILi
        #   Add n to ILi
        #   while ∃ n′ ∈ N where n is an input to n′
        #       Remove n′from N
        #       Add n′to ILi
        #       n ← n′
        #   end
        #   i ++
        # end


if __name__ == "__main__":
    print("Atom.py invoked as script...")
    a = Atom(1, 12, 2)
    print([[node.loc_idx for node in il] for il in a.ILs])
    pass

