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
        self.IL_spikes = self.calculate_spikes()

    def __str__(self):
        s = 'Atom ' + str(self.id) + '\n'
        ILs_idxs = [[node.loc_idx for node in IL] for IL in self.ILs] 
        s += str(ILs_idxs) + '\n'
        s += str(self.IL_spikes) + '\n'
        return s

    def calculate_spikes(self):
        """
        Each interaction list is given a "spike" value, calculated over the length of the attractor cycle
        This is done by adding 1 for each true state and -1 for each false state for each node over the length of 
        the attractor cycle.

        Returns:

        List of integers: A list of the spike values for each interaction list.
        """ 
        # Sum each node's states over the attractor cycle.
        zipped_cycle = list(zip(*self.rbn.attractor_cycle))
        summed_cycle = [sum(item) for item in zipped_cycle]

        # We now make each zero (false) contribute -1 by counting the zeroes, multiplying the count by -1 and 
        # adding it to the current node spike value. The count of the number of zeroes is equal to the length
        # of the attractor minus the summed value for this node since each non-zero value added 1 to the sum.
        node_spikes = list(map(lambda x: (len(self.rbn.attractor_cycle) - x) * -1 + x, summed_cycle))

        # We now have the spike values for each node in the IL. To get the spikes for the whole IL, we have 
        # to sum them up:
        # Get the ILs as lists of local_index values rather than objects.
        ILs_idxs = [[node.loc_idx for node in IL] for IL in self.ILs] 
        # Sum the node spikes referenced by the local_index values of the ILs.
        return list(map(sum, [[node_spikes[i] for i in lst] for lst in ILs_idxs]))

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
    print(a)
    pass

