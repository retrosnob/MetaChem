import rbn
from collections import defaultdict
from operator import add

class Particle:

    """
    Represents an Atom in the Frozen-Node Spiky RBN Chemistry.

    """

    # TODO  Need a way to determine which of the ILs are free to bond and which are occupied.

    def __init__(self):

        """
        The Particle class constructor.

        Parameters:

            id = An id number for the atom.
            n = The number of nodes in the Atom. This is passed directly to the RBN constructor.
            k = The indegree of each node. Passed to the RBN constructor.

        Returns: None

        """
    
    @classmethod
    def compose(cls, particle1, particle2, int_site1, int_site2):
        """
        The Particle factory method.

        Parameters:

        Two particles from which to create the composite.

        Returns: The new composite particle.

        """
        obj = cls()
        obj.rbn = rbn.RBN.compose(rbn1=particle1.rbn, rbn2=particle2.rbn)
        obj.id = "(" + str(particle1.id) + "." + str(particle2.id) + ")"

        # THIS IS WHERE BONDING TAKES PLACE
        # Find the nodes at the odd indices running up to the length of the shorter of the two interaction sites involved in the bond, then switch the edges, e.g in _switch_edges(2, 1, 4, 3):

        # 1 <- 2   becomes    1 <- 4
        # 3 <- 4              3 <- 2
        for index in [i for i in range(1, min(len(int_site1), len(int_site2)))]:
            obj.rbn._switch_edges(int_site1[index], int_site1[index-1], int_site2[index], int_site2[index-1])
        # obj._initialize()
        return obj

    @classmethod
    def new(cls, n, k, id):
        obj = cls()
        obj.rbn = rbn.RBN.new(n, k, rbn.NodeSpace.getInstance())
        obj.id = id
        obj._initialize()
        return obj

    @classmethod
    def fromRBN(cls, rbn, id):
        obj = cls()
        obj.rbn = rbn
        obj.id = id
        obj._initialize()
        return obj

    def _initialize(self):
        """
        Initialization of spikes and interaction sites which is required both for new particles
        and composite particles.
        """
        self.ILs = self._calculate_interaction_sites()
        # self.spike_values = self.krastev_spikes()

        self.spike_values, self.spike_types = self._watson_spikes()
        assert(len(self.ILs) == len(self.spike_values) == len(self.spike_types))

        # The current interaction site method return parallel lists. Here we convert those
        # into a single list of InteractionSite objects. 
        # TODO Decide if parallel lists or single list of objects is best and tidy up. One advantage
        # TODO of the object approach is that the interaction site can be marked as available/unavailable.
        self.interaction_sites = []
        for i, _ in enumerate(self.ILs):
            interaction_site = Interation_Site(self.ILs[i], self.spike_values[i], self.spike_types[i])
            self.interaction_sites.append(interaction_site)

    def __str__(self):
        s = 'Particle ' + str(self.id) + '\n'
        ILs_idxs = [[node.id for node in IL] for IL in self.ILs] 
        s += 'Interaction sites: ' + str(ILs_idxs) + '\n'
        s += 'Spike values: ' +str(self.spike_values) + '\n'
        s += 'Spike types: ' +str(self.spike_types) + '\n'
        return s

    def _krastev_spikes(self):
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
        node_spikes = list(map(lambda x: (len(self.rbn.attractor) - x) * -1 + x, summed_cycle))

        # We now have the spike values for each node in the IL. To get the spikes for the whole IL, we have 
        # to sum them up:
        # Get the ILs as lists of local_index values rather than objects.
        ILs_idxs = [[node.loc_idx for node in IL] for IL in self.ILs] 
        # Sum the node spikes referenced by the local_index values of the ILs.
        return list(map(sum, [[node_spikes[i] for i in lst] for lst in ILs_idxs]))

    def _watson_spikes(self):
        # Create lists such that:
        # frozen_true[node.loc_idx] == True iff node is frozen true
        # frozen_false[node.loc_idx] == True iff node is frozen false
        zipped_cycle = list(zip(*self.rbn.attractor))
        frozen_true = [all(values) for values in zipped_cycle]
        frozen_false = [any(values)==False for values in zipped_cycle]

        spike_values = []
        spike_types = []
        for IL in self.ILs:
            spike_value = 0
            # Set the spike type on the basis of the IL length.
            if len(IL) >= 10:
                spike_types.append(3)
            elif len(IL) >= 5:
                spike_types.append(2)
            else:
                spike_types.append(1)

            # Calculate the spike value for this IL based on whether/how nodes are frozen.
            for node in IL:
                if frozen_true[node.loc_idx]:
                    spike_value += 1
                elif frozen_false[node.loc_idx]:
                    spike_value += -1
            spike_values.append(spike_value)

        return spike_values, spike_types

    def _calculate_interaction_sites(self):
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

class Interation_Site:
    def __init__(self, nodes, spike_value, spike_type):
        self.available = True
        self._nodes = nodes
        self.spike_value = spike_value
        self.spike_type = spike_type
 
    def __str__(self):
        s = ""
        s += "Available: " + str(self.available) + "\n"
        s += "Nodes: " + str(list(map(str, self._nodes))) + "\n"
        s += "Spike value: " + str(self.spike_value) + "\n"
        s += "Spike type: " + str(self.spike_type) + "\n"
        return s

    def __len__(self):
        return len(self._nodes)

    def __getitem__(self, itemnumber):
        return self._nodes[itemnumber]


if __name__ == "__main__":
    print("particle.py invoked as script...")
    rbn.NodeSpace(1000)
    a = Particle.new(12, 2, 0)
    print(a)
