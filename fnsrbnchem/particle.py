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
        # Creates a new particle object.
        obj = cls()
        # This creates a new underlying RBN by performing the edge swaps according to the interaction site nodes list.
        # It also recalculates the attractor.
        obj.rbn = rbn.RBN.compose(rbn1=particle1.rbn, rbn2=particle2.rbn, nodes1=int_site1.nodes, nodes2=int_site2.nodes)
        obj.atoms = particle1.atoms + particle2.atoms # This won't work when we decompose.
        Particle._do_edge_swaps(int_site1, int_site2)

        # We now have a new RBN object with edges between them as defined by the edge swaps specified by the interaction lists.
        # Need to calculate the new cycle and return it to the particle, where new spike details will be calculated.
        obj.rbn.basin, obj.rbn.attractor = rbn._calculate_cycle(obj.rbn)
               
        obj.id = "(" + str(particle1.id) + "." + str(particle2.id) + ")"

        # But what to do we do with InteractionSites? They have to be inherited from the parent particles... ?
        # Concatenate them?


        obj.interaction_sites = particle1.interaction_sites + particle2.interaction_sites

        obj._initialize()

        # First set the bonded interaction sites to unavailable for further bonding.
        int_site1.available = False
        int_site2.available = False

        return obj

    @classmethod
    def new(cls, n, k, id):
        obj = cls()
        obj.rbn = rbn.RBN.new(n, k)
        obj.id = id
        obj.atoms = 1
        # Setting the interaction lists is only ever done once
        # obj.interaction_lists = obj._calculate_interaction_lists()
        obj.interaction_sites = []
        for site in obj._calculate_interaction_lists():
            interaction_site = Interation_Site(site, None, None)
            obj.interaction_sites.append(interaction_site)
        obj._initialize()
        return obj

    """
    @classmethod
    def fromRBN(cls, rbn, id):
        obj = cls()
        obj.rbn = rbn
        obj.id = id
        obj._initialize()
        return obj
    """

    def _initialize(self):
        """
        Initialization of spikes which is required both for new particles
        and composite particles.
        """
        
        # The attractor must have been calculated by this point.
        self.spike_values, self.spike_types = self._watson_spikes()
        assert(len(self.interaction_sites) == len(self.spike_values) == len(self.spike_types))

        # We have the InteractionSite objects but we need to add the spike values and types.
        for i, site in enumerate(self.interaction_sites):
            site.spike_value = self.spike_values[i]
            site.spike_type = self.spike_types[i]

    def __str__(self):
        s = 'Particle ' + str(self.id) + '\n'
        s += f'Atoms: {self.atoms}\n'
        site_idxs = [[node.id for node in site] for site in self.interaction_sites] 
        s += 'Interaction sites: ' + str(site_idxs) + '\n'
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
        site_idxs = [[node.loc_idx for node in site] for site in self.interaction_sites] 
        # Sum the node spikes referenced by the local_index values of the ILs.
        return list(map(sum, [[node_spikes[i] for i in lst] for lst in site_idxs]))

    def _watson_spikes(self):
        # Create lists such that:
        # frozen_true[node.loc_idx] == True iff node is frozen true
        # frozen_false[node.loc_idx] == True iff node is frozen false
        zipped_cycle = list(zip(*self.rbn.attractor))
        frozen_true = [all(values) for values in zipped_cycle]
        frozen_false = [any(values)==False for values in zipped_cycle]

        spike_values = []
        spike_types = []
        for site in self.interaction_sites:
            spike_value = 0
            # Set the spike type on the basis of the IL length.
            if len(site) >= 10:
                spike_types.append(3)
            elif len(site) >= 5:
                spike_types.append(2)
            else:
                spike_types.append(1)

            # Calculate the spike value for this IL based on whether/how nodes are frozen.
            for node in site:
                if frozen_true[node.loc_idx]:
                    spike_value += 1
                elif frozen_false[node.loc_idx]:
                    spike_value += -1
            spike_values.append(spike_value)

        return spike_values, spike_types

    def _calculate_interaction_lists(self):
        """
        Creates and calculates the interaction list for an Atom. This should be done only once, when
        the Atom is created.

        Returns:
          
        ILs (List of lists of RBNNode objects): A list of the interaction lists for this Atom.
        """
        interaction_lists = []
        # Sort nodes in ascending order of influence.
        sortednodes = sorted([node for node in self.rbn.nodes], key = lambda node: len(node.out_edges))
        while sortednodes:
            node = sortednodes.pop(0)
            site = []
            site.append(node)
            nextnodes = [nextnode for nextnode in node.in_edges if nextnode in sortednodes]
            while nextnodes:
                nextnode = nextnodes[0]
                sortednodes.remove(nextnode)
                site.append(nextnode)
                node = nextnode
                nextnodes = [nextnode for nextnode in node.in_edges if nextnode in sortednodes]
            interaction_lists.append(site)
        
        return interaction_lists

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

    @staticmethod
    def _do_edge_swaps(int_site1, int_site2):
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
        for index in [i for i in range(1, min(len(int_site1.nodes), len(int_site2.nodes)))]:
            # Make appropriate changes to the inward edges lists for each node.s
            Particle._switch_edges(int_site1.nodes[index], int_site1.nodes[index-1], int_site2.nodes[index], int_site2.nodes[index-1])
            if (index - 1) % 2 == 1:
                # Swap every other node in the interaction sites themselves. Have to operate on index - 1 node so as not to 
                # affect the edge switching.
                int_site1.nodes[index - 1], int_site2.nodes[index - 1] = int_site2.nodes[index - 1], int_site1.nodes[index - 1] 

        # Need to actually change the interaction sites here but can't change them in the 

    @staticmethod
    def _switch_edges(from_node1, to_node1, from_node2, to_node2):
        """
        In an RBN in which there is an edge from from_node1 to to_node1 and from from_node2 to
        to_node2, this method creates an edge from from_node1 to to_node2 and from from_node2 
        to to_node1, removing the original edges.
        It is important to remember that this method changes in_edges and out_edges, not any interaction site. 
        Interaction sites don't change unless they are part of a (temporary) bond. 
        Spike details should be recalculated after a called to switch_edges.
        Parameters:
        from_node1
        to_node1
        from_node2
        to_node2
        """
       
        # This method preserves k and so doesn't corrupt the RBN and can be called from outside the class.

        # Assert that the nodes are where they should be to start with.
        assert(from_node1 in to_node1.in_edges)
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

class Interation_Site:
    is_counter = 0

    def __init__(self, nodes, spike_value, spike_type):
        self.available = True
        self.nodes = nodes
        self.spike_value = spike_value
        self.spike_type = spike_type

    def __str__(self):
        s = ""
        s += "Available: " + str(self.available) + "\n"
        s += "Nodes: " + str([node.id for node in self.nodes]) + "\n"
        s += "Spike value: " + str(self.spike_value) + "\n"
        s += "Spike type: " + str(self.spike_type) + "\n"
        return s

    def __len__(self):
        return len(self.nodes)

    def __getitem__(self, itemnumber):
        return self.nodes[itemnumber]

if __name__ == "__main__":
    print("particle.py invoked as script...")
    # rbn.NodeSpace(1000)
    a = Particle.new(12, 2, 0)
    print(a)
