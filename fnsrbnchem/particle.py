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
        # ! Try to get by without using this...
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
        
        # This is now in the reaction module.
        # Particle._do_edge_swaps(int_site1, int_site2)

        # We now have a new RBN object with edges between them as defined by the edge swaps specified by the interaction lists.
        # Need to calculate the new cycle and return it to the particle, where new spike details will be calculated.
        # ! Can't do this until the edge swaps have been done
        obj.rbn.basin, obj.rbn.attractor = rbn.RBN.get_cycle(obj.rbn.nodes)
               
        obj.id = "(" + str(particle1.id) + "." + str(particle2.id) + ")"
        obj.interaction_sites = particle1.interaction_sites + particle2.interaction_sites

        # obj._initialize()

        # First set the bonded interaction sites to unavailable for further bonding.
        int_site1.available = False
        int_site2.available = False

        return obj

    @classmethod
    def decompose(cls, particle, int_site1, int_site2):
        pass

    @classmethod
    def new(cls, n, k, id):
        obj = cls()
        obj.rbn = rbn.RBN.new(n, k)
        obj.id = id
        obj.atoms = 1
        obj.parent_composite = None
        obj.interaction_sites = []
        for site in obj._calculate_interaction_lists():
            obj.interaction_sites.append(Interaction_Site(site, parent_atom=obj))
        return obj
     
    def __str__(self):
        s = 'Particle ' + str(self.id) + '\n'
        s += f'Atoms: {self.atoms}\n'
        site_idxs = [[node.id for node in site] for site in self.interaction_sites] 
        s += 'Interaction sites: ' + str(site_idxs) + '\n'
        s += 'Spike values: ' +str([interaction_site.spike_value() for interaction_site in self.interaction_sites]) + '\n'
        s += 'Spike types: ' +str([interaction_site.spike_type() for interaction_site in self.interaction_sites]) + '\n'
        return s

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

class Composite:
    pass

    """
    Represents a composite particle by grouping together bonded atoms.
    """

    def __init__(self, atoms):
        self.atoms = atoms # A list of atoms
        for atom in self.atoms:
            atom.parent_composite = self
    
    def get_nodes(self):
        return [node for nodes in [atom.rbn.nodes for atom in self.atoms] for node in nodes]

    


class Interaction_Site:
    is_counter = 0

    def __init__(self, nodes, parent_atom):
        self.available = True
        self.nodes = nodes
        self.parent_atom = parent_atom

    def spike_type(self):
        if len(self.nodes) >= 10:
            return 3
        elif len(self.nodes) >= 5:
            return 2
        else:
            return 1        

    def spike_value(self):
        value = 0
        parent_nodes = self.parent_atom.parent_composite.get_nodes() if self.parent_atom.parent_composite is not None else self.parent_atom.rbn.nodes
        _, attractor = rbn.RBN.get_cycle(parent_nodes)
        for node in self.nodes:
            i = parent_nodes.index(node) # Get the index into the attractor
            sum_of_states = sum([states[i] for states in attractor])
            if sum_of_states == len(attractor):
                value += 1
            elif sum_of_states == 0:
                value -= 1
        return value

    def __str__(self):
        s = ""
        s += "Available: " + str(self.available) + "\n"
        s += "Nodes: " + str([node.id for node in self.nodes]) + "\n"
        s += "Spike value: " + str(self.spike_value()) + "\n"
        s += "Spike type: " + str(self.spike_type()) + "\n"
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
    for site in a.interaction_sites:
        print(f'{site.spike_value()} {site.spike_type()}')
