import rbn
from collections import defaultdict
from operator import add
import reaction

class Atom:

    """
    Represents an Atom in the Frozen-Node Spiky RBN Chemistry.

    """
    id = 0
    
    def __init__(self):

        """
        The Atom class constructor.

        Parameters:

            id = An id number for the atom.
            n = The number of nodes in the Atom. This is passed directly to the RBN constructor.
            k = The indegree of each node. Passed to the RBN constructor.

        Returns: None

        """
        pass
  
    @classmethod
    def new(cls, n, k):
        obj = cls()
        obj.rbn = rbn.RBN.new(n, k)
        obj.id = Atom.id
        Atom.id += 1
        obj.atoms = 1
        obj.parent_composite = None
        obj.interaction_sites = Atom.calculate_interaction_sites(obj)
        return obj

     
    def __str__(self):
        s = 'Atom ' + str(self.id) + '\n'
        s += f'Atoms: {self.atoms}\n'
        site_idxs = [[node.id for node in site] for site in self.interaction_sites] 
        s += 'Interaction sites: ' + str(site_idxs) + '\n'
        s += 'Spike values: ' +str([interaction_site.spike_value() for interaction_site in self.interaction_sites]) + '\n'
        s += 'Spike types: ' +str([interaction_site.spike_type() for interaction_site in self.interaction_sites]) + '\n'
        return s

    @staticmethod
    def calculate_interaction_sites(parent):
        """
        Calculates and returns the interaction sites for an Atom parent. This should be done only once, when
        the Atom is created.

        Returns:
          
        ILs (List of lists of RBNNode objects): A list of the interaction lists for this Atom.
        """
        interaction_lists = []
        # Sort nodes in ascending order of influence.
        sortednodes = sorted([node for node in parent.rbn.nodes], key = lambda node: len(node.out_edges))
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
        
        return [Interaction_Site(interaction_list, parent) for interaction_list in interaction_lists]

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
    Represents a composite Atom by grouping together bonded atoms.
    """

    def __init__(self, atoms):
        self.atoms = atoms # A list of atoms
        for atom in self.atoms:
            atom.parent_composite = self
    
    def get_nodes(self):
        # Note that this works even if the composite has only one atom.
        return [node for nodes in [atom.rbn.nodes for atom in self.atoms] for node in nodes]

    def is_valid(self):
        """
        Checks that the atoms contained in the composites atoms list is the same as th list of atoms
        discovered by the traversal. If this returns false it implies that the composite comprises
        discrete two or more sets of atoms that are not bonded to each each other. It is easy to 
        create an invalid composite by creating two unbonded atoms and then passing them in a list
        to the Composite constructor.
        """
        atomlist = sorted([atom.id for atom in self.atoms])
        atomtraversal = sorted([atom.id for atom in self.traverse(self.atoms[0])]) 
        return atomlist == atomtraversal

    @classmethod
    # Given any atom, returns a list of atoms in the same parent composite or just the atom if it is not bound.
    def traverse(cls, atom):
        queue = []
        queue.append(atom)
        visited = []
        while queue:
            current = queue.pop(0)
            for site in current.interaction_sites:
                if site.bondedto is not None:
                    next = site.bondedto.parent_atom
                    if next not in visited:	
                        queue.append(next)
            visited.append(current)
        # visited now contains the atoms linked to the first atom
        return visited

    @classmethod
    def decompose(cls, atom):
        """
        
        Parameters: An atom of the composite to be decomposed where the breadth-first traversal begins at the atom supplied.

        Returns: A list of the composite particles into which the original composite has been decomposed.

        The return statements will be:

        decompose(composite 1 atom) + decomposite(composite 2 atom) when the composite belonging to the atom has failed stability somewhere or
        [composite], where the composite has not required decomposition or is a composite of length 1

        But won't this end up with a list of lists? We want it to be a single list. 
        
        """
        if len(atom.parent_composite.atoms) == 1:
            return [atom.parent_composite]
        else:
            queue = []
            queue.append(atom)
            visited = []
            while queue:
                current = queue.pop(0)
                for site in current.interaction_sites:
                    if site.bondedto is not None:
                        # test the bond
                        if reaction.is_stable(site, site.bondedto):
                            next = site.bondedto.parent_atom
                            if next not in visited:	
                                queue.append(next)
                        else:
                            print('Unstable bond found...')
                            site1, site2 = site, site.bondedto
                            reaction.break_bond(site1, site2)
                            # Where do we recalculate the attractor? Need to construct a composite here.
                            # The calculation of the attractor relies on the list of atoms in the parent composite.
                            # Perhaps break bond should return two composites? Why wouldn't it?
                            return Composite.decompose(Composite.traverse(site1.parent_atom)[0]) + \
                            Composite.decompose(Composite.traverse(site2.parent_atom)[0])
                visited.append(current)
            # TODO Will need to go through resulting list of atoms and create new composites or set parent_composite = None
            return [Composite(visited)]




class Interaction_Site:
    is_counter = 0

    def __init__(self, nodes, parent_atom):
        self.bondedto = None
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
        s += "Bonded to: " + str(self.bondedto) + "\n"
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
    a = Atom.new(12, 2)
    print(a)
    for site in a.interaction_sites:
        print(f'{site.spike_value()} {site.spike_type()}')
