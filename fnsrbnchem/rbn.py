from os import linesep
import array
import random

class RBN:

    def __init__(self):
        pass

    @classmethod
    def new(cls, n, k):
        obj = cls()
        obj.k = k
        obj.n = n
        # Create new node objects.
        obj.nodes = [RBNNode(i) for i in list(range(n))]
        for node in obj.nodes:
            # Make a list of all other nodes but this one and select k of them at random from which to add
            # to incoming edges.
            node.in_edges = random.sample(obj._othernodes(node), k)

            # Create the corresponding outgoing edges in the other nodes.
            node.in_edges[0].out_edges.append(node)
            node.in_edges[1].out_edges.append(node)

            # Create a random boolean function for this node.
            node.bool_func = random.randint(0, 2**(2**k)-1) # Nb randint is inclusive on both sides
        
        obj.basin, obj.attractor = obj._calculate_cycle()
        return obj

    @classmethod
    def compose(cls, rbn1, rbn2, nodes1, nodes2):
        # What will the new RBN look like?
        # We know the inward and outward edges of each.
        # There must be a concept of some edge swaps to do.
        # Concatenate the nodes lists and resequence then from 0 to n1+n2.
        # k stays the same; the final n becomes n1+n2.

        assert(rbn1.k == rbn2.k) # An assumption at the moment but could be relaxed later.

        obj = cls()
        obj.k = rbn1.k

        obj.nodes = rbn1.nodes + rbn2.nodes
        for idx, node in enumerate(obj.nodes):
            node.loc_idx = idx
        obj.n = rbn1.n + rbn2.n
        # Edge swaps should happen here?? At the moment they are done in the particle.
        # It's true that RBNs shouldn't know about interaction lists, but it's fine for
        # an RBN to be able to create a new RBN from two separate RBNs and a set of 
        # edges to swap.

        # THIS IS WHERE BONDING TAKES PLACE
        # Find the nodes at the odd indices running up to the length of the shorter of the two interaction sites involved in the bond, then switch the edges, e.g in _switch_edges(2, 1, 4, 3):

        # 1 <- 2   becomes    1 <- 4
        # 3 <- 4              3 <- 2
        for index in [i for i in range(1, min(len(nodes1), len(nodes2)))]:
            obj._switch_edges(nodes1[index], nodes1[index-1], nodes2[index], nodes2[index-1])
        # obj._initialize()

        # We now have a new RBN object with edges between them as defined by the edge swaps specified by the interaction lists.
        # Need to calculate the new cycle and return it to the particle, where new spike details will be calculated.
        obj.basin, obj.attractor = obj._calculate_cycle()
        return obj

    def getnextstate(self, node, currentstates=None):
        """
        Returns a list of the next state of node assuming the current states are as provided in the
        currentstates list argument. 
        If the currentstates argument is not supplied then the actual current states are used.
        """
        if currentstates == None:
            currentstates = self._getcurrentstates()
        # boolean_function_lookup is used to look up the next state from the boolean function
        boolean_function_lookup = 0
        for (i, other_node) in enumerate(node.in_edges):
            # e.g. if k=2 and we have inputs from node 4 (false) and 7 (true) then the boolean
            # function lookup will end up being 10 (i.e. 2, the smaller node value determining the less
            # significant bit). Hence our next state will be bit number 2 from the boolean function.
            boolean_function_lookup = boolean_function_lookup + 2**i * currentstates[other_node.loc_idx]
        return 1 if node.bool_func >> boolean_function_lookup & 1 else 0
    
    def _switch_edges(self, from_node1, to_node1, from_node2, to_node2):
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

    def _getcurrentstate(self, node):
        # return self.nodespace[RBN.offset + node.loc_idx]
        # This is only used to get the initial state, which is set to zero. Needs to be cleaned up.
        # The node should carry its own state, presumably, with an initial value of zero.
        return 0

    def _getcurrentstates(self):
        return [self._getcurrentstate(node) for node in self.nodes]

    def _getnextstates(self, currenstates):
        return [self.getnextstate(node, currenstates) for node in self.nodes]
    
    def _othernodes(self, node):
        return self.nodes[:node.loc_idx] + self.nodes[node.loc_idx+1:]

    def _calculate_cycle(self):
        """
        Calculates the basin and the attractor of the RBN.
         
        Returns:
        basin (list), attractor_cycle (list): A tuple containing the basin and attractor of the RBN as lists of lists of N integers.
        
        """
        cycle = []
        states = self._getcurrentstates()
        while states not in cycle:
            cycle.append(states)
            states = self._getnextstates(states)
        cycle_start = cycle.index(states)
        basin = cycle[0:cycle_start]
        attractor = cycle[cycle_start:]
        return basin, attractor

    def attractorstring(self):
        """ Returns a string representation of the attractor of the RBN. """
        attractorstring = ""
        for count, state in enumerate(self.attractor):
            attractorstring += str(count) + " " + str(state) + linesep
        return attractorstring
    
    def basinstring(self):
        """ Returns a string representation of the basin of the RBN. """
        basinstring = ""
        for count, state in enumerate(self.basin):
            basinstring += str(count) + " " + str(state) + linesep
        return basinstring

    def summarystring(self, verbose=False):
        """
        Returns a string representation of the RBN containing N, K, boolean functions, basin, attractor and
        a "next node state" check to show that the attractor is correct.
        
        Returns:

        summarystring: A string containing the summary of the RBN.
        """
        summarystring = ""
        summarystring += "N=" + str(self.n) + " K=" + str(self.k) + linesep
        for node in self.nodes:
            summarystring += str(node) + linesep
        summarystring += linesep
        summarystring += "Basin:" + linesep
        summarystring += self.basinstring() + linesep
        summarystring += "Attractor:" + linesep
        summarystring += self.attractorstring() + linesep
        summarystring += "Check: next state would be..." + linesep
        summarystring += str(self._getnextstates(self.attractor[-1])) + linesep
        return summarystring

class RBNNode:
    """
    Represents a node within an RBN.
    """
    glob_idx = 0

    def __init__(self, loc_idx, bool_func=0):
        self.loc_idx = loc_idx
        self.id = RBNNode.glob_idx
        RBNNode.glob_idx += 1
        # Need to sort inward edges in ascending order for boolean function lookup to work properly.
        self.in_edges = []
        self.out_edges = []
        self.bool_func = bool_func
        pass

    def __str__(self):
        formatstring = '#0' + str(2**len(self.in_edges) + 2) + 'b'
        return 'id:{:2}'.format(self.id, a=2) + ' bf: ' + format(self.bool_func, formatstring) + ' in:' + str([node.id for node in self.in_edges]) + ' out:' + str([node.id for node in self.out_edges])

if __name__ == "__main__":
    print("Hello rbn...")
    # NodeSpace(100)
    # a = RBN.new(12, 2, NodeSpace.getInstance())
    a = RBN.new(12, 2)
    print(a.summarystring())

