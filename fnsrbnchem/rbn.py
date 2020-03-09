from os import linesep
import array
import random

class RBN:
    N = 12
    K = 2

    def __init__(self):
        pass

    @classmethod
    def new(cls, n=N, k=K):
        obj = cls()
        obj.k = k
        obj.n = n
        # Create new node objects.
        obj.nodes = [RBNNode(i) for i in list(range(n))]
        for node in obj.nodes:
            # Make a list of all other nodes but this one and select k of them at random from which to add
            # to incoming edges.
            node.in_edges = random.sample([n for n in obj.nodes if n is not node], k)

            # Create the corresponding outgoing edges in the other nodes.
            node.in_edges[0].out_edges.append(node)
            node.in_edges[1].out_edges.append(node)

            # Create a random boolean function for this node.
            node.bool_func = random.randint(0, 2**(2**k)-1) # Nb randint is inclusive on both sides
        return obj

    @classmethod
    def compose(cls, rbn1, rbn2, nodes1, nodes2):
        # Replace this method with a get attractor that just returns the attractor for a bunch of nodes
        # supplied by a composite Atom. There is no need to have a composed RBN object. 

        assert(rbn1.k == rbn2.k) # An assumption at the moment but could be relaxed later.

        obj = cls()
        obj.k = rbn1.k

        obj.nodes = rbn1.nodes + rbn2.nodes
        # This loop shouldn't be necessary any more
        for idx, node in enumerate(obj.nodes):
            node.loc_idx = idx
        obj.n = rbn1.n + rbn2.n

        return obj

    @staticmethod
    def _nextstate(node, change_state=False):
        """
        Returns the next state of a node, depending on the current state of the rbn.
        
        Parameters:

        node: The node whose next state is to be calculated.

        change_state: If True then function updates the state of node. If False the new is returned but the node's current state is not changed.
        """
        boolean_function_lookup = 0
        for (i, other_node) in enumerate(node.in_edges):
            # e.g. if k=2 and we have inputs from node 4 (false) and 7 (true) then the boolean
            # function lookup will end up being 10 (i.e. 2, the smaller node value determining the less
            # significant bit). Hence our next state will be bit number 2 from the boolean function.
            boolean_function_lookup = boolean_function_lookup + 2**i * other_node.state
        if change_state:
            node.state = 1 if node.bool_func >> boolean_function_lookup & 1 else 0
        return node.state

    @staticmethod
    def get_cycle(nodes):
        """
        Calculates the basin and the attractor of the RBN.

        Parameters:

        nodes (list): The nodes with which to calculate that attractor cycle.

        Returns:

        basin (list), attractor (list): A tuple containing the basin and attractor of the RBN as lists of lists of len(nodes) integers.
        
        """
        # for node in nodes:
        #    node.state = 0
        cycle = []
        states = [0] * len(nodes)
        while states not in cycle:
            cycle.append(states)
            states = [RBN._nextstate(node, change_state=True) for node in nodes]
        cycle_start = cycle.index(states)
        basin = cycle[0:cycle_start]
        attractor = cycle[cycle_start:]
        return basin, attractor

    def attractorstring(self):
        """ Returns a string representation of the attractor of the RBN. """
        attractorstring = ""
        _, attractor = RBN.get_cycle(self.nodes)
        for count, state in enumerate(attractor):
            attractorstring += str(count) + " " + str(state) + linesep
        return attractorstring
    
    def basinstring(self):
        """ Returns a string representation of the basin of the RBN. """
        basinstring = ""
        basin, _ = RBN.get_cycle(self.nodes)
        for count, state in enumerate(basin):
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
        summarystring += str([RBN._nextstate(node, change_state=False) for node in self.nodes]) + linesep
        return summarystring

class RBNNode:
    """
    Represents a node within an RBN.
    """
    glob_idx = 0

    def __init__(self, loc_idx, bool_func=0):
        self.loc_idx = loc_idx
        self.id = f'{RBNNode.glob_idx // RBN.N}.{RBNNode.glob_idx % RBN.N}'
        RBNNode.glob_idx += 1
        # Need to sort inward edges in ascending order for boolean function lookup to work properly.
        self.in_edges = []
        self.out_edges = []
        self.bool_func = bool_func
        self.state = 0 # All nodes assumed to start at zero
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

