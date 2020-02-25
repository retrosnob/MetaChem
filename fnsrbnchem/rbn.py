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
        
        obj.basin, obj.attractor = _calculate_cycle(obj)
        return obj

    @classmethod
    def compose(cls, rbn1, rbn2, nodes1, nodes2):

        assert(rbn1.k == rbn2.k) # An assumption at the moment but could be relaxed later.

        obj = cls()
        obj.k = rbn1.k

        obj.nodes = rbn1.nodes + rbn2.nodes
        for idx, node in enumerate(obj.nodes):
            node.loc_idx = idx
        obj.n = rbn1.n + rbn2.n

        return obj

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
        summarystring += str([_nextstate(node, change_state=False) for node in self.nodes]) + linesep
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

# Module methods

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

def _calculate_cycle(rbn):
    """
    Calculates the basin and the attractor of the RBN.
        
    Returns:
    basin (list), attractor_cycle (list): A tuple containing the basin and attractor of the RBN as lists of lists of N integers.
    
    """
    cycle = []
    states = [0] * len(rbn.nodes)
    while states not in cycle:
        cycle.append(states)
        states = [_nextstate(node, change_state=True) for node in rbn.nodes]
    cycle_start = cycle.index(states)
    basin = cycle[0:cycle_start]
    attractor = cycle[cycle_start:]
    return basin, attractor

if __name__ == "__main__":
    print("Hello rbn...")
    # NodeSpace(100)
    # a = RBN.new(12, 2, NodeSpace.getInstance())
    a = RBN.new(12, 2)
    print(a.summarystring())

