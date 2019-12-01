from bitarray import bitarray
from os import linesep
import array
import random

"""
The offset + id method of looking up the node value won't work after nodes have been swapped.
"""

class NodeSpace:
    """
    A singleton class for storing all node values in any RBN-based AChem. The RBN is specified
    by an array of indexes into this node space.
    """
    __instance = None

    @staticmethod
    def getInstance():
        if NodeSpace.__instance == None:
            raise Exception("Node space has not been constructed yet.")
        else:
            return NodeSpace.__instance
    
    def __init__(self, size):
        if NodeSpace.__instance != None:
            raise Exception("Cannot instantiate more than one NodeSpace.")
        else:
            self.arr = bitarray(size)
            NodeSpace.__instance = self

    def __getitem__(self, index):
        return 1 if self.arr[index] else 0

    def __len__(self):
        return len(self.arr)

class RBN:
    
    offset = 0 # (static) The offset into the nodespace bitarray.

    def __init__(self, n=None, k=None, offset=None, nodespace=None, rbn1=None, rbn2=None):
        if (rbn1 and rbn2):
            self.composeRBN(rbn1, rbn2)
        elif (n and k):
            self.newRBN(n, k, offset, nodespace)

    def newRBN(self, n, k, offset, nodespace):
        self.n = n
        self.k = k
        self.offset = offset # This doesn't make sense. The nodes should know their place in the nodespace, that's all.
        self.nodespace = nodespace
        # Create new node objects.
        self.nodes = [RBNNode(i, offset+i) for i in list(range(n))]
        # ! Note that the offset+1 method of calculating the index in the nodespace won't work
        # ! after swaps have taken place, so this can only be used in the constructor call.         
        for node in self.nodes:
            # Make a list of all other nodes but this one and select k of them at random from which to add
            # to incoming edges.
            node.in_edges = random.sample(self._othernodes(node), k)

            # Create the corresponding outgoing edges in the other nodes.
            node.in_edges[0].out_edges.append(node)
            node.in_edges[1].out_edges.append(node)

            # Create a random boolean function for this node.
            node.bool_func = random.randint(0, 2**(2**k)-1) # Nb randint is inclusive on both sides
        
        self.basin, self.attractor = self._calculate_cycle()

    
    def composeRBN(self, rbn1, rbn2):
        # What will the new RBN look like?
        # We know the inward and outward edges of each.
        # There must be a concept of some edge swaps to do.
        # Concatenate the nodes lists and resequence then from 0 to n1+n2.
        # k stays the same; the final n becomes n1+n2.
        # We don't need to worry about decomposing them because the binary tree does that.
        # What happens with the boolean functions? Nothing. They don't depend on n. Phew.
        # Recalculate the basin and attractor.
        # The index into the nodespace (if implemented) shouldn't change.
        # Having node objects was an excellent idea!
        assert(rbn1.k == rbn2.k) # An assumption at the moment but could be relaxed later.
        self.k = rbn1.k

        assert(rbn1.nodespace == rbn2.nodespace)
        self.nodespace = rbn1.nodespace

        self.nodes = rbn1.nodes + rbn2.nodes
        for idx, node in enumerate(self.nodes):
            node.loc_idx = idx
        self.n = rbn1.n + rbn2.n
        self.basin, self.attractor = self._calculate_cycle()


    def getnextstate(self, node, currentstates=None):
        """
        Returns a list of the next state of node assuming the current states are as provide in the
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
    
    def _getcurrentstate(self, node):
        return self.nodespace[RBN.offset + node.loc_idx]

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

    def summarystring(self):
        """
        Returns a string representation of the RBN containing N, K, boolean functions, basin, attractor and
        a "next node state" check to show that the attractor is correct.
        
        Returns:

        summarystring: A string containing the summary of the RBN.
        """
        summarystring = ""
        summarystring += "N=" + str(self.n) + " K=" + str(self.k) + linesep
        summarystring += "Boolean functions:" + linesep
        # Ensure that the boolean functions are all 2**k bits long.
        formatstring = '#0' + str(2**self.k + 2) + 'b'
        summarystring += str([format(node.bool_func, formatstring) for node in self.nodes]) + linesep
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
    def __init__(self, loc_idx, glob_idx, bool_func=0):
        self.loc_idx = loc_idx
        self.glob_idx = glob_idx
        # Need to sort inward edges in ascending order for boolean function lookup to work properly.
        self.in_edges = []
        self.out_edges = []
        self.bool_func = bool_func
        pass

    def __str__(self):
        return str(self.loc_idx)

if __name__ == "__main__":
    print("Hello rbn...")
    NodeSpace(100)
    a = RBN(12, 2, 0, NodeSpace.getInstance())
    b = RBN(12, 2, 0, NodeSpace.getInstance())
    print(a.summarystring())
    print(b.summarystring())
    c = RBN(rbn1=a, rbn2=b)
    print(c.summarystring())

