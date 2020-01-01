from os import linesep
from bitarray import bitarray
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

    def __init__(self):
        pass

    @classmethod
    def new(cls, n, k, offset, nodespace):
        obj = cls()
        obj.k = k
        obj.offset = offset # This doesn't make sense. The nodes should know their place in the nodespace, that's all.
        obj.n = n
        obj.nodespace = nodespace
        # Create new node objects.
        obj.nodes = [RBNNode(i, offset+i) for i in list(range(n))]
        # ! Note that the offset+1 method of calculating the index in the nodespace won't work
        # ! after swaps have taken place, so this can only be used in the constructor call.         
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
    def compose(cls, rbn1, rbn2):
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

        obj = cls()
        obj.k = rbn1.k

        # assert(rbn1.nodespace == rbn2.nodespace)
        # obj.nodespace = rbn1.nodespace

        obj.nodes = rbn1.nodes + rbn2.nodes
        for idx, node in enumerate(obj.nodes):
            node.loc_idx = idx
        obj.n = rbn1.n + rbn2.n
        obj.basin, obj.attractor = obj._calculate_cycle()
        return obj

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
    
    def _switch_edge(self, from_node1, to_node1, from_node2, to_node2):
        """
        In an RBN in which there is an edge from from_node1 to to_node1 and from from_node2 to
        to_node2, this method creates an edge from from_node1 to to_node2 and from from_node2 
        to to_node1, removing the original edges.

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
        to_node1.in_edges.insert(to_node1.index(from_node1), from_node2)
        to_node2.in_edges.insert(to_node2.index(from_node2), from_node1)
        from_node1.out_edges.insert(from_node1.out_edges(to_node1), to_node2)
        from_node2.out_edges.insert(from_node2.out_edges(to_node2), to_node1)

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

    def __init__(self, loc_idx, id, bool_func=0):
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
    NodeSpace(100)
    a = RBN.new(12, 2, 0, NodeSpace.getInstance())
    print(a.summarystring())

