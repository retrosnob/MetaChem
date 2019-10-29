import networkx as nx
import os
import networkx.classes.function
from random import choice, randint, sample

class RBN:

    """ The RBN that underlies an Atom in the Chemistry. Any attributes or methods that are specific to RBNs are
     encapsulated here. Each Atom object has an instance of the RBN class.

     """

<<<<<<< HEAD
    def __init__(self, n, k, default=None):
=======
    def __init__(self, n, k):
>>>>>>> 911c8853e66a3e1e5d152a7ffaea81829ff15704

        """
        The RBN constructor. Creates an NK RBN and finds its attractor cycle.

        Parameters:

        n: The number of nodes in the underlying digraph.
        
        k: The indegree of each node.
<<<<<<< HEAD

        default: A default starting node value for all nodes. Can be 1, 0 or None (omitted). If default is None then
        initial node values are random.
=======
>>>>>>> 911c8853e66a3e1e5d152a7ffaea81829ff15704
        """

        # print("Initializing RBN object with n = " + str(n) + ", k = " + str(k))

        self.n = n
        self.k = k

        # make a digraph with n nodes
        nodes = list(range(n))
        nxgraph = nx.DiGraph()
        nxgraph.add_nodes_from(nodes)

        # Setting attributes of nodes involves creating a list of values that is
        # the same length as the nodelist. Here we create lists for the t/f states
        # and the boolean functions of each node. These are populated in the loop
        # and then applied in one statement after the loop.
        node_states = []
        node_boolean_functions = []

        for node in nodes:
            # Update the lists of states and boolean functions
<<<<<<< HEAD
            node_states.append(default if default != None else choice((1, 0)))
=======
            node_states.append(choice((1, 0)))
>>>>>>> 911c8853e66a3e1e5d152a7ffaea81829ff15704
            boolean_function = randint(0, 2**(2**k)-1) # randint is inclusive on both sides
            node_boolean_functions.append(boolean_function)


            # make a list of all other nodes but this one
            other_nodes = nodes[:node] + nodes[node+1:]
            # and select two of them at random
            random_other_k = sample(other_nodes, k)

            # now add an edge from each of the other two nodes
            for i in range(k):
                nxgraph.add_edge(random_other_k[i], node)

        self._initial_node_states = node_states
        self._node_boolean_functions = node_boolean_functions
        self.nxgraph = nxgraph
        self.basin, self.attractor_cycle = self._calculate_cycle(node_states)

    def _get_next_states(self, node_states):
        """
        Calculates the next nodes states from the current node states.

        Parameters:
        node_states (list): A list of N integers describing the current node states of the NK RBN.

        Returns:
        next_node_states (list): A list of N integers describing the current node states of the NK RBN.
        """

        node_next_states = []
        for node in self.nxgraph.nodes:
            # graph.in_edges(node) is a list of k edges expressed as (from, to) tuples.
            in_edges = list(self.nxgraph.in_edges(node))
            # sort the k tuples by the from value
            in_edges.sort(key=lambda tup: tup[0])

            # boolean_function_lookup is used to look up the next state from the boolean function
            boolean_function_lookup = 0
            for (i, edge) in enumerate(in_edges):
                # e.g. if k=2 and we have inputs from node 4 (false) and 7 (true) then the boolean
                # function lookup will end up being 10 (i.e. 2, the smaller node value determining the less
                # significant bit). Hence our next state will be bit number 2 from the boolean function.
                boolean_function_lookup = boolean_function_lookup + 2**i * node_states[edge[0]]

            node_next_states.append(1 if (self._node_boolean_functions[node] >> boolean_function_lookup & 1) else 0)
            # print("ND: " + str(node) + ", BF " + format(self.node_boolean_functions[node], '0' + str(2 ** self.k) + 'b') +
            #      ", ED " + str(in_edges) + ", LK " + str(boolean_function_lookup) + ", NX " + str(node_next_states[node]))
        return node_next_states


    def _calculate_cycle(self, initial_node_states):
        """
        Calculates the basin and the attractor of the RBN.

        Parameters: 
        initial_node_states (list): The initial random node states of the NK RBN as a list of N integers. 
         
        Returns:
        basin (list), attractor_cycle (list): A tuple containing the basin and attractor of the RBN as lists of lists of N integers.
        
        """
        cycle = []
        cycle.append(initial_node_states)
        node_states = self._get_next_states(initial_node_states)
        while not node_states in cycle:
            cycle.append(node_states)
            node_states = self._get_next_states(node_states)
        cycle_start = cycle.index(node_states)
        basin = cycle[0:cycle_start]
        attractor_cycle = cycle[cycle_start:]
        return basin, attractor_cycle

    def attractorstring(self):
        """ Returns a string representation of the attractor of the RBN. """
        attractorstring = ""
        for count, state in enumerate(the_rbn.attractor_cycle):
            attractorstring += str(count) + " " + str(state) + os.linesep
        return attractorstring
    
    def basinstring(self):
        """ Returns a string representation of the basin of the RBN. """
        basinstring = ""
        for count, state in enumerate(the_rbn.basin):
            basinstring += str(count) + " " + str(state) + os.linesep
        return basinstring

    def summarystring(self):
        """
        Returns a string representation of the RBN containing N, K, boolean functions, basin, attractor and
        a "next node state" check to show that the attractor is correct.
        
        Returns:

        summarystring: A string containing the summary of the RBN.
        """
        summarystring = ""
        summarystring += "N=" + str(self.n) + " K=" + str(self.k) + os.linesep
        summarystring += "Boolean functions:" + os.linesep
        # Ensure that the boolean functions are all 2**k bits long.
        formatstring = '#0' + str(2**self.k + 2) + 'b'
        summarystring += str([format(i, formatstring) for i in self._node_boolean_functions]) + os.linesep
        summarystring += "Basin:" + os.linesep
        summarystring += self.basinstring() + os.linesep
        summarystring += "Attractor:" + os.linesep
        summarystring += self.attractorstring() + os.linesep
        summarystring += "Check: next state would be..." + os.linesep
        summarystring += str(the_rbn._get_next_states(self.attractor_cycle[-1])) + os.linesep
        return summarystring

if __name__ == "__main__":
    print("rbn.py running as script...")
    the_rbn = RBN(12, 3)
    print(the_rbn.summarystring())
