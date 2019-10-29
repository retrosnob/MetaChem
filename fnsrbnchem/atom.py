import rbn
from collections import defaultdict
from operator import add

class Atom:

    """ Represents an Atom in the Frozen-Node Spiky RBN Chemistry.

    To do:

    Should bonding be a method in the Atom or in the Chemistry itself?

    """

    def __init__(self, id, n, k):

        """ The Atom class constructor.

        Parameters:

            id = An id number for the atom.
            n = The number of nodes in the Atom. This is passed directly to the RBN constructor.
            k = The indegree of each node. Passed to the RBN constructor.

        Returns: None

        """
        self.id = id
        self.rbn = rbn.RBN(n, k)
        self.interaction_lists = self.calculate_interaction_lists()
        self.interaction_list_spikes = self.calculate_spikes()

    def __str__(self):
        s = 'Atom ' + str(self.id) + '\n'
        s += str(self.interaction_lists) + '\n'
        s += str(self.interaction_list_spikes) + '\n'
        return s

    def calculate_spikes(self):
        """
        Each interaction list is given a "spike" value, calculated over the length of the attractor cycle
        This is done by adding 1 for each true state and -1 for each false state for each node over the length of 
        the attractor cycle.

        Returns:

        interaction_list_spikes (list of integers): A list of the spike values for each interaction list.
        """
        # for node_list in self.rbn.attractor_cycle:
        #     print(node_list)
        zipped_cycle = list(zip(*self.rbn.attractor_cycle))
        spikes = [sum(tpl) for tpl in zipped_cycle]
        # TO DO: Make each 0 contribute -1
        # This can be calculated since the number of zeroes in the spike is equal to
        #  attractor cycle length - current node spike value
        # so if we multiply this value by -1 and add it to the current node spike value then we're done
        node_spikes = list(map(lambda x: (len(self.rbn.attractor_cycle) - x) * -1 + x, spikes))
        # self.spikes now contains the spike value for each node
        # In order to get the spikes for each IL, we have to sum up the spikes for each node on that IL
        total = 0
        interaction_list_spikes = []
        for interaction_list in self.interaction_lists:
            for node in interaction_list:
                total += node_spikes[node]
            interaction_list_spikes.append(total)
            total = 0
        # self.interaction_list_spikes = [reduce(add, self.node_spikes[node]) for interaction_list in self.interaction_lists for node in interaction_list]
        return interaction_list_spikes            

    def calculate_interaction_lists(self):
        """ Creates and calculates the interaction list for an Atom.

          To do: There are a number of things that will need to be done.

          First of all this method was written with a
          view to creating interaction lists just after an Atom has been instantiated. But in fact interaction lists
          have to be recalculated every time a bond is formed or broken so code that is responsible for initialization
          of the interaction lists should be split from code that is responsible for updating interaction lists.
          The code for updating interaction lists will have to take account of existing bonds.

          Secondly the code needs to be optimized. It is certainly not very efficient as it is.

          Returns:
          
          interaction_lists (List of lists of int): A list of the interaction lists for this Atom.

          """

        # Set in_edges_dict
        edges = list(self.rbn.nxgraph.in_edges)

        # Set up incoming and outgoing edges dictionaries where k = node, v = list of edges (either in or out)
        # NOTE: Outgoing edges dict doesn't include nodes which have no outgoing edges
        incoming_edges_by_node = {}
        outgoing_edges_by_node = {}
        for k, v in edges:
            outgoing_edges_by_node.setdefault(k, []).append(v)
        for k, v in edges:
            incoming_edges_by_node.setdefault(v, []).append(k)
        # print("Node out edges dict: " + str(outgoing_edges_by_node))
        # print("Node in edges dict: " + str(incoming_edges_by_node))

        # Sort nodes descending by influence (number of outgoing edges). They will be processed in reverse order by
        # popping them from the end and so actually in ascending order of influence as specified by Krastev.
        nodes_desc_by_inf = sorted(outgoing_edges_by_node, key=lambda j: len(outgoing_edges_by_node[j]), reverse=True)
        # print("Nodes descending on influence: " + str(nodes_desc_by_inf))
        total_nodes_to_assign = len(nodes_desc_by_inf)
        nodes_with_no_edges = [node for node in range(self.rbn.n) if node not in nodes_desc_by_inf]
        # print("Nodes with no edges: " + str(nodes_with_no_edges))

        # Create interaction lists
        # print("Creating interaction lists...")
        interaction_list = []
        interaction_lists = []
        count_nodes_assigned = 0
        node = -1  # Flag value meaning no current node is being processed

        while count_nodes_assigned < total_nodes_to_assign:
            if node == -1:
                # There is no current node so we have to get one from the list by influence
                node = nodes_desc_by_inf.pop()  # popping gets the least influential in the list
                # print("Adding first node of new IL: " + str(node))
                interaction_list.append(node)
                count_nodes_assigned += 1

            # Get nodes that have inputs to the current node and are still in the list by influence
            in_nodes = incoming_edges_by_node.get(node)
            for in_node in in_nodes:
                if in_node in nodes_desc_by_inf:  # while
                    # Add node to interaction list
                    # print(str(in_node) + " found in incoming nodes to " + str(node) + " " + str(nodes_desc_by_inf))
                    # print("Adding next node of current IL: " + str(in_node))
                    interaction_list.append(in_node)
                    count_nodes_assigned += 1
                    # Node was already popped from original list in order of influence so it doesn't need to
                    # be removed here. This must be wrong. It is. This is the second node found so we need to
                    # remove it from nodes_desc_by_inf.
                    nodes_desc_by_inf.remove(in_node)
                    node = in_node
                    break
            else:
                # This only executes if none of the incoming nodes were left in the nodes_desc_by_inf
                # print("New interaction list created: " + str(interaction_list))
                interaction_lists.append(interaction_list)
                interaction_list = []
                node = -1
        if len(interaction_list) > 0:
            interaction_lists.append(interaction_list)

        for remaining_node in nodes_with_no_edges:
            x = [remaining_node]
            interaction_lists.append(x)
        # print("ILs: " + str(interaction_lists))

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


if __name__ == "__main__":
    print("Atom.py invoked as script...")
    a = Atom(1, 12, 2)
    a.calculate_interaction_lists()
    a.calculate_spikes()
    print(a.interaction_lists)
    print(a.interaction_list_spikes)

