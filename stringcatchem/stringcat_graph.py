import stringcatchem.stringcat_nodes as stringcat_nodes
import metachem.container as container
import metachem.control as control
import metachem.graph as graph

# ! Alert comment
# ? Question comment
# TODO To do comment
# * highlight comment

# * Remember: StringCatChem just combines strings and splits them wherever doubled letters occur. The end state must be a population in which all strings have the same letters at both ends, ie they all have the letter 'a' at the start and the finish. In such a population no further combination can take place.

# Declare nodes in graph

# containers
<<<<<<< HEAD
TTanks = container.ListTank() # * TTanks is described as a set of Tanks in the text
=======
TTanks = container.ListTank() # * This is named as TTanks because it's essentially a tank of tanks.
>>>>>>> 911c8853e66a3e1e5d152a7ffaea81829ff15704
Vtime = container.ListEnvironment()
Vtime.add(0)
TTank = stringcat_nodes.StringCatTank()
Scomposite = container.ListSample()
Vreactions = container.ListEnvironment()
Vreactions.add(0)
TLoad = container.ListTank()  # Empty tank used as place holder as no input tank needed

# control nodes
Sload = stringcat_nodes.StringCatLoadSampler(TLoad, TTanks, size=100, tanks=400)
Otime = control.ClockObserver(Vtime, Vtime)
Ssamplertank = control.SimpleSampler(TTanks, TTank) # * Very confusing naming but it's called this because it's a tank that we're choosing.
Ssamplerstring = control.SimpleSampler(TTank, Scomposite) # * Whereas here it's a string that we're choosing.
Ssamplerstringdecomp = control.SimpleSampler(TTank, Scomposite) # * This is exactly the same as Ssamplerstring. The "decomp" is confusing.
Ddecomp = stringcat_nodes.StringCatDecompDecision(2, Scomposite)
Aconcat = stringcat_nodes.StringCatConcatAction(Scomposite, Scomposite)
Asplit = stringcat_nodes.StringCatSplitAction(Scomposite, Scomposite)
Sreturn = control.BruteSampler(Scomposite, TTank)
Scommit = stringcat_nodes.StringCatCommitSampler(TTank, TTanks)
Oreaction = control.ClockObserver(Vreactions, Vreactions)
Dupdate = control.CounterDecision(2, Vreactions, 200)
Stransfers = stringcat_nodes.StringCatTransfersSampler(TTanks, TTanks, gridrows=20, gridcols=20, samplesize=1)

# Implement control edge set to get graph
edges = [[Sload, Otime], [Otime, Ssamplertank], [Ssamplertank, Ssamplerstring], [Ssamplerstring, Ddecomp],
         [Ddecomp, Ssamplerstringdecomp], [Ddecomp, Asplit],   # Ddecomp splits control
         [Ssamplerstringdecomp, Aconcat],
         [Asplit, Sreturn], [Aconcat, Sreturn],  # Control merges again at Sreturn
         [Sreturn, Scommit], [Scommit, Oreaction], [Oreaction, Dupdate],  # Dupdate splits control again
         [Dupdate, Stransfers], [Dupdate, Ssamplertank],  # Loop to the start of reaction
         [Stransfers, Otime]]  # Loop to the start to the generation

system = graph.Graph(edges, verbose=True)

# run graph
system.run_graph(Sload, transitionlimit=200000)
print (TTanks.read())
print (TTank.read())
print (Scomposite.read())
