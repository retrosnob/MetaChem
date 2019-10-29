import string
import random
import metachem.node as node
import metachem.container as container
import metachem.control as control

# List of needed control nodes:
#   s:load - input set of strings
#   o:time - increments time between generations - (ClockObserver)
#   s:samplertank - Chooses a tank to work over that hasn't been chosen in the current generation - (SimpleSampler)
#   s:samplerstring - Choose a string to work with in this reaction - (SimpleSampler)
#   d:decomp - check if string is suitable for decomp if so decomp otherwise link
#   a:concat - join all strings in sample together
#   a:split - seperate string at double letter
#   s:return - return sample to tank - (BruteSampler)
#   s:commit - return tank to set of tanks - (BruteSampler)
#   o:reaction - record tank used and that a reaction has occurred - (ClockObserver)
#   d:update - check if tanks have all been dealt with and enough reactions have been done - (CounterDecision)
#   s:transfers - randomly select a number of pairs of tanks and a number of strings to swap between them

# List of container nodes:
#   T:tanks - set of all tanks in system
#   T:strings - set of initial strings to load
#   V:time - integer generation clock
#   T:tank - current tank of interest
#   S:composite - strings involved in reaction
#   V:reactions - tanks that have been used and total number of reactions


class StringCatLoadSampler(node.Sampler):

    def __init__(self, containersin, containersout, readcontainers=None, size=1, tanks=1):
        # * This is called with size=100, tanks=400
        super(StringCatLoadSampler, self).__init__(containersin, containersout, readcontainers)
        self.size = size
        self.tanks = tanks
        self.sample = None
        pass

    def read(self):
        # * The sample is a list of 400 random lists of one-character strings, each of length 100
        self.sample = [[random.choice(string.ascii_uppercase) for _ in range(0, self.size)] for _ in range(0, self.tanks)]
        pass

    def pull(self):
        self.containersin.remove(self.sample)
        pass

    def push(self):
        self.containersout.add(self.sample)
        pass


class StringCatDecompDecision(node.Decision):

    def __init__(self, options=2, readcontainers=None):
        super(StringCatDecompDecision, self).__init__(options, readcontainers)
        self.samplestring = None
        pass

    def read(self):
        self.samplestring = self.readcontainers.read()

    def process(self):
        doubleindex = [i for i in range(0, len(self.samplestring) - 1)
                       if self.samplestring[i] == self.samplestring[i + 1]]
        if doubleindex:
            return 1
        else:
            return 0


class StringCatConcatAction(node.Action):

    def __init__(self, writesample, readsample, readcontainers=None):
        super(StringCatConcatAction, self).__init__(writesample, readsample, readcontainers)
        self.sample = []
        pass

    def read(self):
        self.sample = self.readsample.read()

    def pull(self):
        # print "pre-remove on pull"
        # print self.sample
        # print self.readsample.read()
        self.readsample.remove(self.sample)

    def check(self):
        return super(StringCatConcatAction, self).check()

    def process(self):
        # print self.sample
        self.sample = "".join(self.sample)
        # print self.sample
        pass

    def push(self):
        self.writesample.add(self.sample)


class StringCatSplitAction(node.Action):

    def __init__(self, writesample, readsample, readcontainers=None):
        super(StringCatSplitAction, self).__init__(writesample, readsample, readcontainers)
        self.sample = None
        pass

    def read(self):
        self.sample = self.readsample.read()

    def pull(self):
        self.readsample.remove(self.sample)

    def check(self):
        return super(StringCatSplitAction, self).check()

    def process(self):
        # * Find all insances of doubled letters
        doubleindices = [i for i in range(0, len(self.sample) - 1) if self.sample[i] == self.sample[i+1]]
        # * Pick a random one
        index = random.choice(doubleindices)
        # ! This next line is supposed to split the sample string at the doubled letter and add the two halves
        # ! back into the writesample, but it contains an error. The rightmost expression should not be 
        # ! self.sample[index:0], it should be self.sample[index+1:] 
        # Removed:
        # self.sample = [self.sample[0:index], self.sample[index:0]]
        # Added:
        self.sample = [self.sample[0:index], self.sample[index+1:]]
        
        # ? Why are the strings written back to the writesample here?
        self.writesample.add(self.sample)
        pass

    def push(self):
        # ? Note that the process method also writes these back to the writesample.
        self.writesample.add(self.sample)


class StringCatTransfersSampler(node.Sampler):

    # ? The grid is set up to be 20x20 so that indices can be used to reference the 400 tanks.

    def __init__(self, containersin, containersout, readcontainers=None,
                 gridrows=1, gridcols=1, samplesize=1):

        super(StringCatTransfersSampler, self).__init__(containersin, containersout, readcontainers)
        self.gridrows = gridrows # 20
        self.gridcols = gridcols # 20
        self.pairs = []
        self.samplesize = samplesize # 1
        pass

    def read(self):
        super(StringCatTransfersSampler, self).read()

    def pull(self):
        sample = self.containersin.read()[0]
<<<<<<< HEAD
=======
        # Set up indices so that cells in the grid can be referenced.
>>>>>>> 911c8853e66a3e1e5d152a7ffaea81829ff15704
        indices = list(range(0, self.gridcols*self.gridrows))
        pairs = []


        # The aim of this loop is to select cells from the grid in pairs of neighboring cells and add them
        # to the pairs list as lists of two indices. Any cells that end up without neighbours are ignored.
        for cell in indices:
            # Set neighbouring cells of this cell. Note that the grid is toroidal. 
             neighbours = [cell+1, cell-1, cell+self.gridrows, cell-self.gridrows]

            for checkcell in neighbours:
                # Remove cell from neighbours if it has already been removed from indices.
                if checkcell not in indices:
                    neighbours.remove(checkcell)
            # If the current cell has no neighbours then remove it from indices and do the next cell.
            if not neighbours:
                indices.remove(cell)
                continue
            # These last lines are only done if the current cell still has neighbours.
            # Surely an else would have been clearer than a continue? 
            othercell = random.choice(neighbours)
<<<<<<< HEAD

            if othercell not in indices:
                raise Exception("Error in pull() method: Trying to remove " + str(othercell) + " from " + str(indices))

=======
            # Remove the cell and a random one of its neighbours from indices.
>>>>>>> 911c8853e66a3e1e5d152a7ffaea81829ff15704
            indices.remove(othercell)
            indices.remove(cell)
            # And append them as a list to pairs
            pairs.append([cell, othercell])


        # sample: TTanks, the list of lists of strings
        # samplesize: 1
        for pair in pairs:
 
            try:
                # Pick a random samplesize = 1 strings from the Nth tank in TTanks where
                # N is the second value in the current pair. 
                sample0 = random.sample(sample[pair[1]], self.samplesize)
            except ValueError:
                # If that raises a ValueError (why should it?) then select the Nth tank as a whole??
                sample0 = sample[pair[1]]
            self.containersin.read()[0][pair[1]].remove(sample0)
            try:
                sample1 = random.sample(sample[pair[0]], self.samplesize)
            except ValueError:
                sample1 = sample[pair[0]]
            self.containersin.read()[0][pair[0]].remove(sample1)
            self.sample.append([sample0, sample1])
        self.pairs = pairs
        pass

    def push(self):
        for i in range(0, len(self.pairs)):
            self.containersout[0, self.pairs[i][0]].append(self.sample[i][0])
            self.containersout[0, self.pairs[i][1]].append(self.sample[i][1])
        pass


class StringCatTank(container.ListTank):

    def __init__(self):
        super(StringCatTank, self).__init__()
        pass

    def read(self):
        return self.list[:]

    def add(self, particles=None):
        if self.list:
            if isinstance(particles[0], list):
                self.list = self.list + particles[0]
            else:
                self.list = self.list + particles
        elif isinstance(particles, list) and isinstance(particles[0], list):
            # print "nested"
            self.list = particles[0]
        elif isinstance(particles, list):
            self.list = particles
        else:
            self.list = [particles]
        # print "List length"
        # print len(self.list)
        return self.list

    def remove(self, particles=None):
        # print self.list
        # print particles
        try:
            [self.list.remove(part) for part in particles]
        except ValueError:
            self.list.remove(particles)
        return self.list


class StringCatCommitSampler(control.BruteSampler):

    def __init__(self, containersin, containersout, readcontainers=None):
        super(StringCatCommitSampler, self).__init__(containersin, containersout, readcontainers)
        pass

    def read(self):
        super(StringCatCommitSampler, self).read()
        pass

    def pull(self):
        # print "pre-pull"
        super(StringCatCommitSampler, self).pull()
        # print self.sample
        pass

    def push(self):
        self.containersout.add([self.sample])
        pass
