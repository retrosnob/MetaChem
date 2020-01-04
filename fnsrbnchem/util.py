import pickle
import rbn
import particle
import datetime

FILEPATH = 'fnsrbnchem/pickle/'

"""
We need to be able to load known particles to test how they are bonding.

Create a particle. Pickle it. Save it's string representation.
We can then examine it and other known particles to see how they bond.

Use time.asctime() as the filename.
"""

def create_pickled_particle(the_particle=None, filename=None):
    """
    Creates a new RBN and pickles it with filename.
    """
    # rbn.NodeSpace(100)
    if the_particle is None:
        the_rbn = rbn.RBN.new(12, 2, 0, rbn.NodeSpace.getInstance())
        the_particle = particle.Particle.fromRBN(the_rbn, 0)
    else:
        the_rbn = the_particle.rbn

    print("Going to pickle: ")
    print(str(the_rbn.summarystring()))
    print(str(the_particle))
    if filename is None:
        filename = str(datetime.datetime.now()).replace(' ', '@')
    with open(FILEPATH + filename + '.pickle', 'wb') as f:
        pickle.dump(the_particle, f, pickle.HIGHEST_PROTOCOL)
    with open(FILEPATH + filename + '.txt', 'w') as f:
        f.write(the_rbn.summarystring())
        f.write(str(the_particle))

def load_pickled_particle(filename):
    with open(FILEPATH + filename + '.pickle', 'rb') as f:
        the_particle = pickle.load(f)
    print('Unpickled:')
    print(the_particle.rbn.summarystring())
    print(the_particle)
    return the_particle


if __name__ == "__main__":
    create_pickled_Particle(filename='00')
    load_pickled_Particle(filename='00')