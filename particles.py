from copy import deepcopy
import random
import matplotlib.pyplot as plt
import scipy.constants as consts
from math import pi, sin, cos, sqrt, fabs, exp
random.seed()

# constants and stuff
hamaker = 3 * consts.physical_constants["electron volt"][0] # gold-gold, as in https://web.science.uu.nl/scm/Articles/2010/SI%2010.1021-nl102705p.pdf
dielectric = 77.46 # dielectric constant of water at 25 degrees centigrade, https://www.sciencedirect.com/science/article/pii/B9780123850133000070
ionicStrength = 0.001 # ionic strength of water
q = consts.physical_constants["elementary charge"][0]
debyeLength = 1/sqrt((consts.k * 298.15 * consts.epsilon_0 * dielectric) / (2*q*q*consts.Avogadro*ionicStrength))
sigma = -2e-3 # surface charge density of gold, https://pubs.acs.org/doi/10.1021/acs.jpcc.5b00568
radius = 50e-9 # radius of gold nanoparticles, https://pubs.acs.org/doi/10.1021/acs.jpcc.5b00568
mass = (19.32e-3) * (4/3 * pi * radius*radius*radius) # calculated from true density https://www.americanelements.com/gold-nanoparticles-7440-57-5

def freeEnergyDL(h):
    return (2*sigma*sigma)/(consts.epsilon_0*dielectric*debyeLength)*exp(-debyeLength*h)

def freeEnergyVDW(h):
    return -(hamaker)/(12*pi*h*h)

def freeEnergy(h):
    return freeEnergyDL(h) + freeEnergyVDW(h)

def force(h):
    return pi*radius*freeEnergy(h)


class Particle:
    def __init__(self, x=None, y=None):
        self.x = x if x else random.random() * 1e-6
        self.y = y if y else random.random() * 1e-6

    def distance(self, p2):
        dx = self.x - p2.x
        dy = self.y - p2.y
        return sqrt(dx*dx + dy*dy)


def forceStep(space):
    timeInterval = 3e-11 # time step; the force is assumed to be linear at steps of this length
    nextSpace = []
    for p1 in space:
        forceSumX, forceSumY = 0, 0
        for p2 in space:
            if p1 == p2:
                continue
            dist = p1.distance(p2)
            (dirx,diry) = ((p1.x-p2.x)/dist,(p1.y-p2.y)/dist)
            f = force(dist)
            forceSumX += f * (p1.x-p2.x) / dist
            forceSumY += f * (p1.y-p2.y) / dist
        #print("Force on {:s} is ({:e},{:e})".format(str(p1),forceSumX,forceSumY))
        dx = forceSumX * timeInterval*timeInterval / (2*mass)
        dy = forceSumY * timeInterval*timeInterval / (2*mass)
        ddist = sqrt(dx*dx+dy*dy)
        print("Particle {:s} moved by ({:e},{:e})".format(str(p1), dx, dy))
        nextSpace.append(Particle(p1.x-dx, p1.y-dy))
    return nextSpace

def initial():
    return [Particle() for i in range(10)]

def simulate(steps):
    space = initial()
    for i in range(steps+1):
        print(i)
        if i % 10 == 0:
            plt.scatter([p.x for p in space], [p.y for p in space])
            plt.title("Step {:d}".format(i))
            plt.xlim(-0.1e-6,1.1e-6)
            plt.ylim(-0.1e-6,1.1e-6)
            plt.savefig("step{:d}.png".format(i))
            plt.close()
        space = forceStep(space)

simulate(100)