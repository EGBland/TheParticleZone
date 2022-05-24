import random
from math import pi, sin, cos, sqrt, fabs
random.seed()

def fuerza(i, shared):
    p = shared.data
    n = shared.count
    pi = p.iloc[i]
    xi = pi.x
    yi = pi.y
    ci = pi.c
    fx, fy = 0, 0
    for j in range(n):
        pj = p.iloc[j]
        cj = pj.c
        dire = (-1)**(1 + (ci * cj < 0))
        dx = xi - pj.x
        dy = yi - pj.y
        factor = dire * fabs(ci - cj) / (sqrt(dx**2 + dy**2) + eps)
        fx -= dx * factor
        fy -= dy * factor
    return (fx, fy)

# make a particle dict with randomised properties
def createParticle():
    pos = { "x": random.random(), "y": random.random() }
    subparticles = []
    for i in range(10):
        r = 0.02 * random.random()
        theta = 2*pi*random.random()
        subparticle = {
            "pos": { "x": pos["x"] + r*cos(theta), "y": pos["y"] + r*sin(theta) },
            "charge": -0.01
        }
        subparticles.append(subparticle)

    return {
        "pos": pos,
        "charge": 1,
        "subparticles": subparticles
    }

# calculate force between two individual (sub)particles
def forcePart(p, q):
    dx = p["pos"]["x"] - q["pos"]["x"]
    dy = p["pos"]["y"] - q["pos"]["y"]
    dist = sqrt(dx*dx + dy*dy)
    direction = (-1) ** (1+(p["charge"]*q["charge"] < 0))
    magnitude = direction * fabs(p["charge"] - q["charge"]) / (1 + dist)
    return (dx*magnitude, dy*magnitude)

# do force between a particle and the other particles in the space
def force(idx, particles):
    p = particles[idx]
    forceSumX, forceSumY = 0,0

    # this is where the forces are calculated
    for q in particles: # for each other particle in the space
        (fx,fy) = forcePart(p, q) # do the force between two big particles
        forceSumX += fx
        forceSumY += fy
        for qSub in q["subparticles"]: # for each of their subparticles
            (fx,fy) = forcePart(p, qSub) # do the force between my big particle and their subparticle
            forceSumX += fx
            forceSumY += fy
        
        for pSub in p["subparticles"]: # for each of my subparticles
            (fx,fy) = forcePart(pSub, q) # do the force between this subparticle and their big particle
            forceSumX += fx
            forceSumY += fy
            for qSub in q["subparticles"]: # for each of their subparticles
                (fx,fy) = forcePart(pSub, qSub) # do the force between these two subparticles
                forceSumX += fx
                forceSumY += fy
    
    return (forceSumX, forceSumY)

particleSpace = []
for i in range(2):
    particleSpace.append(createParticle())

print(force(0,particleSpace)) # these forces are RLY BIG...