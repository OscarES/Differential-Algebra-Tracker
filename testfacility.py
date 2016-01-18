from accelerator import Lattice, Element, LinearElement, Quad, Drift

sigma = 0.001 # the standard deviation that the user will enter
#epsilon = sqrt(sigmax**2*sigmaxp**2-sigmaxxp**2)

quad = Quad('quad', 0.01, 1)
print quad.printInfo()

#quad.evaluate(0,0)


## Non-linear
order = 5 # user sets this

## Lattice construction
drift = Drift('drift', 1)
lattice = Lattice('lattice')
lattice.appendElement(drift)
lattice.appendElement(quad)
lattice.evaluate(0,0)