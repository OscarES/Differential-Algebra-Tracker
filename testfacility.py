from accelerator import LinearElement, Quad

sigma = 0.001 # the standard deviation that the user will enter
#epsilon = sqrt(sigmax**2*sigmaxp**2-sigmaxxp**2)

quad = Quad('quad', 0.01, 1)
print quad.printInfo()




## Non-linear
order = 5 # user sets this