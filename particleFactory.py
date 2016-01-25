def straight(particles):
    x = np.linspace(-1,1,particles) #x and xp for when xp is 0
    xp = np.linspace(0,0,particles)
    y = np.linspace(-1,1,particles) #y and yp for when yp is 0
    yp = np.linspace(0,0,particles)
    return x,xp,y,yp

def scanned(particles):
    x = np.linspace(-0.05,0.05,particles)     # x and xp for when xp is scanned
    xp = np.linspace(-0.0001,0.0001,particles)
    y = np.linspace(-0.05,0.05,particles)     # y and yp for when yp is scanned
    yp = np.linspace(-0.0001,0.0001,particles)
    return x,xp,y,yp
                           
def randomed(particles):
    x = [random.uniform(-0.05, 0.05) for _ in xrange(particles)]
    xp = [random.uniform(-0.00001, 0.00001) for _ in xrange(particles)]
    y = [random.uniform(-0.05, 0.05) for _ in xrange(particles)]
    yp = [random.uniform(-0.00001, 0.00001) for _ in xrange(particles)]
    return x,xp,y,yp

def gaussian(particles):
    x = np.random.normal(0,0.001,particles)
    xp = np.random.normal(0,0.000001,particles)
    y = np.random.normal(0,0.001,particles)
    yp = np.random.normal(0,0.000001,particles)
    return x,xp,y,yp

def gaussiantwiss(particles, alpha, beta, epsilon):
    xi = np.random.normal(0,1,particles)
    xip = np.random.normal(0,1,particles)

    M = np.array([
        [1/sqrt(beta*epsilon), 0],
        [alpha/sqrt(beta*epsilon), sqrt(beta/epsilon)]
        ])

    Minv = np.linalg.inv(M)

    x = np.zeros(particles)
    xp = np.zeros(particles)

    for i in range(particles):
        x[i] = Minv[0,0]*xi[i] + Minv[0,1]*xip[i]
        xp[i] = Minv[1,0]*xi[i] + Minv[1,1]*xip[i]

    return x,xp