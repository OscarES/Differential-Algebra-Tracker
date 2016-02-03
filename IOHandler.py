import numpy as np

def saveMultipart(filename, multipart):
    np.save(filename, multipart)
    return 1

def loadMultipart(filename):
    multipart = np.load(filename)
    return multipart