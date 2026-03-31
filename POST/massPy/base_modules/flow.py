from . import numdiff
import numpy as np

def density(ff):
    return np.sum(ff, axis=1)

def velocity(ff, LX, LY):
    d = density(ff)
    return np.asarray([ (ff.T[1] - ff.T[2] + ff.T[5] - ff.T[6] - ff.T[7] + ff.T[8]) / d,
                        (ff.T[3] - ff.T[4] + ff.T[5] - ff.T[6] + ff.T[7] - ff.T[8]) / d
                      ]).reshape(2, LX, LY)

def vorticity(ff, LX, LY):
    vx, vy = velocity(ff, LX, LY)
    return numdiff.curl2D(vx, vy)

def okubo_weiss(ff, LX, LY):
    # see doi:10.1016/j.dsr2.2004.09.013
    # flow & vorticity field
    vx, vy = velocity(ff, LX, LY)
    # nomral strain
    sn = numdiff.derivX(vx) - numdiff.derivY(vy)
    # shear strain
    ss = numdiff.derivX(vy) + numdiff.derivY(vx)
    # vorticity
    w = numdiff.curl2D(vx, vy)
    return sn**2 + ss**2 - w**2

def strain_rate():
    pass

def rotation_rate():
    pass