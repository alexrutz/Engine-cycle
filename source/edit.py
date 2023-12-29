from math import acos, asin
from side import *
from cp import *
from modules import *

def createSkellet(anchor, beta1, beta2):

    print(beta1, beta2)

    beta1 = beta1 * (pi / 180)
    beta2 = beta2 * (pi / 180)

    print(beta1, beta2)

    r = 0
    dX = sin(beta1) * r

    while r <= dX + 2.1:

        dX = sin(beta1) * r
        r = r + 0.01
        print(r)

    print("initial r found: " + str(r))

    beta2_calc = 3 * pi / 4

    print(beta2_calc, beta2)
    
    while beta2_calc > beta2:

        beta2_calc = asin((sin(beta1) * r + 2) / r)
        r = r + 0.01
        print(beta2_calc, beta2)
        print(r)

    print(beta2_calc * (180 / pi))

    dY = cos(beta1) * r
    
    return circle(int((1.5 * pi + beta1) * (180 / pi)), int((1.5 * pi + beta2) * (180 / pi)), [anchor[0] - dX, anchor[1] - dY], r)

createSkellet([1, 3], 40, 60)