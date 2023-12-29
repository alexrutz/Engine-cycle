def circle(start, end, root, radius):

    start = int(10000 * start)
    end = int(10000 * end)
    lul = [start, end]
    startlul = min(lul)
    endlul = max(lul)
    print("circle start: ", start, end, root ,radius)

    x_coords = []
    y_coords = []

    for degree in range(startlul, endlul + 1):

        rad = degree / 10000

        x = (cos(rad) * radius) + root[0]
        y = (sin(rad) * radius) + root[1]

        x_coords.append(x)
        y_coords.append(y)

    coords = [x_coords, y_coords]

    print("resolution of circle:   ", len(x_coords))

    print("**************************************************")

    return coords

def createSkellet(anchor, beta1, beta2):

    if beta1 == beta2:

        beta1 = beta1 + 0.5

    if beta1 < beta2:

        print("skellet starting values:", beta1 * (180 / pi), beta2 * (180 / pi), anchor)

        r = 0
        dX = sin(beta1) * r

        while r <= dX + 2.001:

            dX = sin(beta1) * r
            r = r + 0.001

        beta2_calc = asin((sin(beta1) * r + 2) / r)

        print("initial r found: ", r)
        print("beta2:           ", beta2 * (180 / pi), beta2)
        print("beta2_calc_start:", beta2_calc * (180 / pi))
        
        while beta2_calc > beta2:

            beta2_calc = asin((sin(beta1) * r + 2) / r)
            r = r + 0.001

        dY = cos(beta1) * r
        dX = sin(beta1) * r

        print("beta2_calc_final:", beta2_calc * (180 / pi))
        print("dY, dX:", dY, dX)
        print("r:     ", r)
        
        return circle((pi / 2) - beta1, (pi / 2) - beta2, [anchor[0] - dX, anchor[1] - dY], r)

    if beta1 > beta2:

        print("skellet starting values:", beta1 * (180 / pi), beta2 * (180 / pi), anchor)

        r = 0
        dX = - sin(beta1) * r

        while r <= dX + 2.001:

            dX = - sin(beta1) * r
            r = r + 0.001

        beta2_calc = acos(- (sin(beta1) * r + 2) / r) - (pi / 2)

        print("initial r found: ", r)
        print("beta2:           ", beta2 * (180 / pi), beta2)
        print("beta2_calc_start:", beta2_calc * (180 / pi), beta2_calc)
        
        while beta2_calc < beta2:

            beta2_calc = acos((- (sin(beta1) * r) + 2) / r) - (pi / 2)
            r = r + 0.001

        dY = - cos(beta1) * r
        dX = - sin(beta1) * r

        print("beta2_calc_final:", beta2_calc * (180 / pi))
        print("dY, dX:", dY, dX)
        print("r:     ", r)
        
        return circle((3 / 2) * pi - beta1, (3 / 2) * pi - beta2, [anchor[0] - dX, anchor[1] - dY], r)