point1 = [224, 453, 547
    ]

sourceRes = [832, 896, 610]
targRes = [ 26625, 28672, 9802 ]

scaleIt = targRes./sourceRes

point2 = point1 .* scaleIt

sprintf('%.0f   %.0f     %.0f ',point2(1),point2(2),point2(3))