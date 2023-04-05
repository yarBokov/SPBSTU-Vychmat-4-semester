import matplotlib.pyplot as plt
from math import sqrt


__author__ = 'Yaroslav Bokov'


# Float range generator
def frange(froom, to, step):
    while froom < to:
        yield froom
        froom += step


# Result checker
def check(x1, x2, p1, p4, p2, p3, p5):
    jacobian = (p1 * p3 * (x1 ** (p3 - 1)) * (1 - p5) /
                (1 + x1 ** p3) ** 2 - (1 + x2)) * (x1 - p4) + x1 * (p2 + x2)
    print ('x1 = ', '%0.2f' % x1, ' x2 = ', '%0.10f' % x2, ' | Jacobian = ', jacobian)


# System solver
def solve_system(p2, p3, p5):
    res = {'x1': [], 'x2': [], 'p1': [], 'p4': []}
    not_turn_dots = {'x1': [], 'x2': [], 'p1': [], 'p4': []}
    for x1 in frange(0.1, 1.5, 0.02):
        # Get coefficients of 3rd equation
        a = 1 - (1 + x1 ** p3) * p3 * (x1 ** p3) * (1 - p5) / \
            (((1 + x1 ** p3) ** 2) * (p5 + x1 ** p3))

        b = ((1 + x1 ** p3) * p3 * (x1 ** (p3 + 1)) * (1 - p5) /
             (((1 + x1 ** p3) ** 2) * (p5 + x1 ** p3))) * (2 - p2) - 2 * x1 + 2 * x1 * p2

        c = (p2 - 1) * (1 + (x1 ** p3)) * p3 * (x1 ** (p3 + 2)) * (1 - p5) / \
            (((1 + x1 ** p3) ** 2) * (p5 + x1 ** p3)) + (x1 ** 2) * (1 - p2)

        # Get discriminant
        D = (b ** 2 - 4 * a * c)
        if D >= 0:
            # Find permissible value p4
            p4_1 = (- b + sqrt(D)) / (2 * a)
            p4_2 = (- b - sqrt(D)) / (2 * a)

            if p4_1 > 0:
                # Find another roots
                x2 = x1 * p2 / (p4_1 - x1)
                if x2 > 0:
                    p1 = x1 * (1 + x2) * (1 + x1 ** p3) / (p5 + x1 ** p3)
                    if p1 > 0:
                        # Check ramification dots
                        if x1 == p4_1 and x1 * x2 == 0:
                            not_turn_dots['x1'].append(x1)
                            not_turn_dots['x2'].append(x2)
                            not_turn_dots['p1'].append(p1)
                            not_turn_dots['p4'].append(p4_1)
                        else:
                            res['x1'].append(x1)
                            res['x2'].append(x2)
                            res['p1'].append(p1)
                            res['p4'].append(p4_1)

            if p4_2 >= 0:
                # Find another roots
                x2 = x1 * p2 / (p4_2 - x1)
                if x2 > 0:
                    p1 = x1 * (1 + x2) * (1 + x1 ** p3) / (p5 + x1 ** p3)
                    if p1 > 0:
                        # Check ramification dots
                        if x1 == p4_2 and x1 * x2 == 0:
                            not_turn_dots['x1'].append(x1)
                            not_turn_dots['x2'].append(x2)
                            not_turn_dots['p1'].append(p1)
                            not_turn_dots['p4'].append(p4_2)
                        else:
                            res['x1'].append(x1)
                            res['x2'].append(x2)
                            res['p1'].append(p1)
                            res['p4'].append(p4_2)

    print (' x1          x2               p1              p4')
    for i in range(0, len(res['x1']), 1):
        print ('%0.2f' % res['x1'][i], '  ', '%0.10f' % res['x2'][i], '  ', \
            '%0.10f' % res['p1'][i], '  ', '%0.10f' % res['p4'][i])

    if len(not_turn_dots['x1']) != 0:
        print
        print ('Not turn dots')
        print ('x1             x2                p1               p4')
        for i in range(0, len(not_turn_dots['x1']), 1):
            print (not_turn_dots['x1'][i], '  ', not_turn_dots['x2'][i], '  ', \
                not_turn_dots['p1'][i], '  ', not_turn_dots['p4'][i])

    print
    print ('Check result:')
    for i in range(0, len(res['x1']), 1):
        check(res['x1'][i], res['x2'][i], res['p1'][i], res['p4'][i], p2, p3, p5)

    return res

# Get result from system solve
result = solve_system(1.5, 3, 0.01)

# Draw bifurcation diagram
plt.title(u'Bifurcation diagram')
plt.xlabel('p4')
plt.ylabel('p1')
plt.plot(result['p4'], result['p1'], '.')
plt.xlim(0, 5)
plt.ylim(0, 22)
plt.show()
