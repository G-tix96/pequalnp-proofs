"""Exact finite searches at x=1 for three sparse cubes (L-gct-237)."""

from fractions import Fraction
from itertools import combinations, product
from math import comb, factorial, gcd, lcm, prod
from time import monotonic


PRIME = 65521


def cube_jet(support, end):
    """Taylor coefficients of the cube through degree end at 1."""
    a, b, c = support
    terms = ((a, c - b), (b, a - c), (c, b - a))
    return taylor_cube(terms, end)


def expand_cube(terms):
    expanded = {}
    for (i, u), (j, v), (k, w) in product(terms, repeat=3):
        n = i + j + k
        expanded[n] = expanded.get(n, 0) + u * v * w
    return expanded


def taylor_cube(terms, end):
    expanded = expand_cube(terms)
    return tuple(sum(v * comb(n, d) for n, v in expanded.items()) for d in range(end + 1))


def primitive(vector):
    """Integer coordinates up to a nonzero real scalar, with canonical sign."""
    scale = gcd(*vector)
    if not scale:
        raise ValueError("coincident points in the projected Taylor coefficients")
    if next(n for n in vector if n) < 0:
        scale = -scale
    return tuple(n // scale for n in vector)


def order_twelve():
    shapes = [(1, 8, 15), (0, 11, 13), (3, 5, 16)]
    jets = [cube_jet(shape, 12) for shape in shapes]
    jet = [sum(Fraction(w * row[d], row[6]) for w, row in zip((-2, 1, 1), jets))
           for d in range(13)]
    if jet != [0] * 12 + [Fraction(198, 7)]:
        raise ArithmeticError(f"Unexpected Taylor expansion: {jet}")
    print(f"Explicit three-cube sum at 1: order 12, leading coefficient {jet[12]}.",
          flush=True)


def double_roots(max_degree, order):
    # Every trinomial with a double root at 1 and this support is a scalar
    # multiple of the one in cube_jet. Its cube has nonzero coefficient at 6.
    supports = list(combinations(range(max_degree + 1), 3))
    jets = [cube_jet(s, order - 1)[6:] for s in supports]
    triples = 0
    for i, a in enumerate(jets):
        directions = {}
        for j in range(i + 1, len(jets)):
            b = jets[j]
            # Normalize the first nonzero coefficient of each point to one.
            # The cross-multiplication clears the two denominators.
            direction = primitive([a[0] * y - b[0] * x for x, y in zip(a[1:], b[1:])])
            triples += directions.get(direction, 0)
            directions[direction] = directions.get(direction, 0) + 1
    print(f"Double-root family: {len(jets)} supports, {comb(len(jets), 3)} distinct "
          f"triples, {triples} collinear triples through order {order - 1}.", flush=True)
    return triples


def moment_kernel(nodes):
    """Two integer basis vectors for sum_j v_j nodes_j^m = 0, m = 1,...,7."""
    pivot = []
    for i, z in enumerate(nodes):
        if z not in (nodes[j] for j in pivot):
            pivot.append(i)
            if len(pivot) == 7:
                break
    basis = []
    for free in range(len(nodes)):
        if free in pivot:
            continue
        z = nodes[free]
        vector = [Fraction(0)] * len(nodes)
        vector[free] = Fraction(1)
        for i in pivot:
            t = nodes[i]
            others = [nodes[j] for j in pivot if j != i]
            vector[i] = -Fraction(z * prod(z - u for u in others),
                                  t * prod(t - u for u in others))
        denom = lcm(*(q.denominator for q in vector))
        column = primitive([int(denom * q) for q in vector])
        if any(sum(v * t**m for v, t in zip(column, nodes)) for m in range(1, 8)):
            raise ArithmeticError("Incorrect moment-kernel basis")
        basis.append(column)
    return basis


def quad(v):
    u, w, b = v
    return 3 * u * b - w * w - 9 * b * b


def det(rows):
    (a, b, c), (d, e, f), (g, h, i) = rows
    return a * (e*i - f*h) - b * (d*i - f*g) + c * (d*h - e*g)


def reciprocal(max_frequency):
    # With x=exp(t), remove the common monomial. Each base is A+B*cosh(k*t).
    # Its cube's coefficients at frequencies k,2k,3k satisfy quad(v) = 0.
    # Testing the even derivatives of orders 2,...,14 gives seven linear
    # equations. Each surviving quadratic uses the same two kernel parameters.
    unresolved = []
    few_nodes = 0
    for ks in combinations(range(1, max_frequency + 1), 3):
        nodes = [(h*k)**2 for k in ks for h in (1, 2, 3)]
        if len(set(nodes)) <= 7:
            # The Vandermonde system forces all nonconstant harmonics to zero.
            few_nodes += 1
            continue
        p, q = moment_kernel(nodes)
        rows = []
        for i in (0, 3, 6):
            a, b = p[i:i+3], q[i:i+3]
            cross = quad([x+y for x, y in zip(a, b)]) - quad(a) - quad(b)
            rows.append((quad(a), cross, quad(b)))
        if det(rows) == 0:
            unresolved.append(ks)
    print(f"Reciprocal family: {comb(max_frequency, 3)} distinct frequency triples, "
          f"{few_nodes} excluded by harmonic count, {len(unresolved)} unresolved by "
          "the determinant test.", flush=True)
    return unresolved


def small_coefficients(degree):
    # Within each base, coefficients are proportional to nonzero integers
    # of magnitude at most two. Allow an independent real scale on each cube.
    jets = []
    for size in (1, 2, 3):
        for support in combinations(range(degree + 1), size):
            for coeffs in product((-2, -1, 1, 2), repeat=size):
                if coeffs[0] > 0 and gcd(*coeffs) == 1:
                    jets.append(taylor_cube(tuple(zip(support, coeffs)), 14))

    # Full rank modulo a prime implies full rank over the reals for integer
    # vectors. Project to four coordinates first to find possible exceptions.
    prime = PRIME
    inverse = [0] + [pow(i, -1, prime) for i in range(1, prime)]
    weights = [[pow(t, d, prime) for d in range(15)] for t in (11, 17, 29, 37)]
    points = []
    for row in jets:
        h, *rest = [sum(x*y for x, y in zip(row, w)) % prime for w in weights]
        if not h:
            raise ValueError("choose another projection with a nonzero first coordinate")
        points.append(tuple(x * inverse[h] % prime for x in rest))
    if len(set(points)) != len(points):
        raise ValueError("choose a projection with distinct images")

    exceptions = 0
    start = last_update = monotonic()
    for i, (x, y, z) in enumerate(points):
        now = monotonic()
        if now - last_update >= 5:
            print(f"Small coefficients: {i}/{len(points)} base points, "
                  f"{now - start:.0f}s elapsed.", flush=True)
            last_update = now
        directions = {}
        for j in range(i + 1, len(points)):
            u, v, w = points[j]
            dx, dy, dz = (u-x) % prime, (v-y) % prime, (w-z) % prime
            if dx:
                key = (dy*inverse[dx] % prime)*prime + dz*inverse[dx] % prime
            elif dy:
                key = prime**2 + dz*inverse[dy] % prime
            else:
                key = prime**2 + prime
            previous = directions.get(key)
            if previous is None:
                directions[key] = j
                continue
            # A projected line may contain many points. Check each prior
            # pair on it against the unprojected integer Taylor coefficients.
            if isinstance(previous, int):
                previous = [previous]
                directions[key] = previous
            for k in previous:
                if dependent(jets[i], jets[j], jets[k]):
                    exceptions += 1
            previous.append(j)
    print(f"Small coefficients: {len(points)} bases, {comb(len(points), 3)} triples, "
          f"{exceptions} linearly dependent triples through order 14.", flush=True)
    return exceptions


def dependent(a, b, c):
    # a and b have distinct projective images, hence a nonzero 2x2 minor.
    i = next(i for i, v in enumerate(a) if v)
    j = next(j for j, v in enumerate(b) if a[i]*v != a[j]*b[i])
    return all(det(((a[i], a[j], x), (b[i], b[j], y), (c[i], c[j], z))) == 0
               for x, y, z in zip(a, b, c))


def hankel_determinant(values):
    # Expand the four-by-four Hankel determinant along its first row.
    first, *lower = [values[i:i+4] for i in range(4)]
    total = 0
    for j, value in enumerate(first):
        minor = [row[:j] + row[j+1:] for row in lower]
        total += (-1)**j * value * det(minor)
    return total % PRIME


def may_have_three_terms(target, inv0, factorials, inverse_factorials, inverse_three):
    # h is the Taylor series of the normalized cube root, with h[0] = 1.
    h = [1]
    square = [1]
    moments = [1]
    for n in range(1, len(target)):
        cross = sum(h[i] * h[n-i] for i in range(1, n))
        rest = cross + sum(h[i] * square[n-i] for i in range(1, n))
        coeff = (target[n] * inverse_factorials[n] * inv0 - rest) * inverse_three % PRIME
        h.append(coeff)
        square.append((2*coeff + cross) % PRIME)
        moments.append(coeff * factorials[n] % PRIME)
        if n >= 6 and hankel_determinant(moments[-7:]):
            return False
    return True


def unrestricted_third():
    # Normalize the proposed third base to value 1 at the candidate root.
    # Its derivative moments in log(x) must be a sum of at most 3 exponentials.
    factorials = [factorial(n) % PRIME for n in range(15)]
    inv_factorials = [pow(f, -1, PRIME) for f in factorials]
    inverse_three = pow(3, -1, PRIME)

    def possible(row):
        return may_have_three_terms(row, pow(row[0], -1, PRIME), factorials,
                                    inv_factorials, inverse_three)

    rows = []
    for support in combinations(range(7), 3):
        for b, c in product((-1, 1), repeat=2):
            poly = expand_cube(tuple(zip(support, (1, b, c))))
            moments = [sum(v * n**j for n, v in poly.items()) for j in range(15)]
            if not possible(moments):
                raise ArithmeticError("Rejected a known trinomial cube")
            rows.append(moments)
    eligible = 0
    for a, b in combinations(rows, 2):
        for u, v in ((1, 1), (1, -1), (1, 2), (1, -2), (2, 1), (2, -1)):
            if u*a[0] + v*b[0] == 0:
                continue
            eligible += 1
            target = [(u*x + v*y) % PRIME for x, y in zip(a, b)]
            if possible(target):
                raise ArithmeticError("A candidate needs exact reconstruction")
    print(f"Unrestricted third base: {eligible} pairs and relative weights excluded.",
          flush=True)


if __name__ == "__main__":
    start = monotonic()
    print("Starting exact multiplicity checks.", flush=True)
    order_twelve()
    if double_roots(16, 12) != 1:
        raise ArithmeticError("Expected the explicit order-12 triple")
    if double_roots(16, 13):
        raise ArithmeticError("Found a dependent double-root triple")
    print(f"Starting reciprocal-family check at {monotonic()-start:.1f}s.", flush=True)
    remaining = reciprocal(20)
    if remaining:
        raise ArithmeticError(f"Unresolved reciprocal frequency triples: {remaining}")
    print(f"Starting bounded-coefficient check at {monotonic()-start:.1f}s.", flush=True)
    if small_coefficients(10):
        raise ArithmeticError("Found a dependent small-coefficient triple")
    print(f"Starting unrestricted-third-base check at {monotonic()-start:.1f}s.", flush=True)
    unrestricted_third()
    print(f"Complete in {monotonic()-start:.1f}s.", flush=True)
