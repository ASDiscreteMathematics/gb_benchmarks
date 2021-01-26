# these two definitions are used to check if
# we have actually found a Groebner basis
def reduce_by_nosig(f, G):
    while any(g.lm().divides(f.lm()) for g in G):
        for g in G:
            if g.lm().divides(f.lm()):
                t = f.lm().quo_rem(g.lm())[0]*f.monomial_coefficient(f.lm())
                f -= t*g
                if f == 0: return f
    raise ValueError("Not a GB")

def check_basis(T, G):
    G = [ g.poly() for g in G ]
    reductions = []
    for i in range(len(G)):
        for j in range(i + 1, len(G)):
            tij = lcm(T[i], T[j])
            ui = tij.quo_rem(T[i])[0]
            uj = tij.quo_rem(T[j])[0]
            s = ui*G[i] - uj*G[j]
            reductions.append(reduce_by_nosig(s, G))
    return reductions

# create a signature
class PolynomialSignature(list):
        # initialize to signature [0], a monomial, and
        # signature[1], a positive integer
    def __init__(self, signature):
        verbose("polysig", level = 3)
        verbose((signature, type(signature)), level = 3)
        self._term = signature[0]
        self._index = signature[1]

    def __getitem__(self, key):
        if key == 0:
            return self._term
        elif key == 1:
            return self._index
        else:
            raise ValueError("Signature Key must be 0 or 1")

    # returns a representation "sigma*e_i" where
    # sigma is signature[0] and i is signature[1]
    def __repr__(self):
        return f"{self._term}*e{self._index}"

    # check for equality
    def __eq__(self, other):
        return (self._term == other._term and self._index == other._index)

    def __ne__(self, other):
        return not self == other

    def term(self):
        return self._term

    def index(self):
        return self._index

    def max_sig(self, other):
        if self._index > other._index:
            return self
        elif self._index == other._index:
            if self._term > other._term:
                return self
            else:
                return other
        else:
            return other

    # multiply the signature by a monomial. This only multiplies
    # sigma by the monomial
    def mult_monomial(self, monomial):
        return PolynomialSignature((monomial*self._term, self._index))

    # returns True if and only if adding t*g to self would
    # corrupt the signature
    def sig_corrupt(self, g, t):
        # if the index of the signature we are trying to
        # reduce by is larger, return True for signature
        # corruption
        if g._index > self._index:
            return True
        # if the indeces are the same, but sigma of the
        # index we are trying to reduce by is larger,
        # return True for signature corruption
        if g._index == self._index and g.mult_monomial(t)._term >= self._term:
            return True
        return False



# initialize the polynomial and its signature
# creates (sig, poly)
class SignedPolynomial(list):
    # initialize to SignedPolynomial[0], a Polynomial Signature,
    # and to SignedPolynomial[1], a polynomial
    def __init__(self, signedpoly):
        self._sig = PolynomialSignature(signedpoly[0])
        self._poly = signedpoly[1]

    def __eq__(self, other):
     return (self._sig == other._sig and self._poly == other._poly)

    def __ne__(self, other):
        return not self == other

    # return the leading monomial
    def lt(self):
        return self._poly.lm()*self._poly.lc()

    # returns the representation
    # (sigma*e_i, poly), where sigma*e_i is the
    # PolynomialSignature of poly
    def __repr__(self):
        return f"({self._sig}, {self._poly})"

    # return the PolynomialSignature
    def sig(self):
        verbose("sig", level = 3)
        return self._sig

    # return the polynomial
    def poly(self):
        return self._poly

    # multiply the polynomial and its signature by monomial
    def multiply(self, monomial):
        sig = self._sig.mult_monomial(monomial)
        poly = monomial*self._poly
        return SignedPolynomial((sig, poly))

    # return the signed S-polynomial of self and poly
    def generate_s_poly(self, poly):
        lcm_of_pair = lcm(self.lt(), poly.lt())
        t = lcm_of_pair.quo_rem(self.lt())[0]
        u = lcm_of_pair.quo_rem(poly.lt())[0]
        s_poly = t*self._poly - u*poly._poly
        # check for signature corruption
        if (self._sig.max_sig(poly.sig())).is_equal(self.sig()):
            sig = self._sig
            return SignedPolynomial((sig, s_poly))
        else:
            raise ValueError("Signature corruption.")

    # return the polynomial reduced by self*mon
    # we only perform this after checking for signature corruption
    def reduce_by(self, poly, mon):
        f = self.poly()
        g = poly.poly()
        t = mon.quo_rem(g.lm())[0]*f.monomial_coefficient(mon)
        f -= t*g
        self._poly = f
        return self

# a class for critical pairs of signed polynomials
class Pair(list):
        # initialize to Pair[0], a PolynomialSignature, Pair[1], the
        # first SignedPolynomial in the critical pair, and to Pair[2],
        # the second SignedPolynomial in the critical pair
    def __init__(self, pair):
        self._sig = PolynomialSignature(pair[0])
        self._poly1 = pair[1]
        self._poly2 = pair[2]

    # represent the critical pair as "(sigma*e_i, f, g)", where
    # sigma*e_i is the PolynomialSignature of the critical pair,
    # f is the first polynomial, and g is the second
    def __repr__(self):
        return f"({self._sig}, {self._poly1}, {self._poly2})"

    # return the signature of the critical pair
    def sig(self):
        verbose("sig", level = 3)
        return self._sig

    # return the first polynomial
    def poly1(self):
        return self._poly1

    # return the second polynomial
    def poly2(self):
        return self._poly2

    def __eq__(self, other):
        if not self._sig == other._sig:
            return False
        if not self._poly1[0] == other._poly1[0]:
            return False
        if not self._poly1[1] == other._poly1[1]:
            return False
        if not self._poly2[0] == other._poly2[0]:
            return False
        if not self._poly2[1] == other._poly2[1]:
            return False
        return True

    def __ne__(self, other):
        return not self == other

# initialize a polynomial ring and return the ring
def initialize_ring(vars, ring):
    if type(vars) == Integer:
        n = vars
        R = PolynomialRing(ring, "x", n)
    else:
        R = PolynomialRing(ring, vars)
    R.inject_variables()
    return R

# initialize the ideal and return the ring R, the ideal I,
# the sorted generators F, the ring generators X, the weight
# vector lp, and Y, a list of variables for the linear program
def initialize_ideal(polys, homogenize = False):
    R = polys[0].parent()
    if homogenize:
        I = R.ideal(polys).homogenize()
        R = I.ring()
        R.inject_variables()
    else:
        I = R.ideal(polys)
    F = list(I.gens())
    F.sort(key = lambda f: f.lm())
    for i in range(len(F)):
        F[i] = SignedPolynomial(((R(1), i), F[i]))
        # F keeps up with the signatures and their corresponding polynomials
        # here, we are just declaring the original signatures f_1 -> (1*e_1, f_1), etc.
    X = I.ring().gens()
    lp, Y = initialize_lp(X, homogenize)
    return R, I, F, X, lp, Y

# initialize the linear program to a vector of ones
# return this weight vector lp and Y, a list of variables for the
# linear program
# to be used in solving the linear program
def initialize_lp(X, homogenize = False):
    lp = MixedIntegerLinearProgram(solver="ppl", maximization=False)
    Y = [ lp[i] for i in range(len(X)) ]
    for y in Y:
        lp.set_integer(y)
        lp.add_constraint(y >= 1)
    lp.set_objective(sum(y for y in Y))
    lp.solve()
    return lp, Y

# update lp based on new lm and return the new weight vector
def update_lp(f, lp, Y, T, X, R, G, S, trace, homogenize = False):
    i = 0
    poly = f
    works = False
    bad_t = False
    lp_temp = copy(lp)
    lp_temp.solve()
    pp_old = 1
    while not works:
        lp_new = copy(lp_temp)
        if bad_t == True:
            U = [ u for u in U if u[0] != t ]
        else: U = compatible_pps(poly, lp_temp, Y)
        # if there are no monomials to compare, exit the loop
        if len(U) == 1: return lp
        # otherwise, use the hilbert heuristic to find the
        # "best" lm
        t = preferred_pp(T, [u[0] for u in U], R)
        a = t.exponents(as_ETuples=False)[0]
        new_expr = sum(a[i]*Y[i] for i in range(len(a)))
        # make it the leading monomial
        for u in [ v[0] for v in U ]:
            if u != t:
                b = u.exponents(as_ETuples=False)[0]
                lp_new.add_constraint(new_expr >= sum(b[i]*Y[i] for i in range(len(b))) + 1)
        try:
            lp_new.solve()
        except:
            # if updating the lp creates an infeasible system,
            # remove t from the options
            bad_t = True
            continue
        # insure that the signatures remained in the same order AND
        # the prior leading monomials remain unchanged
        if verify_signature_order(R, lp_temp, lp_new, Y, S, trace):
            T_temp = copy(T)
            T_temp.append(t)
            works, lp_temp = verify_order(T_temp, G, R, lp_temp, lp_new, Y, X)
            pp_old = t
            bad_t = False
        else:
            bad_t = True
            continue
    lp_temp.solve()
    return lp_temp

# return a matrix order with weight vector v and
# grevlex to break ties
def matrix_order(v, homogenize = False):
    if homogenize:
        vlen = len(v) - 1
    else:
        vlen = len(v)
    A = matrix(ZZ, vlen, vlen)
    for i in range(vlen):
        A[0,i] = v[i]
    for i in range(1,vlen):
        for j in range(vlen-i):
            A[i,j] = 1
    return TermOrder(A)

# returns the monomials of f that are compatible for f
# lp tells us our current weight vector and Y is
# a list of variables for the linear program
def compatible_pps(f, lp, Y):
    result = [(f.lm(),lp.get_values(Y))]
    TOs = list(ray.vector() for ray in lp.polyhedron().rays())
    # test every monomial of f
    for t in f.monomials()[1:]:
        b = vector(t.exponents()[0]) # exponent vector
        for v in TOs:
            success = True # innocent until proven guilty
            # test against every other monomial in f
            for u in result:
                if v * b <= v * vector(u[0].exponents()[0]):
                    # give up on this ray; try another
                    success = False
                    break
            if success:
                # we'll add both t and a "certificate" that t can lead
                result.append((t,v))
                # no point in remaining in the ray loop, since we"ve found a
                # winner for t
                break
    return result

# apply the new monomial order and return the current basis under the
# new weight ordering G, the updated ring R2, and the updated generators F
def redefine_vals(lp, Y, R, G, F):
    TO = matrix_order(lp.get_values(Y))
    R2 = R.change_ring(order=TO)
    # reorder polys in G according to new weight vector
    G = [ SignedPolynomial(((R2(g.sig().term()), g.sig().index()), R2(g.poly()))) for g in G ]
    # reorder polys in F according to new weight vector
    F = [ SignedPolynomial(((R2(f.sig().term()), f.sig().index()), R2(f.poly()))) for f in F ]
    if G[-1].poly() != 0:
        # insure that the lc is 1
        if G[-1].poly().coefficient(G[-1].poly().lm()) != 1:
            G[-1] = SignedPolynomial(((G[-1].sig().term(), G[-1].sig().index()), G[-1].poly().quo_rem(G[-1].poly().coefficient(G[-1].poly().lm()))[0]))
    return G, R2, F

# return the leading monomial that is "best" using the
# hilbert heuristic. T is old basis"s lms, U is list of compatible
# monomials to test
def preferred_pp(T, U, R):
    u0 = U[0]
    V = T + [u0]
    # create ideals with prior lms and the new options
    hp = R.ideal(V).hilbert_polynomial(algorithm="singular")
    hn = R.ideal(V).hilbert_numerator(algorithm="singular")
    result = u0
    for u in U[1:]:
        V = T + [u]
        hp2 = R.ideal(V).hilbert_polynomial(algorithm="singular")
        hn2 = R.ideal(V).hilbert_numerator(algorithm="singular")
        if hilbert_cmp((result,hp,hn),(u,hp2,hn2)) > 0:
            result, hp, hn = u, hp2, hn2
    return result

# return a negative integer if H1 is preferable to H2,
# a positve integer if H2 is preferable, and 0 if H1 and
# H2 are the same. H1 and H2 represent the Hilbert numerator
# and Hilbert polynomial for two choices of leading monomials
def hilbert_cmp(H1, H2):
    hpdiff = H1[1] - H2[1]
    # if the hilbert polynomial is smallest, choose as
    # preferred lm
    if hpdiff != 0:
        result = hpdiff.leading_coefficient()
    # otherwise compare the hilbert numerators
    else:
        hndiff = (H1[2] - H2[2]).coefficients(sparse=False)
        hndiff.reverse()
        while len(hndiff) > 0 and hndiff[0] == 0:
            hndiff = hndiff[1:]
        if len(hndiff) == 0:
                # if the Hilbert data is the same, choose the
                # monomial that is smaller under the monomial order
                # -1 will choose the monomial associated with H1,
                # 1 will choose the monomial associated with H2,
                # 0 means they have the same weighted degree and
                # defaults to the monomial associated with H1
            if H1[0] < H2[0]: result = -1
            elif H1[0] > H2[0]: result = 1
            else: result = 0
        else:
        #choose the monomial associated with the smallest
        #trailing term of the Hilbert numerator
            result = hndiff[0]
    return result

# return a list of "critical pairs" for the generators
def generate_first_pairs(F):
    # create pairs associated with zero for the initial polynomials
    # so we can add f_1, f_2, etc to the basis without an S-poly
    P = []
    S = []
    for i in range(len(F)):
        P.append(Pair((F[i].sig(), (1, i), (0, 0))))
    S.append(P[0])
    P.remove(P[0])
    return P, S

# return a sorted list of the new critical pairs generated by
# adding the most recent element to the basis
def generate_pairs(G, S, R2, lp, Y, T, X, homogenize):
    g = G[-1]
    lm_g = g.poly().lm()
    for i in range(len(G) - 1):
        f = G[i]
        if f.poly() == 0: continue
        least = lcm(lm_g, f.poly().lm())
        t = least.quo_rem(lm_g)[0]
        u = least.quo_rem(f.poly().lm())[0]
        # if they have the same sig, move on
        if f.multiply(u).sig() == g.multiply(t).sig():
            continue
        # assign the max signature to the pair
        sig = f.multiply(u).sig().max_sig(g.multiply(t).sig())
        if not syzygy_check(G, f.multiply(u).sig(), g.multiply(t).sig()):
            if sig == f.multiply(u).sig():
                S.insert(0, Pair((sig, (u, i), (t, len(G) - 1))))
            else:
                S.insert(0, Pair((sig, (t, len(G) - 1), (u, i))))
    return sort_S(S, R2, lp, Y, T, X, G)

# return a sorted list of critical pairs
def sort_S(S, R, lp, Y, T, X, G):
    if len(S) > 1:
        S.sort(key = lambda s:(R(s.sig().term()).degree(), s.sig().term()))
    return S, lp

# compute the s-poly and a signature safe reduction for pair, then
# adds it to the basis G
# returns the number of s-polynomials computed so far (spoly),
# the number of zero polys computed (zero), a list of critical pairs that
# have already been processed (generated), a boolean for whether or not
# we have computed a zero poly (zero_poly), the current basis (G),
# the updated linear program (lp), and the updated list of rewriting
# rules (rule)
def create_poly(spoly, zero, generated, G, F, R, T, X, pair, lp, Y, rule, S, trace, homogenize = False):
    # pair is (sig, (multiple of e_i, e_i #)
    # [Ex. (1, 0) for 1*e_0], (multiple of e_j, e_j #))
    zero_poly = False
    generated.append(pair)
    spoly += 1
    rule.append((pair.sig(), len(G)))
    if pair.poly2()[0] == 0: # <-- if it is one of the original polynomials
        if pair.poly1()[1] != 0:
            s = reduce_all_clo(F[pair.poly1()[1]], G, R)
        else:
            s = F[pair.poly1()[1]]
        t = s.poly();
        # t gives the polynomial associated with e_i
    else:
        # compute the s-poly and reduce by elements with smaller
        # signature
        s = G[pair.poly1()[1]].poly()*(pair.poly1()[0]) - \
        G[pair.poly2()[1]].poly()*(pair.poly2()[0])
        F.append(SignedPolynomial((pair.sig(), R(s))))
        s = reduce_all_clo(F[-1], G, R)
        t = s.poly()
    G.append(s)
    if t == 0:
        zero_poly = True
        zero += 1
    else:
        lp = update_lp(t, lp, Y, T, X, R, G, S, trace, homogenize)
    return spoly, zero, generated, zero_poly, G, lp, rule

# returns True iff sig1 or sig2 is a syzygy signature
def syzygy_check(G, sig1, sig2):
    if check_gs(G, R, sig1) or check_gs(G, R, sig2): return True
    return False

# returns True if there is a g in the basis whose signature has
# a lower index than sig and whose leading monomial divides sig
def check_gs(G, R, sig):
    if any(g.poly() != 0 and R.monomial_divides(g.poly().lm(), sig.term()) for g in G if g.sig().index() < sig.index()):
        return True
    return False

# returns True iff there is an r in rule such that the index of r
# equal the index of p, the signature of r divides the signature of p
# and the polynomial associated with r is not associated with p
# rule should be a list of tuples of the form (sig, poly #), where sig is
# a signature and poly # is the index of a poly in the basis
# p is pair, R is the ring
def rewritten(rule, p, R):
    i = len(rule)-1
    while i > -1 and rule[i][0].index() == p.sig().index():
        if R.monomial_divides(R(rule[i][0].term()), p.sig().term()):
            if rule[i][1] == p.poly1()[1] or rule[i][1] == p.poly2()[1]:
                return False
            else:
                return True
        i=i - 1
    return False

# returns a list of signature redundant SignedPolynomials(sig_red)
# and True iff there exists (tau, g) in the basis such
# that tau divides sig(f) and lm(g) divides lm(f)
def sig_redundant(sig_red, G, f):
    if len(G) == 1:
        return sig_red, False
    if f.poly() == 0:
        return sig_red, False
    if any(g.poly() != 0 and R.monomial_divides(g.sig().term(), f.sig().term()) and R.monomial_divides(g.poly().lm(), f.poly().lm()) for g in G):
        sig_red.append(f.poly())
        return sig_red, True
    return sig_red, False

# returns the signature safe reduction of f by G according to a modified
# version of the division algorithm found in Cox, Little, O"Shea
def reduce_all_clo(f, G, R):
    # division algorithm
    r = 0
    while f.poly() != 0:
        reduced = False
        if any(g.poly() != 0 and R.monomial_divides(g.poly().lm(), f.poly().lm()) and not f.sig().sig_corrupt(g.sig(), f.poly().lm().quo_rem(g.poly().lm())[0]) for g in G):
            for g in G:
                if g.poly() == 0: continue
                if R.monomial_divides(g.poly().lm(), f.poly().lm()):
                    t = f.poly().lm().quo_rem(g.poly().lm())[0]
                    if not f.sig().sig_corrupt(g.sig(), t):
                        f.reduce_by(g, f.poly().lm())
                        reduced = True
                        break
        if not reduced:
            r += f.poly().lm()*(f.poly().coefficient(f.poly().lm()))
            f._poly -= f.poly().lm()*(f.poly().coefficient(f.poly().lm()))
    f._poly = r
    return f

# returns the signature safe reduction of f by the basis elements
def reduced_basis(f, G, R):
    s = SignedPolynomial((f.sig(), f.poly()))
    r = f.poly().lm()
    s._poly -= r
    while s.poly() != 0:
        reduced = False
        if any(R.monomial_divides(g.poly().lm(), s.poly().lm()) and not s.sig().sig_corrupt(g.sig(), s.poly().lm().quo_rem(g.poly().lm())[0]) for g in G):
            for g in G:
                if R.monomial_divides(g.poly().lm(), s.poly().lm()):
                    t = s.poly().lm().quo_rem(g.poly().lm())[0]
                    if not s.sig().sig_corrupt(g.sig(), t):
                        s.reduce_by(g, s.poly().lm())
                        reduced = True
                        break
        if not reduced:
            r += s.poly().lm()*(s.poly().coefficient(s.poly().lm()))
            s._poly -= s.poly().lm()*(s.poly().coefficient(s.poly().lm()))
    s._poly = r
    return s

# returns True iff lm(f) is not divisible by another leading monomial in G
def remove_divisible_lts(f, G, R):
    # pass a polynomial to check (f), remove f from G, then test for divisibility
    lm = f.poly().lm()
    for g in G:
        if g.poly() == 0: continue
        if f != g and R.monomial_divides(g.poly().lm(), lm):
            return True
    return False

# returns True and the updated linear program (lp_temp)
# if the leading monomials of the basis elements are
# unchanged under lp_temp
# returns False and the prior linear program (lp) if updating
# the linear program changes the lms in the basis
def verify_order(T, G, R, lp, lp_temp, Y, X):
    TO = matrix_order(lp_temp.get_values(Y))
    R2 = R.change_ring(order = TO)
    H = [ g.poly() for g in G ]
    for g in G:
        if g.poly() == 0: H.remove(g.poly())
    T_new = [i for i in T if i != 0]
    H_new = [ R2(h).lm() for h in H ]
    if T_new == H_new: return True, lp_temp
    # if the monomial order changed, try adding constraint to the linear
    # program to preserve the old order
    for j in range(len(T_new) - 1):
        if T_new[j] != H_new[j]:
            lp.add_constraint(sum(T_new[j].degree(X[i], std_grading = True)*Y[i] for i in range(len(X))) >= sum(H_new[j].degree(X[i], std_grading = True)*Y[i] for i in range(len(X))) + 1)
    return False, lp

# returns True iff the prior signatures remain in the same
# order under the new linear program AND there are no signatures
# waiting to be computed that have a smaller signature than
# a basis element
def verify_signature_order(R, lp_old, lp_new, Y, S, trace):
    TO = matrix_order(lp_new.get_values(Y))
    R2 = R.change_ring(order = TO)
    trace1 = [ R2(r) for r in trace ]
    trace2 = sorted(trace1)
    # check the prior signature order against their order under the
    # new monomial order
    if trace1 != trace2:
        return False
    # also check that the new monomial order doesn"t cause smaller
    # signatures waiting to be computed
    if len(S) > 0 and R2(trace[-1]) > R2(S[0].sig().term()):
        return False
    return True

# returns a reduced, minimal Groebner basis by removing zero polynomials,
# elements whose leading monomial is divisible by another basis lm,
# and reduces the remaining elements by other basis elements
def clean_g(G, sig_red, R):
    # get rid of polys with same leading term and reduce each poly
    G_new = []
    for g in G:
        if g.poly() != 0 and not (g.poly() in sig_red or remove_divisible_lts(g, G, R)):
            G_new.append(g)
    for g in G_new:
        g = reduced_basis(g, G_new, R)
    return G_new

# returns a signature Groebner basis (G), the corresponding
# weight vector (w), and a list of leading monomials of basis
# elements (T)
# needs as input a list of generators (polys)
# if you want the system to be homogenized, input
# homogenize = True
def f5_dynamic(polys, homogenize = False):
    R, I, F, X, lp, Y = initialize_ideal(polys, homogenize)
    #keep up with the number of spolys
    spoly = 0
    #keep up with signature redundant and zeros
    sig_red = []
    num_red = 0
    zero = 0
    #T keeps up with leading monomials
    T = [];
    #G keeps track of the GB
    G = [];
    #rule keeps up with rewrite rules
    rule = [];
    trace = []
    generated = []
    P, S = generate_first_pairs(F)
    pair = S[0]
    S.remove(pair)
    spoly, zero, generated, zero_poly, G, lp, rule = create_poly(spoly, zero, generated, G, F, R, T, X, pair, lp, Y, rule, S, trace, homogenize)
    G, R2, F = redefine_vals(lp, Y, R, G, F);
    T = [ g.poly().lm() for g in G ];
    S, lp = generate_pairs(G, S, R2, lp, Y, T, X, homogenize);
    G, R2, F = redefine_vals(lp, Y, R, G, F);
    deg = R2(pair.sig().term()).degree()
    while len(P) > 0:
        generated = []
        rule = []
        trace = []
        S.append(P[0])
        P.remove(P[0])
        print(f"working on {S[0]}")
        while len(S) > 0:
            S_empty = False
            pair = S[0]
            S.remove(pair)
            trace.append(pair.sig().term())
            if pair in generated:
                continue
            if not rewritten(rule, pair, R2):
                spoly, zero, generated, zero_poly, G, lp, rule= create_poly(spoly, zero, generated, G, F, R2, T, X, pair, lp, Y, rule, S, trace, homogenize)
                if len(G)%20 == 0:
                    print(f"Size of G: {len(G)}")
                G, R2, F = redefine_vals(lp, Y, R, G, F);
                sig_red, redundant = sig_redundant(sig_red, G[:-1], G[-1])
                if redundant:
                    rule.pop(-1)
                T = [ g.poly().lm() for g in G ];
                if not redundant and not zero_poly:
                    S, lp = generate_pairs(G, S, R2, lp, Y, T, X, homogenize)
                    G, R2, F = redefine_vals(lp, Y, R, G, F);
        G = clean_g(G, sig_red, R2)
        print(f"current size of G:      {len(G)}")
        print(f"current lp:             {[ lp.get_values(y) for y in Y ]}")
        T = [ g.poly().lm() for g in G ];
        num_red += len(sig_red)
        sig_red = []
        print(f"s-polynomials computed: {spoly}")
        print(f"signature redundant:    {num_red}")
        print(f"zero polynomials:       {zero}")
    w = [ lp.get_values(y) for y in Y ]
    G = clean_g(G, sig_red, R2)
    print(f"size of basis: {len(G)}")
    return G, w, [ g.poly().lm() for g in G ]

if __name__ == "__main__":
    R.<x, y, z> = GF(1009)[]
    polys = [x^2+y^2+z^2-1, x*y*z, x*y-1]
    G, w, lms = f5_dynamic(polys)
    print(f"–––––––––––––––––")
    gb = [sig_poly.poly() for sig_poly in G]
    print(f"{gb}")
    print(f"{w}")
    print(f"–––––––––––––––––")
    R = R.change_ring(order="lex")
    gb = [R(poly) for poly in gb]
    print(f"lex       – is GB: {Ideal(gb).basis.is_groebner()}")
    R = R.change_ring(order="degrevlex")
    gb = [R(poly) for poly in gb]
    print(f"degrevlex – is GB: {Ideal(gb).basis.is_groebner()}")
    order = matrix(3,[1,2,3,0,0,-1,0,-1,0])
    R = R.change_ring(order=order)
    gb = [R(poly) for poly in gb]
    print(f"[1,2,3]   – is GB: {Ideal(gb).basis.is_groebner()}")
