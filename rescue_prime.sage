load('CompactFIPS202.sage')

class RescuePrime:
    def __init__(self, p, m, capacity, security_level, N=None):
        assert is_prime(p), f"Rescue Prime is only defined over prime fields."
        self.p = p
        self.Fp = GF(p)
        self.m = m
        self.capacity = capacity
        self.rate = m - capacity
        self.security_level = security_level

        alpha, alphainv = self._get_alphas()
        self.alpha = alpha
        self.alphainv = alphainv

        self.N = N
        if not N:
            self.N = self._get_number_of_rounds()
        self.MDS = self._get_mds_matrix()
        self.round_constants = self._get_round_constants()

    def __call__(self, input_sequence):
        Fp, rate, capacity = self.Fp, self.rate, self.capacity
        input_sequence = [Fp(elem) for elem in input_sequence]
        if len(input_sequence) % rate == 0:
            return self.rescue_prime_hash(input_sequence)
        return self.rescue_prime_wrapper(input_sequence)

    def _get_alphas(self):
        p = self.p
        alpha = 3
        while gcd(alpha, p-1) != 1:
            alpha += 1
        _, alphainv, _ = xgcd(alpha, p-1)
        return alpha, alphainv % (p-1)

    def _get_number_of_rounds(self):
        p, m, rate, capacity, security_level, alpha = self.p, self.m, self.rate, self.capacity, self.security_level, self.alpha
        # get number of rounds for Groebner basis attack
        dcon = lambda N : floor(0.5 * (alpha-1) * m * (N-1) + 2)
        v = lambda N : m*(N-1) + rate
        target = 2^security_level
        for l1 in range(1, 25):
            if binomial(v(l1) + dcon(l1), v(l1))^2 > target:
                break
        # set a minimum value for sanity and add 50%
        return ceil(1.5 * max(5, l1))

    def _get_mds_matrix(self):
        p, Fp, m = self.p, self.Fp, self.m
        # get a primitive element
        g = Fp(2)
        while g.multiplicative_order() != p-1:
            g = g + 1
        # get a systematic generator matrix for the code
        V = matrix([[g^(i*j) for j in range(2*m)] for i in range(m)])
        V_ech = V.echelon_form()
        # the MDS matrix is the transpose of the right half of this matrix
        MDS = V_ech[:, m:].transpose()
        return MDS

    def _get_round_constants(self):
        p, Fp, m, capacity, security_level, N = self.p, self.Fp, self.m, self.capacity, self.security_level, self.N
        # generate pseudorandom bytes
        bytes_per_int = ceil(len(bin(p)[2:]) / 8) + 1
        num_bytes = bytes_per_int * 2 * m * N
        seed_string = f"Rescue-XLIX({p},{m},{capacity},{security_level})"
        byte_string = SHAKE256(bytes(seed_string, "ascii"), num_bytes)
        # process byte string in chunks
        round_constants = []
        for i in range(2*m*N):
            chunk = byte_string[bytes_per_int*i : bytes_per_int*(i+1)]
            integer = sum(256^j * ZZ(chunk[j]) for j in range(len(chunk)))
            round_constants += [Fp(integer)]
        return round_constants

    def rescue_XLIX_permutation(self, state, starting_round=0, num_rounds=0):
        Fp, m, alpha, alphainv, N = self.Fp, self.m, self.alpha, self.alphainv, self.N
        MDS, round_constants = self.MDS, self.round_constants
        state = copy(state) # don't accidentally change inplace
        if num_rounds:
            N = starting_round + num_rounds
        for i in range(starting_round, N):
            # S-box
            for j in range(m):
                state[j,0] = state[j,0]^alpha
            # mds
            state = MDS * state
            # constants
            for j in range(m):
                state[j,0] += round_constants[i*2*m+j]

            # inverse S-box
            for j in range(m):
                state[j,0] = state[j,0]^alphainv
            # mds
            state = MDS * state
            # constants
            for j in range(m):
                state[j,0] += round_constants[i*2*m+m+j]
        return state

    def rescue_prime_wrapper(self, input_sequence):
        Fp, rate, capacity = self.Fp, self.rate, self.capacity
        padded_input = input_sequence + [Fp(1)]
        while len(padded_input) % rate != 0:
            padded_input.append(Fp(0))
        return self.rescue_prime_hash(padded_input)

    def rescue_prime_hash(self, input_sequence):
        Fp, m, rate, capacity = self.Fp, self.m, self.rate, self.capacity
        assert len(input_sequence) % rate == 0
        # initialize state to all zeros
        state = matrix([[Fp(0)]]*m)
        # absorbing
        absorb_index = 0
        while absorb_index < len(input_sequence):
            for i in range(rate):
                state[i,0] += input_sequence[absorb_index]
                absorb_index += 1
            state = self.rescue_XLIX_permutation(state)
        # squeezing
        output_sequence = []
        for i in range(rate):
            output_sequence.append(state[i,0])
        return output_sequence

    def rescue_prime_sponge(self, input_sequence, output_length):
        Fp, m, rate, capacity = self.Fp, self.m, self.rate, self.capacity
        assert len(input_sequence) % rate == 0
        # initialize state to all zeros
        state = matrix([[Fp(0)]]*m)
        # absorbing
        absorb_index = 0
        while absorb_index < len(input_sequence):
            for i in range(rate):
                state[i,0] += input_sequence[absorb_index]
                absorb_index += 1
            state = self.rescue_XLIX_permutation(state)
        # squeezing
        output_sequence = []
        squeeze_index = 0
        while squeeze_index < output_length:
            for i in range(rate):
                output_sequence.append(state[i,0])
                squeeze_index += 1
            if squeeze_index < output_length:
                state = self.rescue_XLIX_permutation(state)
        return output_sequence[:output_length]

class TestRescuePrime:
    def __init__(self):
        if get_verbose() >= 1: print(f"Testing Rescue Prime…")
        p = previous_prime(2^32) # 2^129
        m = 8
        capacity = 4
        security_level = 128

        rp = RescuePrime(p, m, capacity, security_level)
        digest = rp(list(range(52)))

        digest_expected = [3926287221, 2034419010, 1395627778, 657779429]
        assert digest == digest_expected, f"Regression test of Rescue Prime has failed: found {digest}, not {digest_expected}"

        assert test_rescue_prime_last_squeeze_poly_system()
        if get_verbose() >= 1: print(f"Testing of Rescue Prime completed")

def rescue_prime_last_squeeze_poly_system(rp, xs, hash_digest):
    m, rate, cap, alpha, N, MDS, round_constants = rp.m, rp.rate, rp.capacity, rp.alpha, rp.N, rp.MDS, rp.round_constants
    assert len(xs) >= m*N, f"Too few variables ({len(xs)}) for {N} rounds (state size {m}). Need at least {m*N}."
    assert len(hash_digest) == rate, f"The hash digest needs to be of length {rate} (is {len(hash_digest)})."
    MDS_inv = MDS.inverse()
    system = []
    for r in range(N):
        # forward half-round
        if r == 0: # first round
            curr_state = matrix(list(xs[:rate]) + [0]*cap).transpose()
        else:
            curr_state = matrix(xs[(r+0)*m-cap:(r+1)*m-cap]).transpose()
        curr_state = matrix([el[0]^alpha for el in curr_state.rows()]).transpose()
        curr_state = MDS * curr_state
        curr_state = curr_state + matrix(round_constants[(2*r+0)*m:(2*r+1)*m]).transpose()
        # backward half-round
        if r < N - 1:
            next_state = matrix(xs[(r+1)*m-cap:(r+2)*m-cap]).transpose()
        else: # last round
            next_state = matrix(hash_digest + list(xs[(r+1)*m-cap:(r+1)*m])).transpose()
        next_state = next_state - matrix(round_constants[(2*r+1)*m:(2*r+2)*m]).transpose()
        next_state = MDS_inv * next_state
        next_state = matrix([el[0]^alpha for el in next_state.rows()]).transpose()
        system += [n[0] - c[0] for n, c in zip(next_state.rows(), curr_state.rows())]
    return system

def rescue_prime_last_squeeze_trace(rp, preimage):
    m, rate, cap, alpha, alphainv, N, MDS, round_constants = rp.m, rp.rate, rp.capacity, rp.alpha, rp.alphainv, rp.N, rp.MDS, rp.round_constants
    assert len(preimage) == rate, f"The hash preimage needs to be of length {rate} (is {len(preimage)})."
    MDS_inv = MDS.inverse()
    trace = copy(preimage)
    state = matrix(preimage + [0]*cap).transpose()
    for r in range(N):
        # forward half-round
        state = matrix([el[0]^alpha for el in state.rows()]).transpose()
        state = MDS * state + matrix(round_constants[(2*r+0)*m:(2*r+1)*m]).transpose()
        # second half-round, also forward
        state = matrix([el[0]^alphainv for el in state.rows()]).transpose()
        state = MDS * state + matrix(round_constants[(2*r+1)*m:(2*r+2)*m]).transpose()
        if r < N - 1: # not last round
            trace.extend(el[0] for el in state.rows())
        else: # last round
            trace.extend(el[0] for el in state.rows()[rate:m])
    return trace

def test_rescue_prime_last_squeeze_poly_system():
    for p, m, cap, N in [(65519, 2, 1, 1), (101, 7, 4, 2), (1021, 8, 6, 4), (9973, 9, 2, 5)]:
        input_sequence = list(range(m-cap))
        rp = RescuePrime(p, m, cap, 128, N=N)
        ring = PolynomialRing(GF(p), 'z', m*N)
        hash_digest = rp.rescue_prime_hash(input_sequence)
        system = rescue_prime_last_squeeze_poly_system(rp, ring.gens(), hash_digest)
        trace = rescue_prime_last_squeeze_trace(rp, input_sequence)
        for p in system:
            assert p(trace) == 0, f"Polynomial evaluated at trace does not equal zero!"
            print("ok.")
    return True

