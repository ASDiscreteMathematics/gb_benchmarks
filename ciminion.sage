#!/usr/bin/env sage
# coding: utf-8

class Ciminion():
    '''
    A (probably naïve) implementation of arithmetization optimized cipher Ciminion.
    For more information and details, see https://eprint.iacr.org/2021/267
    '''
    def __init__(self, field, constants, master_key, N=None, R=None, IV=1, round_keys=None):
        if field.is_field():
            assert field.is_prime_field(), f"This implementation of Ciminion is only defined over prime fields, not over {field}."
        else:
            assert field.is_ring(), f"The parameter 'field' must be either a field or – for symbolic evaluations – a polynomial ring, not {type(field)}."
            assert field.base_ring().is_prime_field(), f"This implementation of Ciminion is only defined over prime fields, not over {field.base_ring()}."
        # security parameters according to Table 1, ignoring the requirement s ≥ 64
        s = field.characteristic()
        if not N: N = s + 6
        if not R: R = max(ceil((s + 37) / 12), 6)
        assert len(constants) == 4*(N + R), f"Wrong number of constants ({len(constants)}) for desired number of rounds (N = {N}, R = {R})."
        assert len(master_key) == 2, f"Master key must consist of 2 elements, but {len(master_key)} were given."
        assert N > 0, f"pC needs to apply round function f at least 1 time, not {N} times."
        assert R > 0, f"pE needs to apply round function f at least 1 time, not {R} times."
        assert all([constants[4*i + 3] > 1 for i in range(N + R)]), f"The constant 'RC4i' needs to be ≠ 0 and ≠ 1 for all i."
        assert field(IV) >= 0, f"Given initialization vector {IV} cannot be interpreted as a field element."
        if round_keys:
            round_keys = [field(k) for k in round_keys]
        constants = [field(c) for c in constants]
        master_key = [field(mk) for mk in master_key]
        IV = field(IV)
        self.field = field
        self.constants = constants
        self.__master_key = master_key
        self.N = N
        self.R = R
        self.IV = IV
        self.__round_keys = round_keys
        self.__round_key_idx = 0
        self.__key_state = (IV, master_key[0], master_key[1])

    def __call__(self, nonce, plaintext, use_supplied_round_keys=True):
        '''
        Encrypt the given plaintext using supplied nonce.
        If round keys were supplied during initialization, use those, defaulting to the real key schedule
        (picked up at the correct position) after the list of supplied round keys is exhausted.
        '''
        assert is_even(len(plaintext)), f"Ciminion is only defined for an even number of plaintext chunks, not {len(plaintext)}."
        assert len(plaintext) > 0, f"Please call this function only if you want to actually encrypt something :)"
        field = self.field
        pc = self.pc
        pe = self.pe
        rol = self.rol
        next_round_key = lambda : self._next_round_key(use_supplied_round_keys=use_supplied_round_keys)
        nonce = field(nonce)
        plaintext = [field(pt) for pt in plaintext]
        self._reset_key_schedule()
        k_0 = next_round_key()
        k_1 = next_round_key()
        initial_state = vector([nonce, k_0, k_1])
        middle_state = pc(initial_state)
        out_state = pe(middle_state)
        ciphertext = [out_state[0] + plaintext[0], out_state[1] + plaintext[1]]
        for i in range(2, len(plaintext), 2):
            middle_state[1] += next_round_key()
            middle_state[2] += next_round_key()
            middle_state = rol(middle_state)
            out_state = pe(middle_state)
            ciphertext += [out_state[0] + plaintext[i], out_state[1] + plaintext[i + 1]]
        self._reset_key_schedule()
        return ciphertext

    def f(self, i, state):
        '''
        Apply ith round function f_i to state, returning the new state.
        '''
        assert i < self.N + self.R, f"The round function f_{i} is not defined with N = {N} and R = {R}."
        assert len(state) == 3, f"The state must be of size 3, not {len(state)}."
        field = self.field
        constants = self.constants
        a, b, c = [field(s) for s in state]
        new_state = vector([a, b, a*b + c])
        c = constants[4*i + 3]
        mult_matrix = matrix([[0, 0, 1], [1, c, c], [0, 1, 1]])
        round_constants = vector([constants[4*i + 2], constants[4*i + 0], constants[4*i + 1]])
        return mult_matrix*new_state + round_constants

    def __p__(self, state, starting_round, num_rounds):
        '''
        Internal helper function de-duplicating code for pE and pC. Applies round function f_i a number of times.
        '''
        assert len(state) == 3, f"The state must be of size 3, not {len(state)}."
        field = self.field
        new_state = vector([field(s) for s in state])
        for i in range(num_rounds):
            new_state = self.f(starting_round + i, new_state)
        return new_state

    def pc(self, state):
        '''
        Apply the permutation pC to state, returning the new state.
        '''
        assert len(state) == 3, f"The state must be of size 3, not {len(state)}."
        return self.__p__(state, 0, self.N)

    def pe(self, state):
        '''
        Apply the permutation pE to state, returning the new state.
        '''
        assert len(state) == 3, f"The state must be of size 3, not {len(state)}."
        return self.__p__(state, self.N, self.R)

    def rol(self, state):
        '''
        The rolling function applies the Toffoli-gate and returns the new state.
        '''
        assert len(state) == 3, f"The state must be of size 3, not {len(state)}."
        field = self.field
        a, b, c = [field(s) for s in state]
        return vector([a*b + c, a, b])

    def _next_round_key(self, use_supplied_round_keys=True):
        '''
        Returns the next round key in the key schedule.
        '''
        # always advance key state because len(plaintext) is only known at runtime
        self.__key_state = self.pc(self.__key_state)
        round_key = self.__key_state[0]
        if use_supplied_round_keys and self.__round_keys and self.__round_key_idx < len(self.__round_keys):
            round_key = self.__round_keys[self.__round_key_idx]
        self.__round_key_idx += 1
        return round_key

    def _reset_key_schedule(self):
        '''
        Resets the key derivation function.
        '''
        self.__key_state = (self.IV, self.__master_key[0], self.__master_key[1])
        self.__round_key_idx = 0

class TestCiminion():
    def __init__(self):
        self._sym_actual_correspondance(4, 1)
        self._sym_actual_correspondance(6, 5)
        if get_verbose() >= 1: print(f"Testing of Ciminion completed")

    def _sym_actual_correspondance(self, pt_len, nonce):
        constants = [43, 60, 20, 22, 19, 94, 19, 4, 98, 62, 28, 24, 76, 7, 61, 100, 69, 28, 75, 72] # chosen by fair dice roll. guaranteed to be random.
        R = PolynomialRing(GF(101), 'x', pt_len)
        cim = Ciminion(R, constants, (10, 10), 3, 2, round_keys=R.gens())
        sy = cim(nonce, list(range(pt_len)), use_supplied_round_keys=True)
        ct = cim(nonce, list(range(pt_len)), use_supplied_round_keys=False)
        assert [cim._next_round_key(use_supplied_round_keys=True) for _ in range(pt_len)] == list(R.gens()), f"The key schedule does not correctly produce the symbolic keys."
        cim._reset_key_schedule()
        round_keys = [cim._next_round_key(use_supplied_round_keys=False) for _ in range(pt_len)]
        assert [s(round_keys) for s in sy] == list(ct), f"Symbolic and actual evaluation of Ciminion do not correspond."
        return True

if __name__ == "__main__":
    TestCiminion()
