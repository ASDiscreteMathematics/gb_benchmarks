#!/usr/bin/env sage
# coding: utf-8

class Gmimc:
    '''
    A (probably naïve) implementation of arithmetization optimized cipher prime field GMiMC_ERF.
    For more information and details, see https://pure.royalholloway.ac.uk/portal/files/34233554/gmimc_esorics.pdf
    '''
    def __init__(self, field, exponent, master_key, constants, n, round_keys=None):
        if field.is_field():
            assert field.is_prime_field(), f"This implementation of GMiMC is only defined over prime fields, not over {field}."
        else:
            assert field.is_ring(), f"The parameter 'field' must be either a field or – for symbolic evaluations – a polynomial ring, not {type(field)}."
            assert field.base_ring().is_prime_field(), f"This implementation of GMiMC is only defined over prime fields, not over {field.base_ring()}."
        assert gcd(len(field.base_ring()), exponent) == 1, f"For finite field of size {len(field.base_ring())}, exponent {exponent} is not a permutation."
        assert n > 2, f"GMiMC needs at least 2 branches, not {n}."
        if round_keys:
            assert len(round_keys) <= len(constants), f"The supplied round keys ({len(round_keys)}) imply a greater number of round than the constants ({len(constants)})."
            round_keys = [field(k) for k in round_keys]
        master_key = field(master_key)
        constants = [field(c) for c in constants]
        self.field = field
        self.exponent = exponent
        self.__master_key = master_key
        self.constants = constants
        self.n = n
        self.__round_keys = round_keys

    def __call__(self, plaintext, use_supplied_round_keys=True):
        '''
        Encrypt the given plaintext, returning the ciphertext.
        '''
        assert len(plaintext) == self.n, f"Plaintext is of length {len(plaintext)} but should be of length {self.n}."
        ciphertext = plaintext
        for i in range(len(self.constants)):
            ciphertext = self.feistel_erf(i, ciphertext, use_supplied_round_keys=use_supplied_round_keys)
        return ciphertext


    def feistel_erf(self, i, state, use_supplied_round_keys=True):
        '''
        Applies one round of an generalized Feisteln network with expanding round function to state, returning the new state.
        '''
        assert len(state) == self.n, f"ERF input is of length {len(state)} but should be of length {self.n}."
        a = self.f(i, state[0], use_supplied_round_keys=use_supplied_round_keys)
        return [s + a for s in state[1:]] + [state[0]]

    def f(self, i, elem, use_supplied_round_keys=True):
        '''
        Adds up supplied element, ith round key, ith round constant, returning the power of the sum.
        '''
        assert i < len(self.constants), f"Only {len(self.constants)} rounds defined, but round {i} requested."
        assert i >= 0, f"No round function with negative index {i} defined."
        k = self._round_key(i, use_supplied_round_keys=use_supplied_round_keys)
        c = self.constants[i]
        return (elem + c + k)**self.exponent

    def _round_key(self, i, use_supplied_round_keys=True):
        '''
        Returns the ith round key in the key schedule, or optionally the correct supplied round key.
        '''
        round_key = (i + 1) * self.__master_key
        if use_supplied_round_keys and self.__round_keys and i < len(self.__round_keys):
            round_key = self.__round_keys[i]
        return round_key

class TestGmimc:
    def __init__(self):
        self._sym_actual_correspondance(5, [68, 10])
        self._sym_actual_correspondance(5, [12, 28, 74])
        self._sym_actual_correspondance(4, [11, 80, 55, 33])
        if get_verbose() >= 1: print(f"Testing of GMiMC completed")

    def _sym_actual_correspondance(self, n, constants):
        rounds = len(constants)
        R = PolynomialRing(GF(101), 'x', rounds)
        gmimc = Gmimc(R, 3, 10, constants, n, round_keys=R.gens())
        sy = gmimc(list(range(n)), use_supplied_round_keys=True)
        ct = gmimc(list(range(n)), use_supplied_round_keys=False)
        assert [gmimc._round_key(i, use_supplied_round_keys=True) for i in range(rounds)] == list(R.gens()), f"The key schedule does not correctly produce the symbolic keys."
        sub_keys = [gmimc._round_key(i, use_supplied_round_keys=False) for i in range(rounds)]
        assert [s(sub_keys) for s in sy] == list(ct), f"Symbolic and actual evaluation of Gmimc do not correspond."
        return True
