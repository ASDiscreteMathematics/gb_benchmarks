import resource
import fgb_sage
from time import sleep
from concurrent.futures import ThreadPoolExecutor

load('gmimc.sage')
load('poseidon.sage')
load('rescue_prime.sage')
load('ciminion.sage')

#######################################################
################# Questions to answer #################
#######################################################
# • How many rounds can we break?
# • Does this adhere to expectaions?

class MemoryMonitor:
    def __init__(self):
        self.keep_measuring = True

    def measure_usage(self, sleep_time=1, inform_every=None):
        max_usage = 0
        ctr = 0
        while self.keep_measuring:
            cur_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            max_usage = max(max_usage, cur_usage)
            if inform_every and ctr >= inform_every:
                ctr = 0
                print(f"Current memory usage: {cur_usage}")
            ctr += sleep_time
            sleep(sleep_time)
        return max_usage

class ExperimentStarter:
    def __init__(self):
        self.rate = 3
        self.capacity = 4
        self.input_sequence = [1, 2, 3]

    def __call__(self, primitive_name, prime, num_rounds):
        with ThreadPoolExecutor() as executor:
            monitor = MemoryMonitor()
            mem_thread = executor.submit(monitor.measure_usage, 2, 5)
            fn_thread = executor.submit(self.analyze_primitive, primitive_name, prime, num_rounds)
            result = fn_thread.result()
            monitor.keep_measuring = False
            max_usage = mem_thread.result()
        print(f"Resulting Gröbner basis:\n{result}")
        print(f"Peak memory usage: {max_usage} KB")

    def analyze_primitive(self, primitive_name, prime, num_rounds):
        if primitive_name == "poseidon":
            system = self.poseidon_system(prime, num_rounds)
        elif primitive_name == "rescue":
            system = self.rescue_system(prime, num_rounds)
        elif primitive_name == "gmimc":
            system = self.gmimc_system(prime, num_rounds)
        elif primitive_name == "ciminion":
            NotImplementedError(f"Getting the polynomial system for Ciminion is work in progress.")
        else:
            raise ValueError(f"No primitive with name {primitive_name} defined.")
        gb = fgb_sage.groebner_basis(system, threads=8, verbosity=get_verbose())
        gb = list(gb)
        return gb

    def poseidon_system(self, prime, num_rounds):
        R_F, R_P = 2, 1
        t = self.rate + self.capacity
        poseidon = Poseidon(prime=prime, R_F=R_F, R_P=R_P, t=t)
        hash_digest = poseidon(self.input_sequence + [0]*self.capacity)[:self.rate]
        ring = PolynomialRing(GF(prime), 'x', t*(R_F + R_P))
        system = poseidon_last_squeeze_poly_system(poseidon, ring.gens(), hash_digest)
        return system

    def rescue_system(self, prime, num_rounds):
        m = self.rate + self.capacity
        rp = RescuePrime(prime, m, self.capacity, 128, N=num_rounds)
        ring = PolynomialRing(GF(prime), 'x', m*num_rounds)
        hash_digest = rp.rescue_prime_hash(self.input_sequence)
        system = rescue_prime_last_squeeze_poly_system(rp, ring.gens(), hash_digest)
        return system

    def gmimc_system(self, prime, num_rounds):
        constants = [11, 80, 55, 33]
        ring = PolynomialRing(GF(101), 'x', len(constants))
        gmimc = Gmimc(ring, 3, 42, constants, 3, round_keys=ring.gens())
        system = gmimc(self.input_sequence, use_supplied_round_keys=True)
        return system

if __name__ == "__main__":
    set_verbose(1)
    testing = False

    if testing:
        if get_verbose() >=1: print(f"Testing primitives…")
        TestPoseidon()
        TestRescuePrime()
        TestGmimc()
        TestCiminion()
        if get_verbose() >= 1: print(f"Testing of primitives done.")

    prime_small = previous_prime(10^4)
    prime_big = previous_prime(2^128)

    es = ExperimentStarter()
    es(primitive_name, prime, num_rounds)
