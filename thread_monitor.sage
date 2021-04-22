import resource
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

    def measure_usage(self):
        max_usage = 0
        while self.keep_measuring:
            max_usage = max(max_usage, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            sleep(0.1)
        return max_usage

class ExperimentStarter:
    def __init__(self):
        self.rate = 3
        self.capacity = 4
        self.input_sequence = [1, 2, 3]

    def __call__(self, primitive_name, prime):
        with ThreadPoolExecutor() as executor:
            monitor = MemoryMonitor()
            mem_thread = executor.submit(monitor.measure_usage)
            result = self.analyze_primitive(primitive_name, prime)
            monitor.keep_measuring = False
            max_usage = mem_thread.result()
            print(f"Resulting Gröbner basis:\n{result}")
            print(f"Peak memory usage: {max_usage} KB")

    def analyze_primitive(self, primitive_name, prime):
        if primitive_name == "poseidon":
            system = self.poseidon_system(prime)
        elif primitive_name == "rescue":
            system = self.rescue_system(prime)
        elif primitive_name == "gmimc":
            system = self.gmimc_system(prime)
        elif primitive_name == "plookup":
            pass
        else:
            raise ValueError(f"No primitive with name {primitive_name} defined.")
        gb = Ideal(system).groebner_basis()
        return gb

    def poseidon_system(self, prime):
        R_F, R_P = 2, 1
        t = self.rate + self.capacity
        poseidon = Poseidon(prime=prime, R_F=R_F, R_P=R_P, t=t)
        hash_digest = poseidon(self.input_sequence + [0]*self.capacity)[:self.rate]
        ring = PolynomialRing(GF(prime), 'x', t*(R_F + R_P))
        system = poseidon_last_squeeze_poly_system(poseidon, ring.gens(), hash_digest)
        return system

    def rescue_system(self, prime):
        m, cap, N = 7, 4, 1
        rp = RescuePrime(prime, m, cap, 128, N=N)
        ring = PolynomialRing(GF(prime), 'x', m*N)
        hash_digest = rp.rescue_prime_hash(self.input_sequence)
        system = rescue_prime_last_squeeze_poly_system(rp, ring.gens(), hash_digest)
        return system

    def gmimc_system(self, prime):
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
    es("poseidon", 101)
    es("rescue", 101)
    es("gmimc", 101)
