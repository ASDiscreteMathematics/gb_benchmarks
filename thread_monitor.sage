import sys
import pathlib
import resource
import fgb_sage
from time import sleep
from stdout_redirector import stderr_redirector
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

    def measure_usage(self, debug_path, sleep_time=1):
        max_usage = 0
        with open(debug_path + "mem.txt", 'w') as debug_file:
            while self.keep_measuring:
                cur_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                max_usage = max(max_usage, cur_usage)
                debug_file.write(f"{cur_usage}\n")
                debug_file.flush()
                sleep(sleep_time)
        return max_usage

class ExperimentStarter:
    def __init__(self, rate=3, capacity=4, input_sequence=[1,2,3]):
        assert len(input_sequence) == rate, f"Indicated rate and length of input sequence don't correspond ({rate} vs {len(input_sequence)})."
        self.rate = rate
        self.capacity = capacity
        self.input_sequence = input_sequence

    def __call__(self, primitive_name, prime, num_rounds):
        if get_verbose() >= 1: print(f"Starting experiment '{primitive_name}' over F_{prime} with {num_rounds} rounds.")
        result_path = f"./experiments/{primitive_name}_{num_rounds}_"
        with ThreadPoolExecutor() as executor:
            monitor = MemoryMonitor()
            mem_thread = executor.submit(monitor.measure_usage, result_path, 1)
            if get_verbose() >= 2: print(f"Memory measuring thread started.")
            system = self.get_system(primitive_name, prime, num_rounds)
            gb_thread = executor.submit(self.compute_gb, system, result_path)
            gb = gb_thread.result()
            monitor.keep_measuring = False
            max_usage = mem_thread.result()
        with open(result_path + "gb.txt", 'w') as f:
            for p in gb:
                f.write(f"{p}\n")
        with open(result_path + "summary.txt", 'w') as f:
            f.write("A summary will be found here in due time.\n")

    def get_system(self, primitive_name, prime, num_rounds):
        if get_verbose() >= 2: print(f"Retrieving polynomial system for {primitive_name}…")
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
        return system

    def poseidon_system(self, prime, num_rounds):
        R_F, R_P = num_rounds
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

    def compute_gb(self, system, debug_path):
        with open(debug_path + "fgb_debug.txt", 'w+b', buffering=0) as f, stderr_redirector(f):
            if get_verbose() >= 2: print(f"Starting Gröbner basis computation…")
            gb = fgb_sage.groebner_basis(system, threads=8, verbosity=1) # matrix_bound=10**8
        if get_verbose() >= 2: print(f"Finished computing Gröbner basis.")
        return list(gb)

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

    if len(sys.argv) <= 3:
        print("Not enough arguments. Please provide primitive_name, num_rounds, one of {'s', 'b'} for small / big prime.")
        exit()
    primitive_name = sys.argv[1]
    if primitive_name == "poseidon":
        num_rounds = [int(i) for i in sys.argv[2].split(',')]
        assert len(num_rounds) == 2, f"Poseidon takes exactly two round numbers, not {sys.argv[2]}."
    else:
        num_rounds = int(sys.argv[2])
    if sys.argv[3] == 's':
        prime = prime_small
    elif sys.argv[3] == 'b':
        prime = prime_big
    else:
        raise ValueError("Specify either 's' for the small prime or 'b' for the big prime.")

    pathlib.Path("./experiments").mkdir(parents=True, exist_ok=True)

    es = ExperimentStarter()
    es(primitive_name, prime, num_rounds)
