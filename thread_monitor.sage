import os
import sys
import pathlib
import fgb_sage
import subprocess
from time import sleep
from multiprocessing import Process, Pipe
from stdout_redirector import stderr_redirector

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
    def check_pid(self, pid):
        """ Check For the existence of a unix pid. """
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        else:
            return True

    def __init__(self, debug_path, sleep_time=1):
        self.debug_path = debug_path
        self.sleep_time = sleep_time
        self.max_usage = 0

    def measure_usage(self, pid, pipe):
        debug_path, sleep_time = self.debug_path, self.sleep_time
        max_usage = 0
        with open(debug_path + "mem.txt", 'w') as debug_file:
            while self.check_pid(pid):
                ps_res = subprocess.run(["ps", "-p", f"{pid}", "-o", "rss="], capture_output=True, text=True)
                if self.check_pid(pid): # make sure that ps returns parseable output: process to maesure has to live before and after running ps
                    cur_usage = int(ps_res.stdout)
                    max_usage = max(max_usage, cur_usage)
                    debug_file.write(f"{cur_usage}\n")
                    debug_file.flush()
                    sleep(sleep_time)
        self.max_usage = max_usage
        pipe.send(max_usage)
        return max_usage

class ExperimentStarter:
    def __init__(self, result_path, rate=1, capacity=1, input_sequence=[42]):
        assert len(input_sequence) == rate, f"Indicated rate and length of input sequence don't correspond ({rate} vs {len(input_sequence)})."
        self.result_path = result_path
        self.rate = rate
        self.capacity = capacity
        self.input_sequence = input_sequence

    def __call__(self, primitive_name, prime, num_rounds):
        result_path = self.result_path
        if get_verbose() >= 1: print(f"Starting experiment '{primitive_name}' over F_{prime} with {num_rounds} rounds.")
        system = self.get_system(primitive_name, prime, num_rounds)
        gb = self.compute_gb(system, result_path)
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
    prime = fgb_sage.MAX_PRIME

    if testing:
        if get_verbose() >=1: print(f"Testing primitives…")
        TestPoseidon()
        TestRescuePrime()
        TestGmimc()
        TestCiminion()
        if get_verbose() >= 1: print(f"Testing of primitives done.")

    pathlib.Path("./experiments").mkdir(parents=True, exist_ok=True)

    assert len(sys.argv) >= 2, f"Not enough arguments. Please provide primitive_name."
    primitive_name = sys.argv[1]
    if sys.argv[1] == "poseidon":
        assert len(sys.argv) >= 3, f"When running the poseidon experiment, make sure to additionally specifiy the number of partial rounds."
        num_part_rounds = int(sys.argv[2])
        assert num_part_rounds >= 1, f"The number of partial rounds needs to be positive, not {num_part_rounds}."

    r = 1
    memory_exhausted = False
    while not memory_exhausted:
        num_rounds = r
        if primitive_name == "poseidon":
            num_rounds = (2*r, num_part_rounds)

        result_path = f"./experiments/{primitive_name}_{prime}_{num_rounds}_"

        monitor = MemoryMonitor(result_path)
        es = ExperimentStarter(result_path)
        exp_process = Process(target=es, args=(primitive_name, prime, num_rounds))
        exp_process.start()
        if get_verbose() >= 2: print(f"Experiment process started… ({num_rounds} rounds)")

        mem_parent_pipe, mem_child_pipe = Pipe()
        mem_process = Process(target=monitor.measure_usage, args=(exp_process.pid, mem_child_pipe))
        mem_process.start()
        if get_verbose() >= 2: print(f"Memory measuring process started… ({num_rounds} rounds)")

        exp_process.join()
        mem_process.join()
        max_usage = mem_parent_pipe.recv()
        print(f"max memory usage: {max_usage}\n")
        if exp_process.exitcode < 0:
            if get_verbose() >= 2: print(f"Experiment process has terminated with exit code {exp_process.exitcode}")
            memory_exhausted = True
        else:
            r += 1
    if get_verbose() >= 2: print(f"Finished with this line of experiments. Reached round {num_rounds} but couldn't complete it.")