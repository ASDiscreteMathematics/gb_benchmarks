import resource
from time import sleep
from concurrent.futures import ThreadPoolExecutor

## Partly adapted from https://medium.com/survata-engineering-blog/49f027e3d1bau

load('rescue_prime.sage')

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
        list_of_experiments = [
            ("rescue", 101),
        ]
        for experiment in list_of_experiments:
            self.start_experiment(experiment)

    def analyze_experiment(self, experiment):
        exp_name, prime = experiment
        if exp_name == "poseidon":
            pass
        elif exp_name == "rescue":
            pass
        elif exp_name == "gmimc":
            pass
        elif exp_name == "plookup":
            pass
        else:
            raise ValueError(f"No experiment with name {exp_name} defined.")

    def start_experiment(self, experiment):
        with ThreadPoolExecutor() as executor:
            monitor = MemoryMonitor()
            mem_thread = executor.submit(monitor.measure_usage)
            try:
                exp_thread = executor.submit(self.analyze_experiment, experiment)
                result = exp_thread.result()
            finally:
                monitor.keep_measuring = False
                max_usage = mem_thread.result()

            print(f"Peak memory usage: {max_usage} KB")

    def __call__(self):
        for experiment in self.list_of_experiments:
            self.start_experiment(experiment)

if __name__ == "__main__":
    es = ExperimentStarter()
