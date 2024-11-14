import numpy as np
import argparse
import re
import sys
from parser import RefactorParser
import random

class OneSampEngine:

    def __init__(self):
        self.random_flag = None
        self.num_loci = None
        self.input_individuals_count = None
        self.bottleneck_individuals_count = None
        self.bottleneck_length = None
        self.final_individuals_count = None
        self.is_microsats = None
        self.min_allele_frequency = None
        self.omit_threshold = None
        self.theta = None
        self.mutation_rate = None
        self.repetitions = None
        self.program_name = None
        self.syntax_check = False
        self.example = False
        self.raw_stats = False
        self.single_generation = False
        self.example_pop = False
        self.absent_data_extrapolate = False
        self.proportion_missing_data = None

        # Stored arrays for simulation
        self.bottleneck_individuals_count_random_choices = None
        self.bottleneck_length_random_choices = None
        self.theta_random_choices = None
        self.mutation_rate_random_choices = None

        # Data structures for calculations
        self.syntax_results = [0, 0]
        self.double_data = None
        self.number_of_alleles = None
        self.g_type = None
        self.g_count = None
        self.females = None
        self.males = None
        self.final_indivs_data = None


    def reset_arguments(self):
        # Reset the arguments
        self.__init__()

    def report_error(self, message):
        # Prints an error message and exits the program
        print(f"ONESAMP ERROR: {message}")
        exit(1)

    def parse_arguments(self, argv):
        parser = argparse.ArgumentParser(description="OneSamp command line arguments handler.")

        parser.add_argument('-r', type=str, help="Source of random numbers: 'C', 'GFSR', or 'RESET'")
        parser.add_argument('-l', type=int, help="Number of polymorphic loci")
        parser.add_argument('-i', type=int, help="Number of input samples")
        parser.add_argument('-b', type=str, help="Bottleneck size in the form 'min,max'")
        parser.add_argument('-d', type=str, help="Duration of bottleneck in the form 'min,max'")
        parser.add_argument('-m', action='store_true', help="Indicates microsatellite loci")
        parser.add_argument('-s', action='store_true', help="Indicates SNP loci")
        parser.add_argument('-t', type=int, help="Number of repetitions")
        parser.add_argument('-u', type=str, help="Mutation rate in the form 'min,max'")
        parser.add_argument('-v', type=str, help="Theta value in the form 'min,max'")
        parser.add_argument('-x', action='store_true', help="Check syntax")
        parser.add_argument('-e', action='store_true', help="Generate example population")
        parser.add_argument('-w', action='store_true', help="Compute raw statistics")
        parser.add_argument('-g', action='store_true', help="Simulate a single generation")
        parser.add_argument('-p', action='store_true', help="Dump out an example population")
        parser.add_argument('-f', type=float, help="Minimum allele frequency")
        parser.add_argument('-a', action='store_true', help="Extrapolate absent data")
        parser.add_argument('-o', type=float, help="Threshold to omit loci")

        args = parser.parse_args()

        # Process arguments
        self.random_flag = args.r
        self.num_loci = args.l
        self.input_individuals_count = args.i
        if args.b:
            self.bottleneck_individuals_count = [int(x) for x in args.b.split(',')]
        if args.d:
            self.bottleneck_length = [int(x) for x in args.d.split(',')]
        self.is_microsats = True if args.m else (False if args.s else None)
        self.repetitions = args.t
        if args.u:
            self.mutation_rate = [float(x) for x in args.u.split(',')]
        if args.v:
            self.theta = [float(x) for x in args.v.split(',')]
        self.syntax_check = args.x
        self.example = args.e
        self.raw_stats = args.w
        self.single_generation = args.g
        self.example_pop = args.p
        self.min_allele_frequency = args.f
        self.absent_data_extrapolate = args.a
        self.omit_threshold = args.o

        # Argument validation
        if self.num_loci is not None and self.num_loci <= 0:
            self.report_argument_error("Number of loci must be positive.")
        if self.input_individuals_count is not None and self.input_individuals_count <= 0:
            self.report_argument_error("Number of input samples must be positive.")


    def parse_int_pair(self, value, flag_name):
        # Parses a comma-separated pair of integers
        try:
            min_val, max_val = map(int, value.split(','))
            if min_val < 0 or max_val < 0:
                self.report_error(f"Values for {flag_name} must be non-negative integers.")
            return min_val, max_val
        except ValueError:
            self.report_error(f"Invalid input for {flag_name}. Expected a comma-separated pair of integers.")

    def parse_double_pair(self, value, flag_name):
        # Parses a comma-separated pair of doubles
        try:
            min_val, max_val = map(float, value.split(','))
            if min_val < 0 or max_val < 0:
                self.report_error(f"Values for {flag_name} must be non-negative.")
            return min_val, max_val
        except ValueError:
            self.report_error(f"Invalid input for {flag_name}. Expected a comma-separated pair of floats.")

    def flush_arguments(self):
        # Resets or deallocates resources
        self.reset_arguments()


    def allocate_memory(self):
        # Allocate memory for storing results, similar to `allocateOneSampMemory` in C
        self.repetitions = self.repetitions or 10
        self.num_loci = self.num_loci or 10
        self.input_individuals_count = self.input_individuals_count or 100

        # Using numpy to simulate allocated memory
        self.double_data = np.zeros((11, self.repetitions))
        self.number_of_alleles = np.zeros((self.num_loci, self.repetitions))
        self.g_type = np.zeros((self.input_individuals_count, self.repetitions))
        self.g_count = np.zeros((self.input_individuals_count, self.repetitions))

    def allocate_generations(self):
        # Allocate generation data for females and males, similar to `malloc` for `gtype_type` in C
        bottleneck_max = self.bottleneck_max or 50  # Placeholder value
        num_loci = self.num_loci or 10

        self.females = [Gtype(num_loci) for _ in range(bottleneck_max)]
        self.males = [Gtype(num_loci) for _ in range(bottleneck_max)]

    def allocate_sampling_generations(self):
        # Allocate memory for final individuals data over all repetitions
        self.final_indivs_data = [[Gtype(self.num_loci) for _ in range(self.input_individuals_count)]
                                  for _ in range(self.repetitions)]

    def assort(self, count, dest, females, males, bottleneck, iteration, loci):
        # Placeholder for assortative mating logic. Needs implementation based on the original C function.
        pass



    def read_microsatellite_motif_lengths(self, argv):
        motif_lengths = []
        for arg in argv:
            # Check if the argument starts with '-m'
            if not arg.startswith('-m'):
                continue

            # Extract the numbers after '-m'
            motif_str = arg[2:]  # Remove '-m'

            # Use regex to extract digits
            inputs_list = re.findall(r'\d+', motif_str)
            inputs_list = list(map(int, inputs_list))  # Convert to integers

            # Validate that each value is in the set {2, 3, 4, 6}
            if not all(x in {2, 3, 4, 6} for x in inputs_list):
                return -1

            # Append to the motif lengths list
            motif_lengths.extend(inputs_list)

        # Update the class variable or return the result
        if motif_lengths:
            self.motif_lengths = motif_lengths
            return len(motif_lengths)
        else:
            return 0


    def report_error(self, message):
        # Function for error handling similar to reportError in C
        print(f"ERROR: {message}")
        exit(1)

    def run(self, argv):
        syntax_results = [40, 100]
        # Parse arguments and initialize variables
        self.reset_arguments()
        self.parse_arguments(argv)

        # if self.parse_syntax_check():
        #     self.parse_from_file(True, sys.stdin, syntax_results)
        #     print(f"-l{self.num_loci} -i{self.input_individuals_count}")
        #     return

        # Read in motifs if form flag is set
        if self.parse_form_flag():
            if self.read_microsatellite_motif_lengths(argv) != self.num_loci:
                self.report_error("Mangled list of motif lengths passed after -m flag.")

        # Allocate memory for computations
        # self.allocate_memory()
        self.allocate_generations()

        # Parse input if not using an example population
        if not self.example_pop:
            self.parse_from_file(False, sys.stdin)
            if self.num_loci == -1 or self.input_individuals_count == -1:
                self.report_error("Invalid number of loci or individuals provided.")

        if self.single_generation:
            # Allocate memory for final individuals data and simulate a single generation
            self.allocate_sampling_generations()
            self.simulate_single_generation()

        if not self.example_pop:
            # Additional filtering and processing logic here (filter loci, individuals, etc.)
            pass

        self.allocate_sampling_generations()
        self.simulate_generations()

    def simulate_single_generation(self):
        # Placeholder for simulating a single generation and computing statistics
        for i in range(1):
            for j in range(2 * self.bottleneck_max):
                self.final_indivs_data[i][j] = Gtype(self.num_loci)  # Initialize new genotypes
        self.write_output(self.final_indivs_data)

    def simulate_generations(self):
        # Main simulation loop for multiple generations (handling bottlenecks, mating, etc.)
        current = 0
        next = 1
        for i in range(self.repetitions):
            for j in range(self.bottleneck_length[i]):
                # Placeholder for actual mating and bottleneck logic
                pass
            self.assort(self.final_indivs_count, self.final_indivs_data[i], self.females[current], self.males[current],
                        self.bottleneck_max, i, self.num_loci)


    def parse_n_loci(self):
        if self.num_loci is None or self.num_loci <= 0:
            self.report_argument_error("argument -l, num of unlinked polymorphic loci, must be a positive integer.")
        return self.num_loci

    def parse_input_samples(self):
        if self.input_individuals_count is None or self.input_individuals_count <= 0:
            self.report_argument_error("argument -i, num of input samples, must be a positive integer.")
        return self.input_individuals_count


    def write_output(self, samp_data, final_indivs_count):
        num_loci = self.parse_n_loci()
        num_samples = self.parse_iterations()

        print("Auto-generated genotype output.")
        for j in range(num_loci):
            print(j + 1)
        print("Pop")

        for k in range(num_samples):
            for j in range(final_indivs_count):
                print(f"{j + 1}, ", end="")
                for i in range(num_loci):
                    print(f"{samp_data[k][j].mgtype[i]:02d}{samp_data[k][j].pgtype[i]:02d} ", end="")
                print()

# Example usage
if __name__ == "__main__":
    import sys
    engine = OneSampEngine()
    engine.run(sys.argv)


