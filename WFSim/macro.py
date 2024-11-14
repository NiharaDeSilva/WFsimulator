import time
import random
import math
import sys
import parser

parse = parser()
# Constants
MAX_NO_ALLELES = 5 if parse.parse_form_flag() == 0 else 999
ALLELE_TYPE = int  # Equivalent to short in C
DOUBLE_TOLERANCE = 0.00001
TRUE = True
FALSE = False
C_RANDOM_FLAG = parse.parse_r_flag()
ASCII_ZERO = 48
EXTRA_PROPORTION_OF_BUFFER_LOCI = 2

# Struct Equivalent (using dataclass)
from dataclasses import dataclass

@dataclass
class GType:
    pgtype: list  # equivalent to ALLELE_TYPE*
    mgtype: list  # equivalent to ALLELE_TYPE*

    # Global Variables
    initial_indivs_data = None
    final_indivs_data = None

    # Options for GFSR types (equivalent macros)
    GFSR_STYPE = int
    GFSR_STYPE_UNSIGNED = int
    GFSR_STYPE_UNSIGNED_MAX = sys.maxsize

    # Constants for GFSR
    STRUCT_GTYPE_SIZE = 2 * sys.getsizeof(int)
    STRUCT_GTYPE_STAR_SIZE = sys.getsizeof(list)

    # Pseudorandom GFSR definitions
    P = 98
    Q = 27
    GFSR_RESET_PERIOD = 100 * P
    GFSR_FLUSH_ITERATIONS = 5000 * P
    rand_table = [0] * P
    jindic = 0
    filenameGFSR = "INITFILE"

    # # Functions to mimic C macros
    # def gfsr_sbits():
    #     return 8 * sys.getsizeof(GFSR_STYPE)

    # def branch_with_probability(p):
    #     return 1 if p > gfsr4() else 0

    # def uniform_random_distribution(l, t):
    #     return disrand(l, t)
    #
    # def random_quantized_interval_selection(min_val, max_val, step):
    #     return step * disrand(int(min_val / step), int(max_val / step))
    #
    # def reset_gfsr():
    #     assert gfsr_sbits() >= 1
    #     assert GFSR_RESET_PERIOD >= 0
    #     assert GFSR_FLUSH_ITERATIONS >= 0
    #
    #     global rand_table
    #
    #     for i in range(P):
    #         rand_table[i] = 0
    #
    #     for mask_index in range(gfsr_sbits()):
    #         for i in range(P):
    #             rand_table[i] |= 1 << mask_index
    #         for _ in range(GFSR_RESET_PERIOD):
    #             intrand()
    #
    #     for _ in range(GFSR_FLUSH_ITERATIONS):
    #         intrand()
    #
    # def intrand():
    #     global jindic
    #     if C_RANDOM_FLAG:
    #         return random.randint(0, sys.maxsize)
    #     else:
    #         jindic = (jindic + 1) % P
    #         rand_table[jindic] ^= rand_table[(jindic + Q) % P]
    #         return rand_table[jindic]
    #
    # def disrand(l, t):
    #     return (intrand() % (t - l + 1)) + l
    #
    # def gfsr4():
    #     if C_RANDOM_FLAG:
    #         return random.random()
    #     else:
    #         return float(intrand()) / (GFSR_STYPE_UNSIGNED_MAX + 1.0)
    #
    # def opengfsr():
    #     if C_RANDOM_FLAG:
    #         random.seed(time.time())
    #     else:
    #         if len(filenameGFSR) != 0:
    #             try:
    #                 with open(filenameGFSR, "r") as file:
    #                     for i in range(P):
    #                         rand_table[i] = int(file.readline().strip())
    #                 global jindic
    #                 jindic = 0
    #             except FileNotFoundError:
    #                 print(f"I need {filenameGFSR}! Where is {filenameGFSR}?")
    #                 sys.exit(1)
    #         else:
    #             reset_gfsr()
    #
    # def closegfsr():
    #     if not C_RANDOM_FLAG:
    #         with open(filenameGFSR, "w") as file:
    #             for i in range(jindic, P + jindic):
    #                 file.write(f"{rand_table[i % P]}\n")

    # Writing output functions
    def write_output(samp_data, final_indivs_count):
        num_loci = parse.parse_n_loci()
        num_samples = parse.parse_iterations()

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

    def number_of_alleles_dump(number_of_alleles):
        num_loci = parse.parse_n_loci()
        num_samples = parse.parse_iterations()

        print("Auto-generated number of alleles output.")
        for k in range(num_samples):
            for i in range(num_loci):
                print(number_of_alleles[k][i], end=" ")
            print()

    def gtype_dump(g_type, number_of_alleles):
        num_loci = parse.parse_n_loci()
        num_samples = parse.parse_iterations()

        print("Auto-generated gType output.")
        for k in range(num_samples):
            for i in range(num_loci):
                print("(", end="")
                for j in range(number_of_alleles[k][i]):
                    print(f"{g_type[k][i][j]:2d} ", end="")
                print(")", end="")
            print()

    def gcount_dump(g_count, number_of_alleles):
        num_loci = parse.parse_n_loci()
        num_samples = parse.parse_iterations()

        print("Auto-generated gcount output.")
        for k in range(num_samples):
            for i in range(num_loci):
                print("(", end="")
                for j in range(number_of_alleles[k][i]):
                    print(f"{g_count[k][i][j]:2d} ", end="")
                print(")", end="")
            print()


