import random
import math
import numpy as np
import argparse
import re
import sys
# import macro
bottleneck_individuals_count_random_choices = []
bottleneck_length_random_choices = []
mutation_rate_random_choices = []
# Define INT_MAX and GFSR_STYPE_UNSIGNED_MAX based on typical C limits
INT_MAX = sys.maxsize  # This provides the maximum integer size (similar to INT_MAX in C)
GFSR_STYPE_UNSIGNED_MAX = (2**32) - 1  # Assuming a 32-bit unsigned integer max


def initializeMicrosat1(locus):
    return 200

def initializeMicrosat2(locus):
    return 202

# Function to parse the number of loci
def parseNLoci(lociDirectInput):
    return lociDirectInput

# Function to parse theta
def parseTheta(thetaDirectInput):
    return thetaDirectInput

def bottleneck_individual_count(neRange, trials):
    return [random.randint(neRange[0], neRange[1]) for _ in range(trials)]

# Function to parse the form flag (whether microsatellite or SNP)
def parseFormFlag(formFlag):
    return formFlag

# Function to parse minimum allele frequency
def parseMinAlleleFrequency(minAlleleFrequencyDirectInput):
    return minAlleleFrequencyDirectInput

# Function for a falling quotient used in allele probability calculation
def fallingQuotient(s, t1, t2, c):
    if c == 0:
        return s
    else:
        return fallingQuotient(s * t1 / t2, t1 - 1, t2 - 1, c - 1)

# Function to calculate allele probability
def allelePr(val1, val2, theta):
    n = val1 + val2
    result = fallingQuotient(1, n, theta + n - 1, n)
    if val1 != 0:
        result *= theta / val1
    if val2 != 0:
        result *= theta / val2
    if val1 == val2:
        result /= 2
    return result


# Placeholder functions to represent the original C functions
def disrand(l, t):
    return random.randint(l, t)

def rand_bit():
    return random.randint(0, 1)

def branch_with_probability(p):
    return random.random() < p

def gfsr4():
    return (random.randint(0, INT_MAX) / INT_MAX) / (GFSR_STYPE_UNSIGNED_MAX + 1.0)

def mutate_snp(gene):
    gene[0] = 2 * rand_bit() + rand_bit() + 1

def variantOfFinalLocus(sample, finalIndividuals, individualsDirectInput, locusI):
    non_zero_variant = 0
    for i in range(individualsDirectInput):
        geneA1, geneA2 = loadFinalGenotype(sample, finalIndividuals, i, locusI)
        if non_zero_variant == 0:
            non_zero_variant = geneA1
        if non_zero_variant == 0:
            non_zero_variant = geneA2
        if non_zero_variant == 0:
            continue
        if non_zero_variant != geneA1 or non_zero_variant != geneA2:
            return 2  # Two or more variants found
    if non_zero_variant == 0:
        return 0  # No variant found
    return 1


def loadFinalGenotype(sample, finalIndividuals, individual, index):
    genotype1 = finalIndividuals[sample][individual]['pgtype'][index]
    genotype2 = finalIndividuals[sample][individual]['mgtype'][index]
    return genotype1, genotype2

def storeFinalGenotype(sample, finalIndividuals, individual, index, genotype1, genotype2):
    finalIndividuals[sample][individual]['pgtype'][index] = genotype1
    finalIndividuals[sample][individual]['mgtype'][index] = genotype2
    return finalIndividuals

def parse_arguments():
    parser = argparse.ArgumentParser(description="OneSamp command-line argument parser")
    # Randomness options
    parser.add_argument('-r', choices=['C', 'GFSR', 'RESET'], help="Random number generator flag: C for C random, GFSR for GFSR random, RESET to reset GFSR")
    # Loci
    parser.add_argument('-l', type=int, help="Number of loci")
    # Initial individuals
    parser.add_argument('-i', type=int, help="Number of initial individuals")
    # Bottleneck individuals (expects two integers)
    parser.add_argument('-b', nargs=2, type=int, help="Number of individuals in bottleneck generation (must be two positive even integers)")
    # Bottleneck length (expects two integers)
    parser.add_argument('-d', nargs=2, type=int, help="Length of bottleneck period")
    # SNPs versus Microsatellites
    parser.add_argument('-m', action='store_true', help="Use microsatellites")
    parser.add_argument('-s', action='store_true', help="Use SNPs")
    # Repetitions
    parser.add_argument('-t', type=int, help="Number of repetitions")
    # Mutation rate (expects two floats)
    parser.add_argument('-u', nargs=2, type=float, help="Mutation rate")
    # Theta (expects two floats)
    parser.add_argument('-v', nargs=2, type=float, help="Theta parameter")
    # Various modes
    parser.add_argument('-x', action='store_true', help="Syntax check")
    parser.add_argument('-e', action='store_true', help="Example mode")
    parser.add_argument('-w', action='store_true', help="Raw stats mode")
    parser.add_argument('-g', action='store_true', help="Single generation mode")
    parser.add_argument('-p', action='store_true', help="Example population mode")
    # Minimum allele frequency
    parser.add_argument('-f', type=float, help="Minimum allele frequency")
    # Absent data extrapolation
    parser.add_argument('-a', action='store_true', help="Interpolate absent data")
    # Omit loci threshold
    parser.add_argument('-o', type=float, help="Threshold to omit loci")
    args = parser.parse_args()
    # Validate bottleneck individuals input (must be two positive even integers)
    if args.b:
        bottleneck_individuals = args.b
        if any(x <= 0 or x % 2 != 0 for x in bottleneck_individuals):
            print("Error: Bottleneck individuals must be positive even integers.")
            sys.exit(1)

    return args

# Generates the next generation of individuals based on the current one
def assort(next_gen_count, offvec, mothers, fathers, current_gen_count, mutation_rate, num_loci):
    # Read in mutation rate
    p_mutation = mutation_rate
    # Simulate each individual's genotype
    for j in range(next_gen_count):
        # Select random mother from current generation
        m = disrand(0, current_gen_count - 1)
        print(m)
        # Select random father from current generation
        d = disrand(0, current_gen_count - 1)
        print(d)
        # For each allele, select one from parent and possibly mutate
        for i in range(num_loci):
            # Generate a new individual from mother
            offvec[j]['mgtype'][i] = mothers[m]['mgtype'][i] if rand_bit() else mothers[m]['pgtype'][i]
            if branch_with_probability(p_mutation):
                offvec[j]['mgtype'] = mutate_snp(offvec[j]['mgtype'])

            # Generate a new individual from father
            offvec[j]['pgtype'][i] = fathers[d]['mgtype'][i] if rand_bit() else fathers[d]['pgtype'][i]
            if branch_with_probability(p_mutation):
                offvec[j]['pgtype'] = mutate_snp(offvec[j]['pgtype'])


def write_output(finalIndividuals, lociDirectInput, individualsDirectInput):
    # Output simulated data
    print("Auto-generated genotype output.")
    for locus_index in range(lociDirectInput):
        print(f"{locus_index + 1}")
    print("Pop")

    for j in range(individualsDirectInput):
        print(f"{j + 1}, ", end="")
        for locus_index in range(lociDirectInput):
            print(f"{finalIndividuals[0][j]['mgtype'][locus_index]:02d}{finalIndividuals[0][j]['pgtype'][locus_index]:02d} ", end="")
        print()

# Main simulation function
def simulate_population(lociDirectInput, neRange, thetaDirectInput, individualsDirectInput, minAlleleFrequencyDirectInput, formFlag, trials, mutation_rate, durationLength):
    bottleneck_individuals = bottleneck_individual_count(neRange,  trials)
    print(bottleneck_individuals)
    extraProportionOfBufferLoci = 2
    totalLoci = extraProportionOfBufferLoci * lociDirectInput
    for i in range(trials):
        # Initialize simulation variables
        num_genes = 4 * bottleneck_individuals[i]  # total gene copies
        minAlleleCount = math.ceil(parseMinAlleleFrequency(minAlleleFrequencyDirectInput) * num_genes)
        totalAllelePr = 0

        # Coalescent probability memo table for allele frequencies
        coalescentProbabilityMemoTable = [0] * (num_genes // 2 - minAlleleCount + 1)

        # Calculate coalescent probabilities
        for j in range(num_genes // 2 - minAlleleCount):
            coalescentProbabilityMemoTable[j] = allelePr(j + minAlleleCount, num_genes - j - minAlleleCount, thetaDirectInput)
            totalAllelePr += coalescentProbabilityMemoTable[j]

        # Normalize probabilities
        for j in range(num_genes // 2 - minAlleleCount):
            coalescentProbabilityMemoTable[j] /= totalAllelePr

        # # Initialize genotype vectors
        gvec = [0] * num_genes
        # # Store intermediate generations
        # females = [{'pgtype': [0] * totalLoci, 'mgtype': [0] * totalLoci} for _ in range(neRange[1])]
        # males = [{'pgtype': [0] * totalLoci, 'mgtype': [0] * totalLoci} for _ in range(neRange[1])]

        females = [[], []]
        males = [[], []]

        for k in range(2):
            females[k] = []
            males[k] = []
            for j in range(neRange[1]):  # Assuming `parse_bottleneck_max` is defined
                females[i].append({
                    'pgtype': [0] * (totalLoci),  # Creating a list of zeros
                    'mgtype': [0] * (totalLoci)
                })
                males[i].append({
                    'pgtype': [0] * (totalLoci),
                    'mgtype': [0] * (totalLoci)
                })


        #Store Final Generation
        # finalIndividuals = [{'pgtype': [0] * totalLoci, 'mgtype': [0] * totalLoci} for _ in range(individualsDirectInput)]
        finalIndividuals = [[], []]
        for k in range(2):
            finalIndividuals[k] = []
            for j in range(neRange[1]):
                finalIndividuals[i].append({
                    'pgtype': [0] * (totalLoci),  # Creating a list of zeros
                    'mgtype': [0] * (totalLoci) })

        # Iterate over loci and simulate gene frequencies
        for locus_index in range(0, totalLoci):
            # Pick genotypes for SNPs or Microsatellites
            base1 = 2 * random.randint(0, 1) + random.randint(0, 1) + 1 if parseFormFlag(formFlag) == 0 else initializeMicrosat1(locus_index)
            base2 = (base1 + random.randint(0, 2)) % 4 + 1 if parseFormFlag(formFlag) == 0 else initializeMicrosat2(locus_index)

            # Simulate coalescent frequency distribution
            cut = gfsr4()
            for j in range(num_genes // 2 - minAlleleCount):
                cut -= coalescentProbabilityMemoTable[j]
                if cut < 0:
                    break

            # Fill genotype vector based on frequency distribution
            for count2 in range(0, num_genes):
                if count2 < j + minAlleleCount:
                    gvec[count2] = base1
                else:
                    gvec[count2] = base2

            numleft = num_genes

            # Assign alleles to females and males randomly
            for j in range(0, bottleneck_individuals[i]):
                ic = random.randint(0, numleft - 1)
                females[0][j]['pgtype'][locus_index] = gvec[ic]
                gvec[ic] = gvec[numleft - 1]
                numleft -= 1

            for j in range(bottleneck_individuals[i]):
                ic = random.randint(0, numleft - 1)
                females[0][j]['mgtype'][locus_index] = gvec[ic]
                gvec[ic] = gvec[numleft - 1]
                numleft -= 1

            for j in range(bottleneck_individuals[i]):
                ic = random.randint(0, numleft - 1)
                males[0][j]['pgtype'][locus_index] = gvec[ic]
                gvec[ic] = gvec[numleft - 1]
                numleft -= 1

            for j in range(bottleneck_individuals[i]):
                ic = random.randint(0, numleft - 1)
                males[0][j]['mgtype'][locus_index] = gvec[ic]
                gvec[ic] = gvec[numleft - 1]
                numleft -= 1
#text
    #Simulate bottleneck generations from random mating
        # femalesNext = [{'pgtype': [0] * totalLoci, 'mgtype': [0] * totalLoci} for _ in range(bottleneck_individuals[i])]
        # malesNext = [{'pgtype': [0] * totalLoci, 'mgtype': [0] * totalLoci} for _ in range(bottleneck_individuals[i])]
        # parseBottleNeck = random.randint(neRange[0], neRange[1])
        current = 0
        next = 0
        for d in range(durationLength[0], durationLength[1]):
            assort(bottleneck_individuals[i], females[next], females[current], males[current], bottleneck_individuals[i], mutation_rate, totalLoci)
            assort(bottleneck_individuals[i], males[next], females[current], males[current], bottleneck_individuals[i], mutation_rate, totalLoci)
            temp = current
            current = next
            next = temp

    # Then, generate a set of genotypes for final generation
        assort(individualsDirectInput, finalIndividuals[i], females[current], males[current], bottleneck_individuals[i], mutation_rate, totalLoci)

    #Remove monomorphic loci if possible by replacing them with the extra loci
        for l in range(lociDirectInput):
            extraLociIndex = lociDirectInput
            if(variantOfFinalLocus(i, finalIndividuals, individualsDirectInput, l) == 2):
                continue
            while((extraLociIndex < extraProportionOfBufferLoci*lociDirectInput) and variantOfFinalLocus(i, finalIndividuals, individualsDirectInput,extraLociIndex) < 2):
                extraLociIndex += 1
            if((extraLociIndex < extraProportionOfBufferLoci*lociDirectInput)):
                break
            for k in range(individualsDirectInput):
                genotype1, genotype2 = loadFinalGenotype(i, finalIndividuals, k, extraLociIndex)
                storeFinalGenotype(i, finalIndividuals, k, j, genotype1, genotype2)
            extraLociIndex += 1

        write_output(finalIndividuals, lociDirectInput, individualsDirectInput)




lociDirectInput = 10  # Number of loci
thetaDirectInput = 0.0048  # theta
individualsDirectInput = 100  # Number of diploid individuals
minAlleleFrequencyDirectInput = 0.05  # Minimum allele frequency
formFlag = 0  # 0 for SNPs, 1 for microsatellites
trials = 10
durationLength = 2,8
neRange = 150, 250
mutation_rate = 0.000000012


simulate_population(lociDirectInput, neRange, thetaDirectInput, individualsDirectInput, minAlleleFrequencyDirectInput, formFlag, trials, mutation_rate, durationLength)

