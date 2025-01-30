import random
import math
import argparse
import sys
import time
from offspring_generator import assort
from coalescent_simulator import coalescent_population_generation
'''
WFSim generates SNP data in Wright-Fisher populations, forward in time. 
Supports customizable parameters such as Effective population size, number of loci, mutation rate, theta, bottleneck duration etc
Output is generated in GENEPOP format 

'''
start_time = time.time()
'''
CONSTANTS ETC ETC
'''
INT_MAX = sys.maxsize  # This provides the maximum integer size (similar to INT_MAX in C)
GFSR_STYPE_UNSIGNED_MAX = (2**64) - 1  # Maximum value for an unsigned 64-bit integer
random.seed(time.time())
extraProportionOfBufferLoci = 2


parser = argparse.ArgumentParser()
parser.add_argument("--m", type=float, help="Minimum Allele Frequency")
parser.add_argument("--r", type=float, help="Mutation Rate")
parser.add_argument("--lNe", type=int, help="Lower of Ne Range")
parser.add_argument("--uNe", type=int, help="Upper of Ne Range")
parser.add_argument("--lT", type=float, help="Lower of Theta Range")
parser.add_argument("--uT", type=float, help="Upper of Theta Range")
parser.add_argument("--lD", type=float, help="Lower of Duration Range")
parser.add_argument("--uD", type=float, help="Upper of Duration Range")
parser.add_argument("--i", type=int, help="number of individuals")
parser.add_argument("--l", type=int, help="number of loci")
parser.add_argument("--o", type=str, help="The File Name")
parser.add_argument("--n", type=bool, help="whether to filter the monomorphic loci", default=False)
args = parser.parse_args()



'''
INITIALIZING PARAMETERS
'''
minAlleleFreq = 0.05
if (args.m):
    minAlleleFreq = float(args.m)

mutationRate = 0.000000012
if (args.r):
    mutationRate = float(args.r)

lowerNe = 150
if (args.lNe):
    lowerNe = int(args.lNe)

upperNe = 250
if (args.uNe):
    upperNe = int(args.uNe)

if (int(lowerNe) > int(upperNe)):
    print("ERROR:main:lowerNe > upperNe. Fatal Error")
    exit()

if (int(lowerNe) < 1):
    print("ERROR:main:lowerNe must be a positive value. Fatal Error")
    exit()

if (int(upperNe) < 1):
    print("ERROR:main:upperNe must be a positive value. Fatal Error")
    exit()

rangeNe = (lowerNe, upperNe)

lowerTheta = 0.000048
if (args.lT):
    lowerTheta = float(args.lT)

upperTheta = 0.0048
if (args.uT):
    upperTheta = float(args.uT)

rangeTheta = (lowerTheta, upperTheta)

lowerDuration = 2
if (args.lD):
    lowerDuration = float(args.lD)

upperDuration = 8
if (args.uD):
    upperDuration = float(args.uD)

rangeDuration = (lowerDuration, upperDuration)

lociDirectInput = 40  # Number of loci
individualsDirectInput = 10  # Number of diploid individuals


def parseNLoci(lociDirectInput):
    return lociDirectInput

def parseTheta(thetaRange):
    return random.uniform(thetaRange[0], thetaRange[1])

def bottleneck_individual_count(neRange):
    return random.randint(neRange[0], neRange[1])

def bottleneck_length_random_choices(duration):
    return random.randint(duration[0], duration[1])

def parseMinAlleleFrequency(minAlleleFrequencyDirectInput):
    return minAlleleFrequencyDirectInput

def disrand(l, t):
    return random.randint(l, t)

def rand_bit():
    return random.randint(0, 1)


def variantOfFinalLocus(finalIndividuals, individualsDirectInput, locusI):
    non_zero_variant = 0
    for i in range(individualsDirectInput):
        geneA1, geneA2 = loadFinalGenotype(finalIndividuals, i, locusI)
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


def loadFinalGenotype(finalIndividuals, individual, index):
    genotype1 = finalIndividuals[individual]['pgtype'][index]
    genotype2 = finalIndividuals[individual]['mgtype'][index]
    return genotype1, genotype2

def storeFinalGenotype(finalIndividuals, individual, index, genotype1, genotype2):
    finalIndividuals[individual]['pgtype'][index] = genotype1
    finalIndividuals[individual]['mgtype'][index] = genotype2
    return finalIndividuals


def filter_monomorphic_loci(finalIndividuals, lociDirectInput, individualsDirectInput):
    #Remove monomorphic loci if possible by replacing them with the extra loci
    for l in range(lociDirectInput):
        extraLociIndex = lociDirectInput
        if(variantOfFinalLocus(finalIndividuals, individualsDirectInput, l) == 2):
            continue
        while((extraLociIndex < extraProportionOfBufferLoci*lociDirectInput) and variantOfFinalLocus(finalIndividuals, individualsDirectInput,extraLociIndex) < 2):
            extraLociIndex += 1
        if((extraLociIndex < extraProportionOfBufferLoci*lociDirectInput)):
            break
        for k in range(individualsDirectInput):
            genotype1, genotype2 = loadFinalGenotype(finalIndividuals, k, extraLociIndex)
            storeFinalGenotype(finalIndividuals, k, j, genotype1, genotype2)
        extraLociIndex += 1
    return finalIndividuals


def write_output(finalIndividuals, lociDirectInput, individualsDirectInput):
    # Output simulated data
    print("Auto-generated genotype output.")
    for locus_index in range(lociDirectInput):
        print(f"{locus_index + 1}")
    print("Pop")

    for j in range(individualsDirectInput):
        print(f"{j + 1}, ", end="")
        for locus_index in range(lociDirectInput):
            print(f"{finalIndividuals[j]['mgtype'][locus_index]:02d}{finalIndividuals[j]['pgtype'][locus_index]:02d} ", end="")
        print()



'''
Main Simulation function
'''
def simulate_population(lociDirectInput, neRange, rangeTheta, individualsDirectInput, minAlleleFrequencyDirectInput, mutation_rate, durationLength):
    bottleneck_individuals = bottleneck_individual_count(neRange)
    bottleneck_length_choices = bottleneck_length_random_choices(durationLength)
    theta = parseTheta(rangeTheta)
    num_genes = 4 * bottleneck_individuals  # total gene copies
    minAlleleCount = math.ceil(parseMinAlleleFrequency(minAlleleFrequencyDirectInput) * num_genes)
    totalLoci = extraProportionOfBufferLoci * lociDirectInput
    minAlleleCount = math.ceil(parseMinAlleleFrequency(minAlleleFrequencyDirectInput) * num_genes)
    females, males = coalescent_population_generation(neRange, individualsDirectInput, bottleneck_individuals, bottleneck_length_choices, minAlleleCount, theta, totalLoci, num_genes)

    #Simulate bottleneck generations from random mating
    femalesNext = [{'pgtype': [0] * totalLoci, 'mgtype': [0] * totalLoci} for _ in range(bottleneck_individuals)]
    malesNext = [{'pgtype': [0] * totalLoci, 'mgtype': [0] * totalLoci} for _ in range(bottleneck_individuals)]
    #Store Final Generation
    finalIndividuals = [{'pgtype': [0] * totalLoci, 'mgtype': [0] * totalLoci} for _ in range(individualsDirectInput)]

    for d in range(bottleneck_length_choices):
        assort(bottleneck_individuals, femalesNext, females, males, bottleneck_individuals, mutation_rate, totalLoci)
        assort(bottleneck_individuals, malesNext, females, males, bottleneck_individuals, mutation_rate, totalLoci)
        females = femalesNext
        males = malesNext
    # Then, generate a set of genotypes for final generation
    assort(individualsDirectInput, finalIndividuals, femalesNext, malesNext, bottleneck_individuals, mutation_rate, totalLoci)
    finalIndividuals = filter_monomorphic_loci(finalIndividuals, lociDirectInput, individualsDirectInput)
    write_output(finalIndividuals, lociDirectInput, individualsDirectInput)



# simulate_population(lociDirectInput, rangeNe, rangeTheta, individualsDirectInput, minAlleleFreq, mutationRate, rangeDuration)
simulate_population(args.l, (args.lNe, args.uNe), (args.lT, args.uT), args.i, args.m, args.r, (args.lD, args.uD))

print("----- %s seconds -----" % (time.time() - start_time))
