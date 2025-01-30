import random
import sys
import time

'''
Initial population generation using coalescent model
'''

INT_MAX = sys.maxsize  # This provides the maximum integer size (similar to INT_MAX in C)
GFSR_STYPE_UNSIGNED_MAX = (2**64) - 1  # Maximum value for an unsigned 64-bit integer
random.seed(time.time())
def disrand(l, t):
    return random.randint(l, t)

def rand_bit():
    return random.randint(0, 1)

def gfsr4():
    return (random.randint(0, INT_MAX) / INT_MAX) / (GFSR_STYPE_UNSIGNED_MAX + 1.0)

def fallingQuotient(s, t1, t2, c):
    while c > 0:
        s *= t1 / t2
        t1 -= 1
        t2 -= 1
        c -= 1
    return s

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


def coalescent_population_generation(neRange, individualsDirectInput, bottleneck_individuals, bottleneck_length_choices, minAlleleCount, theta, totalLoci, num_genes):

    totalAllelePr = 0
    # Coalescent probability memo table for allele frequencies
    coalescentProbabilityMemoTable = [0] * (num_genes // 2 - minAlleleCount + 1)

    # Calculate coalescent probabilities
    for j in range(num_genes // 2 - minAlleleCount):
        coalescentProbabilityMemoTable[j] = allelePr(j + minAlleleCount, num_genes - j - minAlleleCount, theta)
        totalAllelePr += coalescentProbabilityMemoTable[j]

    # Normalize probabilities
    for j in range(num_genes // 2 - minAlleleCount):
        coalescentProbabilityMemoTable[j] /= totalAllelePr

    # # Initialize genotype vectors
    gvec = [0] * num_genes
    # # Store intermediate generations
    females = [{'pgtype': [0] * totalLoci, 'mgtype': [0] * totalLoci} for _ in range(neRange[1])]
    males = [{'pgtype': [0] * totalLoci, 'mgtype': [0] * totalLoci} for _ in range(neRange[1])]


    # Iterate over loci and simulate gene frequencies
    for locus_index in range(0, totalLoci):
        base1 = 2 * random.randint(0, 1) + random.randint(0, 1) + 1
        base2 = ((base1 + random.randint(0, 2)) % 4) + 1

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
        for j in range(bottleneck_individuals):
            ic = random.randint(0, numleft - 1)
            females[j]['pgtype'][locus_index] = gvec[ic]
            gvec[ic] = gvec[numleft - 1]
            numleft -= 1

        for j in range(bottleneck_individuals):
            ic = random.randint(0, numleft - 1)
            females[j]['mgtype'][locus_index] = gvec[ic]
            gvec[ic] = gvec[numleft - 1]
            numleft -= 1

        for j in range(bottleneck_individuals):
            ic = random.randint(0, numleft - 1)
            males[j]['pgtype'][locus_index] = gvec[ic]
            gvec[ic] = gvec[numleft - 1]
            numleft -= 1

        for j in range(bottleneck_individuals):
            ic = random.randint(0, numleft - 1)
            males[j]['mgtype'][locus_index] = gvec[ic]
            gvec[ic] = gvec[numleft - 1]
            numleft -= 1

    return females, males
