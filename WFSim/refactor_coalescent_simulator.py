import random
import math

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

def parseBottleneck(individualsDirectInput):
    return individualsDirectInput // 2

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

# Main simulation function
def simulate_population(lociDirectInput, thetaDirectInput, individualsDirectInput, minAlleleFrequencyDirectInput, formFlag):
    # Initialize simulation variables
    num_genes = 4 * parseBottleneck(individualsDirectInput)  # total gene copies
    minAlleleCount = math.ceil(parseMinAlleleFrequency(minAlleleFrequencyDirectInput) * num_genes)
    totalAllelePr = 0

    # Coalescent probability memo table for allele frequencies
    coalescentProbabilityMemoTable = [0] * (num_genes // 2 - minAlleleCount + 1)

    # Calculate coalescent probabilities
    for j in range(num_genes // 2 - minAlleleCount + 1):
        coalescentProbabilityMemoTable[j] = allelePr(j + minAlleleCount, num_genes - j - minAlleleCount, thetaDirectInput)
        totalAllelePr += coalescentProbabilityMemoTable[j]

    # Normalize probabilities
    for j in range(num_genes // 2 - minAlleleCount + 1):
        coalescentProbabilityMemoTable[j] /= totalAllelePr

    # Initialize genotype vectors
    gvec = [0] * num_genes
    females = [{'pgtype': [0] * lociDirectInput, 'mgtype': [0] * lociDirectInput} for _ in range(parseBottleneck(individualsDirectInput))]
    males = [{'pgtype': [0] * lociDirectInput, 'mgtype': [0] * lociDirectInput} for _ in range(parseBottleneck(individualsDirectInput))]

    # Iterate over loci and simulate gene frequencies
    for locus_index in range(lociDirectInput):
        # Pick genotypes for SNPs or Microsatellites
        base1 = 2 * random.randint(0, 1) + random.randint(0, 1) + 1 if parseFormFlag(formFlag) == 0 else initializeMicrosat1(locus_index)
        base2 = (base1 + random.randint(0, 2)) % 4 + 1 if parseFormFlag(formFlag) == 0 else initializeMicrosat2(locus_index)

        # Simulate coalescent frequency distribution
        cut = random.random()
        for j in range(num_genes // 2 - minAlleleCount + 1):
            cut -= coalescentProbabilityMemoTable[j]
            if cut < 0:
                break

        # Fill genotype vector based on frequency distribution
        for count2 in range(num_genes):
            if count2 < j + minAlleleCount:
                gvec[count2] = base1
            else:
                gvec[count2] = base2

        numleft = num_genes

        # Assign alleles to females and males randomly
        for j in range(parseBottleneck(individualsDirectInput)):
            ic = random.randint(0, numleft - 1)
            females[j]['pgtype'][locus_index] = gvec[ic]
            gvec[ic] = gvec[numleft - 1]
            numleft -= 1

        for j in range(parseBottleneck(individualsDirectInput)):
            ic = random.randint(0, numleft - 1)
            females[j]['mgtype'][locus_index] = gvec[ic]
            gvec[ic] = gvec[numleft - 1]
            numleft -= 1

        for j in range(parseBottleneck(individualsDirectInput)):
            ic = random.randint(0, numleft - 1)
            males[j]['pgtype'][locus_index] = gvec[ic]
            gvec[ic] = gvec[numleft - 1]
            numleft -= 1

        for j in range(parseBottleneck(individualsDirectInput)):
            ic = random.randint(0, numleft - 1)
            males[j]['mgtype'][locus_index] = gvec[ic]
            gvec[ic] = gvec[numleft - 1]
            numleft -= 1

#text

    # Output simulated data
    print("Auto-generated genotype output.")
    for locus_index in range(lociDirectInput):
        print(f"{locus_index + 1}")
    print("Pop")

    for j in range(parseBottleneck(individualsDirectInput)):
        print(f"{j + 1}, ", end="")
        for locus_index in range(lociDirectInput):
            print(f"{females[j]['mgtype'][locus_index]:02d}{females[j]['pgtype'][locus_index]:02d} ", end="")
        print()

    for j in range(parseBottleneck(individualsDirectInput)):
        print(f"{parseBottleneck(individualsDirectInput) + j + 1}, ", end="")
        for locus_index in range(lociDirectInput):
            print(f"{males[j]['mgtype'][locus_index]:02d}{males[j]['pgtype'][locus_index]:02d} ", end="")
        print()





lociDirectInput = 10  # Number of loci
thetaDirectInput = 100  # theta
individualsDirectInput = 50  # Number of diploid individuals
minAlleleFrequencyDirectInput = 0.05  # Minimum allele frequency
formFlag = 0  # 0 for SNPs, 1 for microsatellites

simulate_population(lociDirectInput, thetaDirectInput, individualsDirectInput, minAlleleFrequencyDirectInput, formFlag)

