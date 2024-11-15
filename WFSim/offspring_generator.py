import random

'''
Simulation of offspring by random mating
Generates the next generation of individuals based on the current one
'''
def disrand(l, t):
    return random.randint(l, t)

def rand_bit():
    return random.randint(0, 1)

def branch_with_probability(p):
    return random.uniform(0.00000000, 0.0000001) < p

def mutate_snp(gene):
    gene[0] = 2 * rand_bit() + rand_bit() + 1
    return gene[0]


# Generates the next generation of individuals based on the current one
def assort(next_gen_count, offvec, mothers, fathers, current_gen_count, mutation_rate, num_loci):
    # p_mutation = mutation_rate
    # Simulate each individual's genotype
    for j in range(next_gen_count):
        m = disrand(0, current_gen_count - 1)
        d = disrand(0, current_gen_count - 1)
        # For each allele, select one from parent and possibly mutate
        for i in range(num_loci):
            offvec[j]['mgtype'][i] = mothers[m]['mgtype'][i] if rand_bit() else mothers[m]['pgtype'][i]
            if branch_with_probability(mutation_rate):
                offvec[j]['mgtype'][i] = mutate_snp(offvec[j]['mgtype'])

            offvec[j]['pgtype'][i] = fathers[d]['mgtype'][i] if rand_bit() else fathers[d]['pgtype'][i]
            if branch_with_probability(mutation_rate):
                offvec[j]['pgtype'][i] = mutate_snp(offvec[j]['pgtype'])
    return offvec
