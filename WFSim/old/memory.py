 import numpy as np

class OneSampMemoryManager:
    def __init__(self):
        self.initial_indivs_data = None
        self.final_indivs_data = None
        self.number_of_alleles = None
        self.double_data = None
        self.gtype = None
        self.gcount = None

    def allocate_struct1(self, num_samples, num_loci_allocation):
        """Allocates arrays to store number of alleles."""
        self.number_of_alleles = np.zeros((num_samples, num_loci_allocation), dtype=int)

    def allocate_struct2(self, num_samples, num_loci_allocation):
        """Allocates arrays to track number of alleles at each locus."""
        self.gtype = [[[np.zeros(MAX_NO_ALLELES, dtype=int) for _ in range(num_loci_allocation)]
                       for _ in range(num_samples)]]

    def allocate_struct3(self, num_samples, num_loci_allocation):
        """Allocates arrays to track frequency of alleles at each locus."""
        self.gcount = [[[np.zeros(MAX_NO_ALLELES, dtype=int) for _ in range(num_loci_allocation)]
                        for _ in range(num_samples)]]

    def allocate_struct4(self, initial_indivs_count_allocation, num_loci_allocation):
        """Allocates arrays to store initial data."""
        self.initial_indivs_data = [{
            'pgtype': np.zeros(num_loci_allocation, dtype=int),
            'mgtype': np.zeros(num_loci_allocation, dtype=int)
        } for _ in range(initial_indivs_count_allocation)]

    def allocate_struct6(self, num_samples):
        """Allocates arrays to store generation of statistics."""
        self.double_data = np.zeros(11 * num_samples, dtype=float)

    def allocate_one_samp_memory(self, initial_indivs_count_allocation, bottleneck_indivs_count, final_indivs_count,
                                 num_samples, num_loci_allocation):
        """Invokes calls to allocate all global data structures."""
        self.allocate_struct1(num_samples, num_loci_allocation)
        self.allocate_struct2(num_samples, num_loci_allocation)
        self.allocate_struct3(num_samples, num_loci_allocation)
        self.allocate_struct4(initial_indivs_count_allocation, num_loci_allocation)
        self.allocate_struct6(num_samples)

    def deallocate_one_samp_memory(self):
        """Deallocate all data structures."""
        # In Python, memory deallocation is automatic through garbage collection,
        # but we can delete the references to these arrays to free up memory.
        del self.double_data
        del self.initial_indivs_data
        del self.gcount
        del self.gtype
        del self.number_of_alleles

    def load_initial_genotype(self, individual, index):
        """Loads an initial genotype from memory."""
        genotype1 = self.initial_indivs_data[individual]['pgtype'][index]
        genotype2 = self.initial_indivs_data[individual]['mgtype'][index]
        return genotype1, genotype2

    def store_initial_genotype(self, individual, index, genotype1, genotype2):
        """Stores an initial genotype to memory."""
        self.initial_indivs_data[individual]['pgtype'][index] = genotype1
        self.initial_indivs_data[individual]['mgtype'][index] = genotype2

    def load_final_genotype(self, sample, individual, index):
        """Loads a genotype from the final sample."""
        genotype1 = self.final_indivs_data[sample][individual]['pgtype'][index]
        genotype2 = self.final_indivs_data[sample][individual]['mgtype'][index]
        return genotype1, genotype2

    def store_final_genotype(self, sample, individual, index, genotype1, genotype2):
        """Stores a genotype to the final sample."""
        self.final_indivs_data[sample][individual]['pgtype'][index] = genotype1
        self.final_indivs_data[sample][individual]['mgtype'][index] = genotype2


# Constants (These values would need to be defined properly for the actual use case)
MAX_NO_ALLELES = 10
STRUCT_GTYPE_SIZE = 10  # Placeholder for the actual size; depends on the implementation
STRUCT_GTYPE_STAR_SIZE = 10  # Placeholder for the actual size; depends on the implementation

# Example usage
manager = OneSampMemoryManager()
manager.allocate_one_samp_memory(initial_indivs_count_allocation=5, bottleneck_indivs_count=10,
                                 final_indivs_count=15, num_samples=20, num_loci_allocation=30)

# Example of loading and storing genotypes
genotype1, genotype2 = manager.load_initial_genotype(individual=0, index=0)
manager.store_initial_genotype(individual=0, index=0, genotype1=genotype1, genotype2=genotype2)

# Deallocate memory
manager.deallocate_one_samp_memory()
