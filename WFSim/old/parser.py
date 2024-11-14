import sys
import re

class RefactorParser:
    def __init__(self):
        self.eof_error_string = "Unexpected EOF"
        self.comma_error_string = "Unexpected comma"
        self.genePop_error_string = "Genotype data does not conform to GenePop 4.0 standard: malformed number."
        self.blankLine_error_string = "GenePop 4.0 format does not allow blank lines."
        self.unexpected_input_string = "Unexpected data in input file: incorrect number of individuals or loci specified to program, or inconsistent number of alleles specified for each individual."
        self.invalid_SNP_string = "Invalid SNP: must be a number that is 0, 1, 2, 3, or 4."

        self.line_num = 1   # Current line in file
        self.char_num = 0   # Current character past last line
        self.next_char = None
        self.queue = [None] * 8  # Data stored in queue (length is QUEUE_LENGTH)
        self.queue_index = 0
        self.assert_eof = False
        self.parse_num_genes = 0
        self.individual = 0
        self.data_file = None
        self.mode_read = 'r'

    def parse(self, syntax_check_flag, data_file_name, syntax_results):
        with open(data_file_name, self.mode_read) as self.data_file:
            self.parse_from_file(syntax_check_flag, syntax_results)

    def parse_from_file(self, syntax_check_flag, syntax_results):
        self.initialize_parser_queue()
        self.next_char = self.data_file.read(1)

        syntax_results[0] = 0  # Number of loci per sample
        syntax_results[1] = 0  # Number of individuals

        # Parse the first line (skip)
        while self.next_char not in ['\n', '\r', '']:
            self.next_char = self.data_file.read(1)

        assert self.line_num == 1

        self.line_num += 1
        self.char_num = 0

        # Check for unexpected EOF.
        if self.next_char == '':
            self.report_parse_error(self.eof_error_string)

        assert self.line_num == 2

        # Parse the second line, count commas as proxy for number of genes until "POP"
        self.parse_num_genes = 0
        while not self.match_pop():
            if self.next_char == '':
                self.report_parse_error(self.eof_error_string)

            if self.next_char in [',', '\n', '\r']:
                self.parse_num_genes += 1

            self.accept_character()

        self.parse_num_genes -= 1

        # Ensure we have a reasonable number of genes.
        assert self.parse_num_genes > 0

        # Process genotype information
        locus_counter = 0
        while self.next_char != '':
            while self.next_char not in ['\n', '\r', '']:
                self.skip_past_comma()
                for locus_counter in range(self.parse_num_genes):
                    if self.next_char in [' ', '\t']:
                        self.accept_character()

                    if self.match_end_of_line():
                        break

                    while self.match_soft_delimiter() and self.next_char != '':
                        self.accept_character()

                    self.accept_character()
                    return_gene1, return_gene2 = self.parse_genotype_token()

                    if (return_gene1 > 4 or return_gene2 > 4) and self.parse_form_flag() == 0:
                        self.report_parse_error(self.invalid_SNP_string)

                    if self.next_char == '':
                        break

                    if not syntax_check_flag:
                        self.store_initial_genotype(self.individual, locus_counter, return_gene1, return_gene2)

                self.individual += 1

                while self.match_whitespace():
                    self.accept_character()

            syntax_results[1] = self.individual
            self.individual = 0
            self.accept_character()

        syntax_results[0] = self.parse_num_genes

    def initialize_parser_queue(self):
        self.queue_index = 0
        self.assert_eof = False
        self.parse_num_genes = 0

        # Initialize queue with EOF tokens
        self.queue = [None] * 8

        # Track file position
        self.line_num = 1
        self.char_num = 0

    def accept_character(self):
        self.next_char = self.data_file.read(1)
        self.enqueue_parser_token(self.next_char)

    def enqueue_parser_token(self, token):
        self.queue_index = (self.queue_index + 1) % len(self.queue)
        self.queue[self.queue_index] = token
        self.char_num += 1

        if self.assert_eof and token != '':
            self.report_parse_error(self.blankLine_error_string)

    def parse_genotype_token(self):
        # Using the last 8 tokens in the queue to extract the genotype
        prior_tokens = [self.access_prior_element(i) for i in range(8)]

        # Check for legal 4-digit number
        if (self.match_soft_delimiter(prior_tokens[0]) and all(self.match_digit(prior_tokens[i]) for i in range(1, 5))
                and self.match_soft_delimiter(prior_tokens[5])):
            return_gene1 = 10 * (int(prior_tokens[4]) - ord('0')) + (int(prior_tokens[3]) - ord('0'))
            return_gene2 = 10 * (int(prior_tokens[2]) - ord('0')) + (int(prior_tokens[1]) - ord('0'))
            return return_gene1, return_gene2

        # Check for legal 6-digit number
        if (self.match_soft_delimiter(prior_tokens[0]) and all(self.match_digit(prior_tokens[i]) for i in range(1, 7))
                and self.match_soft_delimiter(prior_tokens[7])):
            return_gene1 = 100 * (int(prior_tokens[6]) - ord('0')) + 10 * (int(prior_tokens[5]) - ord('0')) + (int(prior_tokens[4]) - ord('0'))
            return_gene2 = 100 * (int(prior_tokens[3]) - ord('0')) + 10 * (int(prior_tokens[2]) - ord('0')) + (int(prior_tokens[1]) - ord('0'))
            return return_gene1, return_gene2

        # Otherwise, data is malformed
        self.report_parse_error(self.genePop_error_string)

    def access_prior_element(self, ago):
        assert ago >= 0 and ago < len(self.queue)
        return self.queue[(self.queue_index + (len(self.queue) - ago)) % len(self.queue)]

    def match_pop(self):
        return (self.match_end_of_line(self.access_prior_element(0)) and
                self.match_char(self.access_prior_element(1), 'P') and
                self.match_char(self.access_prior_element(2), 'O') and
                self.match_char(self.access_prior_element(3), 'P') and
                self.match_end_of_line(self.access_prior_element(4)))

    def match_end_of_line(self, c=None):
        c = c if c is not None else self.next_char
        return c in ['\n', '\r']

    def match_whitespace(self, c=None):
        c = c if c is not None else self.next_char
        return c in [' ', '\t', '\n', '\r']

    def match_soft_delimiter(self, c=None):
        return self.match_whitespace(c) or c in [',']

    def match_digit(self, c):
        return c in '0123456789'

    def match_char(self, c, target):
        return c.lower() == target.lower()

    def skip_past_comma(self):
        while self.next_char != ',' and self.next_char != '':
            self.accept_character()

    def report_parse_error(self, message):
        print(f"Parse Error: {message} at line {self.line_num}, character {self.char_num}")
        sys.exit(1)

    def parse_form_flag(self):
        # Placeholder for `parseFormFlag()` implementation
        return 0

    def store_initial_genotype(self, individual, locus_counter, gene1, gene2):
        # Placeholder for storing genotype
        pass

    def pair_distance(self, geneA1, geneA2, geneB1, geneB2):
        if geneA1 > geneA2:
            return self.pair_distance(geneA2, geneA1, geneB1, geneB2)
        if geneB1 > geneB2:
            return self.pair_distance(geneA1, geneA2, geneB2, geneB1)
        if geneA2 == 0 or geneB2 == 0:
            return 2
        if geneA1 == 0:
            return 1 if geneA2 in [geneB1, geneB2] else 2
        if geneB1 == 0:
            return 1 if geneB2 in [geneA1, geneA2] else 2
        return (0 if geneA1 == geneB1 else 1) + (0 if geneA2 == geneB2 else 1)

    def hamming_distance(self, indivI, indivJ):
        result = 0
        for k in range(self.num_loci):
            geneA1, geneA2 = self.load_initial_genotype(indivI, k)
            geneB1, geneB2 = self.load_initial_genotype(indivJ, k)
            result += self.pair_distance(geneA1, geneA2, geneB1, geneB2)
        return result

    def fill_in_missing_data(self):
        distances = [[self.hamming_distance(i, j) for j in range(self.num_individuals)]
                     for i in range(self.num_individuals)]

        for i in range(self.num_individuals):
            for j in range(self.num_loci):
                geneA1, geneA2 = self.load_initial_genotype(i, j)
                if geneA1 == 0 or geneA2 == 0:
                    smallest_distance = 2 * self.num_loci
                    geneB1, geneB2 = 0, 0

                    # Find smallest distance
                    for k in range(self.num_individuals):
                        if i == k:
                            continue
                        geneB1, geneB2 = self.load_initial_genotype(k, j)
                        if geneB1 == 0 or geneB2 == 0:
                            continue
                        if distances[i][k] < smallest_distance:
                            smallest_distance = distances[i][k]

                    # Select most frequently occurring nearest neighbor
                    frequencies = [0] * self.num_individuals
                    for k in range(self.num_individuals):
                        if distances[i][k] == smallest_distance:
                            if self.pair_distance(geneA1, geneA2, geneB1, geneB2) == 0:
                                frequencies[k] += 1

                    # Find the most frequent neighbor
                    max_index = frequencies.index(max(frequencies))

                    # Fill in missing data
                    geneB1, geneB2 = self.load_initial_genotype(max_index, j)
                    self.store_initial_genotype(i, j, geneB1, geneB2)

    def delete_locus(self, locusI):
        for i in range(self.num_individuals):
            geneA1, geneA2 = self.load_initial_genotype(i, self.num_loci - 1)
            self.store_initial_genotype(i, locusI, geneA1, geneA2)

        if self.parse_form_flag():
            self.motif_lengths[locusI] = self.motif_lengths[self.num_loci - 1]

    def delete_individual(self, indivI):
        for i in range(self.num_loci):
            geneA1, geneA2 = self.load_initial_genotype(self.num_individuals - 1, i)
            self.store_initial_genotype(indivI, i, geneA1, geneA2)

    def locus_coverage(self, locusI):
        total_individuals_without_missing_data = sum(
            1 for i in range(self.num_individuals)
            if self.load_initial_genotype(i, locusI)[0] != 0 and self.load_initial_genotype(i, locusI)[1] != 0
        )
        return total_individuals_without_missing_data

    def variants_of_locus(self, locusI):
        non_zero_variant = 0
        for i in range(self.num_individuals):
            geneA1, geneA2 = self.load_initial_genotype(i, locusI)
            if non_zero_variant == 0:
                non_zero_variant = geneA1 or geneA2
            if non_zero_variant == 0:
                continue
            if non_zero_variant != geneA1 or non_zero_variant != geneA2:
                return 2
        return 1 if non_zero_variant != 0 else 0

    def indiv_coverage(self, indivI):
        total_loci_without_missing_data = sum(
            1 for i in range(self.num_loci)
            if self.load_initial_genotype(indivI, i)[0] != 0 and self.load_initial_genotype(indivI, i)[1] != 0
        )
        return total_loci_without_missing_data

    def filter_low_coverage_loci(self, omit_threshold):
        i = 0
        while i < self.num_loci:
            if self.locus_coverage(i) < omit_threshold * self.num_individuals:
                self.delete_locus(i)
                self.num_loci -= 1
            else:
                i += 1

    def filter_low_coverage_individuals(self, omit_threshold):
        i = 0
        while i < self.num_individuals:
            if self.indiv_coverage(i) < omit_threshold * self.num_loci:
                self.delete_individual(i)
                self.num_individuals -= 1
            else:
                i += 1

    def filter_monomorphic_loci(self):
        i = 0
        while i < self.num_loci:
            if self.variants_of_locus(i) < 2:
                self.delete_locus(i)
                self.num_loci -= 1
            else:
                i += 1

    def obtain_missing_proportion(self):
        total_coverage = sum(self.indiv_coverage(j) for j in range(self.num_individuals))
        proportion_missing_data = total_coverage / (self.num_loci * self.num_individuals)
        self.set_proportion_missing_data(proportion_missing_data)

    def load_initial_genotype(self, individual, locus):
        # Placeholder for function that loads initial genotype
        return (0, 0)

    def store_initial_genotype(self, individual, locus, gene1, gene2):
        # Placeholder for function that stores the initial genotype
        pass

    def parse_form_flag(self):
        # Placeholder for function that determines the data format
        return False

    def set_proportion_missing_data(self, value):
        # Placeholder for setting proportion of missing data
        pass

