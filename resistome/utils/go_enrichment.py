import networkx
import scipy.stats
from collections import defaultdict

__author__ = 'jdwinkler'


class GONetwork:

    def __init__(self, parent_child_tuples):

        """
        
        Generates a GONetwork object representing a gene ontology network.
        
        It is not actually specific to GO data, so other similar DAGs can be represented using this class as well
        
        :param parent_child_tuples: (parent, child) iterable representing an adjacency list
        """

        edges = []

        for (parent, child) in parent_child_tuples:
            edges.append((parent, child))

        self.G = networkx.DiGraph()
        self.G.add_edges_from(edges)

        # find roots for bio/molecular function/cellular component
        in_degree_nodes = self.G.in_degree()

        roots = []

        for node in in_degree_nodes:

            if in_degree_nodes[node] == 0:
                roots.append(node)

        self.roots = roots
        self.tag_root_distance = dict()

        for tag in self.G.nodes():
            self.tag_root_distance[tag] = self._compute_tag_root_distance(tag)

    def _compute_tag_root_distance(self, tag):

        """
        
        Computes tag distance from root nodes. Will only work for DAGs.
        
        :param tag: tag of interest
        :return: minimum distance between tag and root(s)
        """

        distances = []

        stack = []
        current_tag = tag
        stack.extend([(x, 1) for x in self.G.predecessors(current_tag)])

        while len(stack) != 0:

            temp_stack = []
            for (tag, d) in stack:

                predecessors = self.G.predecessors(tag)

                if len(predecessors) == 0:
                    distances.append(d)
                else:
                    temp_stack.extend([(x, d + 1) for x in predecessors])

            stack = temp_stack

        if len(distances) == 0:
            distances.append(0)

        return min(distances)

    def filter_by_minimum_root_distance(self, tags, threshold=0):

        if tags is None:
            return []

        remaining_tags = filter(lambda x: x[1] > threshold,
                                [(x, self.tag_root_distance[x]) for x in tags])

        return [x[0] for x in remaining_tags]


class EnrichmentCalculator:

    def __init__(self, gene_term_association):

        """
        
        Wrapper for performing fisher's exact tests using some underlying ontology.
        
        :param gene_term_association: (gene, label) iterable
        """

        self.go_to_gene = defaultdict(set)
        self.gene_to_go = defaultdict(set)

        self.go_tag_occurrences = defaultdict(int)

        for (gene, label) in gene_term_association:

            self.go_to_gene[label].add(gene)
            self.gene_to_go[gene].add(label)
            self.go_tag_occurrences[label] += 1

        self.total_tag_count = sum(self.go_tag_occurrences.values())

    def compute_enrichment(self, gene_list, p_value_filter=1.0, apply_multih_correction=True, feature_conversion=True):

        """
        
        Returns a list of (tag, p-values) that are considered significant. The p-value threshold will be divided by the
        length of gene_list as a simple Bonferroni-style multi-hypothesis correction.
        
        :param gene_list: list of tags (genes, etc) to test
        :param p_value_filter: minimum p-value to consider significant, default = 1.0/len(gene_list)
        :return: 
        """

        tag_counter = defaultdict(int)

        for gene in gene_list:

            if feature_conversion:
                go_tags = self.gene_to_go[gene]
            else:
                go_tags = [gene]

            for tag in go_tags:
                tag_counter[tag] += 1

        total_go_tags = sum(tag_counter.values())

        tag_pvalues = []

        counter = 0

        total_tag_count = self.total_tag_count
        go_tag_occurrences = self.go_tag_occurrences

        for tag, tag_value in tag_counter.items():

            contingency_table = [[tag_value, go_tag_occurrences[tag]],
                                 [total_go_tags - tag_value,
                                  total_tag_count - go_tag_occurrences[tag]]]

            values = [tag_value, go_tag_occurrences[tag], total_go_tags, total_tag_count]
            if all(x > 10 for x in values):
                # approx but good for large numbers. otherwise this is far too slow to use. SciPy docs suggest a min
                # of 5 counts in all cells.
                (odds_ratio, pvalue, dof, poptd) = scipy.stats.chi2_contingency(contingency_table)
            else:
                # use Fisher's exact at low counts.
                (odds_ratio, pvalue) = scipy.stats.fisher_exact(contingency_table)

            tag_pvalues.append((tag, pvalue, odds_ratio))

            counter += 1

        multi_h_correction = 1.0 if not apply_multih_correction else float(len(tag_pvalues))

        return list(filter(lambda x: x[1] <= p_value_filter / multi_h_correction, tag_pvalues))

    def compute_enrichment_specify_counts(self, tag_counter, p_value_filter=1.0):

        """
        
        Same as compute_enrichment, but with pre-counted tags.
        
        :param tag_counter: dict (str : int) representing tag_counts
        :param p_value_filter: minimum p-value to consider significant, default = 1.0/len(gene_list)
        :return: 
        """

        total_go_tags = sum(tag_counter.values())

        tag_pvalues = []

        for tag in tag_counter:
            contingency_table = [[tag_counter[tag], self.go_tag_occurrences[tag]],
                                 [total_go_tags - tag_counter[tag],
                                  self.total_tag_count - self.go_tag_occurrences[tag]]]

            values = [tag_counter[tag], self.go_tag_occurrences[tag], total_go_tags, self.total_tag_count]
            if all(x > 10 for x in values):
                # approx but good for large numbers. otherwise this is far too slow to use. SciPy docs suggest a min
                # of 5 counts in all cells.
                (odds_ratio, pvalue, dof, poptd) = scipy.stats.chi2_contingency(contingency_table)
            else:
                # use Fisher's exact at low counts.
                (odds_ratio, pvalue) = scipy.stats.fisher_exact(contingency_table)

            tag_pvalues.append((tag, pvalue, odds_ratio))

        return list(filter(lambda x: x[1] <= p_value_filter / float(len(tag_pvalues)), tag_pvalues))
