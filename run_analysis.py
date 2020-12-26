import AminoAcidAnalysis
import DescriptiveStatistics
import GenomeAnalysis
import ProteinAnalysis
import FluxAnalysis
from resistome.graphics import circos

"""

This script is meant to automate the routine analyses you might wish to run using the Resistome. Many of these 
figures have appeared in previous publications.

"""


def main():

    AminoAcidAnalysis.run_analysis()
    DescriptiveStatistics.run_analysis()
    GenomeAnalysis.run_analysis()
    ProteinAnalysis.run_analysis()
    FluxAnalysis.run_analysis()
    circos.generate_ecoli_config(output_prefix='gene_')


if __name__ == '__main__':
    main()
