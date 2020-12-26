### RegulonDB Data

Database files for RegulonDB 10.6. Both the regulatory network and genomic features are extracted by the 
regulondb_parser.py script.

Relevant DB tables are interactions, genomic_features, and genomic_feature_association under each strain specific 
schema.

### Table Summary

As of December 2020, these files are:

* attenuator_terminator.txt: Anti-terminators
* conformation.txt: regulator active conformations (need to go from active regulator => protein)
* gene.txt: gene information
* gene_product_link.txt: gene to protein (protein/RNA)
* genetic_network.txt: regulatory relationships
* object_external_db_link.txt: cross-references to other DBs
* object_synonym.txt: synonyms for ECK accession
* object_properties.txt: a catch-all file for a lot of object properties, used to generate display pages; includes 
genes, sRNAs, operons, promoters, binding sites, etc
* product.txt: product (protein/RNA) information
* product_tf_link.txt: product => transcription factor relationships
* promoter.txt: promoter details
* promoter_feature.txt: promoter properties, included TATA box locations
* regulatory_interaction.txt: interactions between regulators and associated genetic elements
* rfam.txt: riboswitches
* shine_dalgarno.txt: SD boxes
* sigma_tmp.txt: Sigma factor details
* srna_interaction: short RNA => gene regulatory interactions
* site.txt: DNA binding site details
* terminators.txt: terminator properties, including rho (in)dependent
* tf_gene_interaction: transcription factor => gene interactions
* tu_gene_link.txt: transcription unit memberships

These relationships don't appear to be documented anywhere, so these files may change or disappear as RegulonDB is 
updated.