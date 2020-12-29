### *E. coli* Resistome Database Repository

Welcome to the Resistome database repository! The code here is meant to minimally generate the 
Resistome database if the required inputs are present, as well as run basic analyses that you 
might be interested in using in your own research. It is not really meant to be used as a library 
but might be useful nonetheless.

The public Resistome website is [here](https://resistome-web-interface.herokuapp.com/) if you want to quickly 
search for high-level data concerning genotype-phenotype relationships.

### About

The Resistome contains standardized representations of *E. coli* mutants resistant or sensitized to over 
500 hundred types of inhibition (solvents, environmental stresses, antibiotics, and many others). This unique 
data source is intended to help synthetic biologists and evolutionary engineers to identify loci likely to affect 
their phenotype of interest for forward genetic engineering. Ultimately, we will hopefully be able to design explainable 
machine learning approaches to predict genotypes required to confer desired phenotypes using these curated data, or 
to improve the effectiveness of (targeted) adaptive laboratory evolution. Our curated dataset contains information on 
\>10,000 *E. coli* mutants, including the specific variants detected in the strains from the library or selection 
experiments carried out in more than 440 studies (as of December 2020).

This database does not include information concerning phenotypes endowed by foreign genetic elements (either libraries 
or mobilizable elements/plasmids), as other databases such as CARD will cover those in more detail. A related database 
focused solely on ALE was recently published by the Palsson laboratory at [aledb.org](aledb.org), in case it is also 
useful for your research.

### Quick Setup

Here is the quick setup procedure:

0. Setup a python >=3.8 (or equivalent) virtual environment. 
1. Install the python requirements using `pip install -r requirements.txt`.
2. Install and setup Postgres (if needed).
2. Download the latest database dump on the public Resistome website 
([dump](https://resistome-web-interface.herokuapp.com/static/resistome-202012212206.sql), 
[summary page](https://resistome-web-interface.herokuapp.com/summary)).
3. Restore the custom format SQL dump using `pg_restore` into your target database.
4. Adjust the credentials in /inputs/db_credentials/credentials.txt to match your target database.
5. Resist away!

You can see below for more details concerning the expected inputs and how to manually build the database. This repo 
will usually have a more updated database by virtue of being easier to deploy.

### Required Inputs/Infrastructure

#### Database Server
The Resistome assumes a Postgres server is running on the localhost at port 5432. The default DB name 
is `resistome` with username and password defined in [here](inputs/db_credentials/credentials.txt). 
These constants can be changed if desired.

#### Database Schema

The annotation, Resistome, and "extra" schemas are stored under inputs/sql. An ER diagram can be generated using most 
available database management tools, such as DBeaver or pgAdmin. The supporting annotation databases are structured 
to maximize referential integrity but there are no cross schema (public-resistome) constraints.

#### Curated Data

The curated study input files are stored under `inputs/database_store`. The data are stored in a custom format that is 
hard to read/modify, but it should work well enough for analysis. Ideally, we will switch over to a combination of XML 
typed records with a validating schema after repeating the resequencing analysis for the relevant studies, but 
Resistome development is currently unfunded.

These files are generated using a parser like `resistome/examples/record_generator.py`. Every study essentially 
requires a custom parser, so there is little value in code re-use. In workflow management systems, you should try 
to automate the read preparation => variant calling => output generation using a more standard format like VCF, 
GenomeDiff, etc using a common reference.

#### NCBI/Uniprot Inputs

NCBI and Uniprot annotations are used to populate the annotation tables used by the Resistome to disambiguated curated 
datasets. These strains are:

* W3110: GCF_000010245.2 (RefSeq) + the Genbank file for the Genbank assembly (GCA_000010245.1)
* MDS42: GCA_000350185.1 (Genbank)
* BW25113: GCA_000750555.1 (GenBank)
* BL21(DE3): GCA_000022665.2 (GenBank)
* BL21: GCA_013166975.1 (GenBank)
* MG1655: GCA_000005845.2 (Genbank)
* REL606: GCA_000017985.1 (Genbank)
* W: GCA_000184185.1 (Genbank)

The Genbank feature files (GBFF), RNA (rna_from_genomic), protein, CDS, genome, and feature tables are required for 
each strain. See [this readme file for more informaton](inputs/biocyc/readme.md). The Uniprot table defining the K-12
proteome is also used for K-12 derivatives strains help disambiguate strain gene names further.

We previously utilized Biocyc databases to provide this same information for *E. coli* MG!655, REL606, and W strains.

#### Standardization

Gene names, phenotypes, and compound names are standardized to enable qualitative cross comparison between experiments. 
The files defining these mappings are included under inputs/standardization and inputs/settings. If you are adding more 
data with new phenotypes/compounds mentioned, you will need to update the files mentioned in the standardization 
[READ ME](inputs/standardization/readme.md). An error will be thrown if you attempt add a paper to the database that 
contains unknown phenotypes.

Gene name mappings can be generated using bidirectional best hits between *E. coli* strains if mappings are not 
publicly available. For *E. coli* strains, it is generally possible to determine the mapping from pre-computed annotations 
as strains will be tagged with the equivalent MG1655 gene (if present) by PGAP or other public annotation pipelines. A 
future update will examine adding a pipeline to handle *de novo* mapping generation using BLAST. 

#### Protein-Protein Interactions

We use the [EcoliNet](https://www.inetbio.org/ecolinet/) functional gene-gene interaction network to add some additional 
predictive power for understanding how genotypes translate into phenotypes. However, this approach has not been extensively 
tested so you may wish to look at more detailed datasets.

#### Protein Change Effect Prediction Inputs

SNAP2, INPS, and DeMaSk are used to predict the effect of amino acid substitutions on protein function. These inputs 
were provided by external collaborators but are too large to distribute in this repository. You can download the original 
SNAP2 datasets generated in 2016 [here](https://zenodo.org/record/4394374); a future update will re-run these analyses 
for the strains currently represented in the Resistome. DeMaSk predictions for all strains can be found 
[here](10.5281/zenodo.4399936). See [the README](inputs/protein_stability/readme.md) for more details.

Once (if) AlphaFold2 becomes generally available, it may be possible to include assessments of structural impacts 
directly.

Citations:

* SNAP2: Better prediction of functional effects for sequence variants", BMC Genomics (2015) 
[DOI](10.1093/bioinformatics/btw192)
* INPS: INPS-MD: a web server to predict stability of protein variants from sequence and structure, Bioinformatics 
(2016) [DOI](10.1186/1471-2164-16-S8-S1)
* DeMaSk: DeMaSk: a deep mutational scanning substitution matrix and its use for variant impact prediction, 
Binformatics (2020) [DOI](10.1093/bioinformatics/btaa1030)

#### RegulonDB Extraction

Currently, regulatory interactions are extracted from the provided NCBI databases and 
[RegulonDB](http://regulondb.ccg.unam.mx/). See [the README](inputs/regulondb/readme.md) for more information. 
Genome annotations for promoters, operons, terminators, and DNA binding sites are also extracted from the database tables
provided by the maintainers.

#### Metabolic Models

Both iJO1366 and iML1515 are included for use in simulating the metabolic impacts of mutations.

### Full Setup

Assuming you have all the required data and have setup your postgres server, you are now ready to build the 
database locally. From the repository root, you should be able to run the following command:

```python
python3 install_db.py
```
to start the database construction process. This script will first construct the NCBI support tables, followed by 
the Resistome table. See the `resistome/sql/ncbi_data_parser.py` and `resistome/sql/resistome_builder.py` for the 
actual work of constructing the public and Resistome tables respectively. Name mappings are constructed using the 
`resistome/utils/name_helper.py` script to help disambiguate gene names in curated data. A basic validation of uploaded 
data is performed by the `resistome/sql/validator.py` script.

You can build each database (support/Resistome) separately, but if you rebuild the NCBI tables, you should rebuild 
the Resistome as well. The build process should require no more than 10-15 minutes on a modern laptop. The design of 
the database tries to minimize inconsistencies arising from curation and reporting errors, but if you notice any 
problems or oddities in the data, please contact us.

### Curation Accuracy

As of December 2020, ~95.5% of variant calls in the Resistome pass validation (e.g. WT bases or residues are correct, 
genomic locations are valid, genes are matched to accessions, etc). See [resistome/sql/validator.py](resistome/sql/validator.py) 
for more information. If genes involved in large (deletion, amplification, inversion) mutations are included, then 
\>99% of annotations appear to be correct.

### Adding New Data

It is usually pretty straightforward to add data to the Resistome by using the `resistome/examples/record_generator.py` 
script. You will need to alter this script to extract genetic or transcriptional data from the study or 
studies of interest depending on how they format their data. Mutation types must be only those found in 
[inputs/settings/FK_InternalFields.txt](inputs/settings/FK_InternalFields.txt) entries with a "mutation" tag.

Unfortunately the annotation style can get pretty complicated, and there is no *de jure* definition of the annotation 
grammar. However, I suggest searching through the existing database_store files to see formats for each mutation type in 
the meantime. Some common ones:

* nuc_snps: nucleotide SNPs, (REF base)(location)(ALT base): A123456T
* aa_snps: amino acid residue changes, (REF residue)(protein location)(ALT residue): M1V
* indel: insertion or deletion, (location)|(indel size)|(type = absolute, relative): 123456|-200|absolute
* is_insertion: name of IS or transposon insertion source into sequence, e.g. IS1/Mariner/Tn5/etc
* frameshift: frameshift mutation, no annotation required
* del: complete gene deletion, no annotation required
* rep: repression, free form text comment (e.g. CRISPRi)
* oe: overexpression, free form text comment (usually promoter name like tac, T7, etc)
* plasmid: plasmid expression, free form text comment (usually copy number)
* large_\[deletion/inversion\]: first gene|last gene
* large_amplification: first gene|last gene|fold amplification

The `resistome_builder.py` script will do as much QA checking as possible to make sure your entries match the expected 
format, but this checking will not be foolproof. You should also try to use absolute genomic coordinates if possible to 
simplify data analysis; at the moment, the positions in the following strains are checked: BW25113, BL21, BL21(DE3), 
MG1655, REL606, W, W3110, and MDS42. Other strains will be added as they become used for ALE or library experiments 
in the future.

Note: large mutations are automatically expanded into the complete mutation set: all genes between gene 1...gene N 
are included if they can be found in the supporting databases. Otherwise, only the first/last genes are associated 
with the mutations.

### Web Interface 

The code for the web interface hosted on the public Resistome website [here](https://resistome-web-interface.herokuapp.com/)
can be found on [Github](https://github.com/jdwinkler/resistome_web_interface). It will work with the database 
dump present in the repo above but may require updates to work with this repo directly as the database schema 
changes over time.

### Basic Analysis

You can run the `run_analysis.py` script to generate many of the figures that have previously appeared in Resistome 
publications. All databases accesses go through the "resistome/sql/sql_interface.py" class if you need ideas on 
how to access the data easily from Python without worrying about writing raw SQL to do so.

### Help/Collaborations

Please file an issue describing your problem with any error output you are getting. The Resistome is very much a 
research product, so expect rough edges unfortunately. You can also email [James Winkler](mailto:james.winkler@gmail.com) 
if you need assistance or would like to collaborate on a research project.

### Resistome Publications

1. Winkler, JD et al. "The Resistome: A Comprehensive Database of Escherichia coli Resistance Phenotypes", 
ACS Synthetic Biology (2016) [DOI](https://doi.org/10.1021/acssynbio.6b00150).
2. Erickson, KE, Winkler, JD et al. "The Tolerome: A Database of Transcriptome-Level Contributions to Diverse Escherichia coli Resistance and Tolerance Phenotypes
", ACS Synthetic Biology (2017)  [DOI](https://doi.org/10.1021/acssynbio.7b00235).
3. Winkler, JD. "The Resistome: updating a standardized resource for analyzing resistance phenotypes", BiorXiv (2018) 
[DOI](https://doi.org/10.1101/418814)

### License

This work is licensed under CC BY-NC 4.0: see [this webpage](https://creativecommons.org/licenses/by-nc/4.0/) 
for additional information. RegulonDB has its own separate license and can only be used for 
academic/non-commercial research. Please see [here](http://regulondb.ccg.unam.mx/menu/download/full_version/terms_and_conditions.jsp) 
for the full license terms.



