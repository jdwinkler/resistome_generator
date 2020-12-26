create schema resistome

-- note: in many places I reference the term explanation table using a foreign key constraint.
-- you usually want to match int/bigint => primary key (int/bigint), not string => string.
-- postgres allows you to do it this way, but it is slower.
-- however, since the Resistome is a small DB, speed is not so important compared to enforcing referential integrity.
-- the terms in term_explanation are now defined in the FK_InternalFields.txt file under settings, so if you want
-- to add any mutation types/phenotypes/etc, you will have to add them there.

-- However, given the table design some of the constraints that should be referential (phenotype names, etc) are
-- enforced by the resistome uploader. it will take some effort to fix this unfortunately. For another day...

create table resistome.term_explanation (

    term_id serial primary key,
	term_type	text not null,
	internal_name	text not null,
	explanation	text not null,

	constraint unique_term_pairs UNIQUE (term_type, internal_name),
	constraint unique_internal_name UNIQUE (internal_name)

);

create table resistome.phenotype_standardization (

    phenotype_id    serial primary key,
	-- name	text not null,
	standard_name	text unique not null,
	phenotype_type	text not null,
	root_class	text not null,
	specific_classes	text[] not null,

	constraint fk_phenotype_type FOREIGN KEY (phenotype_type) REFERENCES resistome.term_explanation (internal_name)

);

create table resistome.phenotype_mapping (

    phenotype_map_id    serial primary key,
    standardized_id integer not null,
    original_name   text not null,

    constraint fk_standard_id FOREIGN KEY (standardized_id) REFERENCES resistome.phenotype_standardization (phenotype_id)

);

CREATE INDEX ps_std_idx on resistome.phenotype_standardization (standard_name);


create table resistome.papers (
	
	paper_id	serial primary key,
	title	text not null,
	doi     text not null,
	year	smallint not null,
	research_group	text not null,
	journal	text not null,
	methods text[] not null,
	score	smallint not null,
	reason	text[] not null,
	genome_reference	text,
	designs	bigint not null,
	comments	text,

	constraint positive_designs CHECK (designs >= 0),
	constraint positive_year CHECK (year > 0),
	constraint unique_doi UNIQUE (doi)
	
);

create index papers_idx on resistome.papers (paper_id);

create table resistome.paper_tags (

    id  serial primary key,
	paper_id integer not null,
	tag	text not null,

	CONSTRAINT fk_paper_id FOREIGN KEY ( paper_id ) REFERENCES resistome.papers ( paper_id ),
	CONSTRAINT fk_tags FOREIGN KEY (tag) REFERENCES resistome.term_explanation (internal_name)

);

create index tag_index on resistome.paper_tags (tag);
create index pid_index on resistome.paper_tags (paper_id);

create table resistome.mutants (

	mutant_id	serial primary key,
	paper_id	integer not null,
	name	text not null,
	species	text not null,
	original_strain text not null, --Added Dec2020 to keep the original strain used in the paper so users can contextualize better
	strain	text not null,
	oxygenation	text,
	medium	text,
	carbon_source	text[],
	media_supplements	text[],
	ph	real,
	vessel_type	text,
	culture_volume	real,
	vessel_volume	real,
	temperature	real,
	rotation	real,
	fold_improvement	real,
	initial_fitness	real,
	final_fitness	real,
	fitness_unit	text,
	comments	text,
	
	CONSTRAINT fk_paper_id FOREIGN KEY ( paper_id ) REFERENCES resistome.papers ( paper_id ),
	CONSTRAINT valid_numbers CHECK ( (ph IS NULL or ph > 0) AND
	                                 (culture_volume IS NULL or culture_volume >= 0) AND
	                                 (vessel_volume IS NULL or vessel_volume >= 0) AND
	                                 (temperature IS NULL or temperature > -273.15) AND
	                                 (rotation IS NULL or rotation >= 0) AND
	                                 (initial_fitness IS NULL or initial_fitness >= 0) AND
	                                 (final_fitness IS NULL or final_fitness >= 0) ),
	CONSTRAINT chk_paper_mutant_name UNIQUE (paper_id, name)
);

create index mutant_id_idx on resistome.mutants (mutant_id);

create table resistome.mutant_methods (

    id serial primary key,
	mutant_id	integer not null,
	method text not null,
	
	CONSTRAINT fk_source_id FOREIGN KEY ( mutant_id ) REFERENCES resistome.mutants(mutant_id),
	CONSTRAINT fk_methodology FOREIGN KEY (method) REFERENCES resistome.term_explanation (internal_name)

);

create table resistome.phenotypes (

    id  serial primary key,
	mutant_id	integer not null,
	phenotype	text not null,
	phenotype_class	text not null,
	phenotype_type	varchar(1) not null,
	ontology_root	text not null,
	resistance_level text,
	resistance_unit text,
	
	CONSTRAINT fk_source_id FOREIGN KEY ( mutant_id ) REFERENCES resistome.mutants(mutant_id),
	constraint chk_phenotype_type check (phenotype_type = 'R' or phenotype_type = 'S'),
	constraint fk_phenotype_std FOREIGN KEY (phenotype) REFERENCES resistome.phenotype_standardization (standard_name),
	constraint fk_phenotype_class FOREIGN KEY (phenotype_class) REFERENCES resistome.term_explanation (internal_name)
	
);

create index combined_phenotypes_idx on resistome.phenotypes (mutant_id);


create table resistome.mutations (

    gene_id	serial primary key,
	paper_id	integer not null,
	mutant_id	integer not null,
	name	text not null,
	species	text not null,
	strain	text not null,
	effects text[] not null,
	original	boolean not null,
	
	CONSTRAINT fk_source_id FOREIGN KEY ( mutant_id ) REFERENCES resistome.mutants(mutant_id),
	CONSTRAINT fk_source2_id FOREIGN KEY ( paper_id ) REFERENCES resistome.papers(paper_id),
	CONSTRAINT chk_gene_length CHECK (length(name) > 0),
	CONSTRAINT chk_species_length CHECK (length(species) > 0),
	constraint unique_gene_name_pairing unique (gene_id, name)

);

create index combined_mutation_idx on resistome.mutations (mutant_id);


create table resistome.expression_studies (

    study_id	serial primary key,
	paper_id	integer not null,
	mutant_id	integer not null,
	accession	text,
	measurement_method	text,
	statistical_method	text,
	duration	text,
	growth_phase	text,
	stressor_level	text,
	stressor_units	text,
	
	CONSTRAINT fk_source_id FOREIGN KEY ( mutant_id ) REFERENCES resistome.mutants(mutant_id),
	CONSTRAINT fk_source2_id FOREIGN KEY ( paper_id ) REFERENCES resistome.papers(paper_id)
	
);

create table resistome.expressions (

	paper_id	integer not null,
	mutant_id	integer not null,
	study_id	integer not null,
	gene_id	serial primary key,
	name	text not null,
	species	text not null,
	strain	text not null,
	status varchar(1) not null,
	pvalue	double precision,
	fold_change	double precision,
	
	CONSTRAINT fk_source_id FOREIGN KEY ( mutant_id ) REFERENCES resistome.mutants(mutant_id),
	CONSTRAINT fk_source2_id FOREIGN KEY ( paper_id ) REFERENCES resistome.papers(paper_id),
	CONSTRAINT fk_source3_id FOREIGN KEY ( study_id ) REFERENCES resistome.expression_studies(study_id),
	constraint transdir_chk check (status = '-' or status = '+')
	
);

create index expression_idx on resistome.expressions (mutant_id);


create table resistome.annotations (

	annotation_id	serial primary key,
	gene_id	integer not null,
	mutation_type text not null,
	annotation	jsonb not null,
	--until I add in the validation code, just assume all mutations are valid.
	valid   boolean DEFAULT TRUE,
	valid_reason    text,
	
	constraint fk_gene_id foreign key (gene_id) references resistome.mutations (gene_id),
	constraint fk_mutation_type foreign key (mutation_type) references resistome.term_explanation (internal_name)

);

create index annotations_idx on resistome.annotations (gene_id);

create table resistome.gene_standardization (

    gs_id   serial primary key,
	strain	text	not null,
	species_accession	text not null,
	mg1655_accession	text not null,
	
	constraint unique_tuple_gs unique (strain, species_accession)

);

CREATE INDEX gs_name_idx ON resistome.gene_standardization (species_accession);

create table resistome.gene_metadata (

    gm_id   serial primary key,
	name	text	not null,
	species	text	not null,
	strain	text	not null,
	species_accession	text not null,
	mg1655_accession	text not null,
	
	constraint unique_tuple_gmd unique (name, species, strain)

);

CREATE INDEX gmd_name_idx ON resistome.gene_metadata (name);

create table resistome.abbreviations (

    abbrev_id   serial primary key,
	entry	text not null,
	entry_type text not null,
	converted_entry	text not null,
	
	constraint entry_type_unique unique (entry, entry_type)

);

--create table resistome.metabolomics (
--
--    ion_id  serial primary key,
--	accession	text not null,
--	metabolite	text not null,
--
--	constraint fk_metabolite FOREIGN KEY (metabolite) references public.ms_ions (ion_id)
--
--);
--
--create index metab_idx on resistome.metabolomics (accession);
--
--create table resistome.gene_ontology (
--
--    go_id   serial primary key,
--	accession	text not null,
--	go_term	text not null,
--
--	constraint fk_go_term FOREIGN KEY (go_term) references public.go_table (go_term)
--
--);
--
--create index go_idx on resistome.gene_ontology (accession);