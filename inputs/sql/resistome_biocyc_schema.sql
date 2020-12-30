create table IF NOT EXISTS  strain (

    -- table IF NOT EXISTS  to link against for all strains.

    strain_id  serial primary key,
    species text not null,
    strain  text not null,

    constraint unique_species_strain UNIQUE (species, strain)

);

create index IF NOT EXISTS strain_idx on strain (strain_id);

create table IF NOT EXISTS  genes (

	gene_id	serial primary key,
	strain_id   integer not null,
	biocyc_id	text not null,
	name	text not null,
	synonyms	text[] default array[]::text[],
	products	text[] default array[]::text[],
	accession	text not null,
	description text,
	pseudogene  boolean not null,
	essential	boolean not null,

    constraint fk_strain_genes FOREIGN KEY (strain_id) REFERENCES strain (strain_id),
	constraint unique_accession UNIQUE (accession),
	constraint unique_biocyc_id UNIQUE (biocyc_id),
	constraint name_accession_length CHECK (length(name) > 0 and length(accession) > 0)
);

create index IF NOT EXISTS name_idx on genes (name);
create index IF NOT EXISTS accession_strain_idx on genes (strain_id, accession);
create index IF NOT EXISTS gene_id_source_idx on genes (gene_id);

create table IF NOT EXISTS  gene_locations (

	loc_id	serial primary key,
	gene_id integer not null,
	start	integer not null,
	stop	integer not null,
	direction	varchar(1) not null,

    constraint transdir_chk check (direction = '-' or direction = '+'),
    constraint gene_fk FOREIGN KEY (gene_id) references genes(gene_id)

);

create index IF NOT EXISTS location_gene_idx on gene_locations (gene_id);

create table IF NOT EXISTS  interactions (

    interaction_id  serial primary key,
    strain_id   integer not null,
	interaction_type	text not null,
	regulator_name  text not null,
	target_name text not null,
	regulator	integer,
	target	integer,
	direction varchar(2) not null,
	db_source	text not null,

    constraint fk_strain_int FOREIGN KEY (strain_id) REFERENCES strain (strain_id),
    constraint fk_regulator FOREIGN KEY (regulator) REFERENCES genes (gene_id),
    constraint fk_target FOREIGN KEY (target) REFERENCES genes (gene_id),
	constraint transdir_chk check (direction = '-' or direction = '+' or direction = '?' or direction = '+-'
	                               or direction = '+?' or direction = '-?')

);

create index IF NOT EXISTS  regulator_idx on interactions (regulator);
create index IF NOT EXISTS  target_idx on interactions (target);

create table IF NOT EXISTS  genome (

    chromosome_id   serial primary key,
    strain_id   integer not null,
	chromosome	text not null,
	sequence	text not null,

    constraint fk_strain_int FOREIGN KEY (strain_id) REFERENCES strain (strain_id),
	constraint size_minimum CHECK ( length(chromosome) > 0 and length(sequence) > 0)

);

create table IF NOT EXISTS metabolomics (

    ion_id  serial primary key,
	gene_id	integer not null,
	metabolite_id	integer not null,

	constraint fk_gene_link FOREIGN KEY (gene_id) REFERENCES genes (gene_id),
	constraint fk_metabolite_id FOREIGN KEY (metabolite_id) REFERENCES ms_ions (ion_pk)
);

create index IF NOT EXISTS metabolomics_gene_idx on metabolomics (gene_id);

create table IF NOT EXISTS protein_stability (

	stable_id	serial primary key,
	gene_id integer not null,
	position	smallint	not null,
	wt_aa	varchar(1) not null,
	mutant_aa	varchar(1) not null,
	score	float not null,
	method	text not null,

	constraint accession_fk foreign key (gene_id) references genes (gene_id)

);

create index IF NOT EXISTS protein_stability_gene_idx on protein_stability (gene_id);
create index IF NOT EXISTS stability_tuple_idx on protein_stability (gene_id, position, wt_aa, mutant_aa, method);
create index IF NOT EXISTS stability_tuple_idx on protein_stability (gene_id, position, mutant_aa, method);

create table IF NOT EXISTS go_terms (

    species_go_id   serial primary key,
	gene_id integer not null,
	go_id	integer not null,
	
	constraint accession_fk foreign key (gene_id) references genes (gene_id),
	constraint go_term_fk foreign key (go_id) references public.go_table (go_id)
);

create index IF NOT EXISTS go_term_gene_idx on go_terms (gene_id);

create table IF NOT EXISTS uniprot (

    uniprot_id  serial primary key,
	gene_id integer not null,
	region	text not null,
	start	integer not null,
	stop	integer not null,
	note	text,

	constraint accession_fk foreign key (gene_id) references genes (gene_id)

);

create index IF NOT EXISTS uniprot_gene_idx on uniprot (gene_id);

create table IF NOT EXISTS aa_sequences (

    aa_id   serial primary key,
    gene_id integer not null,
	aa_seq	text not null,
	
	constraint accession_fk foreign key (gene_id) references genes (gene_id)

);

create index IF NOT EXISTS aa_seq_gene_idx on aa_sequences (gene_id);


create table IF NOT EXISTS  nt_sequences (

    nt_id   serial primary key,
	gene_id integer not null,
	nt_seq	text not null,
	
	constraint accession_fk foreign key (gene_id) references genes (gene_id)
);

create index IF NOT EXISTS nt_seq_gene_idx on nt_sequences (gene_id);


create table IF NOT EXISTS  genomic_features (

	feature_id serial primary key,
	strain_id   integer not null,
	biocyc_id	text unique not null,
    feature_type	text not null,
    start	integer not null,
    stop	integer not null,
    associated_data jsonb,
    source  text not null,

    constraint fk_strain_int FOREIGN KEY (strain_id) REFERENCES strain (strain_id)
    
);

create index IF NOT EXISTS genomic_features_gene_idx on genomic_features (strain_id);

create table IF NOT EXISTS  genomic_feature_association (

    -- Biocyc Regulator Genes

	interaction_id serial primary key,
	gene_id integer,
	feature_id integer not null,
	relationship    text not null,
	name	text not null,
	
	constraint feature_id_fk foreign key (feature_id) references genomic_features(feature_id),
	constraint accession_fk foreign key (gene_id) references genes (gene_id)
	-- constraint regulator_fk foreign key (name) references genes (accession)
   
);

create index IF NOT EXISTS gfa_gene_idx on genomic_feature_association (gene_id);
create index IF NOT EXISTS gfa_feature_idx on genomic_feature_association (feature_id);


