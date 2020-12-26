create schema public

create table go_table (

    go_id   serial primary key,
	go_term	text not null,
	ancestors	text[] not null,
	name	text not null,

	constraint unique_go UNIQUE (go_term)

);

create index go_tag_idx on go_table (go_term);

create table go_adj_list (

    edge_id serial primary key,
	parent	text not null,
	child	text not null,
	parent_id  integer not null,
	child_id    integer not null,

	constraint unique_edge UNIQUE (parent, child),
	constraint fk_parent FOREIGN KEY (parent_id) REFERENCES go_table (go_id),
	constraint fk_child FOREIGN KEY (child_id) REFERENCES go_table (go_id)

);

create index go_adj_source_idx on go_adj_list ( parent );
create index go_adj_sink_idx on go_adj_list ( child );

create table ms_ions (

    ion_pk  serial primary key,
	ion_id	text not null,
	kegg_id	text,
	name	text,

	constraint unique_ion UNIQUE (ion_id)
);

create index ms_idx on ms_ions (ion_pk);

create table genetic_code (

    code_id serial primary key,
    code    integer not null,
	codon	varchar(3) not null,
	aa		varchar(1) not null,
	frequency	real not null

);
