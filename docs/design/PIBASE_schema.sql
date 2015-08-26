
CREATE TABLE scop_cla (
    scop_id char(7) not null,
    pdb_id char(4) not null,
    chain_id char(1) not null,
    start_resno char(10) not null,
    end_resno char(10) not null,
    class_id char(1) not null,
    fold_id smallint unsigned not null,
    superfam_id smallint unsigned not null,
    fam_id smallint unsigned not null,
    cl_id mediumint unsigned not null,
    cf_id mediumint unsigned not null,
    sf_id mediumint unsigned not null,
    fa_id mediumint unsigned not null,
    dm_id mediumint unsigned not null,
    sp_id mediumint unsigned not null,
    px_id mediumint unsigned not null,
    PRIMARY KEY (scop_id,chain_id,start_resno,end_resno)
);

CREATE TABLE interface_contacts_prototype (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    chain_no_1 integer unsigned not null,
    chain_id_1 char(1) not null,
    resno_1 char(10) not null,
    resna_1 char(3) not null,
    chain_no_2 integer unsigned not null,
    chain_id_2 char(1) not null,
    resno_2 char(10) not null,
    resna_2 char(3) not null,
    min_dist float not null,
    num_contacts integer unsigned not null,
    num_contacts_4 integer unsigned not null,
    num_contacts_4p5 integer unsigned not null,
    num_contacts_5 integer unsigned not null,
    num_contacts_5p5 integer unsigned not null,
    PRIMARY KEY (bdp_id,subset_id_1,subset_id_2,chain_no_1,resno_1,chain_no_2,resno_2)
);

CREATE TABLE pdb_entry_type (
    pdb_id char(4) not null,
    entry_type char(20) not null,
    experiment_type char(11) not null,
    PRIMARY KEY (pdb_id)
);

CREATE TABLE cath_domain_list (
    domain_name char(7) binary not null,
    pdb_id char(4) not null,
    chain_id char(1) binary not null,
    domain_no tinyint unsigned not null,
    class smallint unsigned not null,
    arch smallint unsigned not null,
    topol smallint unsigned not null,
    homol smallint unsigned not null,
    s35no smallint unsigned not null,
    s95no smallint unsigned not null,
    s100no smallint unsigned not null,
    domain_length smallint unsigned not null,
    PRIMARY KEY (domain_name)
);

CREATE TABLE interface_secstrx_basic_contacts_tables (
    bdp_id integer unsigned not null,
    table_name char(40),
    source_file char(255),
    PRIMARY KEY (bdp_id)
);

CREATE TABLE pdb_release (
    pdb_id char(4) not null,
    release_date date,
    PRIMARY KEY (pdb_id)
);

CREATE TABLE pilig_pep_lig_overlap (
    subset_id char(50) not null,
    chain_subset_id char(50) not null,
    chain_length integer not null,
    ligbs_id char(50) not null,
    ligbs_subset_id char(50) not null,
    ligcode char(3) not null,
    numres_p integer not null,
    numres_p_ident integer not null,
    numres_l integer not null,
    numres_l_ident integer not null,
    numres_l_and_p integer not null,
    numres_l_and_p_ident integer not null,
    numres_domain_nongap integer not null,
    numres_domain_nongap_ident integer not null,
    bdp_id integer not null,
    PRIMARY KEY (subset_id,chain_subset_id,ligbs_id,ligbs_subset_id)
);

CREATE TABLE external_data_sources (
    data_source char(50) not null,
    table_name char(100) not null,
    upload_date date not null,
    data_date date not null,
    version char(50),
    url char(100),
    comments char(100),
    PRIMARY KEY (data_source,table_name)
);

CREATE TABLE bdp_numdom (
    bdp_id integer unsigned not null,
    subset_source_id integer not null,
    num_dom integer not null,
    PRIMARY KEY (bdp_id,subset_source_id)
);

CREATE TABLE bdp_files (
    bdp_id integer unsigned auto_increment,
    file_path char(255) not null,
    file_base char(30) not null,
    pdb_id char(4) not null,
    raw_pdb bool,
    PRIMARY KEY (bdp_id)
);

CREATE TABLE subsets (
    subset_id char(50) binary not null,
    bdp_id integer,
    pdb_id char(4),
    description char(250),
    subset_source_id integer not null,
    class char(70) not null,
    PRIMARY KEY (subset_id)
);

CREATE TABLE pilig_lig_resno (
    ligbs_id char(50) not null,
    subset_id char(50) not null,
    cluster_num integer not null,
    scop_level char(5) not null,
    resno_string text not null,
    PRIMARY KEY (ligbs_id,subset_id,scop_level)
);

CREATE TABLE pilig_pep_lig_overlap_summary (
    subset_id char(50) not null,
    chain_subset_id char(50) not null,
    chain_length integer not null,
    numres_p integer not null,
    cum_numres_l_and_p integer not null,
    max_l_and_p integer not null,
    lig_max_l_and_p char(20) not null,
    cum_l_and_p_perseqid_20 integer not null,
    max_l_and_p_perseqid_20 integer not null,
    cum_l_and_p_perseqid_50 integer not null,
    max_l_and_p_perseqid_50 integer not null,
    cum_l_and_p_perseqid_90 integer not null,
    max_l_and_p_perseqid_90 integer not null,
    PRIMARY KEY (subset_id,chain_subset_id)
);

CREATE TABLE pilig_pep (
    subset_id char(50) not null,
    chain_subset_id char(50) not null,
    cluster_num integer not null,
    scop_level char(5) not null,
    class char(70) not null,
    chain_length integer not null,
    alnpos_string text not null,
    numres_bs integer not null,
    aln_length integer not null,
    bdp_id integer not null,
    PRIMARY KEY (subset_id,chain_subset_id,scop_level)
);

CREATE TABLE pdb_obsolete (
    old_pdb_id char(4) not null,
    new_pdb_id char(4),
    obsolete_date date,
    PRIMARY KEY (old_pdb_id)
);

CREATE TABLE interface_secstrx_contacts_tables (
    bdp_id integer unsigned not null,
    table_name char(40),
    source_file char(255),
    PRIMARY KEY (bdp_id)
);

CREATE TABLE bdp_secstrx_tables (
    bdp_id integer unsigned not null,
    table_name char(40),
    source_file char(255),
    PRIMARY KEY (bdp_id)
);

CREATE TABLE subsets_files (
    subset_id char(50) binary,
    file_path char(200) not null,
    PRIMARY KEY (subset_id)
);

CREATE TABLE bdp_xpack (
    bdp_id integer unsigned not null,
    xpack tinyint(1) not null,
    PRIMARY KEY (bdp_id)
);

CREATE TABLE bindingsite_contacts_tables (
    bdp_id integer unsigned not null,
    cutoff float not null,
    table_name char(40),
    source_file char(255),
    PRIMARY KEY (bdp_id)
);

CREATE TABLE pisa_index (
    pisa_id char(10) not null,
    pdb_id char(4) not null,
    split_no smallint unsigned,
    entry_type char(30),
    homo_hetero enum('HETERO', 'HOMO','OVERLAP','IMPSYM','TOO_BIG') not null,
    num_chain int unsigned,
    num_resid int unsigned,
    total_sasa float,
    buried_sasa float,
    num_hbonds int unsigned,
    num_salt_bridges int unsigned,
    num_disulfide int unsigned,
    diss_energy float,
    assembly_formula char(75),
    PRIMARY KEY (pisa_id)
);

CREATE TABLE bdp_residues_tables (
    bdp_id integer unsigned not null,
    table_name char(40),
    source_file char(255),
    PRIMARY KEY (bdp_id)
);

CREATE TABLE bindingsite_secstrx_basic_contacts_prototype (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    subset_id char(50) not null,
    sse_id_1 integer unsigned not null,
    sse_id_2 integer unsigned not null,
    num_res_1 integer unsigned not null,
    num_res_2 integer unsigned not null,
    num_res_12 integer unsigned not null,
    sse_1 enum('H', 'B', 'T', ' ') not null,
    sse_2 enum('H', 'B', 'T', ' ') not null,
    PRIMARY KEY (bdp_id,subset_id_1,subset_id_2,subset_id)
);

CREATE TABLE bindingsite_secstrx_basic_contacts_tables (
    bdp_id integer unsigned not null,
    table_name char(40),
    source_file char(255),
    PRIMARY KEY (bdp_id)
);

CREATE TABLE pilig_pep_lig_overlap_seqaln (
    subset_id char(50) not null,
    chain_subset_id char(50) not null,
    ligbs_id char(50) not null,
    ligbs_subset_id char(50) not null,
    seqalnstring_p text not null,
    seqalnstring_l text not null,
    PRIMARY KEY (subset_id,chain_subset_id,ligbs_id,ligbs_subset_id)
);

CREATE TABLE interface_sasa (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    dsasa_all float not null,
    dsasa_sc float not null,
    dsasa_mc float not null,
    dsasa_polar float not null,
    dsasa_nonpolar float not null,
    PRIMARY KEY (subset_id_1,subset_id_2)
);

CREATE TABLE pilig_pi (
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    subset_id char(50) not null,
    cluster_num integer not null,
    scop_level char(5) not null,
    class char(70) not null,
    alnpos_string text not null,
    numres_bs integer not null,
    aln_length integer not null,
    bdp_id integer not null,
    PRIMARY KEY (subset_id_1,subset_id_2,subset_id,scop_level)
);

CREATE TABLE patch_residues_tables (
    bdp_id integer unsigned not null,
    cutoff float not null,
    table_name char(40),
    source_file char(255),
    PRIMARY KEY (bdp_id)
);

CREATE TABLE interface_size (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) binary not null,
    subset_id_2 char(50) binary not null,
    num_res_1 integer unsigned not null,
    num_res_2 integer unsigned not null,
    subset_size_1 integer unsigned not null,
    subset_size_2 integer unsigned not null,
    interface_res_1 text,
    interface_res_2 text,
    PRIMARY KEY (subset_id_1,subset_id_2)
);

CREATE TABLE pilig_pi_lig_overlap_seqaln (
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    subset_id char(50) not null,
    ligbs_id char(50) not null,
    ligbs_subset_id char(50) not null,
    seqalnstring_p text not null,
    seqalnstring_l text not null,
    PRIMARY KEY (subset_id_1,subset_id_2,subset_id,ligbs_id,ligbs_subset_id)
);

CREATE TABLE interface_secstrx_profile (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    num_sse_1 integer unsigned not null,
    num_sse_2 integer unsigned not null,
    num_sse_12 integer unsigned not null,
    sse_contact_vector char(36) not null,
    sse_vector_1 char(8) not null,
    sse_vector_2 char(8) not null,
    sse_contact_vector_norm integer unsigned not null,
    sse_vector_1_norm integer unsigned not null,
    sse_vector_2_norm integer unsigned not null,
    num_sse_1_H integer unsigned not null,
    num_sse_1_B integer unsigned not null,
    num_sse_1_E integer unsigned not null,
    num_sse_1_G integer unsigned not null,
    num_sse_1_I integer unsigned not null,
    num_sse_1_T integer unsigned not null,
    num_sse_1_S integer unsigned not null,
    num_sse_1_u integer unsigned not null,
    num_sse_2_H integer unsigned not null,
    num_sse_2_B integer unsigned not null,
    num_sse_2_E integer unsigned not null,
    num_sse_2_G integer unsigned not null,
    num_sse_2_I integer unsigned not null,
    num_sse_2_T integer unsigned not null,
    num_sse_2_S integer unsigned not null,
    num_sse_2_u integer unsigned not null,
    num_sse_12_HH integer unsigned not null,
    num_sse_12_HB integer unsigned not null,
    num_sse_12_HE integer unsigned not null,
    num_sse_12_HG integer unsigned not null,
    num_sse_12_HI integer unsigned not null,
    num_sse_12_HT integer unsigned not null,
    num_sse_12_HS integer unsigned not null,
    num_sse_12_Hu integer unsigned not null,
    num_sse_12_BB integer unsigned not null,
    num_sse_12_BE integer unsigned not null,
    num_sse_12_BG integer unsigned not null,
    num_sse_12_BI integer unsigned not null,
    num_sse_12_BT integer unsigned not null,
    num_sse_12_BS integer unsigned not null,
    num_sse_12_Bu integer unsigned not null,
    num_sse_12_EE integer unsigned not null,
    num_sse_12_EG integer unsigned not null,
    num_sse_12_EI integer unsigned not null,
    num_sse_12_ET integer unsigned not null,
    num_sse_12_ES integer unsigned not null,
    num_sse_12_Eu integer unsigned not null,
    num_sse_12_GG integer unsigned not null,
    num_sse_12_GI integer unsigned not null,
    num_sse_12_GT integer unsigned not null,
    num_sse_12_GS integer unsigned not null,
    num_sse_12_Gu integer unsigned not null,
    num_sse_12_II integer unsigned not null,
    num_sse_12_IT integer unsigned not null,
    num_sse_12_IS integer unsigned not null,
    num_sse_12_Iu integer unsigned not null,
    num_sse_12_TT integer unsigned not null,
    num_sse_12_TS integer unsigned not null,
    num_sse_12_Tu integer unsigned not null,
    num_sse_12_SS integer unsigned not null,
    num_sse_12_Su integer unsigned not null,
    num_sse_12_uu integer unsigned not null,
    numres_sse_1_H integer unsigned not null,
    numres_sse_1_B integer unsigned not null,
    numres_sse_1_E integer unsigned not null,
    numres_sse_1_G integer unsigned not null,
    numres_sse_1_I integer unsigned not null,
    numres_sse_1_T integer unsigned not null,
    numres_sse_1_S integer unsigned not null,
    numres_sse_1_u integer unsigned not null,
    numres_sse_2_H integer unsigned not null,
    numres_sse_2_B integer unsigned not null,
    numres_sse_2_E integer unsigned not null,
    numres_sse_2_G integer unsigned not null,
    numres_sse_2_I integer unsigned not null,
    numres_sse_2_T integer unsigned not null,
    numres_sse_2_S integer unsigned not null,
    numres_sse_2_u integer unsigned not null,
    numres_sse_12_HH integer unsigned not null,
    numres_sse_12_HB integer unsigned not null,
    numres_sse_12_HE integer unsigned not null,
    numres_sse_12_HG integer unsigned not null,
    numres_sse_12_HI integer unsigned not null,
    numres_sse_12_HT integer unsigned not null,
    numres_sse_12_HS integer unsigned not null,
    numres_sse_12_Hu integer unsigned not null,
    numres_sse_12_BB integer unsigned not null,
    numres_sse_12_BE integer unsigned not null,
    numres_sse_12_BG integer unsigned not null,
    numres_sse_12_BI integer unsigned not null,
    numres_sse_12_BT integer unsigned not null,
    numres_sse_12_BS integer unsigned not null,
    numres_sse_12_Bu integer unsigned not null,
    numres_sse_12_EE integer unsigned not null,
    numres_sse_12_EG integer unsigned not null,
    numres_sse_12_EI integer unsigned not null,
    numres_sse_12_ET integer unsigned not null,
    numres_sse_12_ES integer unsigned not null,
    numres_sse_12_Eu integer unsigned not null,
    numres_sse_12_GG integer unsigned not null,
    numres_sse_12_GI integer unsigned not null,
    numres_sse_12_GT integer unsigned not null,
    numres_sse_12_GS integer unsigned not null,
    numres_sse_12_Gu integer unsigned not null,
    numres_sse_12_II integer unsigned not null,
    numres_sse_12_IT integer unsigned not null,
    numres_sse_12_IS integer unsigned not null,
    numres_sse_12_Iu integer unsigned not null,
    numres_sse_12_TT integer unsigned not null,
    numres_sse_12_TS integer unsigned not null,
    numres_sse_12_Tu integer unsigned not null,
    numres_sse_12_SS integer unsigned not null,
    numres_sse_12_Su integer unsigned not null,
    numres_sse_12_uu integer unsigned not null,
    PRIMARY KEY (subset_id_1,subset_id_2)
);

CREATE TABLE scop_interface_clusters (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    scopclass_pair char(150) not null,
    cluster_level enum('fam','sf') not null,
    interface_class char(150) not null,
    cluster_no integer not null,
    member_no integer not null,
    PRIMARY KEY (bdp_id,subset_id_1,subset_id_2)
);

CREATE TABLE interface_contacts_tables (
    bdp_id integer unsigned not null,
    cutoff float not null,
    table_name char(40),
    source_file char(255),
    PRIMARY KEY (bdp_id)
);

CREATE TABLE bdp_chains (
    bdp_id integer unsigned not null,
    real_chain_no integer unsigned not null,
    real_chain_id char(1) not null,
    chain_type enum('p', 'n') not null,
    pdb_chain_no integer unsigned,
    pdb_chain_id char(1),
    start_resno char(10) not null,
    start_resno_int integer not null,
    end_resno char(10) not null,
    end_resno_int integer not null,
    num_res integer unsigned not null,
    num_atoms integer unsigned not null,
    num_hetatm integer unsigned not null,
    sequence text not null,
    PRIMARY KEY (bdp_id,real_chain_no)
);

CREATE TABLE subsets_sequence (
    bdp_id integer not null,
    subset_id char(50) binary not null,
    num_chains integer not null,
    num_res integer not null,
    sequence text not null,
    PRIMARY KEY (subset_id)
);

CREATE TABLE patch_residues_prototype (
    bdp_id integer unsigned not null,
    subset_id char(50) not null,
    chain_no integer unsigned not null,
    chain_id char(1) not null,
    resno char(10) not null,
    resna char(3) not null,
    num_contacts integer unsigned,
    PRIMARY KEY (bdp_id,subset_id,chain_no,resno)
);

CREATE TABLE bindingsite_sse_topology (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    subset_id char(50) not null,
    num_domains integer unsigned not null,
    num_edges integer unsigned not null,
    nodelist text not null,
    nodelist_sse text not null,
    edgelist text not null,
    edgelist_sse text not null,
    PRIMARY KEY (bdp_id,subset_id_1,subset_id_2,subset_id)
);

CREATE TABLE bdp_residues_prototype (
    bdp_id integer unsigned not null,
    chain_no integer unsigned not null,
    chain_id char(1) not null,
    resno_serial integer unsigned not null,
    resno char(10) not null,
    resno_int integer not null,
    resna char(3) not null,
    chain_type enum('p', 'n') not null,
    PRIMARY KEY (bdp_id,chain_no,resno_serial)
);

CREATE TABLE pdb_patchres_bychain (
    pdb_id char(4) not null,
    chain_id char(1) not null,
    patch_no integer not null,
    resnos text not null,
    PRIMARY KEY (pdb_id,chain_id,patch_no)
);

CREATE TABLE interface_continuity (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    seq_segments_1 integer unsigned not null,
    strx_patches_1 integer unsigned not null,
    seq_segments_2 integer unsigned not null,
    strx_patches_2 integer unsigned not null,
    cutoff float not null,
    PRIMARY KEY (subset_id_1,subset_id_2)
);

CREATE TABLE pilig_pi_resno (
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    subset_id char(50) not null,
    cluster_num integer not null,
    scop_level char(5) not null,
    resno_string text not null,
    PRIMARY KEY (subset_id_1,subset_id_2,subset_id,scop_level)
);

CREATE TABLE interface_planarity (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    planarity float not null,
    PRIMARY KEY (subset_id_1,subset_id_2)
);

CREATE TABLE subsets_residues_tables (
    bdp_id integer unsigned not null,
    table_name char(40),
    source_file char(255),
    PRIMARY KEY (bdp_id)
);

CREATE TABLE interface_resvector (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    contact_vector char(210) not null,
    res_vector_1 char(20) not null,
    res_vector_2 char(20) not null,
    contact_vector_norm integer not null,
    res_vector_1_norm integer not null,
    res_vector_2_norm integer not null,
    PRIMARY KEY (subset_id_1,subset_id_2)
);

CREATE TABLE bdp_interaction_topology_graph (
    bdp_id integer unsigned not null,
    subset_source_id integer unsigned not null,
    file_path char(200) not null,
    PRIMARY KEY (bdp_id,subset_source_id)
);

CREATE TABLE cath_domain_description (
    domain_name char(7) binary not null,
    pdb_id char(4) not null,
    chain_id char(1) binary not null,
    domain_no tinyint unsigned not null,
    class_id smallint unsigned not null,
    arch_id smallint unsigned not null,
    topol_id smallint unsigned not null,
    homol_id smallint unsigned not null,
    class char(30) not null,
    arch char(50) not null,
    topol char(150) not null,
    homol char(150) not null,
    domain_length smallint unsigned not null,
    domain_sequence_header char(25) not null,
    domain_sequence text not null,
    segment_id tinyint not null,
    start_resno char(10) not null,
    end_resno char(10) not null,
    segment_length smallint unsigned not null,
    segment_sequence_header char(25) not null,
    segment_sequence text not null,
    PRIMARY KEY (domain_name,segment_id)
);

CREATE TABLE interface_secstrx_contacts_prototype (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    sse_id_1 integer unsigned not null,
    sse_id_2 integer unsigned not null,
    num_res_1 integer unsigned not null,
    num_res_2 integer unsigned not null,
    num_res_12 integer unsigned not null,
    sse_1 enum('H', 'B', 'E', 'G', 'I', 'T', 'S', ' ') not null,
    sse_2 enum('H', 'B', 'E', 'G', 'I', 'T', 'S', ' ') not null,
    PRIMARY KEY (bdp_id,subset_id_1,subset_id_2)
);

CREATE TABLE pilig_ligand_info (
    ligcode char(3) not null,
    name text,
    formula char(30) not null,
    mol_weight float,
    num_atoms integer,
    num_atoms_nonh integer,
    PRIMARY KEY (ligcode)
);

CREATE TABLE interface_sc (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    sc float,
    median_dist float,
    PRIMARY KEY (subset_id_1,subset_id_2)
);

CREATE TABLE pilig_pep_resno (
    subset_id char(50) not null,
    chain_subset_id char(50) not null,
    cluster_num integer not null,
    scop_level char(5) not null,
    resno_string text not null,
    PRIMARY KEY (subset_id,chain_subset_id,scop_level)
);

CREATE TABLE interface_contacts_special_prototype (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    chain_no_1 integer unsigned not null,
    chain_id_1 char(1) not null,
    resno_1 char(10) not null,
    resna_1 char(3) not null,
    chain_no_2 integer unsigned not null,
    chain_id_2 char(1) not null,
    resno_2 char(10) not null,
    resna_2 char(3) not null,
    type enum('salt', 'hbond', 'ssbond'),
    PRIMARY KEY (bdp_id,subset_id_1,subset_id_2,chain_no_1,resno_1,chain_no_2,resno_2,type)
);

CREATE TABLE cath_names (
    class smallint unsigned not null,
    arch smallint unsigned not null,
    topol smallint unsigned not null,
    homol smallint unsigned not null,
    representative_dom char(7) not null,
    description tinytext,
    PRIMARY KEY (class,arch,topol,homol)
);

CREATE TABLE interface_contacts_special_tables (
    bdp_id integer unsigned not null,
    cutoff float not null,
    table_name char(40),
    source_file char(255),
    PRIMARY KEY (bdp_id)
);

CREATE TABLE interface_secstrx_prototype (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    subset_id char(50) not null,
    sse_id integer unsigned not null,
    num_res integer unsigned not null,
    sse enum('H', 'B', 'E', 'G', 'I', 'T', 'S', ' ') not null,
    PRIMARY KEY (bdp_id,subset_id_1,subset_id_2,subset_id)
);

CREATE TABLE scop_hie (
    self_sun_id mediumint unsigned not null,
    parent_sun_id mediumint unsigned,
    kids_sun_id text,
    PRIMARY KEY (self_sun_id)
);

CREATE TABLE pilig_pi_lig_overlap_summary (
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    numres_p_1 integer not null,
    cum_l_and_p_1 integer not null,
    max_l_and_p_1 integer not null,
    lig_max_l_and_p_1 char(20) not null,
    cum_l_and_p_1_perseqid_20 integer not null,
    max_l_and_p_1_perseqid_20 integer not null,
    cum_l_and_p_1_perseqid_50 integer not null,
    max_l_and_p_1_perseqid_50 integer not null,
    cum_l_and_p_1_perseqid_90 integer not null,
    max_l_and_p_1_perseqid_90 integer not null,
    numres_p_2 integer not null,
    cum_l_and_p_2 integer not null,
    max_l_and_p_2 integer not null,
    lig_max_l_and_p_2 char(20) not null,
    cum_l_and_p_2_perseqid_20 integer not null,
    max_l_and_p_2_perseqid_20 integer not null,
    cum_l_and_p_2_perseqid_50 integer not null,
    max_l_and_p_2_perseqid_50 integer not null,
    cum_l_and_p_2_perseqid_90 integer not null,
    max_l_and_p_2_perseqid_90 integer not null,
    PRIMARY KEY (subset_id_1,subset_id_2)
);

CREATE TABLE subsets_sasa (
    bdp_id integer unsigned not null,
    subset_id char(50) not null,
    sasa_all float not null,
    sasa_sc float not null,
    sasa_mc float not null,
    sasa_polar float not null,
    sasa_nonpolar float not null,
    PRIMARY KEY (subset_id)
);

CREATE TABLE interface_vol (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    vol_12 float not null,
    delta_vol float not null,
    PRIMARY KEY (subset_id_1,subset_id_2)
);

CREATE TABLE bdp_interaction_topology (
    bdp_id integer unsigned not null,
    subset_source_id integer unsigned not null,
    num_domains integer unsigned not null,
    num_edges integer unsigned not null,
    num_domain_classes integer unsigned not null,
    nodestring text not null,
    edgestring text not null,
    PRIMARY KEY (bdp_id,subset_source_id)
);

CREATE TABLE bindingsite_contacts_prototype (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    subset_id char(50) not null,
    chain_id_1 char(1) not null,
    resno_1 char(10) not null,
    resna_1 char(3) not null,
    chain_id_2 char(1) not null,
    resno_2 char(10) not null,
    resna_2 char(3) not null,
    min_dist float not null,
    num_contacts integer unsigned not null,
    num_contacts_4 integer unsigned not null,
    num_contacts_4p5 integer unsigned not null,
    num_contacts_5 integer unsigned not null,
    num_contacts_5p5 integer unsigned not null,
    PRIMARY KEY (bdp_id,subset_id_1,subset_id,resno_1,resno_2)
);

CREATE TABLE interface_sse_topology (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    num_domains integer unsigned not null,
    num_edges integer unsigned not null,
    nodelist text not null,
    nodelist_sse text not null,
    edgelist text not null,
    edgelist_sse text not null,
    PRIMARY KEY (bdp_id,subset_id_1,subset_id_2)
);

CREATE TABLE subsets_source (
    subset_source_id integer unsigned auto_increment not null,
    subset_source char(50) not null,
    version char(20) not null,
    PRIMARY KEY (subset_source_id,subset_source)
);

CREATE TABLE subsets_class (
    subset_source_id integer unsigned not null,
    class char(70) not null,
    description char(250),
    PRIMARY KEY (subset_source_id,class)
);

CREATE TABLE bdp_numres (
    bdp_id integer unsigned not null,
    num_res integer not null,
    PRIMARY KEY (bdp_id)
);

CREATE TABLE scop_des (
    sun_id mediumint unsigned not null,
    entry_type enum('cl', 'cf', 'sf', 'fa', 'dm' ,'sp', 'px') not null,
    class_id char(1) not null,
    fold_id smallint unsigned not null,
    superfam_id smallint unsigned not null,
    fam_id smallint unsigned not null,
    scop_id char(7),
    description char(250),
    PRIMARY KEY (sun_id)
);

CREATE TABLE pilig_pi_lig_overlap_resno (
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    subset_id char(50) not null,
    ligbs_id char(50) not null,
    ligbs_subset_id char(50) not null,
    resno_l_onto_pbs text not null,
    resno_lp_onto_pbs text not null,
    resno_p_onto_lbs text not null,
    resno_lp_onto_lbs text not null,
    PRIMARY KEY (subset_id_1,subset_id_2,subset_id,ligbs_id,ligbs_subset_id)
);

CREATE TABLE interface_secstrx_basic_contacts_prototype (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    sse_id_1 integer unsigned not null,
    sse_id_2 integer unsigned not null,
    num_res_1 integer unsigned not null,
    num_res_2 integer unsigned not null,
    num_res_12 integer unsigned not null,
    sse_1 enum('H', 'B', 'T', ' ') not null,
    sse_2 enum('H', 'B', 'T', ' ') not null,
    PRIMARY KEY (bdp_id,subset_id_1,subset_id_2)
);

CREATE TABLE pilig_pi_lig_overlap (
    subset_id_1 char(50) not null,
    subset_id_2 char(50) not null,
    subset_id char(50) not null,
    ligbs_id char(50) not null,
    ligbs_subset_id char(50) not null,
    ligcode char(3) not null,
    numres_p_side integer not null,
    numres_p_ident integer not null,
    numres_l integer not null,
    numres_l_ident integer not null,
    numres_l_and_p_side integer not null,
    numres_l_and_p_ident integer not null,
    numres_domain_nongap integer not null,
    numres_domain_nongap_ident integer not null,
    bdp_id integer not null,
    PRIMARY KEY (subset_id_1,subset_id_2,subset_id,ligbs_id,ligbs_subset_id)
);

CREATE TABLE subsets_residues_prototype (
    bdp_id integer unsigned not null,
    chain_no integer unsigned not null,
    chain_id char(1) not null,
    resno_serial integer unsigned not null,
    resno char(10) not null,
    resno_int integer not null,
    subset_id char(50) not null,
    PRIMARY KEY (bdp_id,chain_no,resno_serial,subset_id)
);

CREATE TABLE intersubset_contacts (
    bdp_id integer unsigned not null,
    subset_id_1 char(50) binary not null,
    class_1 char(70) not null,
    subset_id_2 char(50) binary not null,
    class_2 char(70) not null,
    num_contacts integer unsigned not null,
    cutoff float not null,
    num_contacts_4 integer unsigned not null,
    num_contacts_4p5 integer unsigned not null,
    num_contacts_5 integer unsigned not null,
    num_contacts_5p5 integer unsigned not null,
    hbond integer unsigned not null,
    salt integer unsigned not null,
    ssbond integer unsigned not null,
    chains enum('same', 'diff', 'both') not null,
    PRIMARY KEY (subset_id_1,subset_id_2)
);

CREATE TABLE subsets_vol (
    subset_id char(50) binary not null,
    volume float not null,
    PRIMARY KEY (subset_id)
);

CREATE TABLE pilig_lig (
    ligbs_id char(50) not null,
    subset_id char(50) not null,
    cluster_num integer not null,
    scop_level char(5) not null,
    class char(70) not null,
    alnpos_string text not null,
    numres_bs integer not null,
    aln_length integer not null,
    ligcode char(3) not null,
    bdp_id integer not null,
    PRIMARY KEY (ligbs_id,subset_id,scop_level)
);

CREATE TABLE subsets_details (
    subset_id char(50) binary not null,
    segment_id char(50) not null,
    chain_no integer unsigned,
    chain_id char(1) not null,
    start_resno char(10) not null,
    start_resno_int integer,
    end_resno char(10) not null,
    end_resno_int integer,
    PRIMARY KEY (subset_id,segment_id)
);

CREATE TABLE subsets_env (
    bdp_id integer not null,
    subset_id char(50) binary not null,
    num_domains_bdp integer not null,
    num_interactions integer not null,
    num_domains_chain integer not null,
    PRIMARY KEY (subset_id)
);

CREATE TABLE bdp_secstrx_prototype (
    bdp_id integer unsigned not null,
    chain_id char(1) not null,
    resno char(10) not null,
    sse enum('H', 'B', 'E', 'G', 'I', 'T', 'S', ' ') not null,
    sse_basic enum('H', 'B', 'T', ' ') not null,
    sse_id integer unsigned not null,
    sse_basic_id integer unsigned not null,
    PRIMARY KEY (bdp_id)
);

CREATE TABLE pdb_entries (
    pdb_id char(4) not null,
    header char(100),
    accession_date date,
    compound text,
    source char(150),
    author text,
    resolution float,
    experiment_type char(90),
    PRIMARY KEY (pdb_id)
);

CREATE TABLE interface_secstrx_tables (
    bdp_id integer unsigned not null,
    table_name char(40),
    source_file char(255),
    PRIMARY KEY (bdp_id)
);

CREATE TABLE pilig_pep_lig_overlap_resno (
    subset_id char(50) not null,
    chain_subset_id char(50) not null,
    ligbs_id char(50) not null,
    ligbs_subset_id char(50) not null,
    resno_l_onto_pbs text not null,
    resno_lp_onto_pbs text not null,
    resno_p_onto_lbs text not null,
    resno_lp_onto_lbs text not null,
    PRIMARY KEY (subset_id,chain_subset_id,ligbs_id,ligbs_subset_id)
);

