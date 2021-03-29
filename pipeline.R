#putting it all together
options(java.parameters = '-Xmx16g')
dtexbuilder::syn_login()

#Chembl
chembl <- dtexbuilder::process_chembl(activity_file_id = "syn25151631",
                                      structures_file_id = "syn25153278",
                                      names_file_id = "syn25151647")

#DGIDB
dgidb_df <- dtexbuilder::prepare_dgidb_data("https://www.dgidb.org/data/monthly_tsvs/2021-Jan/interactions.tsv")

dtexbuilder::prepare_dgidb_names_for_pubchem(dgidb_df)

dgidb <- dtexbuilder::process_dgidb(dgidb_df, 
                                    path_to_inchi_tsv = "~/Downloads/ids.txt",
                                    path_to_inchikey_tsv = "~/Downloads/inchikey.txt",
                                    path_to_smiles_tsv = "~/Downloads/smiles.txt")

#ChemicalProbes
chemicalprobes <- dtexbuilder::process_chemicalprobes('syn25253561')

#DrugBank
drugbank <- dtexbuilder::process_drugbank(db_synapse_id = "syn25323646",
                                          structures_synapse_id = "syn25323608")

structures <- dtexbuilder::process_structures(chembl$structures, 
                                              dgidb$structures,
                                              chemicalprobes$structures,
                                              drugbank$structures)

activities <- dtexbuilder::process_activities(processed_structures = structures,
                                              chembl$activities, 
                                              dgidb$activities,
                                              chemicalprobes$activities,
                                              drugbank$activities)

names <- dtexbuilder::process_names(processed_structures = structures,
                                    processed_activities = activities,
                                    chembl$names,
                                    dgidb$names,
                                    chemicalprobes$names)

external_ids <- dtexbuilder::process_external_ids(processed_structures = structures, 
                                                  processed_activities = activities)

distinct_structures <- dtexbuilder::filter_structures(processed_structures = structures, processed_activities = activities)

#save.image("imports")

fingerprints_circular <- generate_fingerprints(structure_df = distinct_structures, type = "circular")
fingerprints_maccs <- generate_fingerprints(structure_df = distinct_structures, type = "maccs")
fingerprints_extended <- generate_fingerprints(structure_df = distinct_structures, type = "extended")


saveRDS(activities, 'compound_target_associations_v4.rds')
readr::write_csv(activities, 'compound_target_associations_v4.csv')
saveRDS(names$names, 'distinct_compound_names_v4.rds')
readr::write_csv(names$names, 'distinct_compound_names_v4.csv')
saveRDS(names$synonyms, 'distinct_compound_synonyms_v4.rds')
readr::write_csv(names$synonyms, 'distinct_compound_synonyms_v4.csv')
saveRDS(distinct_structures, 'compound_structures_v4.rds')
readr::write_csv(distinct_structures, 'compound_structures_v4.csv')


saveRDS(fingerprints_circular, 'db_fingerprints_circular_v4.rds')
saveRDS(fingerprints_extended, 'db_fingerprints_extended_v4.rds')
saveRDS(fingerprints_maccs, 'db_fingerprints_maccs_v4.rds')

synapse <- reticulate::import("synapseclient")
syn <- syn$Synapse()

syn$store(synapse$File("compound_target_associations_v4.rds",parentId = "syn25328061"))
syn$store(synapse$File("compound_target_associations_v4.csv",parentId = "syn25328061"))
syn$store(synapse$File("distinct_compound_names_v4.rds",parentId = "syn25328061"))
syn$store(synapse$File("distinct_compound_names_v4.csv",parentId = "syn25328061"))
syn$store(synapse$File("distinct_compound_synonyms_v4.csv",parentId = "syn25328061"))
syn$store(synapse$File("distinct_compound_synonyms_v4.rds",parentId = "syn25328061"))
syn$store(synapse$File("compound_structures_v4.rds",parentId = "syn25328061"))
syn$store(synapse$File("compound_structures_v4.csv",parentId = "syn25328061"))
syn$store(synapse$File("db_fingerprints_circular_v4.rds",parentId = "syn25328061"))
syn$store(synapse$File("db_fingerprints_extended_v4.rds",parentId = "syn25328061"))
syn$store(synapse$File("db_fingerprints_maccs_v4.rds",parentId = "syn25328061"))
