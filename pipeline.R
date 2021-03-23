##putting it all together
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
drugbank <- dtexbuilder::process_drugbank("syn12685088", 
                                          "syn12973257")

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

fingerprints_circular <- dtexbuilder::generate_fingerprints(structure_df = structures, type = "circular")
fingerprints_maccs <- dtexbuilder::generate_fingerprints(structure_df = structures, type = "maccs")
fingerprints_extended <- dtexbuilder::generate_fingerprints(structure_df = structures, type = "extended")


saveRDS(activities, 'drug_target_associations_v4.rds')
saveRDS(names$names, 'distinct_compound_names_v4.rds')

saveRDS(fingerprints_circular, 'db_fingerprints_circular_v4.rds')
saveRDS(fingerprints_extended, 'db_fingerprints_extended_v4.rds')
saveRDS(fingerprints_maccs, 'db_fingerprints_maccs_v4.rds')

