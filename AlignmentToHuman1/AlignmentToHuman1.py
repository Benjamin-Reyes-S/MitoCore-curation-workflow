import cobra
import re
import subprocess
import pandas as pd

from src.core.parsing import get_model_ids, map_ids_to_db, add_ids_to_model, add_ids_to_model_notes, get_model_ids_for_gpr, add_gpr_to_model, clean_invalid_annotations
from src.core.processingDbFiles import mnx_processing_file_mnx, mnx_processing_file_eccode
from src.core.APIRequests import get_hsa_from_uniprot, get_kegg_from_hsa

mitocore= cobra.io.read_sbml_model("/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_MitoMammal.xml")
human1= cobra.io.read_sbml_model('/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Input_Models/Human-GEM.xml')

print("Human1 update started.")
def main():
    
    print("starting alignmento to Human1")

    # Process mnx file (MNX processing function 1.)
    mnx_reactions= mnx_processing_file_mnx('/Users/benjaminreyes/Desktop/Masterarbeit/neutrophil_modeling/neutrophil-modeling/curation/Files_Databases/MNX_reac_xref.tsv','Recon2','reaction',("bigg.reaction:", "biggR:"))
    print('processes MNX reaction file')
    print(mnx_reactions)

    # get kegg and bigg ids from model
    print('geting ids for all possible databases ')
    df_react = get_model_ids(mitocore, 'reactions', 'notes', ['KEGG', 'Recon2'])
    print(df_react)

    #  map model ids (bigg) to metanetx processed  csv file
    print('Mapping  metabolites to metanetx.')
    df_mnx_react = map_ids_to_db(df_react, mnx_reactions,'Recon2', 'Recon2', 'metanetx.reaction', ',')
    print('MNX metabolites \n',df_mnx_react)

    # Add mnx ids to the model
    print('adding MNX ids to the model')
    mnx_react = add_ids_to_model(mitocore, df_mnx_react, 'reactions','metanetx.reaction', )


    #-----------BiGG to Human1 mapping and adding Human1 ids to the model-----------------
    print('geting ids for all bigg')
    df_react = get_model_ids(mitocore, 'reactions', '_annotation', ['bigg'])
    print(df_react)

    # map bigg ids to Human1 mapping file
    human1_react= pd.read_csv('/Users/benjaminreyes/Desktop/Masterarbeit/neutrophil_modeling/neutrophil-modeling/curation/Files_Databases/Human1_reactions.tsv', sep='\t')
    df_human1_react = map_ids_to_db(df_react, human1_react,'bigg', 'rxnBiGGID', 'rxns', '\t')
    print('human1 reactions \n',df_human1_react)
    df_human1_react = df_human1_react.rename(columns={'model_id': 'model_id', 'rxns': 'Human1'})
    print(df_human1_react)

    # add Human1 ids to model
    print('adding human1 ids to the model notes')
    mitocore = add_ids_to_model_notes(mitocore, df_human1_react, 'reactions','Human1')

    #-----------KEGG to Human1 mapping and adding Human1 ids to the model-----------------
    print('geting ids for all kegg')
    df_react = get_model_ids(mitocore, 'reactions', '_annotation', ['kegg'])
    print(df_react)
    df_react = df_react.explode('kegg')
    print(df_react)

    #  map bigg ids to Human1 mapping file
    print('Mapping kegg reactions to human1')
    df_human1_react = map_ids_to_db(df_react, human1_react,'kegg', 'rxnKEGGID', 'rxns', '\t')
    print('human1 reactions \n',df_human1_react)
    df_human1_react = df_human1_react.rename(columns={'model_id': 'model_id', 'rxns': 'Human1'})
    print(df_human1_react)

    # add Human1 ids to model
    print('adding human1 ids to the model')
    mitocore = add_ids_to_model_notes(mitocore, df_human1_react, 'reactions','Human1')

    #-----------MetaNetX to Human1 mapping and adding Human1 ids to the model-----------------
    print('geting ids for all MNX')
    df_react = get_model_ids(mitocore, 'reactions', '_annotation', ['metanetx'])
    print(df_react)
    # metanetx ids to Human1 mapping file
    print('Mapping mnx reactions to human1')
    df_human1_react = map_ids_to_db(df_react, human1_react,'metanetx', 'rxnMetaNetXID', 'rxns', '\t')
    print('human1 reactions \n',df_human1_react)
    df_human1_react = df_human1_react.rename(columns={'model_id': 'model_id', 'rxns': 'Human1'})
    print(df_human1_react)
    # add Human1 ids to model
    print('adding human1 ids to the model')
    mitocore = add_ids_to_model_notes(mitocore, df_human1_react, 'reactions','Human1')



    # extract Human1 ids in mitocore
    human1_mitocore=get_model_ids(mitocore, 'reactions', 'notes', ['Human1'])
    #print(human1_mitocore)
    # extract GPR rules from models
    gpr_mitocore= get_model_ids_for_gpr(mitocore, 'reactions', '_gpr')
    #print(gpr_mitocore)
    gpr_human1= get_model_ids_for_gpr(human1, 'reactions', '_gpr')

    print(gpr_human1)
    # merge mitocore dfs to have Human1 ids and GPR rules.
    gpr_mitocore = gpr_mitocore.merge(human1_mitocore, on='model_id')
    print(gpr_mitocore)
    print('df with _gpr and model ids in Mitocore', gpr_mitocore)
    # df with gpr from human1 mapped to mitocore
    gpr_parsed=map_ids_to_db(gpr_mitocore, gpr_human1, 'Human1', 'model_id', '_gpr', ',')
    print(gpr_parsed)

    # add gpr from human1 to mitocore and Memote report after adding GPRs
    gpr_add=add_gpr_to_model(mitocore, gpr_parsed, 'reactions', '_gpr')

    #--------------------GPR duplicates processing----------------------------------
    #create df with uncleaned GPRs 
    gpr_dict={'model_id':[], 'gpr':[]}

    for reaction in mitocore.reactions:
        gpr=reaction.gene_reaction_rule
        if gpr:
            gpr_dict['model_id'].append(reaction.id)
            gpr_dict['gpr'].append(gpr)
        else:
            gpr_dict['model_id'].append(reaction.id)
            gpr_dict['gpr'].append('')

    gpr_df = pd.DataFrame(gpr_dict)


    # 1. Split GPR on 'or'
    for index, row in gpr_df.iterrows():
        gpr = row['gpr']
        elements_in_gpr = str(gpr).split(' or ')

        # 2. remove duplicates for each reaction with function set()
        unique_parentheses = set() # colect single genes in parenthesis separated by 'and' (ENSG00000115361 and ENSG00000171503) --> ('ENSG00000115361', 'ENSG00000115361')
        unique_single_genes = set()
        cleaned_gprs = [] 
        pattern = r"\((.*?)\)"  # matches content inside parentheses

        # 3. iterate through each element in the GPR: (parenthesis or single gene) 
        for element in elements_in_gpr:
            element = element.strip()

            # genes in  a parehteses such as: (ENSG00000115361 and ENSG00000171503)
            if element.startswith('(') and 'and' in element:
                match = re.search(pattern, element)
                if match:
                    genes = match.group(1).split(' and ')
                    #randomize order with tuple: order-independent for comparisson
                    genes = tuple(sorted(g.strip() for g in genes)) 
                    if genes not in unique_parentheses:
                        # colect single genes in parenthesis separated by 'and' 
                        # (ENSG00000115361 and ENSG00000171503) --> ('ENSG00000115361', 'ENSG00000115361')
                        unique_parentheses.add(genes)
                        #join single genes back with separator 'and' in the cleaned_gprs list
                        cleaned_gprs.append(f"({' and '.join(genes)})")

            # Single gene (not in parentheses)
            elif element not in unique_single_genes:
                unique_single_genes.add(element)
                cleaned_gprs.append(element)

        # 4. update the GPR back to the DataFrame with ' or ' as separator
        gpr_df.at[index, 'gpr'] = ' or '.join(cleaned_gprs)

    # update new GPRs from df to the model, save and test model
    for reaction in mitocore.reactions:
        id = reaction.id
        model_gpr = reaction.gene_reaction_rule

        for index, row in gpr_df.iterrows():
            if row['model_id'] == id:  # Use an if statement to check the condition
                new_gpr = row['gpr']  # Access the '_gpr' column of the current row
                reaction.gene_reaction_rule = new_gpr
                break
    #--------------------ec-codes processing----------------------------------
    # 1.extract ec-codes and corresponding mnx id from metanetx file and create a csv file 
    ec_code = mnx_processing_file_eccode('/Users/benjaminreyes/Desktop/Masterarbeit/neutrophil_modeling/neutrophil-modeling/curation/Files_Databases/reac_prop.tsv', 'metanetx.reaction')

    # 3. extract metanetx ids from model
    print('geting ids for all mnx reactions')
    df_model_id = get_model_ids(mitocore, 'reactions', '_annotation', ['metanetx.reaction'])
    print(df_model_id)

    # 4. map metanetx ids to ec-codes
    print('Mapping  reactions')
    df_react_map = map_ids_to_db(df_model_id, ec_code ,'metanetx.reaction', 'metanetx.reaction', 'ec-code', ',')
    print('ec_code reactions \n',df_react_map)

    # 5. add ec-codes to model and Memote test
    print('adding ec-codes to the model')
    mnx_react = add_ids_to_model(mitocore, df_react_map, 'reactions','ec-code')

     #--------------------overtake Human1 reaction annotations----------------------------------
    for reaction_hum in human1.reactions:
        human1_id = reaction_hum.id

        for reaction_mit in mitocore.reactions:
            notes = reaction_mit.notes

            # check if Human1 mapping exists
            if "Human1" not in notes:
                continue

            if notes["Human1"] != human1_id:
                continue

            # merge annotations: Human1 to MitoCore
            for ann_key, ann_val in reaction_hum.annotation.items():
                # only add if MitoCore does NOT already contain the annotation
                if ann_key not in reaction_mit.annotation:
                    reaction_mit.annotation[ann_key] = ann_val


     #--------------------Genes----------------------------------
    #Human1 ids --> notes 
    Human1_ids=[]

    #loop over genes in mitocore
    for gene in mitocore.genes:
        #get Human1 id if existing otherwise 'None' entry
        human1_id= gene.notes.get('Human1', None)
        # use gene id as Human1 id if no Human1 id existing
        id=gene.id
        Human1_ids.append(id)

        if human1_id is None:
            gene.notes['Human1'] = id

    #ensembl ids --> the annotations
    ensembl_ids=[]
    #loop over genes in mitocore
    for gene in mitocore.genes:
        #get Human1 id if existing otherwise 'None' entry
        ensembl_id= gene.annotation.get('ensembl', None)
        # use gene id as Human1 id if no Human1 id existing
        id=gene.id
        ensembl_ids.append(id)

        if human1_id is None:
            gene.annotation['ensembl'] = id

#--------------------overtake Human1 gene annotations----------------------------------
    for gene_hum in human1.genes:
        human1_id = gene_hum.id

        for gene_mit in mitocore.genes:
            notes = gene_mit.notes

            # check if Human1 mapping exists
            if "Human1" not in notes:
                continue

            if notes["Human1"] != human1_id:
                continue

            # merge annotations: Human1 to MitoCore
            for ann_key, ann_val in gene_hum.annotation.items():
                # only add if MitoCore does NOT already contain the annotation
                if ann_key not in gene_mit.annotation:
                    gene_mit.annotation[ann_key] = ann_val

# API for hsa and kegg genes
    uniprot_ids= get_model_ids(mitocore, 'genes', '_annotation', ['uniprot'])
    # API request function call
    hsa_ids_Genes = get_hsa_from_uniprot(uniprot_ids)

    mapping_df= get_kegg_from_hsa()

        # Match MIRIAM pattern for kegg.genes
    for index, row in mapping_df.iterrows():
        kegg_id = row["kegg.genes"]
        if not kegg_id.startswith("hsa:"):  # Check if prefix is missing
            kegg_id = f"hsa:{kegg_id}"  # Add prefix
            mapping_df.at[index, 'kegg.genes'] = kegg_id
        hsa_id = row["hsa"] 
        if kegg_id.startswith("hsa:"):
            hsa_id = row["hsa"].split(":")[1]
            mapping_df.at[index, 'hsa'] = hsa_id


#-------hsa--------
    print('geting ids for all uniprot')
    df_genes = get_model_ids(mitocore, 'genes', '_annotation', ['uniprot'])
    print(df_genes)

    # 3. map uniprot ids to hsa ids 
    print('Mapping  metabolites to hsa')
    df_hsa_genes = map_ids_to_db(df_genes, hsa_ids_Genes,'uniprot', 'uniprot', 'hsa', ',')
    print('hsa metabolites \n',df_hsa_genes)

    # 4. add hsa ids to model
    print('adding hsa ids to the model')
    hsa_react = add_ids_to_model(mitocore, df_hsa_genes, 'genes','hsa')

#---------kegg--------

    print('geting ids for all uniprot')
    df_genes = get_model_ids(mitocore, 'genes', '_annotation', ['hsa'])
    print(df_genes)

    # 3. map hsa orthology to kegg ids 
    print('Mapping  metabolites to hsa')
    df_kegg_genes = map_ids_to_db(df_genes, mapping_df,'hsa', 'hsa', 'kegg.genes', ',')
    print('hsa metabolites \n',df_kegg_genes)

    # 4. add hsa ids to model
    print('adding kegg ids to the model')
    kegg_genes = add_ids_to_model(mitocore, df_kegg_genes, 'genes','kegg.genes')

#---------MNX-------------
    # load MNX database file for metabolites
    mapping_kegg= mnx_processing_file_mnx('/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Files_Databases/MNX_chem_xref.tsv','kegg.compound','metabolites',( "kegg.compound:", "keggC:"))
    mapping_bigg= mnx_processing_file_mnx('/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Files_Databases/MNX_chem_xref.tsv','bigg.metabolite','metabolites',("bigg.metabolite:", "biggM:"))

    # merge files into one mapping file
    metanetx_mapp_metabolites = pd.merge(mapping_bigg, mapping_kegg, on='metanetx.metabolites', how='outer')
    print(metanetx_mapp_metabolites)
    #rename columns in dataframe to match MIRIAM pattern
    metanetx_mapp_metabolites = metanetx_mapp_metabolites.rename(columns={'bigg.metabolite': 'bigg.metabolite', 'metanetx.metabolites': 'metanetx.chemical', 'kegg.compound': 'kegg.compound'})
    print(metanetx_mapp_metabolites)

    # extract kegg and bigg ids from model into dataframe
    df_metab = get_model_ids(mitocore, 'metabolites', '_annotation', ['kegg', 'bigg'])
    print(df_metab)
    # map kegg ids to Human1 metabolites mapping file
    print('Mapping metabolites to Human1 ids.')
    human1_mets= pd.read_csv('/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Files_Databases/Human1_metabolites.tsv', sep='\t')
    df_Human1_metab = map_ids_to_db(df_metab, human1_mets,'kegg', 'metKEGGID', 'mets','\t')
    print(df_Human1_metab)
    #rename columns to match naming convention
    df_Human1_metab = df_Human1_metab.rename(columns={'model_id': 'model_id', 'mets': 'Human1'})
    print(df_Human1_metab)
    # add Human1 ids to model
    print('adding Human1 ids to the model')
    add_ids_to_model_notes(mitocore, df_Human1_metab,  'metabolites','Human1' )

    # 2. extract kegg and bigg ids from model into dataframe
    print('geting ids for all possible databases ')
    df_metab = get_model_ids(mitocore, 'metabolites', '_annotation', ['kegg.compound', 'bigg.metabolite'])
    print(df_metab)

    # 3. map kegg ids to metanetx ids metabolites mapping file
    print('Mapping  metabolites to metanetx.')
    df_MNX_metab = map_ids_to_db(df_metab, metanetx_mapp_metabolites,'kegg.compound', 'kegg.compound', 'metanetx.chemical', ',')
    print('MNX metabolites \n',df_MNX_metab)

    # add metanetx ids to model
    print('adding MNX ids to the model')
    MNX_metab = add_ids_to_model(mitocore, df_MNX_metab,  'metabolites','metanetx.chemical' )

#--------------------overtake Human1 metabolite annotations----------------------------------
for met_hum in human1.metabolites:
    human1_id = met_hum.id

    for met_mit in mitocore.metabolites:
        notes = met_mit.notes

        # check if Human1 mapping exists
        if "Human1" not in notes:
            continue

        if notes["Human1"] != human1_id:
            continue

        # merge annotations: Human1 to MitoCore
        for ann_key, ann_val in met_hum.annotation.items():
            # only add if MitoCore does NOT already contain the annotation
            if ann_key not in met_mit.annotation:
                met_mit.annotation[ann_key] = ann_val


        # clean model from N/A values, save and run memote test
        model_clean = clean_invalid_annotations(mitocore)
        print(type(model_clean))
        # save cleaned model as new SBML file
        cobra.io.write_sbml_model(model_clean, "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_Human1.xml")
        subprocess.run(
        [
            "memote", "report", "snapshot",
            "--filename", "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_Human1.html",
            "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_Human1.xml",
        ],
        check=True )     

if __name__ == "__main__":
    main()

print("Human1 update finished.")
