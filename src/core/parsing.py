#specific dependencies for the parsing module
import re
import pandas as pd
import os
import ast
import pandas as pd
import cobra

#  extract ids from model  
def get_model_ids(model, component_type, dict, databases, pattern):
    """
    Extract ids from model annotations for various databases across a specified component type.
    
    Parameters:
    - model: COBRApy model object (mitocore)
    - component_type: str, one of 'reactions', 'metabolites', 'genes'
    - dict: dictionary where database ids are located ('annotations', 'notes')
    - databases: list of str, e.g., ['kegg', 'bigg', 'metanetx', 'uniprot'].
    - pattern: regex pattern to extract ids from annotation strings. e.g for KEGG: r"R\d+"
    
    Returns:
    - DataFrame with component ID and extracted database ids.
    """
    # Initialize dictionary
    ids_dict = {'model_id': []}
    for db in databases: 
        ids_dict[db] = []

    # Access model core-component (reactions, genes, metabolites)
    components = getattr(model, component_type)

    for comp in components:
        ids_dict['model_id'].append(comp.id)

        for db in databases:
            ids= None
            value = None
            annotation = getattr(comp, dict, {})

            # look for different naming conventions 
            for key in [db, 
                        f"{db}.reaction", 
                        f"{db}.compound", 
                        f"{db}.metabolite", 
                        f"{db}.chemical", 
                        f"{db}.genes",
                        f"{db}.gene",
                        #reactions:'Recon2 id': 'r1447', 'KEGG id': 'R03857'
                        #metabolites: 'RECON2': '2oxoadp', 'KEGG ID': 'C00322'                        
                        f"{db} id",
                        f"{db.upper()} id",
                        f"{db.upper()} ID",
                        f"{db.upper()}"
                        ]:
                # search for ids and extract them as lis with re.findall 
                if key in annotation:
                    ids = re.findall(pattern, annotation[key])
                    if len(ids) == 0:
                        continue
                    break

            ids_dict[db].append(ids)
            os.makedirs("Test_output", exist_ok=True)
            pd.DataFrame(ids_dict).to_csv(f"Test_output/{component_type}_ids.csv", index=False)

    return pd.DataFrame(ids_dict)

#  map them to a target database in a mapping file
def map_ids_to_db(id_df, input_df, id_col_input, id_col_map, target_db_col, sep_map_file):
    """
    Map ids from a column in a dataframe to another column of a dataframe and add the target identifier.
    
    Parameters:
    - id_df: input Data frame with model_ids and ids to a specific database to be mapped.
    -input_df: dataframe for mapping ids to the same database as id_df and target database
    -id_col_input: column in id_df that contains ids to be mapped.
    -id_col_map: the column in mapping_file that contains the ids to be mapped.
    -target_db_col: the column in mapping_file that contains the target database ids to be extracted.
    -sep_map_file: the separator used in the mapping file (default is ',').
    
    Returns:
    - DataFrame with model ids and extracted target database ids.
    """
    import pandas as pd
    ids_dict = {'model_id': [],target_db_col: []}

    #load csv mapping file
    mapping_df = input_df

    # iterate over ids 
    for index, row in id_df.iterrows():
        id = row[id_col_input]

        # check for model ids are present in mapping file 
        match= mapping_df[mapping_df[id_col_map] == id]
        # if id model existing extract target id (target database)
        if not match.empty:
            target_id= match[target_db_col].values[0]
            print(f"Match found: {id} -> {target_id}")
        # if id model not found extract 'None' value
        else:
            target_id= None
            print(f"No match found for {id}, using default value None")

        #collect target ids into ids dictionary
        ids_dict['model_id'].append(row['model_id'])
        ids_dict[target_db_col].append(target_id)
    
    #return dictionary as dataframe
    ids_dict= pd.DataFrame(ids_dict)
    return ids_dict

    # 3. add them to current model 
    # compare which add function is robuster (this or Human1 alignment defined one)
def add_ids_to_model(model, df, component_type, database):
        """
        Add ids from a DataFrame to the model annotations for a specified component type.

        Parameters:
        - model: COBRApy model object (e.g., mitocore)
        - df: DataFrame with model_id and database id columns.
        - component_type: str, one of 'reactions', 'metabolites', 'genes'
        - database: str, the name of the database column in df to map (e.g., 'kegg.reaction')
        
        Returns:
        - Updated model with new annotations added.
        """

        components = getattr(model, component_type)

        # iterate over dataframe rows
        for _, row in df.iterrows():
            model_id = row["model_id"]
            db_value = row[database]

            # normalize db_value 
            # skip missing values
            if db_value is None or (isinstance(db_value, float) and pd.isna(db_value)):
                continue

            # parse stringified lists (e.g. "['R00209', 'R01699']")
            if isinstance(db_value, str) and db_value.startswith("["):
                try:
                    db_value = ast.literal_eval(db_value)
                except Exception as e:
                    print(f"Error parsing {db_value}: {e}")
                    continue

            # enforce list semantics
            if not isinstance(db_value, list):
                db_value = [db_value]

            if len(db_value) == 0:
                continue

            # update model 
            for comp in components:
                if comp.id == model_id:
                    existing_annotation = comp.annotation.get(database)

                    # normalize existing annotation to list
                    if existing_annotation is None:
                        merged_values = db_value
                    else:
                        # if model contain local ids, it eliminate duplicates
                        if not isinstance(existing_annotation, list):
                            existing_annotation = [existing_annotation]
                        merged_values = list(
                            dict.fromkeys(existing_annotation + db_value)
                        )

                    comp.annotation[database] = merged_values
                    print(f"Added {merged_values} ({type(merged_values)}) to {model_id}")
                    break

        return model

#  get attribute (not inside of attribute 'annotation', e.g.GPR rule)
def get_model_ids_for_attr(model, component_type, attr_name):
    """
    Extract attribute in the model in the target core-component.
    
    Parameters:
    - model: COBRApy model object (mitocore)
    - component_type: Core component of GEMs, str, one of 'reactions', 'metabolites', 'genes'
    - attr_name: name of the dictionary where target object is located ('notes', '_gpr')

    
    Returns:
    - DataFrame with component ID and extracted database ids.
    """
    # dictionary to collect ids and attributes
    ids_dict = {'model_id': [], attr_name:[]}

    # Access model component (e.g., reactions, genes)
    components = getattr(model, component_type)

    # iterate over components and extract attribute
    for comp in components:
        #add model id and attribute to dictionary
        ids_dict['model_id'].append(comp.id)
        value = comp.__dict__.get(attr_name, None) #None in case no attribute present
        ids_dict[attr_name].append(value) 

    #returns dataframe with attributes
    return pd.DataFrame(ids_dict)

# add ids to model notes (not inside of 'annotation' attribute)
def add_ids_to_model_notes(model, df, component_type, database):
    """
    Add ids from a DataFrame to the model annotations for a specified component type.
    new defined in automation pipeline to handle not stringtified lists in dataframes insted of .csv stringtified lists 

    Parameters:
    - model: COBRApy model object (e.g., mitocore)
    - df: DataFrame with model_id and database id columns.
    - component_type: str, one of 'reactions', 'metabolites', 'genes'
    - database: str, the name of the database column in df to map (e.g., 'kegg', 'bigg')
 
    
    Returns:
    - Updated model with new annotations added (only if not already present).
    """
    import ast
    components = getattr(model, component_type)

    for index, row in df.iterrows():
        model_id = row['model_id']
        db_value = row[database]

        # skip missing or None values 
        if db_value is None or (isinstance(db_value, float) and pd.isna(db_value)):
            continue

        # convert stringified lists (from CSVs) into real lists
        if isinstance(db_value, str) and db_value.startswith("[") and db_value.endswith("]"):
            try:
                db_value = ast.literal_eval(db_value)
            except Exception as e:
                print(f"Error parsing {db_value}: {e}")
                continue

        # now db_value is either a string or a real list
        for comp in components:
            if comp.id == model_id:
                if isinstance(db_value, list) and len(db_value) == 1:
                    db_value = db_value[0]
                comp.notes[database] = db_value
                print(f"Added {db_value} (type {type(db_value)}) to {model_id}")
                break

    return model


def get_model_ids_for_gpr(model, component_type, attr_name):

    ids_dict = {'model_id': [], attr_name: []}

    components = getattr(model, component_type)

    for comp in components:
        ids_dict['model_id'].append(comp.id)

        # SPECIAL CASE FOR GPR
        if attr_name in ['gpr', 'gene_reaction_rule', '_gpr']:
            value = comp.gene_reaction_rule
        else:
            value = getattr(comp, attr_name, None)

        ids_dict[attr_name].append(value)

    return pd.DataFrame(ids_dict)

# add GPR rule to the model and complete if already present
def add_gpr_to_model(model, df, component_type, gpr_header_in_df):
    """
    check if attribute GPR is empty in the model or already present,
    Add ids from a dataFrame to the gpr attribute if empty or append to existing GPR rule.

    Parameters:
    - model: COBRApy model object (e.g., mitocore)
    - df: DataFrame with model_id and formulas column.
    - component_type: str, one of 'reactions', 'metabolites', 'genes'
    - gpr_header_in_df: str, the name of the database column in df to map (e.g., '_gpr')
    
    Returns:
    - Updated model with new annotations added (only if not already present).
    """
    # Access model component (e.g., reactions, genes)
    components = getattr(model, component_type )

    # iterate over ids (rows)
    for index, row in df.iterrows():
        model_id = row['model_id']
        new_rule = str(row[gpr_header_in_df]).strip() if pd.notnull(row[gpr_header_in_df]) else ''

        # find matching component based on model_id
        for comp in components:
            if comp.id == model_id and model_id not in ['CBPS', 'ASPCT', 'DHORTS', 'DHORD9', 'DM_orot_c', 'CI_MitoCore', 'CIV_MitoCore']:  # Exclude newly added reactions
                existing_rule = comp.gene_reaction_rule.strip()
                #check if GPR is empty)
                if existing_rule == '' and new_rule != '':
                    comp.gene_reaction_rule = new_rule
                    print(f"Added GPR: {comp.id} -> {new_rule}")
                #check if GPR already exists and append new rule with 'or' (does not avoid duplicates, this issue is corrected later)
                elif existing_rule != '' and new_rule != '':
                    comp.gene_reaction_rule = existing_rule + ' or ' + new_rule
                    print(f"GPR already exists for {comp.id}: {existing_rule}. New rule {new_rule} added.")
                break  # Stop once matched
    return model

#  splits word separated by "and", "or", "/", " " into strings in a list
def split_gene_rule(gene_string):
        """Strips all brackets and "and" and "or" from gene rule and splits it into a list of gene codes"""
        space_split = gene_string.split(" ")

        if len(space_split) == 1:
            return space_split

        space_split = [re.sub(r'[\(\),]', '', string) for string in space_split]

        split = [string for string in space_split if not string == "or" and not string == "and"]

        return split

 # 5 parse reaction notes to get gene annotations 
def parse_reaction(reaction, gene_hgnc_mapping, reaction_stats):

            notes = reaction.notes
        
            ensebml_annotation_string = "" #to store GENE_ASSOCIATION ids
            gene_hgnc_string = "" #to store reaction id with hgnc gene code
            gene_name_string = "" # to store GENE_LIST

            non_gene_annotations = ["Unknown", "N/A", "Non-Enzymatic", "Non-enzymatic"]

            #parse notes 
            #all stored in "reaction_stats" is for statistical counting 
            for key in notes.keys():
                reaction_stats["annotations"].add(key)
                # count reactions with hgnc gene code
                if key == "HGNC" or key == "(HGNC":
                    reaction_stats["reactions_with_hgnc"].add(reaction.id)
                    gene_hgnc_string = notes[key]                

                # count reactions with gene associations
                # handle case of non_gene_annotations where "GENE_ASSOCIATION" key exists but is N/A or unknown
                if key == "GENE_ASSOCIATION" and not (notes[key] in non_gene_annotations):
                    reaction_stats["reactions_with_genes"].add(reaction.id)

                    ensebml_annotation_string = notes[key]

                if key == "GENE_LIST":
                    gene_name_string = notes[key]
            
            # check how many reactions have gene and hgnc code (pure statistic quantification)
            if gene_hgnc_string == "":
                reaction_stats["reaction_without_hgnc"].add(reaction.id)
            if gene_hgnc_string == "" and ensebml_annotation_string == "":
                reaction_stats["reaction_without_genes"].add(reaction.id)

            ensembl_code = []
            hgnc_ids = []
            hgnc_symbols = []

            #if not emoty string in following cases, process the IDs with split_gene_rule() function
            if not ensebml_annotation_string == "":
                ensembl_code = split_gene_rule(ensebml_annotation_string)
                
            if not gene_hgnc_string == "":
                hgnc_ids = split_gene_rule("HGNC:" + gene_hgnc_string)

            if not gene_name_string == "":
                hgnc_symbols = split_gene_rule(gene_name_string)

            #dictionary to store parameters obtained from hgnc_ids
            for ensembl, hgnc_id, hgnc_symbol in zip(ensembl_code, hgnc_ids, hgnc_symbols):
                gene_hgnc_mapping[ensembl] = {
                    "hgnc": hgnc_id,
                    "hgnc.symbol": hgnc_symbol,
                    "ensembl": ensembl
                }

#  Main function to get gene annotations from reactions
def get_gene_annotations_from_reactions(reactions):

        '''Parses notes in reaction objects, retrieves uniprot ids for each gene id in the gpr rules and returns dict containing mapping of gene id to gene name, hgnc id and uniprot id'''
        # gene code in gpr -> hgnc, gene name
        gene_hgnc_mapping = {}

        reaction_stats = {
            "annotations": set(),           # all annotations in reactions
            "reactions_with_genes": set(),  # reactions with gene rules
            "reactions_with_hgnc": set(),   # reactions with hgnc gene codes
            "reaction_without_hgnc": set(), # reactions without hgnc gene codes
            "reaction_without_genes": set() # reactions without gene rules
        }
        #
        reaction_index = 0
        for reaction in reactions:
            print("                                                                       ", end="\r")
            print("Parsing reaction " + str(reaction_index) + "/" + str(len(reactions)) + "(reaction id: " + reaction.id + ")", sep=" ", end='\r')
            parse_reaction(reaction, gene_hgnc_mapping, reaction_stats)
            reaction_index += 1

        print("occuring annotations: " + str(reaction_stats["annotations"]))

        print("reactions with gene rule / reactions with hgnc ids: " + str(len(reaction_stats["reactions_with_genes"])) + "/" + str(len(reaction_stats["reactions_with_hgnc"])))

        return gene_hgnc_mapping

# clean model from N/A values function
def clean_invalid_annotations(model, invalid_values={'N/A', 'None', 'nan', 'NaN','n/a',''}):
    """
    Remove invalid placeholder values from model annotations across reactions, metabolites, and genes to avoid false positives in MEMOTE test.
    
    -model: COBRApy model object (already read with COBRApy, e.g. Mitocore)
    -invalid_values: dict, keys to be removed
    """
    #acces model's objects core-components 
    for component in model.reactions + model.metabolites + model.genes:
        #collect keys to remove
        keys_to_remove = []
        for key, val in component.annotation.items():
            if isinstance(val, (str, type(None))):
                #check if value inside key is in invalid dict
                if str(val).strip() in invalid_values:
                    keys_to_remove.append(key)
        for key in keys_to_remove:
            del component.annotation[key]
    return model

def clean_not_human_genes(model):
    '''Remove non-human genes from the model'''
    non_human_genes = [gene for gene in model.genes if not gene.id.startswith("ENSG")]
    cobra.manipulation.delete.remove_genes(model, non_human_genes, remove_reactions=False)
    if non_human_genes:
        print(f"Removed non-human genes: {[gene.id for gene in non_human_genes]}")    
    return model

