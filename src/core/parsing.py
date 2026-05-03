#specific dependencies for the parsing module
import re
import pandas as pd
import os

# 1. extract ids from model  
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
