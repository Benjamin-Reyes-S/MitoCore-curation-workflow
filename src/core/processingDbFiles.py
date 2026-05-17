import re 
import pandas as pd 



# 1.MetaNetX identifier
def mnx_processing_file_mnx(mnx_file,db_map,component,start_patterns):
    """Process the MNX mapping file to extract relevant columns and save it as a CSV.

    Parameters:
    - mnx_file: metanetx _xref.tsv file to process
    - db_map: reference ids to map for ['bigg', 'kegg'...]
    - component: either 'reactions' or 'metabolites' can be tracked in metanetx database.
    - start_patterns: Tupple of patterns to look for in file for id processing:
        e.g. ( "kegg.compound:", "keggC:") or ("bigg.metabolite:", "biggM:")
    
    Returns:
    - DataFrame with component id and extracted database ids.

    """ 
    #dictionary to collect database ids and respective metanetx ids
    mapping = {db_map: [], f"metanetx.{component}": []}

    #open and read metanetx mapping file
    with open(mnx_file, "r", encoding="utf-8") as file:
        # iterate over lines in file
        for line in file:
            #check for comments (#) and empty lines
            if line.startswith('#') or not line.strip():
                continue
            
            # split line by tab and check if has at least 3 columns
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue  # Skip malformed lines
            
            # mnx id is in 2°nd part out of 3 parts
            source = parts[0]
            target_mnx_id = parts[1]
            rest = parts[2]

            #generalize starting patterns to look for   KEGG or BIGG
            if source.startswith(start_patterns):
                #recon_id is a generic name for database id
                recon_id = source.split(":", 1)[1].strip()
                mnx_id = target_mnx_id.strip()
                mapping[db_map].append(recon_id)
                mapping[f"metanetx.{component}"].append(mnx_id)

    #return dataframe of collecting dictionary
    metanetx_react = pd.DataFrame.from_dict(mapping, orient="columns")
    return metanetx_react
# 2.ec-code
def mnx_processing_file_eccode(mnx_file, db_map):
    """Process the MNX mapping file to extract ec-codes and database ids.

    Parameters:
    - mnx_file: Path to the MetaNetX file with reaction or metabolite data
    - db_map: Database name (e.g., 'bigg', 'kegg') to store as a column

    Returns:
    - DataFrame with db_map and EC code columns (EC is a list if >1, else string)
    """
    # dictionary to collect database ids and respective ec-code
    mapping = {db_map: [], 'ec-code': []}
    # Regular expression pattern to match ec-code
    ec_pattern = re.compile(r"\b\d+\.\d+\.\d+(?:\.\d+|\.n)?\b")
    
    #open and read metanetx mapping file
    with open(mnx_file, "r", encoding="utf-8") as file:
        # iterate over lines in file
        for line in file:
            if line.startswith('#') or not line.strip():
                continue

            # split line by tab and check if has at least 4 columns
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue  # Skip malformed lines

            # ec-code is in 4°nd part out of 4 parts
            mnx_id = parts[0].strip()
            ec_codes = parts[3].strip()


            if ec_pattern.search(ec_codes):
                codes_patter= re.compile(r'^\d+\.\d+\.\d+\.\d+$')
                codes = [code.strip() for code in ec_codes.split(';') if code.strip()]
                codes = [code for code in codes if codes_patter.match(code)]
                # Convert to string if only one ec-code
                formatted_ec = codes[0] if len(codes) == 1 else codes
                mapping[db_map].append(mnx_id)
                mapping['ec-code'].append(formatted_ec)
                print(f"{formatted_ec} is  {type(formatted_ec)}")

    # return dataframe of collecting dictionary for ec-codes and respective mnx ids
    df = pd.DataFrame(mapping)
    return df