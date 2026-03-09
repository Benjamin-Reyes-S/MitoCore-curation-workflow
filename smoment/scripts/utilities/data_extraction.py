import time
import pandas
import json
import requests
import cobra

def read_data(path, sheet_names, columns):
    """
    Reads data from an Excel file and merges the data from different sheets into a single DataFrame.
"""
    data = {}

    for sheet_name in sheet_names:
        df = pandas.read_excel(path, sheet_name=sheet_name, usecols=columns)
        df['category'] = sheet_name
        data[sheet_name] = df

    merged_data = pandas.concat(
        [sheet for sheet in data.values()],
        ignore_index=True,
        keys=data.keys(),
        names=columns
    )
    
    return merged_data

def request_mol_weights(protein_list: list[str], output_file_name: str, write_file: bool = True):

    """
    Function to request the molecular weights of a list of proteins from the Uniprot API.
    """

    from pathlib import Path

    response_dict = {}

    if Path(output_file_name).is_file():
        mol_masses_df = pandas.read_excel(output_file_name)

        for index, protein in mol_masses_df.iterrows():
            response_dict[protein['uniprot_accession']] = protein['molar_mass']
        
        return response_dict

    url = 'https://rest.uniprot.org/uniprotkb/search'

    for protein in protein_list:

        params = {
            'query': 'accession:' + protein,
            'fields': 'mass'
        }

        response = json.loads(requests.get(url, params).text)
        try:
            response_dict[protein] = response['results'][0]['sequence']['molWeight']
        except KeyError:
            print(protein)
    
    if write_file:
        response_df = pandas.DataFrame.from_dict(response_dict, orient='index').reset_index().rename(columns={'index': 'uniprot_accession', 0: 'molar_mass'})
        response_df.to_excel(output_file_name, index=False)

    return response_dict

def request_mol_weights_batch(accessions, output_file_name: str, batch_size: int = 100, write_file: bool = True):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    final_results = {}

    # 1. Split the list into chunks of 'batch_size'
    batches = [accessions[i:i + batch_size] for i in range(0, len(accessions), batch_size)]
    
    print(f"Starting retrieval of {len(accessions)} IDs in {len(batches)} batches...")

    for index, batch in enumerate(batches):
        # 2. Construct the OR query for this specific batch
        query = " OR ".join([f"accession:{acc}" for acc in batch])
        
        params = {
            "query": query,
            "fields": "accession,mass",
            "format": "tsv",
            "size": batch_size  # Ensure we get all results in the batch
        }

        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            
            # 3. Parse the TSV lines
            lines = response.text.strip().splitlines()
            if len(lines) > 1:
                # Skip header, map Accession (col 0) to Mass (col 1)
                batch_dict = {}
                for line in lines[1:]:
                    if len(line.split('\t')) >= 2 and not line.split('\t')[1] == '':  # Ensure there are at least 2 columns and the second column is a number
                        batch_dict[line.split('\t')[0]] = int(line.split('\t')[1])
                    else:
                        print(f"Skipping malformed line in batch {index + 1}: {line}")
                        
                final_results.update(batch_dict)
            
            print(f"Batch {index + 1}/{len(batches)} completed.")
            
            # Optional: Short sleep to be polite to the API
            time.sleep(0.1) 

        except requests.exceptions.RequestException as e:
            print(f"Error fetching batch {index + 1}: {e}")
    
    if write_file:
        response_df = pandas.DataFrame.from_dict(final_results, orient='index').reset_index().rename(columns={'index': 'uniprot_accession', 0: 'molar_mass'})
        response_df.to_excel(output_file_name, index=False)
    
    return final_results

def write_mol_weights_json(model: cobra.Model, mol_weights_dict, output_path):
    gene_to_uniprot_dict = gene_id_to_uniprot_and_pw(model)

    gene_to_mol_weights_dict = {}

    for gene in gene_to_uniprot_dict.keys():
        uniprot_accession = gene_to_uniprot_dict[gene]['uniprot']
        
        if uniprot_accession is None:
            print(f'No uniprot accession found for gene {gene}, skipping...')
            continue
        
        mol_weight = mol_weights_dict[uniprot_accession]

        gene_to_mol_weights_dict[gene] = int(mol_weight)
    
    with open(output_path, 'w') as fp:
        json.dump(gene_to_mol_weights_dict, fp)

def gene_id_to_uniprot_and_pw(model: cobra.Model):

    pathway_groups = model.groups

    gene_to_uniprot_and_pw = {}

    for gene in model.genes:
       
        pathways = []

        # for every reation the gene connects to, get the pathways the reaction participates in
        for reaction in gene.reactions:
            for group in pathway_groups:
                if reaction in group.members:
                    pathways.append(group.name)
        
         # enter information into dict
        gene_to_uniprot_and_pw[gene.id] = {
            'uniprot': gene._annotation['uniprot'] if 'uniprot' in gene._annotation else None,
            'pathway_groups': pathways
        }

    return gene_to_uniprot_and_pw
        
def write_smoment_parameters_excel(
    abundances: pandas.Series,
    total_protein_fraction: float,
    model_protein: float,
    output_path: str,
    average_saturation: float = 1.0):
    # smoment global protein parameters
    total_protein_data_dict = {
        "Total protein content [g/gDW]:":total_protein_fraction,
        "Fraction of masses of model-included enzymes in comparison to all enzymes (0.0 to 1.0):": model_protein,
        "Average saturation level (0.0 to 1.0):": average_saturation,
        }
            
    # smoment individual protein parameters
    single_protein_data_dict = {
        "Protein ID (as in SBML model)": list(abundances.index),
        "Protein concentration [mmol/gDW]": list(abundances.values)
    }

    # write to xlsx
    with pandas.ExcelWriter(output_path + '.xlsx') as writer:
        pandas.DataFrame.from_dict(total_protein_data_dict, orient='index').to_excel(writer, sheet_name='Total protein data', header=False)
        pandas.DataFrame.from_dict(single_protein_data_dict).to_excel(writer, sheet_name='Single protein data', index=False)

def write_protein_data_for_smoment(
    transformed_abundances: pandas.Series,
    smoment_parameters: pandas.DataFrame,
    output_path):

    """
    Writes the protein data in the format for the '_protein_data.xlsx' files generated by sMOMENTs 'get_initial_spreadsheets' function. """

    total_proteins_mass = smoment_parameters['total_protein_mass'].values[0]
    model_enzymes_mass = smoment_parameters['model_included_enzymes_mass'].values[0]
    
    write_smoment_parameters_excel(
        transformed_abundances.dropna(),
        smoment_parameters['fraction_of_protein_biomass'].values[0],
        model_enzymes_mass / total_proteins_mass,
        output_path,
        average_saturation=1.0,
        )

def calculate_smoment_data_inputs(
        dataset_df_mean,
        molar_masses,
        ens_to_uniprot_df,
        dw_per_cell,
        output_path
        ):
    """
    Calculates all information from protein data required as input for sMOMENT: individual protein abundances in mmol/gDW,
    protein fraction in biomass, mass fraction of model included enzymes compared to total protein mass fraction.
    """
    
    avogadro = 6.02214076e23
    
    # output data structures
    smoment_parameters = {
        'fraction_of_protein_biomass': [],
        'fraction_of_model_included_enzymes': [],
        'total_protein_mass': [],
        'model_included_enzymes_mass': []
        }
    

    # output df for individual protein data for control
    result_dict = {'ens_id': [], 'uniprot_accession': [], 'copies_cell': [], 'mol': [], 'mmol/gDW': []}
            
    # parameters
    total_protein_mass = 0
    model_included_enzymes_mass = 0

    for index, row in dataset_df_mean.iterrows():        
        uniprot_accession = row['uniprot_accession']
        copy_number = row['copies_cell']

        protein_mol = copy_number / avogadro
        
        try:
            mol_weight = molar_masses[uniprot_accession]
        except KeyError:
            print(f'No molecular weight found for {uniprot_accession}')
            continue

        protein_g = protein_mol * mol_weight

        # if the protein is in the model, put it on the output list
        if uniprot_accession in ens_to_uniprot_df['uniprot'].values:
            
            # transfer copies/cell to mmol/gDW
            protein_mmol_gdw = protein_mol / dw_per_cell * 1e3 # 1e3 to convert mol to mmol

            ens_gene = str(ens_to_uniprot_df[ens_to_uniprot_df['uniprot'] == row['uniprot_accession']].index[0])

            result_dict['ens_id'].append(ens_gene)
            result_dict['uniprot_accession'].append(uniprot_accession)
            result_dict['copies_cell'].append(copy_number)
            result_dict['mol'].append(protein_mol)
            result_dict['mmol/gDW'].append(protein_mmol_gdw)
            
            # sum weights of proteins mapped to model genes, divide by total protein weight later
            model_included_enzymes_mass += protein_g
        
        # calculate total protein weight/cell for protein fraction
        total_protein_mass += protein_g

    smoment_parameters['fraction_of_protein_biomass'].append(total_protein_mass / dw_per_cell)
    smoment_parameters['fraction_of_model_included_enzymes'].append(model_included_enzymes_mass / total_protein_mass)
    smoment_parameters['total_protein_mass'].append(total_protein_mass)
    smoment_parameters['model_included_enzymes_mass'].append(model_included_enzymes_mass)
        
    
    # create dataframe with all protein abundances from all datasets for sMOMENT
    result_sMOMENT_df = pandas.DataFrame(result_dict)

    # write smoment parameters
    smoment_parameters_df = pandas.DataFrame(smoment_parameters)
    
    transformed_df = result_sMOMENT_df.set_index('ens_id')
    transformed_df = transformed_df.loc[:,'mmol/gDW']

    # transfer data to format for sMOMENT/combine data from the two possible dataset types according to dataset_combinations
    write_protein_data_for_smoment(transformed_df, smoment_parameters_df, output_path)
