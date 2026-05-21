import cobra
import subprocess
import pandas as pd

from src.core.parsing import get_model_ids, clean_invalid_annotations, add_ids_to_model
from src.core.APIRequests import get_uniprot_from_hgnc

mitocore = cobra.io.read_sbml_model("/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_Original.xml")


print("Starting preliminary curation.")
def main ():

    #convert units to standard mmol/h gDW (from umol/min gDW) by multiplying with factor 0.06
    for reaction in mitocore.reactions:
        ub = reaction.upper_bound
        lb = reaction.lower_bound
        if not ub == 1000: # keep default constraints
            reaction.upper_bound = ub * 0.06 # factor 0.6 for unit convertion umol/min gDW -> mmol/h gDW  
        if not lb == -1000:
            reaction.lower_bound = lb * 0.06
        #why keeping default bounds and changing them just if not == 1000?

        # check if any ubs are -1000 or lbs are 1000 (just to be save)
        if reaction.upper_bound == -1000 or reaction.lower_bound == 1000:
            print(reaction.id)
    

    # load BIGG metabolites model
    bigg_metabolites = pd.read_csv(
        "/Users/benjaminreyes/Desktop/Masterarbeit/MitoCore_for_Disease_Modelling/Files_Databases/bigg_models_metabolites.txt",
        sep="\t"
    )
    print("BIGG metabolites txt file loaded.")

    # mapping dictionary to store bigg id to kegg id {bigg:[kegg1, kegg2....]}
    bigg_id_mapping = {}

    for index, row in bigg_metabolites.iterrows():

        db_links = row["database_links"]
        bigg = row["universal_bigg_id"]

        if pd.isna(db_links):
            bigg_id_mapping[bigg] = None
            continue

        db_entries = db_links.split("; ")
        kegg_ids = []

        for entry in db_entries:
            if "KEGG Compound" not in entry:
                continue

            identifier = entry.split("identifiers.org/")[1]
            kegg_id = identifier.split("kegg.compound/")[1]

            kegg_ids.append(kegg_id)

        #no kegg ids
        if len(kegg_ids) == 0:
            continue

        # manage multiple kegg ids
        for kegg in kegg_ids:
            bigg_id_mapping[kegg] = bigg

    print(bigg_id_mapping)


    for met in mitocore.metabolites:
        annotations = met.annotation

        for key, val in list(annotations.items()):
            if key == 'kegg.compound':
                kegg_list = val
                for kegg in kegg_list: 
                    if kegg in bigg_id_mapping.keys():
                        met.annotation['bigg.metabolite'] = bigg_id_mapping[kegg]

    # Manual curation of metabolites with mutiple bigg ids
    met.annotation['bigg.metabolite'] = "2hb"

    met = mitocore.metabolites.get_by_id("2hb_e")  # used in recon 3d
    met.annotation['bigg.metabolite'] = "2hb"

    met = mitocore.metabolites.get_by_id("3hbcoa_m")  # used in recon 3d
    met.annotation['bigg.metabolite'] = "3hbcoa"

    met = mitocore.metabolites.get_by_id("3hibutcoa_m")  # used in recon 3d
    met.annotation['bigg.metabolite'] = "3hibutcoa"

    met = mitocore.metabolites.get_by_id("dd2coa_m")  # used in recon 3d
    met.annotation['bigg.metabolite'] = "dd2coa"

    # API request to get uniprot ids from hgnc ids
    df_hgnc_genes = get_model_ids(mitocore, 'genes', 'annotation', [ 'hgnc'], r"HGNC:\d+")
    df_uniprot_genes= get_uniprot_from_hgnc(df_hgnc_genes)
    print(df_uniprot_genes) 
    kegg_met = add_ids_to_model(mitocore, df_uniprot_genes, 'genes','uniprot')

    #count genes with uniprot annotation
    counter=0
    for gene in mitocore.genes:
        if 'uniprot' in gene.annotation:
            counter += 1
    print(f"Total genes with uniprot ids: {counter}")

    # clean model from N/A values, save and run memote test
    model_clean = clean_invalid_annotations(mitocore)
    print(type(model_clean))
    # save cleaned model as new SBML file
    cobra.io.write_sbml_model(model_clean, "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_Preliminary.xml")
    subprocess.run(
    [
        "memote", "report", "snapshot",
        "--filename", "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_Preliminary.html",
        "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_Preliminary.xml",
    ],
    check=True )


if __name__ == "__main__":
    main()

print("Preliminary curation finished.")