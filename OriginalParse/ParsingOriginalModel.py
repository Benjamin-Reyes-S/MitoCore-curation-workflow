import cobra
import subprocess
from src.core.parsing import get_gene_annotations_from_reactions, get_model_ids,add_ids_to_model, clean_invalid_annotations

mitocore = cobra.io.read_sbml_model("/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Input_Models/MitoCore_Original_2017.xml")


def main():
    print("starting original model parsing")

# ----------------------------------------------------------------------------------------------------------------------------
    # KEGG ids from notes for possible databases and collect in dataframe
    print('geting ids for all possible databases ')
    df_kegg_react = get_model_ids(mitocore, 'reactions', 'notes', ['KEGG'],r"R\d+")
    print(df_kegg_react)
    df_react_kegg = df_kegg_react.rename(columns={'model_id': 'model_id','KEGG':'kegg.reaction', 'Recon2': 'bigg.reaction'})
    print(df_react_kegg)
    print('adding kegg ids to the model')
    kegg_react = add_ids_to_model(mitocore, df_react_kegg, 'reactions','kegg.reaction')

    # ----------------------------------------------------------------------------------------------------------------------------
    #  ec-codes from notes for possible databases and collect in dataframe
    print('geting ec-codes from notes ')
    df_ec_react = get_model_ids(mitocore, 'reactions', 'notes', ['EC Number'], r"\d+\.\d+\.\d+\.\d+")
    print(df_ec_react)  
    df_react_ec = df_ec_react.rename(columns={'model_id': 'model_id','EC Number':'ec-code'})
    print(df_react_ec) 
    print('adding ec-codes to the model')
    ec_react = add_ids_to_model(mitocore, df_react_ec, 'reactions','ec-code') 

    #Drop malformed genes
    mitocore.reactions.HtmB_MitoCore.notes['GENE_LIST'] = "UCP2 or UCP3"
    cobra.manipulation.delete.remove_genes(
        mitocore, 
        ["Non-enzymatic", "Non-Enzymatic", "Unknown", "N/A"], remove_reactions = False)
    
    #rewrite variable subunits are added at the end of every complex
    CV_MitoCore_gpr = "(ENSG00000152234 and ENSG00000110955 and ENSG00000165629 and ENSG00000099624 and ENSG00000124172 and ENSG00000116459 and ENSG00000167863 and ENSG00000169020 and ENSG00000154723 and ENSG00000241468 and ENSG00000167283 and ENSG00000249222 and ENSG00000241837 and ENSG00000198899 and ENSG00000228253 and ENSG00000159199) or (ENSG00000152234 and ENSG00000110955 and ENSG00000165629 and ENSG00000099624 and ENSG00000124172 and ENSG00000116459 and ENSG00000167863 and ENSG00000169020 and ENSG00000154723 and ENSG00000241468 and ENSG00000167283 and ENSG00000249222 and ENSG00000241837 and ENSG00000198899 and ENSG00000228253 and ENSG00000135390) or (ENSG00000152234 and ENSG00000110955 and ENSG00000165629 and ENSG00000099624 and ENSG00000124172 and ENSG00000116459 and ENSG00000167863 and ENSG00000169020 and ENSG00000154723 and ENSG00000241468 and ENSG00000167283 and ENSG00000249222 and ENSG00000241837 and ENSG00000198899 and ENSG00000228253 and ENSG00000154518)"

    # iterate through all combinations of subunits
    first_subunits = ["ENSG00000124406", "ENSG00000143515", "ENSG00000081923", "ENSG00000104043", "ENSG00000054793", "ENSG00000166377", "ENSG00000206190", "ENSG00000145246", "ENSG00000068650", "ENSG00000058063", "ENSG00000101974"]
    second_subunits = ["ENSG00000112697", "ENSG00000182107"]

    PCFLOPm_PSFLIPm_PEFLIPm_gpr = ""

    #Discuss with Emanuel the output! shouldn't it work with a While conditional??
    count = 1

    for first_subunit in first_subunits:
        for second_subunit in second_subunits:
            PCFLOPm_PSFLIPm_PEFLIPm_gpr += "(" + first_subunit + " and " + second_subunit + ")"
            #If not the last combination, " or " is appended to separate subunit pairs.
            if count < len(first_subunits) * len(second_subunits):
                PCFLOPm_PSFLIPm_PEFLIPm_gpr += " or "
            count += 1
            print(PCFLOPm_PSFLIPm_PEFLIPm_gpr)
    # replace subunits with new GPR
    mitocore.reactions.CV_MitoCore.gene_reaction_rule = CV_MitoCore_gpr
    mitocore.reactions.PCFLOPm.gene_reaction_rule = PCFLOPm_PSFLIPm_PEFLIPm_gpr
    mitocore.reactions.PSFLIPm.gene_reaction_rule = PCFLOPm_PSFLIPm_PEFLIPm_gpr
    mitocore.reactions.PEFLIPm.gene_reaction_rule = PCFLOPm_PSFLIPm_PEFLIPm_gpr


    
    # ----------------------------------------------------------------------------------------------------------------------------
    # Transfer KEGG metabolite ids
    print('geting ids for all possible databases ')
    df_kegg_met = get_model_ids(mitocore, 'metabolites', 'notes', [ 'KEGG ID'], r"C\d+")
    print(df_kegg_met)
    df_met_kegg = df_kegg_met.rename(columns={'model_id': 'model_id', 'KEGG ID': 'kegg.compound'})
    print(df_met_kegg)
    print('adding bigg ids to the model')
    kegg_met = add_ids_to_model(mitocore, df_met_kegg, 'metabolites','kegg.compound')
    

    # MIRIAM annotaions were only included in the notes of the reactions
    # Mitocore -> MIRIAM
    # HGNC -> hgnc (id)
    # GENE_LIST -> hgnc.symbol
    # Ensembl -> ensembl

    # Added annotations
    gene_id_annotation_dict = get_gene_annotations_from_reactions(mitocore.reactions)
    for gene in mitocore.genes:
        if gene.id in gene_id_annotation_dict.keys():
            gene._annotation = gene_id_annotation_dict[gene.id]

    # load and clean model from N/A values
    model_clean = clean_invalid_annotations(mitocore)
    print(type(model_clean))
    # save cleaned model as new SBML file
    cobra.io.write_sbml_model(model_clean, "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_Original.xml")
    subprocess.run(
    [
        "memote", "report", "snapshot",
        "--filename", "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_Original.html",
        "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_Original.xml",
    ],
    check=True )
    print(" Original model parsing finished")
if __name__== "__main__":
    main()