import cobra 
from cobra import Model, Reaction, Metabolite
import pandas as pd 
import subprocess

from src.core.parsing import clean_invalid_annotations


mitocore= cobra.io.read_sbml_model("/Users/benjaminreyes/Desktop/Masterarbeit/MitoCore_for_Disease_Modelling/Output_Models_MitoCore/Mitocore_Preliminary.xml")
mitomammal= cobra.io.read_sbml_model("/Users/benjaminreyes/Desktop/Masterarbeit/MitoCore_for_Disease_Modelling/Input_Models/MitoMAMMAL_08.25.xml")


def main():

# classification of reactions in MitoMammal and MitoCore
    mitomammal_in_mitocore= {'mitocore_id':[], 'mitomammal_id':[], 'mitocore_GPR':[], 'mitomammal_GPR':[]}
    mitomammal_not_in_mitocore= { 'mitomammal_id':[], 'mitomammal_GPR':[]}


    for reaction in mitomammal.reactions:
        if reaction.id in mitocore.reactions:
                #mitomammal
                mitomammal_reaction = mitocore.reactions.get_by_id(reaction.id)
                mitmammal_GPR = reaction.gene_reaction_rule
                #mitocore
                mitocore_reaction = mitocore.reactions.get_by_id(reaction.id)
                mitocore_GPR = mitocore_reaction.gene_reaction_rule

                mitomammal_in_mitocore['mitocore_id'].append(mitocore_reaction.id)
                mitomammal_in_mitocore['mitomammal_id'].append(reaction.id)
                mitomammal_in_mitocore['mitocore_GPR'].append(mitocore_GPR)
                mitomammal_in_mitocore['mitomammal_GPR'].append(mitmammal_GPR)

        elif reaction.id not in mitocore.reactions:
                mitomammal_not_in_mitocore['mitomammal_id'].append(reaction.id)
                mitomammal_not_in_mitocore['mitomammal_GPR'].append(reaction.gene_reaction_rule)
                


    mitomammal_in_mitocore_df= pd.DataFrame(mitomammal_in_mitocore)
    mitomammal_not_in_mitocore_df= pd.DataFrame(mitomammal_not_in_mitocore)




    mitocore_not_in_mitomammal= { 'mitomammal_id':[], 'mitomammal_GPR':[]}
    for reaction in mitocore.reactions:
        if reaction not in mitomammal.reactions:
                mitocore_reaction = mitocore.reactions.get_by_id(reaction.id)
                mitocore_GPR = mitocore_reaction.gene_reaction_rule
                mitocore_not_in_mitomammal['mitomammal_id'].append(reaction.id)
                mitocore_not_in_mitomammal['mitomammal_GPR'].append(mitocore_GPR)
    print(f"Number of reactions in mitocore but not in mitomammal: {len(mitocore_not_in_mitomammal)}")
    mitocore_not_in_mitomammal_df= pd.DataFrame(mitocore_not_in_mitomammal)


    reactions=['CBPS', 'ASPCT', 'DHORTS', 'DHORD9', 'DM_orot_c'] #list of  5 new reactions in MitoMAMMAL
    genes_mitomammal= set([])
    metabolites_mitomammal= set([])

    #
    for reaction in mitomammal.reactions:
        if reaction.id in reactions:
            annotations = reaction.annotation
            notes= reaction.notes
            metabolites= reaction.metabolites
            genes= reaction.genes
            
            print("--------------------------------------------------")
            print(reaction.id)
            print("--------------------------------------------------")
            for metabolite in metabolites:
                metabolites_mitomammal.add(metabolite)
                #print(f"for reaction{reaction.id} metabolite.id: {metabolite.id} annotation {metabolite.annotation}, formula {metabolite.formula}")
            print("--------------------------------------------------")
            for gene in genes:
                genes_mitomammal.add(gene)
                #print(f"for reaction{reaction.id} gene.id: {gene.id}", gene.annotation)


    model = mitocore 
    # 1. add reactions from mitomammal to mitocore
    for reaction in mitomammal.reactions:
        if reaction.id in ['CBPS', 'ASPCT', 'DHORTS', 'DHORD9', 'DM_orot_c']:
            print(reaction.id)
            model.add_reactions([reaction])

    # 2. add genes and metabolites from mitomammal to mitocore
    for gene in genes_mitomammal:
        if gene not in mitocore.genes:
            print(gene.id)
            model.genes.append(gene)

    # 3. add metabolites from mitomammal to mitocore
    for metabolite in metabolites_mitomammal:
        if metabolite not in mitocore.metabolites:
            print(metabolite.id)
            model.metabolites.append(metabolite)


    for reaction in mitocore.reactions:
        # 1. update GPR for Complex I
        if reaction.id =='CI_MitoCore':
            gpr_old= reaction.gene_reaction_rule
            gpr_new= mitomammal.reactions.get_by_id('CI_mitoMap').gene_reaction_rule
            print(f"Old GPR: {gpr_old}")
            print(f"New GPR: {gpr_new}")
            reaction.gene_reaction_rule = gpr_new
            print(f"Updated GPR for reaction {reaction.id}: {reaction.gene_reaction_rule}")
        # 2. update GPR for Complex IV
        if reaction.id =='CIV_MitoCore':
            gpr_old= reaction.gene_reaction_rule
            gpr_new= mitomammal.reactions.get_by_id('CIV_mitoMap').gene_reaction_rule
            print(f"Old GPR: {gpr_old}")
            print(f"New GPR: {gpr_new}")
            reaction.gene_reaction_rule = gpr_new
            print(f"Updated GPR for reaction {reaction.id}: {reaction.gene_reaction_rule}")

    # check added reactions
    for reaction in mitocore.reactions:
        if reaction.id in ['CBPS', 'ASPCT', 'DHORTS', 'DHORD9', 'DM_orot_c']:
            annotations = reaction.annotation
            notes= reaction.notes
            metabolites= reaction.metabolites
            bounds= (reaction.lower_bound, reaction.upper_bound)
            print(f"reaction: {reaction.id}, GPR: {reaction.gene_reaction_rule}, annotations: {annotations}")
            print("--------------------------------------------------")

    #check after adding genes and metabolites from MitoMAMMAL

    for metabolite in metabolites_mitomammal:
        #check if metabolites are not in MitoCore (they should be there now after adding them)
        if metabolite not in mitocore.metabolites:
            print(f"metabolite.id: {metabolite.id} annotation {metabolite.annotation}, formula {metabolite.formula}")


    for gene in genes_mitomammal:
        #check if genes are not in MitoCore (they should be there now after adding them)
        if gene not in mitocore.genes:
            print(f"gene.id: {gene.id}", gene.annotation)



    def clean_not_human_genes(model):
        '''Remove non-human genes from the model'''
        non_human_genes = [gene for gene in model.genes if not gene.id.startswith("ENSG")]
        cobra.manipulation.delete.remove_genes(model, non_human_genes, remove_reactions=False)
        if non_human_genes:
            print(f"Removed non-human genes: {[gene.id for gene in non_human_genes]}")    
        return model

   # clean model from N/A values, save and run memote test
    model_clean = clean_invalid_annotations(mitocore)
    print(type(model_clean))
    # save cleaned model as new SBML file
    cobra.io.write_sbml_model(model_clean, "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_MitoMammal.xml")
    subprocess.run(
    [
        "memote", "report", "snapshot",
        "--filename", "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_MitoMammal.html",
        "/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Output_Models_MitoCore/Mitocore_MitoMammal.xml",
    ],
    check=True )

if __name__ == "__main__":
    main()  