import cobra
from src.core.parsing import get_model_ids
#Load SBML model
mitocore = cobra.io.read_sbml_model("/Users/benjaminreyes/Desktop/Projects/MitoCore_Modular_Curation/Input_Models/MitoCore_Original_2017.xml")

def main():

    #Test 
    #reactions:'Recon2 id': 'r1447', 'KEGG id': 'R03857'
    # BiGG ids from notes for possible databases and collect in dataframe
    print('geting ids for all possible databases ')
    df_bigg_react = get_model_ids(mitocore, 'reactions', 'notes', [ 'Recon2'], r"^[A-Za-z0-9_]+$")
    print(df_bigg_react)

if __name__== "__main__":
    main()