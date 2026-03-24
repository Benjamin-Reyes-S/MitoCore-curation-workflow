import re
from typing import Dict

import numpy
import pandas
import cobra
import statistics
import math
import matplotlib.pyplot as plt
from cobra.util import linear_reaction_coefficients
from optlang.symbolics import Zero


from utilities.kcat_sensitivity import get_split_reactions_in_direction, scale_pseudo_metabolite_stoichiometry

def get_default_kcat(reactions_kcat_mapping_database, type_of_default_kcat_selection):
    
    # Get all kcats which are not math.nan and calculate the median of them, which will be used as default kcat
    all_kcats = [x["forward"] for x in reactions_kcat_mapping_database.values()] + [
        x["reverse"] for x in reactions_kcat_mapping_database.values()
    ]
    all_kcats = [x for x in all_kcats if not math.isnan(x)]

    if type_of_default_kcat_selection == "median":
        default_kcat = statistics.median(all_kcats)
    elif type_of_default_kcat_selection == "mean":
        default_kcat = statistics.mean(all_kcats)
    elif type_of_default_kcat_selection == "random":
        default_kcat = math.random.choice(all_kcats)
    
    return default_kcat

def get_original_reaction_id(smoment_reaction_id: str):
    
    # print(smoment_reaction_id)
    
    reverse = False
    
    if '_TG_' in smoment_reaction_id:
        reverse = smoment_reaction_id.split('_TG_')[1] == 'reverse'
    
    isoenzyme_split_index = -1
    
    if '_GPRSPLIT_' in smoment_reaction_id.split('_TG_')[0]:
        isoenzyme_split_index = smoment_reaction_id.split('_TG_')[0].split('_GPRSPLIT_')[1]
    
    return smoment_reaction_id.split('_TG_')[0].split('_GPRSPLIT_')[0], reverse, isoenzyme_split_index

def get_deliveries_for_reaction_metabolites(reaction: cobra.Reaction):
    
    enzymes = []
    
    for metabolite in reaction.metabolites.keys():
        if metabolite.id.startswith('ENZYME_'):
            gene = metabolite.id.split('ENZYME_')[1]
            delivery_reaction = f'ENZYME_DELIVERY_{gene}'
            enzymes.append(delivery_reaction)

    return enzymes

def get_barplot_for_fva(fba_results, fva_results, reactions, limits=[-1100, 1100]):
    
    fva_and_fba = fva_results.loc[reactions, :]
    fva_and_fba['fba'] = fba_results.fluxes.loc[reactions]
    fva_and_fba.reset_index(inplace=True)
    fva_and_fba.columns = ['reaction', 'minimum', 'maximum', 'fba'] # rename columns
    fva_and_fba.sort_values(by='reaction')
    
    fig, ax = plt.subplots()

    ax.barh(
        y = fva_and_fba['reaction'],
        width = fva_and_fba['maximum'] - fva_and_fba['minimum'],
        left = fva_and_fba['minimum'],
        align='center',
        color = 'lightblue',
        edgecolor = 'black',
        alpha = 0.5,)
    ax.scatter(
        x = fva_and_fba['fba'],
        y = fva_and_fba['reaction'],
        color = 'black'
    )
    ax.set_xlim(limits)
    ax.set_xlabel('Flux (mmol/(gdw h))')
    ax.set_ylabel('Reactions')
    ax.set_title('Flux variability analysis (FVA) and flux balance analysis (FBA)')

    plt.gca().invert_yaxis()

    plt.show()
    
def apply_kcats_to_model(model_id, model_path, protein_pool, kcat_factor_dict) -> cobra.Model:

    model =  cobra.io.read_sbml_model(model_path)

    print('\n\nScaling ' + model_id)
    
    # step 1) insert calibrated prot pool as stoichiometric coefficient

    prot_pool_delivery = model.reactions.get_by_id('ER_pool_TG_')
    prot_pool_met = list(prot_pool_delivery.metabolites.keys())[0]

    new_met_dict = {prot_pool_met: protein_pool}

    prot_pool_delivery.add_metabolites(new_met_dict, combine = False )

    print('Scaled protein pool to ' + str(model.reactions.get_by_id('ER_pool_TG_').metabolites))

    # step 2) scale reactions with calibrated kcats

    for reaction, factor in kcat_factor_dict.items():

        reaction_id = re.sub('^R_','', reaction)
        reaction_obj = model.reactions.get_by_id(reaction_id)
            
        print('\nScaled ' + reaction_id + ' from ')

        unscaled_reaction_equation = str(reaction_obj)
        print(unscaled_reaction_equation + ' to ')

        # scaling factor in this functions refers to kcat which is in the denominator
        # the matlab script scales the numerator
        # for correct application, the reciprokal value needs to be used
        scale_pseudo_metabolite_stoichiometry(reaction_obj, factor)

        print(str(model.reactions.get_by_id(reaction_id)))
                
    # step 3) Perform a test simulation

    current_solution = model.optimize()
    current_objective = current_solution.fluxes['OF_ATP_MitoCore']

    print('\nobjective value for ' + model_id + ': ' + str(current_objective))

    return model

def ko_reverse_reactions(model, reaction, reaction_id_pattern):
    # knocks out all split reverse reactions
    if "_GPRSPLIT_" in reaction.id and reaction.id.endswith("_forward"):
        for reac in model.reactions:
            # if pattern.match(reac.id) and not reac.id == reaction.id
            if reaction_id_pattern.match(reac.id) and reac.id.endswith("_reverse"):
                model.reactions.get_by_id(reac.id).knock_out()
    # knocks out all split forward reactions
    elif "_GPRSPLIT_" in reaction.id and reaction.id.endswith("_reverse"):
        for reac in model.reactions:
            if reaction_id_pattern.match(reac.id) and reac.id.endswith("_forward"):
                model.reactions.get_by_id(reac.id).knock_out()
    # if reactions were not split up:
    elif reaction.id.endswith("reverse"):
        forward_reaction_name = reaction.id.replace(
            "reverse", "forward")
        model.reactions.get_by_id(forward_reaction_name).knock_out()
    elif reaction.id.endswith("forward"):
        reverse_reaction_name = reaction.id.replace(
            "forward", "reverse")
        model.reactions.get_by_id(reverse_reaction_name).knock_out()
            
def fix_objective_reactions(model, fraction_of_optimum):
    sol = model.optimize()
    if sol.status != 'optimal':
        return # Handle error
    
    target_value = sol.objective_value * fraction_of_optimum
    
    # Create the constraint expression (sum of coeff * flux)
    obj_expression = model.objective.expression
    
    # Define a new constraint that the current objective must be at least target_value
    constraint = model.problem.Constraint(
        obj_expression,
        lb=target_value,
        name="fva_objective_constraint"
    )
    model.add_cons_vars(constraint)
    
    model.objective = Zero
    
def get_fva_statistics(model, fraction_of_optimum = 1.0):
    """adapted from autopacmen

    :param model:
    :return:
    """
    
    fva_dict = {}
    
    with model:
        # set lower bound of formr objective
        fix_objective_reactions(model, fraction_of_optimum)
        
        for reaction in model.reactions:

            orig_id = reaction.id.split("_GPRSPLIT_")[0]
            pattern = re.compile(orig_id + "_GPRSPLIT_\d+")

            # Minimum
            with model:                
                # knocks out all split reverse reactions
                ko_reverse_reactions(model, reaction, pattern)

                objective_expression = -1 * model.reactions.get_by_id(reaction.id).flux_expression       
                
                model.objective = model.problem.Objective(
                    objective_expression
                    )
                
                minimum = -1 *  model.slim_optimize()

            # Maximum
            with model:
                ko_reverse_reactions(model, reaction, pattern)
                
                objective_expression = 1 * model.reactions.get_by_id(reaction.id).flux_expression

                model.objective = model.problem.Objective(
                    objective_expression)
                
                maximum = model.slim_optimize()

            variability = abs(maximum - minimum)
            
            if variability > 1000:
                print(f"High variability detected in reaction {reaction.id}: {variability}")

            fva_dict[reaction.id] = {}
            fva_dict[reaction.id]["minimum"] = minimum
            fva_dict[reaction.id]["maximum"] = maximum
            fva_dict[reaction.id]["variability"] = variability

    return fva_dict

def ec_to_og_fva_mapping(fva_df, original_reactions, ec_model):
    mapped_results_df = pandas.DataFrame(index=original_reactions, columns=['range'])

    for reaction in ec_model.reactions:
        reaction_og = ''
        reaction_ec = ''
        
        # omit: ENZYME delivery, port_pool, armm reactions
        if reaction.id.startswith("ENZYME_DELIVERY_") or reaction.id == "ER_pool_TG_" or reaction.id.startswith("armr_"):
            continue
        
        if "_GPRSPLIT_" in reaction.id:
            split_reactions_in_direction = get_split_reactions_in_direction(reaction, ec_model)

            # for reversible reactions
            direction = ''
            if '_TG_' in reaction.id:
                direction = "_TG_" + reaction.id.split('_TG_')[1]
            
            original_reaction_id = reaction.id.split('_TG_')[0].split('_GPRSPLIT_')[0]
            
            if len(split_reactions_in_direction) == 1:
                # no parallel reactions: _GPRSPLIT_ is removed
                reaction_og = original_reaction_id + direction
                reaction_ec = reaction.id

            else:
                # get armm reaction in the same direction and enter min/max directly
                
                # no reaction direction
                arm_reaction_id = "armr_" + original_reaction_id
                
                # ... if reaction has a direction
                if not direction.split('_TG_')[-1] == '':
                    # arm reactions dont include the seperator _TG_, split reactions do
                    arm_reaction_id = f"{arm_reaction_id}_{direction.split('_TG_')[-1]}"
                
                reaction_og = original_reaction_id + direction
                reaction_ec = arm_reaction_id    
                
        else:
            reaction_og = reaction.id
            reaction_ec = reaction.id

        mapped_results_df.loc[reaction_og] = abs(fva_df.loc[reaction_ec, 'maximum'] - fva_df.loc[reaction_ec, 'minimum'])

    return mapped_results_df

def set_axis_style(ax, group_labels, x_label):
    ax.set_xticks(numpy.arange(1, len(group_labels) + 1), labels=group_labels,)
    ax.set_xlim(0.25, len(group_labels) + 0.75)
    ax.set_xlabel(x_label)

def fva_scatter_plot(fva_results: Dict[str, pandas.DataFrame]):
    
    rng = numpy.random.default_rng(seed=42)

    fig, ax = plt.subplots()

    for index, key in enumerate(fva_results.keys(), start = 1):
        print(key)
        vals = fva_results[key]
        x_positions = rng.normal(index, 0.05, size=len(vals))
        ax.scatter(x_positions, vals, s=30, alpha=0.3,
                label="Reaction ranges" if index == 1 else "")

        ax.set_ylabel("Flux Range in mmol/(gDW h)")
        ax.set_title("Flux Ranges per Model")
        set_axis_style(ax, list(fva_results.keys()), 'Models')
        # ax.set_ylim([0, 1])

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

def statistics_and_barplots(fva_results_mapped_dict):
    # get cumulated variability and average variability
    cumulated_variabilies = {}
    mean_variabilities = {}

    for model, ranges in fva_results_mapped_dict.items():
        cumulated_variabilies[model] = ranges.sum()[0]
        mean_variabilities[model] = ranges.mean()[0]
    
    print("Cumulated variabilities:")
    print(cumulated_variabilies)
    print("Mean variabilities:")
    print(mean_variabilities)
    
    # barplots
    models = list(cumulated_variabilies.keys())
    values = [value*100/list(cumulated_variabilies.values())[0] for value in list(cumulated_variabilies.values())]

    fig, ax = plt.subplots()

    bars = ax.bar(models, values)
    plt.xticks(rotation=45, ha='right')
    ax.set_ylim((0.0, 110.0))
    ax.bar_label(bars, padding=3, fmt="%.1f")
    ax.set_ylabel('relative cumulated variabilities')

    plt.show()
    
    models = list(mean_variabilities.keys())
    values = [value for value in list(mean_variabilities.values())]

    fig, ax = plt.subplots()

    bars = ax.bar(models, values)
    plt.xticks(rotation=45, ha='right')
    ax.set_ylim((0, 690))
    ax.bar_label(bars, padding=3, fmt="%.1f")
    ax.set_ylabel('mean variabilities in mmol/(gdw h)')

    plt.show()

def cumulative_flux_variability_graph(
        variability_list_dict: dict, savefig_path: str):
    fig, ax = plt.subplots(figsize=(8, 4))
    n_bins = 2000
    linewidth = 1.2

    index = 0

    for model in variability_list_dict.keys():
        major_variabilities = [x for x in variability_list_dict[model]['range']] # if x > 1e-3]
        variabilities_above_500 = [x for x in variability_list_dict[model]['range'] if x >= 500]
        variabilites_equal_1000 = [x for x in variability_list_dict[model]['range'] if x == 1000]

        print('------------------' + model + '-----------------')
        print('# vars above 500: ' + str(len(variabilities_above_500)))
        print('# vars above 1000: ' + str(len(variabilites_equal_1000)))

        if index == 0:
            n, bins, patches = ax.hist(major_variabilities, n_bins, density=True, histtype='step', linewidth=linewidth,
                                       cumulative=True, label=f'{model} (n={len(major_variabilities)})')
            ax.set_ylim((0.3, 1.02))
        else:
            # Overlay a reversed cumulative histogram.
            ax.hist(major_variabilities, bins=bins, density=True, histtype='step', cumulative=True, linewidth=linewidth,
                    label=f'{model} (n={len(major_variabilities)})')

        index += 1

    # Add titles
    ax.grid(True)
    ax.legend(loc='lower center')
    ax.set_xlabel('Variabilität in $\mathregular{mmol(g_{DW}h)^{-1}}$')
    ax.set_ylabel('Kumulative Flussvariabilitäten')
    fig1 = plt.gcf()
    plt.show()

    if savefig_path != "":
        fig1.savefig(savefig_path, format="svg")
        
def check_orphan_deliveries(model):
    
    orphan_deliveries = []
    
    for reaction in model.reactions:
        if reaction.id.startswith('ENZYME_DELIVERY_'):
            metabolites = list(reaction.metabolites.keys())
            if len(metabolites) > 1:
                continue
            metabolite = metabolites[0]
            
            if len(metabolite.reactions) == 1:
                orphan_deliveries.append(reaction.id)
    
    print(f'model contains {len(orphan_deliveries)} orphan deliveries.')
    print(orphan_deliveries)
    
    return orphan_deliveries