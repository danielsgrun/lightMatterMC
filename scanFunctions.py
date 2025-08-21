# Suite of functions for running scans on the simulation #

import numpy as np

def run_scan(
    scan_vars, scan_values, param_template, static_params,
    filename_template, scan_labels, runner_args):

    """
    General-purpose scan function for SimulationRunner.

    Parameters
    ----------
    scan_vars : list[str]
        Names of variables to scan (e.g., ['P'] or ['s0', 'deltas']).
    scan_values : list[list]
        Lists of values corresponding to each variable in scan_vars.
    param_template : dict
        Base parameter set where scan variable values will be inserted.
    static_params : dict
        Additional fixed parameters for the simulation (e.g., T, titf).
    filename_template : str
        Template string for output file naming.
    scan_labels : list[str]
        Human-readable labels for the variables being scanned.
    runner_args : dict
        Additional arguments passed directly to SimulationRunner (e.g., n_samples, nAtoms).
    """
    
    from copy import deepcopy
    from wrapperFunctions import simClass

    
    is_double_scan = len(scan_vars) == 2
    results = []

    # Ensure templates aren't mutated
    base_template = deepcopy(param_template)

    if is_double_scan:
        var1, var2 = scan_vars
        vals1, vals2 = scan_values

        for j, val2 in enumerate(vals2):
            row = []
            for i, val1 in enumerate(vals1):
                param_set = deepcopy(base_template)
                param_set[var1] = val1
                param_set[var2] = val2

                # Create runner instance
                runner = simClass(
                    **{**param_set, **static_params, **runner_args}
                )
                row.append(runner.survival_prob())
            results.append(row)

    else:  # Single-variable scan
        var = scan_vars[0]
        vals = scan_values[0]

        for val in vals:
            param_set = deepcopy(base_template)
            param_set[var] = val

            runner = simClass(
                **{**param_set, **static_params, **runner_args}
            )
            results.append(runner.survival_prob())

    # Build fileName with available values
    # Safely pull out only variables needed for formatting
    format_dict = {**param_template, **static_params}
    try:
        fileName = filename_template.format(**format_dict)
    except KeyError:
        fileName = "result_unknown_format"

    # Build paramsDict to record metadata
    paramsDict = deepcopy(param_template)
    for key, val in static_params.items():
        paramsDict[key] = val
    paramsDict['scanvar'] = scan_labels if is_double_scan else scan_labels[0]

    return results, fileName, paramsDict


def getSurvivalPerAtom(results):
    nAtoms = len(results[0][0])
    survivals = [[] for j in range(nAtoms)]
        
    for j in range(nAtoms):
        survivals[j].append([results[i][0][j] for i in range(len(results))])
        
    phScatt = [results[i][1] for i in range(len(results))]
        
    survivals = np.array(survivals)[:,0]
    
    return [survivals, phScatt]

