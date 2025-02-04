"""
Program: FERE (Fitted Elemental Reference Energies) calculator
Author: Shuji Kanayama, Soungmin Bae*, and Hannes Raebiger*
E-mail: soungminbae@gmail.com; hannes@ynu.ac.jp

[ Description ]
This program calculates the Formation Enthalpy of compounds using the FERE (Fitted Elemental Reference Energies) method from a single input CSV file.
The input file contains both compound data and element reference energy data in one file.

Input file format (e.g., input.csv, default):
--------------------------------------------------------------------------------
Compounds,Total energy (eV),"Formation enthalpy (experiment, eV)",Element,Total energy (eV)
Ag2O,-67.79822767,-0.11,Ag,-30.52927095
Ag2O2,-74.12923815,-0.06,Al,-7.73443733
Ag2S,-71.09878532,-0.11,As,-20.19399038
Ag2Se,-81.5628843,-0.15,Au,-67.95332567
AlAs,-29.13136633,-0.61,Ba,-39.53105312
...
--------------------------------------------------------------------------------

Note:
  - The first "Total energy (eV)" column represents the compound's total energy (in eV).
  - The second "Total energy (eV)" column (automatically renamed by pandas as "Total energy (eV).1") represents the reference energy for the element.
  - The column "Formation enthalpy (experiment, eV)" contains the experimental formation enthalpy (in eV).

Processing:
  1. For each compound, the program parses the chemical formula (from the "Compounds" column) to determine its elemental composition.
  2. It then builds a reference energy (E_ref) for each element using the "Element" column and the second "Total energy (eV)" column.
     If an element appears multiple times, the average is used.
  3. For each compound, the following are computed:
       - N_total: Total number of atoms (parsed from the chemical formula)
       - S_ref: The sum over elements (n_i * E_ref(i))
       - y = (E_compound – S_ref) – (N_total × ΔH_exp),
         where E_compound is the compound's total energy (from the first "Total energy (eV)" column)
         and ΔH_exp is the experimental formation enthalpy.
  4. A composition matrix M (rows: compounds, columns: elements) is built from the atomic counts.
  5. The linear system M · δ = y is solved (in the least-squares sense) to obtain the FERE correction values δ (in eV) for each element.
  6. Finally, for each compound, the program calculates:
       - Uncorrected (calculated) formation enthalpy: ΔH_calc = (E_compound – S_ref) / N_total
       - FERE corrected formation enthalpy: ΔH_FERE = (E_compound – Σ_i[n_i*(E_ref(i)+δ(i))]) / N_total

Output:
  Two CSV files are generated:
    1. output_formation_enthalpy.csv: For each compound, includes:
           Chemical Formula, Experimental Enthalpy (eV), Calculated Enthalpy (uncorrected, eV),
           FERE Corrected Enthalpy (eV)
    2. output_FERE_correction.csv: For each element, includes:
           Element, FERE Correction (eV)

Additionally, the program prints to the screen:
  - The FERE Correction value for each element.
  - The Mean Absolute Error (MAE) between the experimental enthalpy and both the uncorrected and FERE corrected formation enthalpies.

Usage:
  To run the program (default input file is "./input.csv"):
      python FERE_calculator.py --input ./input.csv
  Use the -h or --help option for detailed usage and file format information.
"""

import os
import argparse
import numpy as np
import pandas as pd
import re

# --------------------------
# Chemical Formula Parsing Functions
# --------------------------
def parse_formula(formula, pattern=r'([A-Z][a-z]*)(\d*)'):
    """
    Parse a chemical formula string and return a dictionary with element counts.
    Example: "Fe2O3" -> {"Fe": 2, "O": 3}
    """
    matches = re.findall(pattern, formula)
    counts = {}
    for element, count_str in matches:
        count = int(count_str) if count_str else 1
        counts[element] = counts.get(element, 0) + count
    return counts

def get_total_atoms(formula):
    """
    Return the total number of atoms in the chemical formula.
    """
    counts = parse_formula(formula)
    return sum(counts.values())

# --------------------------
# Input Data Reading Function
# --------------------------
def read_input(input_csv_path):
    """
    Read the input CSV file (e.g., input.csv) and return a pandas DataFrame.
    The file is expected to be comma-separated.
    """
    return pd.read_csv(input_csv_path, sep=",")

# --------------------------
# FERE Correction Calculation Functions
# --------------------------
def compute_FERE_corrections(df):
    """
    Using the DataFrame (df), compute for each compound:
      - N_total: Total number of atoms (parsed from the chemical formula)
      - S_ref: Sum over elements (n_i * E_ref(i)), where E_ref(i) is the reference energy for element i.
      - y = (E_compound – S_ref) – (N_total × ΔH_exp),
        where E_compound is the compound's total energy (first "Total energy (eV)" column)
        and ΔH_exp is the experimental formation enthalpy.
    
    The reference energy for each element is built using the "Element" column and the second
    "Total energy (eV)" column (renamed by pandas as "Total energy (eV).1"). If an element appears
    multiple times, the average is used.
    
    A composition matrix M (rows: compounds, columns: elements) is constructed from the atomic counts.
    The linear system M · δ = y is solved in the least-squares sense to obtain the FERE correction
    values δ (in eV) for each element.
    
    Returns:
      - corrections: A dictionary { element: δ }
      - comp_data: A list of dictionaries for each compound:
            { 'formula': ..., 'counts': {element: count, ...}, 'N_total': ..., 'S_ref': ...,
              'E_comp': ..., 'ΔH_exp': ..., 'y': ... }
    """
    compounds = df["Compounds"].tolist()
    # Compound energy from the first "Total energy (eV)" column
    E_compound_arr = df["Total energy (eV)"].values.astype(float)
    # Experimental formation enthalpy from the "Formation enthalpy (experiment, eV)" column
    exp_enthalpy_arr = df['Formation enthalpy (experiment, eV)'].values.astype(float)
    
    # The second "Total energy (eV)" column (element energy) is automatically renamed to "Total energy (eV).1"
    if "Total energy (eV).1" in df.columns:
        E_element_arr = df["Total energy (eV).1"].values.astype(float)
    else:
        raise ValueError("The input file must include the second 'Total energy (eV)' column (element energy).")
    
    # Get element names from the "Element" column (ignoring blanks or NaN)
    element_list_input = df["Element"].tolist()
    
    # Build reference energy: if an element appears multiple times, use the average.
    ref_energy = {}
    counts_ref = {}
    for el, e_val in zip(element_list_input, E_element_arr):
        if pd.isna(el):
            continue
        ref_energy[el] = ref_energy.get(el, 0) + e_val
        counts_ref[el] = counts_ref.get(el, 0) + 1
    for el in ref_energy:
        ref_energy[el] /= counts_ref[el]
    
    # Build compound data
    comp_data = []
    # Collect all elements present in the compounds (for columns of the composition matrix)
    elements_in_compounds = set()
    for formula in compounds:
        comp_counts = parse_formula(formula)
        elements_in_compounds.update(comp_counts.keys())
        N_total = sum(comp_counts.values())
        comp_data.append({
            'formula': formula,
            'counts': comp_counts,
            'N_total': N_total
        })
    elements_in_compounds = sorted(list(elements_in_compounds))
    
    # For each compound, compute S_ref = sum_{i in formula} [n_i * E_ref(i)]
    # If an element is missing from ref_energy, raise an error.
    y_list = []
    for idx, comp in enumerate(comp_data):
        counts = comp['counts']
        S_ref = 0
        for el, n in counts.items():
            if el not in ref_energy:
                raise ValueError(f"In compound {comp['formula']}, reference energy for element {el} is missing in the input file.")
            S_ref += n * ref_energy[el]
        comp['S_ref'] = S_ref
        comp['E_comp'] = E_compound_arr[idx]
        comp['ΔH_exp'] = exp_enthalpy_arr[idx]
        # y = (E_compound - S_ref) - (N_total * ΔH_exp)
        y = (comp['E_comp'] - S_ref) - (comp['N_total'] * comp['ΔH_exp'])
        comp['y'] = y
        y_list.append(y)
    
    # Build the composition matrix M (rows: compounds, columns: elements)
    M = np.zeros((len(comp_data), len(elements_in_compounds)))
    for i, comp in enumerate(comp_data):
        for j, el in enumerate(elements_in_compounds):
            M[i, j] = comp['counts'].get(el, 0)
    
    y_arr = np.array(y_list)
    
    # Solve the linear system M * δ = y in the least-squares sense
    δ, residuals, rank, s = np.linalg.lstsq(M, y_arr, rcond=None)
    
    corrections = {el: δ_val for el, δ_val in zip(elements_in_compounds, δ)}
    return corrections, comp_data

def compute_formation_enthalpy(comp, corrections, ref_energy):
    """
    For a given compound (dictionary comp), compute:
      - Uncorrected (calculated) formation enthalpy: ΔH_calc = (E_compound - S_ref) / N_total
      - FERE corrected formation enthalpy: ΔH_FERE = (E_compound - Σ_i[n_i*(E_ref(i)+δ(i))]) / N_total
    Return both values.
    """
    counts = comp['counts']
    S_ref_corr = 0
    for el, n in counts.items():
        S_ref_corr += n * (ref_energy[el] + corrections[el])
    ΔH_FERE = (comp['E_comp'] - S_ref_corr) / comp['N_total']
    ΔH_calc = (comp['E_comp'] - comp['S_ref']) / comp['N_total']
    return ΔH_calc, ΔH_FERE

# --------------------------
# Main Processing Function
# --------------------------
def main():
    parser = argparse.ArgumentParser(
        description="""FERE Formation Enthalpy Calculation Script (Independent of U-values)

Input file (e.g., input.csv) format (comma-separated):
--------------------------------------------------------------------------------
Compounds,Total energy (eV),"Formation enthalpy (experiment, eV)",Element,Total energy (eV)
Ag2O,-67.79822767,-0.11,Ag,-30.52927095
Ag2O2,-74.12923815,-0.06,Al,-7.73443733
Ag2S,-71.09878532,-0.11,As,-20.19399038
Ag2Se,-81.5628843,-0.15,Au,-67.95332567
AlAs,-29.13136633,-0.61,Ba,-39.53105312
...
--------------------------------------------------------------------------------
Note: The duplicate header "Total energy (eV)" will have its second occurrence automatically renamed to "Total energy (eV).1".

Output files:
  1. output_formation_enthalpy.csv: For each compound, includes:
         Chemical Formula, Experimental Enthalpy (eV), Calculated Enthalpy (uncorrected, eV), FERE Corrected Enthalpy (eV)
  2. output_FERE_correction.csv: For each element, includes:
         Element, FERE Correction (eV)
        """,
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "-i", "--input",
        default="./input.csv",
        help="Input CSV file path (default: ./input.csv)"
    )
    args = parser.parse_args()
    
    input_csv = args.input
    if not os.path.exists(input_csv):
        print(f"Error: The input file '{input_csv}' does not exist.")
        return
    
    # Read the input file (comma-separated)
    df = read_input(input_csv)
    
    # Compute FERE corrections and build compound data
    corrections, comp_data = compute_FERE_corrections(df)
    
    # Build the reference energy dictionary (average if an element appears multiple times)
    if "Total energy (eV).1" in df.columns:
        E_element_series = df["Total energy (eV).1"]
    else:
        raise ValueError("The input file must include the second 'Total energy (eV)' column (element energy).")
    ref_energy = {}
    counts_ref = {}
    for el, e_val in zip(df["Element"].tolist(), E_element_series.tolist()):
        if pd.isna(el):
            continue
        ref_energy[el] = ref_energy.get(el, 0) + e_val
        counts_ref[el] = counts_ref.get(el, 0) + 1
    for el in ref_energy:
        ref_energy[el] /= counts_ref[el]
    
    # Print each element's FERE Correction value on the screen
    print("\nFERE Correction values for each element:")
    for element, corr in sorted(corrections.items()):
        print(f"{element}: FERE Correction: {corr:.4f} eV")
    
    # Compute uncorrected and FERE corrected formation enthalpies for each compound
    calc_enthalpies = []
    fere_enthalpies = []
    formulas = []
    exp_enthalpies = []
    
    for comp in comp_data:
        ΔH_calc, ΔH_FERE = compute_formation_enthalpy(comp, corrections, ref_energy)
        calc_enthalpies.append(ΔH_calc)
        fere_enthalpies.append(ΔH_FERE)
        formulas.append(comp['formula'])
        exp_enthalpies.append(comp['ΔH_exp'])
    
    # Compute Mean Absolute Error (MAE) for both uncorrected and FERE corrected formation enthalpies vs. experiment
    mae_corrected = np.mean(np.abs(np.array(fere_enthalpies) - np.array(exp_enthalpies)))
    mae_uncorrected = np.mean(np.abs(np.array(calc_enthalpies) - np.array(exp_enthalpies)))
    
    print(f"\nFERE Corrected Formation Enthalpy MAE (eV): {mae_corrected:.4f}")
    print(f"FERE Uncorrected Formation Enthalpy MAE (eV): {mae_uncorrected:.4f}")
    
    # Create the compounds result DataFrame and save to output_formation_enthalpy.csv
    final_df = pd.DataFrame({
        "Chemical Formula": formulas,
        "Experimental Enthalpy (eV)": exp_enthalpies,
        "Calculated Enthalpy (uncorrected, eV)": calc_enthalpies,
        "FERE Corrected Enthalpy (eV)": fere_enthalpies
    })
    final_df.to_csv("./output_formation_enthalpy.csv", index=False)
    print("\nCompound results have been saved to 'output_formation_enthalpy.csv'.")
    
    # Create the element corrections DataFrame and save to output_FERE_correction.csv
    correction_df = pd.DataFrame({
        "Element": list(corrections.keys()),
        "FERE Correction (eV)": list(corrections.values())
    })
    correction_df.to_csv("./output_FERE_correction.csv", index=False)
    print("Element correction values have been saved to 'output_FERE_correction.csv'.")

if __name__ == "__main__":
    main()

