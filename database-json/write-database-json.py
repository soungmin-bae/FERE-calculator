import os
import json
import re
import pandas as pd
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from monty.json import MontyEncoder

# Create a log file to track errors
error_log_file = "error_log.txt"

def log_error(message):
    """Logs error messages to a file and prints them."""
    with open(error_log_file, "a") as log_file:
        log_file.write(message + "\n")
    print(message)

# Load FERE Correction data
fere_file = "output_FERE_correction.csv"
fere_data = pd.read_csv(fere_file, index_col=0).to_dict()["FERE Correction (eV)"]

# Load Formation Enthalpy data
enthalpy_file = "output_formation_enthalpy.csv"
enthalpy_data = pd.read_csv(enthalpy_file, index_col=0).to_dict()

# Get base directory and subdirectories
base_dir = os.getcwd()
element_dir = os.path.join(base_dir, "01-ELEMENT")
compound_dir = os.path.join(base_dir, "02-COMPOUNDS")

# Dictionary to store all data
all_data = {
    "elements": {},
    "compounds": {}
}

# Function to read POSCAR safely
def read_poscar_safe(poscar_file, compound_name):
    """Attempts to read POSCAR and logs an error if empty."""
    try:
        structure = Structure.from_file(poscar_file)
        return structure
    except ValueError as e:
        if "Empty POSCAR" in str(e):
            error_msg = f"❌ ERROR: Empty POSCAR file in directory: {compound_name} ({poscar_file})"
            log_error(error_msg)
            return None
        else:
            raise e  # Re-raise other ValueErrors

# Function to read outcar_summary.dat
def read_outcar_summary(filename):
    outcar_data = {}
    magmom_list = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("NION"):
                outcar_data["NION"] = int(line.split("=")[1].strip())
            elif line.startswith("Total_energy"):
                outcar_data["Total_energy"] = float(line.split("=")[1].strip())
            elif line.startswith("ISPIN"):
                outcar_data["ISPIN"] = int(line.split("=")[1].strip())
            elif line.startswith("MAGMOM"):
                magmom_list.extend([float(x) for x in line.split("=")[1].strip().split()])
            elif line.startswith("mangetic phase"):
                outcar_data["magnetic_phase"] = line.split(":")[1].strip()

    if magmom_list:
        outcar_data["MAGMOM"] = magmom_list

    return outcar_data

# Function to read INCAR file
def read_incar(filename):
    incar_data = {}
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            key_value = re.split(r"\s*=\s*", line, maxsplit=1)
            if len(key_value) == 2:
                key, value = key_value
                if re.match(r"^\d+$", value):  # Integer
                    value = int(value)
                elif re.match(r"^\d+\.\d+$", value):  # Float
                    value = float(value)
                elif value.lower() in ["true", "false"]:  # Boolean
                    value = value.lower() == "true"
                elif " " in value:  # List
                    value = [float(v) if "." in v else int(v) for v in value.split()]
                incar_data[key] = value
    return incar_data

# Extract data for **02-COMPOUNDS** (Compounds)
if os.path.exists(compound_dir):
    compound_folders = [d for d in os.listdir(compound_dir) if os.path.isdir(os.path.join(compound_dir, d))]

    for compound in compound_folders:
        compound_path = os.path.join(compound_dir, compound)

        # Ensure required files exist
        required_files = ["POSCAR", "INCAR", "outcar_summary.dat"]
        if not all(os.path.exists(os.path.join(compound_path, f)) for f in required_files):
            continue  # Skip if files are missing

        # Read POSCAR safely
        poscar_file = os.path.join(compound_path, "POSCAR")
        structure = read_poscar_safe(poscar_file, compound)
        if structure is None:
            continue  # Skip this compound if POSCAR is empty

        # Read INCAR and OUTCAR
        incar_file = os.path.join(compound_path, "INCAR")
        outcar_file = os.path.join(compound_path, "outcar_summary.dat")
        incar_info = read_incar(incar_file)
        outcar_info = read_outcar_summary(outcar_file)

        # Extract symmetry information
        sga = SpacegroupAnalyzer(structure)
        symmetry_data = {
            "space_group_symbol": sga.get_space_group_symbol(),
            "space_group_number": sga.get_space_group_number(),
            "crystal_system": sga.get_crystal_system(),
            "hall_symbol": sga.get_hall(),
            "point_group": sga.get_point_group_symbol(),
            "symmetry_operations": [op.as_dict() for op in sga.get_symmetry_operations()]
        }

        # Extract formation enthalpy values
        enthalpy_info = {}
        if compound in enthalpy_data["Experimental Enthalpy (eV)"]:
            enthalpy_info["DFT_formation_enthalpy"] = enthalpy_data["Calculated Enthalpy (uncorrected, eV)"][compound]
            enthalpy_info["Experimental_formation_enthalpy"] = enthalpy_data["Experimental Enthalpy (eV)"][compound]
            enthalpy_info["FERE_corrected_enthalpy"] = enthalpy_data["FERE Corrected Enthalpy (eV)"][compound]

        # Extract FERE correction values for compound elements
        fere_corrections = {}
        elements_in_compound = structure.composition.elements
        for element in elements_in_compound:
            symbol = element.symbol
            if symbol in fere_data:
                fere_corrections[symbol] = fere_data[symbol]

        # Store compound data
        all_data["compounds"][compound] = {
            "directory": os.path.relpath(compound_path, base_dir),
            "structure": structure.as_dict(),
            "symmetry": symmetry_data,
            "incar": incar_info,
            "outcar_summary": outcar_info,
            "formation_enthalpy": enthalpy_info,
            "FERE_corrections": fere_corrections
        }

# Save everything to a JSON file
json_filename = "vasp_all_data.json"
with open(json_filename, "w") as f:
    json.dump(all_data, f, cls=MontyEncoder, indent=4)

print(f"Data successfully saved to {json_filename}")
print(f"⚠️ Any errors encountered have been logged to {error_log_file}")
