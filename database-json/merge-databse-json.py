import os
import json

# Get base directory (current working directory)
base_dir = os.getcwd()

# Define U values and corresponding directories
U_values = ["U=0", "U=1", "U=2", "U=2.5"]

# Final merged data structure
merged_data = {
    "database": "R2SCAN+U",
    "entries": {}
}

# Loop over each U directory
for U in U_values:
    U_dir = os.path.join(base_dir, U)
    json_file = os.path.join(U_dir, "vasp_all_data.json")

    # Check if vasp_all_data.json exists
    if os.path.exists(json_file):
        with open(json_file, "r") as f:
            data = json.load(f)

        # Label each entry with the corresponding Hubbard U value
        for key, value in data["compounds"].items():
            if key not in merged_data["entries"]:
                merged_data["entries"][key] = {}  # Initialize sub-dictionary

            merged_data["entries"][key][U] = value  # Store data under U label

print("✅ Successfully merged data from all U values.")

# Save the merged data into a single JSON file
output_json = "merged_vasp_all_data_R2SCANU.json"
with open(output_json, "w") as f:
    json.dump(merged_data, f, indent=4)

print(f"📁 Final merged data saved as {output_json}")
