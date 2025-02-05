#!/usr/bin/env python
"""
=======================================================================
Periodic Table FERE Plotter
=======================================================================
This script generates a periodic table plot where each element is
colored based on its "FERE Correction (eV)" value from an input CSV file.
Each element's box displays the element symbol and its corresponding
FERE value. A colorbar is also included to indicate the mapping of
FERE values to colors.

Usage:
------
Run the script from the command line with the following options:

    -i, --input   : (Optional) Input CSV file containing FERE data.
                    Default is "output_FERE_correction.csv".
    -o, --output  : (Optional) Output PDF file for the generated plot.
                    Default is "periodic-FERE.pdf".
    -s, --show    : (Optional) If specified, displays the plot window.
                    By default, the plot window is not shown.

Example:
    python plotter-periodic-FERE.py -i mydata.csv -o myplot.pdf -s

Input CSV File Format:
----------------------
The input CSV file must contain at least the following columns:
    - "Element": The chemical symbol of the element (e.g., H, He, Li, etc.).
    - "FERE Correction (eV)": A numeric value representing the FERE
      correction in electron volts.

The script uses a predefined periodic table layout to map each element
to its respective position in the plot.

Output:
-------
The script generates a PDF file (default: periodic-FERE.pdf) containing
the periodic table plot with colored boxes and a colorbar.
Additionally, the script prints the input file used, the output file name,
and a confirmation message once the PDF is saved.

Error Handling:
----------------
If the specified input CSV file does not exist in the current folder,
the script will print an error message along with the usage instructions
and exit.

=======================================================================
"""

import argparse
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase
import matplotlib as mpl

def main(args):
    # Print starting messages with the input and output file names
    print(f"Using input file: {args.input}")
    print(f"Output will be saved as: {args.output}")

    # Set the default font for the entire plot
    mpl.rcParams['font.family'] = 'Arial'

    # Adjustable parameters for font sizes, font weights, and figure size
    element_fontsize = 6       # Font size for the element symbol
    fere_fontsize = 6          # Font size for the FERE value

    element_fontweight = 'bold'    # Font weight for the element symbol
    fere_fontweight = 'normal'     # Font weight for the FERE value

    linewidth = 1              # Line width for the element box borders
    fig_width_cm = 14          # Figure width in centimeters
    fig_height_cm = 7.0        # Figure height in centimeters

    # Adjustable parameter for the XC functional type.
    # This variable will be displayed in the plot title.
    xc_functional = "SCAN+U (U= 1eV)"

    # Adjustable colormap parameters:
    # Available colormaps include (but are not limited to):
    # 'coolwarm' (default), 'viridis', 'plasma', 'inferno', 'magma', 'cividis' (Perceptually Uniform),
    # 'bwr', 'seismic' (Diverging),
    # 'Greys', 'Blues', 'Reds', etc. (Sequential)
    colormap_name = "coolwarm"  # Choose the colormap (default: 'coolwarm')
    cmap = getattr(plt.cm, colormap_name)   # Obtain the colormap object from plt.cm

    # Adjustable range for colormap normalization
    colormap_range_min = -1.0
    colormap_range_max = 1.0

    # Load the specified input CSV file
    data = pd.read_csv(args.input)

    # Define the periodic table layout (row, column positions for each element)
    layout = {
        'H': (0, 0), 'He': (0, 17),
        'Li': (1, 0), 'Be': (1, 1), 'B': (1, 12), 'C': (1, 13), 'N': (1, 14), 'O': (1, 15), 'F': (1, 16), 'Ne': (1, 17),
        'Na': (2, 0), 'Mg': (2, 1), 'Al': (2, 12), 'Si': (2, 13), 'P': (2, 14), 'S': (2, 15), 'Cl': (2, 16), 'Ar': (2, 17),
        'K': (3, 0), 'Ca': (3, 1), 'Sc': (3, 2), 'Ti': (3, 3), 'V': (3, 4), 'Cr': (3, 5), 'Mn': (3, 6), 'Fe': (3, 7),
        'Co': (3, 8), 'Ni': (3, 9), 'Cu': (3, 10), 'Zn': (3, 11), 'Ga': (3, 12), 'Ge': (3, 13), 'As': (3, 14), 'Se': (3, 15),
        'Br': (3, 16), 'Kr': (3, 17),
        'Rb': (4, 0), 'Sr': (4, 1), 'Y': (4, 2), 'Zr': (4, 3), 'Nb': (4, 4), 'Mo': (4, 5), 'Tc': (4, 6), 'Ru': (4, 7),
        'Rh': (4, 8), 'Pd': (4, 9), 'Ag': (4, 10), 'Cd': (4, 11), 'In': (4, 12), 'Sn': (4, 13), 'Sb': (4, 14), 'Te': (4, 15),
        'I': (4, 16), 'Xe': (4, 17),
        'Cs': (5, 0), 'Ba': (5, 1), 'La': (5, 2), 'Hf': (5, 3), 'Ta': (5, 4), 'W': (5, 5), 'Re': (5, 6), 'Os': (5, 7),
        'Ir': (5, 8), 'Pt': (5, 9), 'Au': (5, 10), 'Hg': (5, 11), 'Tl': (5, 12), 'Pb': (5, 13), 'Bi': (5, 14), 'Po': (5, 15),
        'At': (5, 16), 'Rn': (5, 17),
        'Fr': (6, 0), 'Ra': (6, 1), 'Ac': (6, 2)
    }

    # Map the layout information (row, column) into the data for plotting
    data['Row'] = data['Element'].map(lambda x: layout[x][0] if x in layout else np.nan)
    data['Column'] = data['Element'].map(lambda x: layout[x][1] if x in layout else np.nan)

    # Define the colormap normalization based on the adjustable range
    norm = mcolors.Normalize(vmin=colormap_range_min, vmax=colormap_range_max)
    colormap = cmap  # Use the selected colormap

    # Map colors to elements based on their FERE Correction (eV) value
    color_mapping = {row['Element']: colormap(norm(row['FERE Correction (eV)']))
                     for _, row in data.iterrows()}

    # Create the periodic table plot
    fig, ax = plt.subplots(figsize=(fig_width_cm / 2.54, fig_height_cm / 2.54))  # Convert cm to inches

    # Offsets for placing the element symbol and FERE value inside each box
    element_text_offset_y = 0.65  # Vertical offset for the element symbol (top of the box)
    fere_text_offset_y = 0.35     # Vertical offset for the FERE value (bottom of the box)

    # Add boxes and text for each element in the layout
    for element, (r, c) in layout.items():
        rect_color = "#ffffff"  # Default color is white if no data is available
        if element in color_mapping:
            rect_color = color_mapping[element]
            value = data.loc[data['Element'] == element, 'FERE Correction (eV)'].values[0]
            # Place the element symbol (upper part of the box)
            ax.text(c + 0.5, -r + element_text_offset_y, element,
                    ha='center', va='center', fontsize=element_fontsize, fontweight=element_fontweight)
            # Place the FERE value (lower part of the box)
            ax.text(c + 0.5, -r + fere_text_offset_y, f"{value:.2f}",
                    ha='center', va='center', fontsize=fere_fontsize, fontweight=fere_fontweight)
        else:
            # If data is not available for the element, just display the symbol in the center
            ax.text(c + 0.5, -r + 0.5, element,
                    ha='center', va='center', fontsize=element_fontsize, fontweight=element_fontweight)
        # Draw the rectangle for the element box
        rect = plt.Rectangle((c, -r), 1, 1, color=rect_color, ec='black', lw=linewidth)
        ax.add_patch(rect)

    # Set axis limits and turn off the axis display
    ax.set_xlim(-0.5, 18)
    ax.set_ylim(-7, 1)
    ax.axis('off')

    # Set the plot title with the XC functional type
    plt.title(f"FERE correction of elements, XC functional: {xc_functional}",
              fontsize=12, weight='bold')

    # Add a colorbar corresponding to the colormap
    cbar_ax = fig.add_axes([0.93, 0.2, 0.03, 0.6])  # [left, bottom, width, height]
    ColorbarBase(cbar_ax, cmap=colormap, norm=norm, orientation='vertical', label='FERE Correction (eV)')

    # Save the figure to the specified output PDF file
    plt.savefig(args.output, format="pdf", bbox_inches='tight')
    print(f"Output saved to: {args.output}")

    # If the --show / -s option is specified, display the plot window
    if args.show:
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to accommodate the colorbar
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a periodic table plot with FERE values."
    )
    parser.add_argument("-i", "--input", type=str, default="output_FERE_correction.csv",
                        help="Input CSV file containing FERE data (default: output_FERE_correction.csv)")
    parser.add_argument("-o", "--output", type=str, default="periodic-FERE.pdf",
                        help="Output PDF file for the generated plot (default: periodic-FERE.pdf)")
    parser.add_argument("-s", "--show", action="store_true", default=False,
                        help="Flag to display the plot window (default: do not show)")
    args = parser.parse_args()

    # If the input file does not exist in the current folder, print an error message, display help, and exit
    if not os.path.exists(args.input):
        print(f"Error: The input file '{args.input}' does not exist in the current folder.")
        parser.print_help()
        sys.exit(1)

    main(args)
