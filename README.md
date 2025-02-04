# FERE-calculator
The program calculates the Formation Enthalpy of compounds using the FERE (Fitted Elemental Reference Energies) correction from DFT-calculated energies and experimental formation enthalpies

Author: Shuji Kanayama, Soungmin Bae*, and Hannes Raebiger*
E-mail: soungminbae@gmail.com; hannes@ynu.ac.jp

# Reference
Correcting density functional theory for accurate predictions of compound enthalpies of formation: Fitted elemental-phase reference energies, Vladan Stevanović, Stephan Lany, Xiuwen Zhang, Alex Zunger, Phys. Rev. B 85, 115104 (2012).
doi: https://doi.org/10.1103/PhysRevB.85.115104

# FERE database in json file (database-FERE-json.tar.gz)
Contains information of 251 compounds and 50 elements including input files of VASP (Vienna Ab initio Simulation Package) package, calculated and experimental formation enthalpies, structure, and symmetry (pymatgen format).
- database-01-r2scan+u.json
- database-02-scan+u.json
- database-03-scan+u_Artrith.json
- database-04-pbe+u_Moore.json

# Usage
python FERE-calculator.py -i inputfile.csv
python -h FERE-calculator.py

## r2SCAN+U
r2SCAN+U calculations with the effecive U parameters (U= 0, 1, 2, 2.5 eV)
- input-01-r2scan-u=0.csv
- input-01-r2scan-u=1.csv
- input-01-r2scan-u=2.csv
- input-01-r2scan-u=2.5.csv

## SCAN+U
SCAN+U calculations with the effecive U parameters (U= 0, 1, 2, 2.5 eV)
- input-02-scan-u=0.csv
- input-02-scan-u=1.csv
- input-02-scan-u=2.csv
- input-02-scan-u=2.5.csv

## SCAN+U, U parameter of Artrith et al.,
Reference: Data-driven approach to parameterize SCAN +𝑈 for an accurate description of 3⁢𝑑 transition metal oxide thermochemistry, Nongnuch Artrith, José Antonio Garrido Torres, Alexander Urban, and Mark S. Hybertsen, Phys. Rev. Materials 6, 035003 (2022).
doi: https://doi.org/10.1103/PhysRevMaterials.6.035003
- input-03-scan-u_Artrith.csv
- u_parameter_set-Artrith.yaml

## PBE+U, U parameter of Artrith et al.,
Reference: High-throughput determination of Hubbard 𝑈 and Hund 𝐽 values for transition metal oxides via the linear response formalism, Guy C. Moore, Matthew K. Horton, Edward Linscott, Alexander M. Ganose, Martin Siron, David D. O'Regan, and Kristin A. Persson, Phys. Rev. Materials 8, 014409 (2024).
doi: https://doi.org/10.1103/PhysRevMaterials.8.014409
- input-04-pbe-u_Moore.csv
- u_parameter_set-Moore.yaml

