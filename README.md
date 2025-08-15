# Aqueous-Homopolymer-Circuit-Topology
Code to perform a multichain Circuit Topology (CT) analysis given trajectory data of entangled homopolymers. The original trajectory data was generated with Gromacs and will be made publicly available on final publication. 

The circuit_topology_analysis.c performs the CT analysis on an ensemble of equilibrated multipolymer systems, where monomer coordinates are recorded in files with the general structure "%dK%03d", with first argument the temperature and second argument the simulation number (index). Compile using the associated Makefile, which yields circuit_topology_analysis.exe.

The circuit_topology_analysis.exe program takes four arguments (separated by spaces): #polymers #monomers-per-polymer cutoff temperature. The cutoff field defines the threshold (in units of monomer spacing) for which the proximity of two monomers can be considered a contact (excluding neighbours). Running the script produces files containing the multichain motif fractions, the contact map for a given simulation index (here chosen as i=50), and a file containing the loops between contacts and their sizes. 

The bash script run_circuit_analysis.sh is included to perform the CT analysis over a range of temperature and cutoff values, depending on the available trajectory data. Sample trajectory data is included in the data.zip file.
