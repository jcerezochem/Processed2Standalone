# Processed2Standalone
Generate stand alone GROMACS topology from processed one (generated from gompp). The script generates a minimal standalone topology explicitely including all parameters instead of using [ Xtype ] databases, taking a preprocessed topology generated from `grompp`

## Quick guide

To get a minimal topology from a topolgy with #include calls (e.g. topol.top) one needs to follow these steps:

1. Run `grompp` to generate a prepocessed topology file. This assumes that you have a structure file (e.g. conf.gro) that is consistent with your topology (topol.top). You can normally use a blank mdp file (e.g. null.mdp):<br>
`gmx grompp -f null.mdp -p topol.top -c conf.gro -pp processed.top`

2. Then, apply the script to processed.top:<br>
`procesed2standalone.py -f processed.top > minimal.top`

In the process, a file with the FF in Gaussian MM format is generated (`ff_gau.prm`)
