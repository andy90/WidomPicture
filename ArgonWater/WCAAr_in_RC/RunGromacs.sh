gmx grompp -f SPCE.mdp -p water.top -c start.gro -o md.tpr
gmx mdrun -v -deffnm md -table table.xvg
