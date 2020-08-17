gmx grompp -f pull.mdp -c start.gro -p water.top -o md.tpr
gmx mdrun -v -deffnm md -px pullx.xvg -pf pullf.xvg -table table.xvg
