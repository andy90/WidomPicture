gmx grompp -f equil.mdp -c start_0.gro -p water.top -o md.tpr -maxwarn 1
gmx mdrun  -deffnm md -table table.xvg
