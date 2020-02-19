# Carrega as variaveis de ambiente do AMBER aqui.
source ~/software/amber18/amber.sh

# Calculo quantico semi-empirico (AM1-BCC) para as cargas do ligante
antechamber -i Etileno_glicol.mol2 -fi mol2 -o ligand.mol2 -fo mol2 -c bcc -s 2 -eq 2 -pf y -rn MOL -at gaff2

# Designa os parametros do campo de forças GAFF2 para o ligante
parmchk2 -i ligand.mol2 -o ligand.frcmod -pf 2 -w yes -a yes -f mol2

