""" Definition des parametres du probleme : """

nvar = 3

# Parametre temporels :
tFinal = 0.002914094
timeFunction = 'SSP2'
cfl = 0.33

# Degre de l'EF :
degre = 1

# Parametres du maillage :
xDebut = 0.0
xFin = 1.
nCells = 100
deltax = (xFin - xDebut) / nCells