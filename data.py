""" Definition des parametres du probleme : """

nvar = 2
g = 9.8

# Parametre temporels :
tFinal = 0.000000000001
timeFunction = 'SSP1'
cfl = 0.5

# Degre de l'EF :
degre = 0

eps = 0.000001

# Parametres du maillage :
xDebut = 0.0
xFin = 4.
nCells = 400
deltax = (xFin - xDebut) / nCells