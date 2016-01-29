""" Definition des parametres du probleme : """

nvar = 2
g = 9.8
eps = 0.000001

# Parametre temporels :
tFinal = 1.0

degree = 1
timeFunction = 'SSP2'
cfl = 0.2

# Boundary condition :
boundaryType = {
    1 : 'wall',
    2 : 'wall'
}




