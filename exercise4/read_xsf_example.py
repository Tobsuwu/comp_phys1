from numpy import *
from scipy.integrate import simps

def read_example_xsf_density(filename):
    lattice=[]
    density=[]
    grid=[]
    shift=[]
    i=0
    start_reading = False
    with open(filename, 'r') as f:
        for line in f:
            if "END_DATAGRID_3D" in line:
                start_reading = False
            if start_reading and i==1:
                grid=array(line.split(),dtype=int)
            if start_reading and i==2:
                shift.append(array(line.split(),dtype=float))
            if start_reading and i==3:
                lattice.append(array(line.split(),dtype=float))
            if start_reading and i==4:
                lattice.append(array(line.split(),dtype=float))
            if start_reading and i==5:
                lattice.append(array(line.split(),dtype=float))
            if start_reading and i>5:            
                density.extend(array(line.split(),dtype=float))
            if start_reading and i>0:
                i=i+1
            if "DATAGRID_3D_UNKNOWN" in line:
                start_reading = True
                i=1
    
    rho=zeros((grid[0],grid[1],grid[2]))
    ii=0
    for k in range(grid[2]):
        for j in range(grid[1]):        
            for i in range(grid[0]):
                rho[i,j,k]=density[ii]
                ii+=1

    # convert density to 1/Angstrom**3 from 1/Bohr**3
    a0=0.52917721067
    a03=a0*a0*a0
    rho/=a03
    return rho, array(lattice), grid, shift

def number_of_electrons(filename1, filename2):
    # If rho is the numerical density, we can calculate the number of electrons inside the unit cell by integrating
    # it over the lattice vectors

    rho1, lattice1, grid1, shift1 = read_example_xsf_density(filename1)
    rho2, lattice2, grid2, shift2 = read_example_xsf_density(filename2)

    x1 = linspace(0, lattice1[0][0], grid1[0])
    y1 = linspace(0, lattice1[1][1], grid1[1])
    z1 = linspace(0, lattice1[2][2], grid1[2])

    x2 = linspace(0, lattice2[0][0], grid2[0])
    y2 = linspace(0, lattice2[1][1], grid2[1])
    z2 = linspace(lattice2[2][0], lattice2[2][2], grid2[2])

    num_elec1 = simps(simps(simps(rho1, z1), y1), x1)
    num_elec2 = simps(simps(simps(rho2, z2), y2), x2)

    print('Number of electrons: ')
    print('dft_chargedensity1.xsf: ', round(num_elec1))
    print('dft_chargedensity2.xsf: ', round(num_elec2))

def reciprocal_vectors(filename1, filename2):
    # Reciprocal lattice vectors can be calculated from the real lattice vectors by equation B^T*A = 2pi*I (lec4 slide 4)

    rho1, lattice1, grid1, shift1 = read_example_xsf_density(filename1)
    rho2, lattice2, grid2, shift2 = read_example_xsf_density(filename2)

    B1 = transpose(2*pi*linalg.inv(transpose(lattice1)))
    B2 = transpose(2*pi*linalg.inv(lattice2.T))

    print('Reciprocal lattices: ')
    print('dft_chargedensity1.xsf:')
    print(B1)
    print()
    print('dft_chargedensity2.xsf: ')
    print(B2)

def main():

    filename1 = 'dft_chargedensity1.xsf'
    filename2 = 'dft_chargedensity2.xsf'
    number_of_electrons(filename1, filename2)
    print()
    reciprocal_vectors(filename1, filename2)


if __name__=="__main__":
    main()



