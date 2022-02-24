# hRotor

## About

**_hRotor_** is a program to compute the rotational partition function for a molecular system at a temperature T or a range of temperatures.

## Platforms

hRotor has been testen on the following platforms:

- Ubuntu 20.04 LTS (Focal Fossa)

### Aditional Linux erquirements

- g++ 9.3.0 or later
- Lapack 3.9.0 or later

## Citation

The **_hRotor_** program is based on the following article:

- **"The 1-D hindered rotor approximation"**, J. Pfaendtner, X. Yu and L. J. Broadbelt,
  [Theor Chem Account 118, 881 (2007)](https://doi.org/10.1007/s00214-007-0376-5)

The solution of Schödinger equation is based on the article (FGH):

- **"The Fourier grid Hamiltonian method for bound state eigenvalues and eigenfunctions"**,
  C. C. Marston and G. G. Balint-Kurti,
  [J. Chem. Phys. 91, 3571 (1989)](https://doi.org/10.1063/1.456888)

## Installation

Let's start by copying or cloning the project in your system. Move to the directory
**src** containing a file **Makefile**.
Now run the file by typing **make** inside the **src** directory.

```
$ make
```

After the compilation you have a executable file called **_hRotor.x_**, you can run by typing.

```
$./hrotor.x -g GEOMETRY -p POTENTIAL -t ROTATING_TOPS -o OutputName
```

There are mandatory flags such as:

| Flag | Type   | Description                              |
| ---- | ------ | ---------------------------------------- |
| `-g` | string | Assign a file with the geometry          |
| `-p` | string | Indicate a file with the potential       |
| `-t` | string | Contain the information of rotation tops |

There are also optional flags such as:

| flag        | Type            | Description                                     |
| ----------- | --------------- | ----------------------------------------------- |
| `-o`        | string          | Asign a name for the outpu                      |
| `-T`        | float           | Indicate a Temperature T different to 298.15 K  |
| `-r`        | float float int | Indicate a range of temperature: Ti, Tf, stepsT |
| `-s`        | int             | Indicate the symmetry number                    |
| `-I`        | int             | Indicate the n in I(2,n) for the Inertia moment |
| `-n`        | int             | Indicate the size in FGH method                 |
| `-l`        |                 | Display LAPACK version                          |
| `-v` o `-V` |                 | Display version                                 |
| `-h` o `-H` |                 | Display help                                    |

## Code contributors

- Raymundo Hernandez-Esparza (rayhe88@gmail.com)
- Julio Manuel Hernandez-Perez ()

![hrotor](https://drive.google.com/file/d/1MrTIBqFVbqaSYdoWa4wc_SawPjhznI2A/view?usp=sharing)
