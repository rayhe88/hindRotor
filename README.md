![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/rayhe88/hindRotor?label=version)
![GitHub last commit](https://img.shields.io/github/last-commit/rayhe88/hindRotor?label=last%20modified)
![GitHub top language](https://img.shields.io/github/languages/top/rayhe88/hindRotor?color=green)
# hRotor

![rotorCompleto](https://user-images.githubusercontent.com/30420297/155459487-9683adcc-9a37-4995-8d6b-9393bb1dcec1.jpeg)

## About

**_hRotor_** is a program to compute the rotational partition function for a molecular system at a temperature T or a range of temperatures.

## Platforms

hRotor has been tested on the following platforms:

- Ubuntu 20.04 LTS (Focal Fossa)

### Additional Linux requirements

- g++ 9.3.0 or later
- cmake 3.8 or later
- Lapack 3.9.0 or later

## Citation

The **_hRotor.x_** program is based on the following article:

- **"The 1-D hindered rotor approximation"**, J. Pfaendtner, X. Yu and L. J. Broadbelt,
  [*Theor Chem Account*, 118, 881 (2007)](https://doi.org/10.1007/s00214-007-0376-5)

The solution of Schödinger equation is based on the article (FGH):

- **"The Fourier grid Hamiltonian method for bound state eigenvalues and eigenfunctions"**,
  C. C. Marston and G. G. Balint-Kurti,
  [*J. Chem. Phys.*, 91, 3571 (1989)](https://doi.org/10.1063/1.456888)

## Installation

Let's start by copying or cloning the project in your system. Move or create a **build** directory.
The file tree should look like the following:

```
-> hrotor/$ tree
.
├── bin
│   └── hRotor.x
├── build
├── example
│   ├── methanol-opt.wfn
│   ├── methanolScanT.txt
│   └── topsFile.txt
├── include
│   ├── Atom.h
│   ├── FourierH.h
│   ├── HPot.h
│   ├── Inertia.h
│   ├── Matrix.h
│   ├── Molecule.h
│   ├── runCommands.h
│   ├── Rvector.h
│   ├── screen.h
│   ├── thermo.h
│   └── version.h.in
├── README.md
└── src
    ├── Atom.cpp
    ├── CMakeLists.txt
    ├── FourierH.cpp
    ├── HPot.cpp
    ├── Inertia.cpp
    ├── main.cpp
    ├── Makefile
    ├── Matrix.cpp
    ├── Molecule.cpp
    ├── runCommands.cpp
    ├── Rvector.cpp
    └── thermo.cpp

```

Go to the build directory and run the following commands:

```
 -> hrotor/build/$ cmake ../src
```

And then run

```
-> hrotor/build/$ cmake --build .
```

The executable **_hRotor.x_** will be created and moved to the **hrotor/bin** directory, you can run by typing.

```
$./hRotor.x -g GEOMETRY -p POTENTIAL -t ROTATING_TOPS -o OutputName
```

There are mandatory flags such as:

| Flag | Type   | Description                              |
| ---- | ------ | ---------------------------------------- |
| `-g` | string | Assign a file with the geometry          |
| `-p` | string | Indicate a file with the potential       |
| `-t` | string | Contain the information of rotation tops |

There are also optional flags such as:

| flag        | Type            | Description                                              |
| ----------- | --------------- | -------------------------------------------------------- |
| `-o`        | string          | Asign a name for the outpu                               |
| `-T`        | float           | Indicate a Temperature T different to 298.15 K           |
| `-r`        | float float int | Indicate a range of temperature: Ti, Tf and nT           |
| `-s`        | int             | Indicate the symmetry number                             |
| `-I`        | int             | Indicate the n in I(2,n) for the Inertia moment          |
| `-m`        | float           | Assign a numerical value for the Inertia moment          |
| `-n`        | int             | Indicate the size in FGH method                          |
| `-l`        |                 | Display LAPACK version                                   |
| `-v` o `-V` |                 | Display the code version, compilation date, and git info |
| `-h` o `-H` |                 | Display help                                             |

### Outputs

The output from the terminal is as follows:

```
$./hRotor.x -g ../example/methanol-opt.wfn -p ../example/methanolScanT.txt -o output_ -t ../example/topsFile.txt  -T 1500. -s 3
[ INFO  ]  Geometry  name : ../example/methanol-opt.wfn
[ INFO  ]  Potential name : ../example/methanolScanT.txt
[ INFO  ]  Tops      name : ../example/topsFile.txt
[ INFO  ]  Output    name : output_
[ INFO  ]  Single    Temp :  true
[ INFO  ]  Range     Temp :  false
[WARNING]  No Inertia moment assigned
           The Inertia moment will be taken by default I (2,1)
[WARNING]  No Hamiltonian size assigned
           The Hamiltonian size will be taken by default H size = 501
[ INFO  ]   Hamilton size : 501
[ INFO  ]    Kinetic size : 250
[ INFO  ]     delta   x   : 0.0125413
[ INFO  ]     delta   k   : 1
[ INFO  ]     Data in pot : 35
=================================================================
 Print the values of thermo chemical properties
      Temperature  :    1500.00000000  K
          kB T     :      12.47169393  kJ / mol
 Highest potential :       5.57918738  kJ / mol
      I_red(2,1)   :       0.67176249  Da A^2
=================================================================
                  Free Rotor    Hindered Rotor
=================================================================
    Qr :          3.80813558      3.08477950
    S  :          0.01527483      0.01518053   kJ / (mol K )
    E  :          6.23584696      8.72168005   kJ / mol
    Cv :          0.00415723      0.00433721   kJ / (mol K )
=================================================================

```

The internal energy (E), entropy (S), head capacity at constant volume (Cv), and
the partition function are reported at the end ofthe run.

On the other hand, the following files are also generated:

| File              | Data it contains                                          |
| ----------------- | --------------------------------------------------------- |
| output_Energy.dat | eigenvalues of FGH                                        |
| output_Pot.dat    | Potential using in the FGH evaluate in te grid            |
| outPut_Rho.dat    | Desity for each eigenstate evaluate in the grid           |
| outPut_Wf.dat     | Wavefunction for each eigenstate evaluated in the grid    |
| outPut_chem.dat   | The information when range Temperature has been activated |

## Code contributors

- Raymundo Hernandez-Esparza (rayhe88@gmail.com)
- Julio Manuel Hernandez-Perez
- Sebastian García Pineda
