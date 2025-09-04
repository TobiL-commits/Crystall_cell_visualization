# Crystall cell visualization

This project is inspired by code I have written for an exercise sheet in university. I want to make the code easy to integrade into other projects, like for exmple adding some form of physics to it. The first stage of the project will be to give the lattice atoms and cells precise positions and use plotly to achieve different ways of simulating the resulting cells and atoms.

The long term prospects are, that I will write some form of physics engine, that is able to integrate with various projects, that simulate differend kinds of behaviour. Maybe shifting the calculations into C++ at some point, to recieve better performance. All these prospects are still written in the stars.

## Code structure idea

Create two different moduldes for making it easier to integrate it into other programs or applications. Therefore there should be

1. Positional module: A module saving positions and making calculations or recieving an APIs information
2. Visualization module: A module taking the input from the previous module ans drawing graphs.

To accomplish this, the two modules are to be written somewhat simultaneously, where the visualization is coded in a testing oriented manner, that first does not concern itself with creating classes or lager structures. It is supposed to preordain the following code structure, by creating simple functionality, that can be later integrated into the actual module.

## A similar more comprehensive project

-   https://github.com/manny405/mcse

## APIs to consider including

-   https://ase-lib.org/ase/optimize.html
-   https://pymatgen.org/
-   https://github.com/spglib/spglib
