# aspp-project 2022 - Lucas Baldo Mesa Casa
Project for the ASPP course

## Note

During this past week I was able to implement several parts of my plan, but there is much of what I had originally planned left to be done. Since this code is part of my reasearch I will slowly implement the remaining ideas, related to the techniques I learned in this course.

## Plan

In this project I intend to use the techniques I'm learning in the course to improve a code repository for tightbinding calculations that I have developed and have been using in my research. I feel the code has some room for improvement in terms of speed optimization and I would like to implement some MPI features to it, so I can make use of the Rackham cluseter I have access to. I also need to properly document the project.

A better outline of the project still needs to be developed, but the main points I want the project to cover are:
1. Documenting:
  So far I have been working alone on this code, so many things were developed undocumented. I plan to fix this so that others, as well as myself, can understand my code more quickly and effortlessly.
2. Implementing MPI capabilities:
  I want to be able to grant at least some of the routines capability to run calculations in parallel, so that I can take advantage of clusters where I get access to several processing cores and nodes.
3. Create package:
  Transform the code collection into a package by creating a __init__.py file. There is a method called reload_modules() that is a custom method that is very useful for me when I'm implementing new features. I need to find a way to keep this method in the initializing file.

Bonus:
  In case there is enough time I would also like to complete the following tasks.

a. Speed optimization:
  I think the way I currently build matrices could be optimized. Ideally I want to spend most of the time of the code diagonalizing the matrices, and the least amount of time doing other things (like constructing the matrices). This is because the part of the code that diagonalizes these matrices comes from standard libraries that are probably almost completely optimized.

b. Improve workflow and versatility:
  The python library is currently accompanied by jupyter notebook files which I use to input configurations, and manipulate and plot data. With time these interfaces became very bloated and convoluted. If possible I would like to overhaul them to a simpler system and add some command-line interface that will allow me to run the code remotely in clusters, for example.

c. GPU acceleration:
  Although I currently don't have access to a GPU-based computing centre, it would be interesting to learn the packages related to GPU acceleration in python as well as its strengths and weaknesses.

## Completion list:

- [ ] Documented the code:
    - [x] bands.py
    - [x] constants.py
    - [ ] graphics.py
    - [x] hamiltonians.py
    - [ ] multiterminal.py
    - [ ] nanowires.py
    - [ ] sharpness.py
    - [ ] solvers.py
    - [ ] supercurrent.py
    - [x] temperature.py
    - [ ] utilities.py
    - [ ] wavefunctions.py
- [ ] MPI:
    - [x] Create MPI_solver.py
    - [x] Implement code to separate parameter tuples and ask each process to solve them in parallel to other processes
    - [x] First round of debugging and testing
    - [x] Added a tester file "test_mpi.sh"
    - [ ] More rounds of debugging and testing
    - [ ] Organizing the code
    - [ ] Documenting
- [x] Transform code into a package by creating __init__.py file with necessary imports and my custom reload function
