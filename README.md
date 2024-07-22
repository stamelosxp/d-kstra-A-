# Dijkstra and A* Algorithm Implementation

## Overview

This repository contains an implementation of the Dijkstra-SP algorithm and the A* algorithm for finding the shortest path between two nodes in a graph. The project was developed as a part of a programming exercise to compare the performance of these algorithms using the BOOST library.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Algorithms Implemented](#algorithms-implemented)
- [Heuristic Function for A*][#heuristic-function-for-a*]
- [Test Cases][#test-cases]
- [Requirements][#requirements]
### Introduction 

Dijkstra's Algorithm is a classic algorithm used in graph theory to determine the shortest path between nodes in a weighted graph. This repository offers a Python implementation that can be easily integrated into other projects or used for educational purposes.
### Installation

To use this project, you need to have the LEDA library installed on your system. Follow the instructions below to set up the project:

1. Clone the repository:
    ```sh
    git clone https://github.com/stamelosxp/dijkstra-implementation.git
    ```

2. Navigate to the `bin` directory:
    ```sh
    cd dijkstra-implementation/bin
    ```

3. Compile the project and build the main excecutable:
    ```sh
    make main
    ```

## Usage

Run the compiled `main` executable:
```sh
./main
```


### Algorithms Implemented

1. **Dijkstra-SP Algorithm**: A variant of Dijkstra's algorithm that stops as soon as the target node `t` is reached.
2. **A* Algorithm**: An extension of Dijkstra's algorithm that uses a heuristic function `ht` to guide the search.

### Heuristic Function for A*

1. **Euclidean Distance**: For grid graphs, using the Euclidean distance between nodes.
2. **Landmark-based Lower Bound**: For random graphs, using distances to and from strategically chosen landmark nodes.

### Test Cases

1. **Grid Graphs**
2. **Random Graphs**

## Requirements

- Boost Library