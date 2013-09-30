ar
==

Build:

```
$ make
```

Example usage:

```
./seq 200 200 50 50 5000 500
```

This program will use sequential method to calculate voltage inside 200 by 200 grid.
Wire inside is 50 by 50 with constant voltage 5.0 V.
Iteration limit is 500.

```
mpirun -np 9 ./par 200 200 50 50 5000 500
```

This program will use parallel method to calculate voltage inside 200 by 200 grid.
Wire inside is 50 by 50 with constant voltage 5.0 V.
Iteration limit for each node is 500.
Program will utilise 9 nodes.
