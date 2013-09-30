ar
==

Build
===

```
$ make
```

Examples
===

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

Graphs
====

![sequential graph](http://dl.dropbox.com/u/709944/Screenshots/6wp9.png)

![parallel graph](http://dl.dropbox.com/u/709944/Screenshots/8932.png)
