# sparse_pauli
Implementation of large, sparse Pauli operators using pairs of sets. Contains absolutely minimal functionality, and is an example of extreme "physicist code".

## Rough Description
Sometimes, with a Pauli, you only care about where it has X support and Z support (like when you're checking cnot-based circuits that measure stabilisers in topological codes). For this, I'm going to try storing these supports as sets. This whole project should use IntSet's in Julia, but I want to get this done ASAP. Anyone who wants to port, feel free.