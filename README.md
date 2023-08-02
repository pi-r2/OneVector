This repository contains the code implementing the attacks against UOV described in the paper https://eprint.iacr.org/2023/1131.

To run the test, execute the following code on a machine with sagemath 9.5 avalaible:

```
$ sage Attack.sage
```

The functions of interest are:

- attack(G,x) which returns a basis of the secret subspace of the public key G if x belongs to the secret subspace.
- in_secret_subspace(G,x) which returns True if x belongs to the secret subspace of the public key G and False otherwise.


