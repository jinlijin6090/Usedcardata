# ReaxFF sanity checks

force field: `ff/ffield.reax.FC`

| quantity | ReaxFF | reference | unit |
|---|---|---|---|
| Ni fcc lattice constant | 3.612 | 3.524 (exp) | A |
| Ni cohesive energy | 4.497 | 4.44 (exp) | eV/atom |
| graphene C-C bond | 1.433 | 1.42 (exp) | A |
| graphene cohesive energy | 7.849 | 7.9 (exp, ~7.4-7.9) | eV/atom |
| CO bond length | 1.126 | 1.128 (exp) | A |
| CO dissociation energy | 11.530 | 11.1 (exp) | eV |
| O2 bond length | 1.263 | 1.207 (exp) | A |
| O2 dissociation energy | 5.262 | 5.1-5.2 (exp) | eV |
| CO2 C=O bond length | 1.253 | 1.160 (exp) | A |
| CO2 atomization energy | 16.935 | 16.6 (exp) | eV |
| CO2 -> CO + O(atom) energy | 5.405 | 5.5 (exp) | eV |
| C on Ni(111) fcc (vs gas atom) | -6.729 | -6.5 to -7.0 (GGA) | eV |
| C on Ni(111) hcp (vs gas atom) | -6.721 | -6.5 to -7.0 (GGA) | eV |
| O on Ni(111) fcc (vs gas atom) | -5.886 | -5.0 to -5.6 (GGA) | eV |
| O on Ni(111) hcp (vs gas atom) | -5.758 | -5.0 to -5.6 (GGA) | eV |
| CO on Ni(111) top (vs gas CO) | -1.740 | -1.3 (exp) / -1.5 to -1.9 (GGA) | eV |
| CO on Ni(111) fcc (vs gas CO) | -1.740 | -1.3 (exp) / -1.5 to -1.9 (GGA) | eV |
| C* + O* -> CO* (surface reaction energy) | -0.663 | ~ -1 to +0.5 (GGA, coverage dep.) | eV |
| C* + O* -> CO(g) | 1.077 | ~ +0.3 to +1.0 (GGA) | eV |

Reference values are approximate; replace with your own DFT numbers (same functional for all entries) before drawing quantitative conclusions.
