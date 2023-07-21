This repository contains code used in the backpropagation processing of radio occultations.

Code log:

- Original: Code used in the AMT publication.

- 2023/07/21 (branch): Antenna pattern correction included at the step where the instrument noise to be added to the signal is calculated.
Previously, there was extra noise added to MPS signal, which yielded a bias in S4 index.