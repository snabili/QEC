This repo is intended to simulate QEC codes

Import the following packages:

```
conda create --name qecenv python=3.11 # qiskit_aer does not work with newer python versions
conda activate qecenv
conda install pip

pip install stim # to simulate circuits using surface code using Google Quantum AI library, for high speed simulation of error correction
pip install pymatching # this will import numpy, matplotlib, scipy; used to decode the error
pip install qiskit
pip install qiskit_aer
pip install qiskit_experiments
pip install pylatexenc # to visualize circuits
pip install pandas
pip install seaborn
```

The CSV files dowloaded from IBM contain these error information per qubit:

$T_1$ (Decoherence): The time it takes for a qubit to relax from |1> to $|0>

$T_2$ (Dephasing): The time it takes for the qubit to lose its quantum phase.

Readout Error: The probability that a |0> is measured as a |1> (or vice versa).

Gate Errors: The specific error rate for operations like sx, x, and cx.

