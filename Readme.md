# Quantum Error Correction (QEC) Study
## Purpose
Quantum computers are significantly more susceptible to noise than classical computers due to short qubit coherence times and imperfect gate operations. This repository contains implementations designed to study and benchmark Quantum Error Correction algorithms, specifically targeting IBM’s superconducting architectures, using calibration data and noise profiles from IBM Quantum's online backends. For more details about this study please read [these slides](QEC_IBM.pdf) .

## Technical Stack
This project uses a mix of:

- Simulation: Stim (High-speed Clifford circuit simulation)
- Decoding: PyMatching (Minimum Weight Perfect Matching decoder)
- Hardware Interface: Qiskit (Accessing IBM Quantum backend properties)

## Clone & Setting Environment
To clone this repository to your area:

```
git clone https://github.com/snabili/QEC.git
```

To setup the proper environment create a conda environment:

```
conda create --name qecenv python=3.10 # qiskit_aer does not work with newer python versions
conda activate qecenv
conda install pip
```

Import the following packages:

```
pip install stim # to simulate circuits using Google Quantum AI library, for high speed simulation of error correction
pip install pymatching # this will import numpy, matplotlib, scipy; used to decode the error
pip install qiskit
pip install qiskit_aer
pip install qiskit_experiments
pip install pylatexenc # to visualize circuits
pip install pandas
pip install seaborn
```

The above packages are added to qec_env.yml, to skip the above commands and to buid the environment in one run:

```
conda env create -f qec_env.yml
```

This is a first time setup and for later use run the following command to set the environment every time you want to run the code:

```
conda activate qecenv
```

## Datasets
To access the IBM calibrated dataset:
- Open an account from [ibm_webpage](https://cloud.ibm.com/login)
- Download the api key

The calibration datasets are categorized in three main classes:
- Coherence times in $\mu s$ including the decoherence and dephasing times
- Gate properties in terms of average errors for single-qubit, double-qubits (entangled qubits), and gate duration in $ns$ 
- Measurement errors including the readout error 


## Code Structure:

```
├── qec --> the source codes
│   ├── __init__.py
│   ├── config.py
│   ├── data.py
│   ├── functions.py
│   └── utils.py
├── qec_env.yml
├── Readme.md
├── test
│   ├── comp_threshold.py --> to plot $P_L$ vs $P_{phys}$ to extract threshold
│   ├── heavyhex_lattice.py --> diagnozes surface code errors
│   ├── ibm_api.py --> keep it hidden in .gitignore file
│   ├── leakage.py --> Computes leakage error based on pure simulation
│   ├── plotting.py
│   ├── simulate_noise_model.py --> Plots errors from downloaded CSV file
│   ├── stabilizer.py --> running stabilizer on surface code 
│   ├── xt_leakage.py --> crosstalk effect on EPC using fully simulated backend
│   └── xtalk_error.py --> crosstalk effect on EPC using online IBM backend; open-instance account have limited time usage
├── files
├── plots
```

## Surface Code
The gold standard of error correction code for QEC in industry (Google, IBM), because of its high fault-tolerance threshold (∼ 1%) & physical qubits space connectivity. It maps physical complex noises (e.g. leakage, decoherence) into discrete bit-flip
(X) and phase-flip (Z) errors. The elements of the surface code are: 

### Logical Qubits
Software-designed qubits grouped from many physical qubits to act as a single, reliable unit to detect and collect errors. To compute logical error for various number of qubits:
```
python test/comp_threshold.py
```
The plot produced by this code:
![My Figure](p_threshold.png)

### Stabilizers:
Operators that represent the health condition of the system.
First step in applying stabilizer in QEC code is to define the right patch with a desired distance. IBM heron heavy-hex machines are Heron class featuring 156 qubits.

To select the proper number of qubits and maitain their connectivity:

```
python test/heavyhex_lattice.py X
```

Where `X` is the central qubit in the patch. In `test/heavyhex_lattice.py` the patch distance is `d=3`, thus there are 17 qubits in total: $tot_{qubits} = 2d^2-1=17$, with $data_{qubit}=d^2=9$ and $ancilla_{qubit}=d^2-1=8$. In an ideal patch X and Z stabilizers (ancillas that measure Z and X errors are the same).

The code above produces this plot:

![My Figure](logicpatch_syndromes.png)

Stabilizers role are to find qubit flips' errors via certain qubits called Syndrome or Ancillas. Based on the type of the errors two separate stabilizers are designed

### X_Stabilizers:
Their roles is to identify the $\textbf{Phase-Flip}$ errors, the kind of errors that occurs due to Z-Pauli operation on a qubit. To identify this error, X-Syndromes are used. The steps to identify the phase-flip error on quantum circuits:

- Reset data/ancilla qubits' state to $\lvert 0 \rangle$
- Change qubit basis to X-Pauli basis by applying Hadamart gate; from $\lvert 0 \rangle \Rightarrow \lvert + \rangle$
- Apply Z (Phase-Flip) error to data qubits
- Apply CNOT gate by setting ancillas as control and data as target qubit: CNOT(A,D). This will cause the $\textbf{Phase Kickback}$ effect.
- Apply Hadamart gate to ancillas to change basis to Z-basis
- Measure ancillas

<img src="x_stabilizer.png" width="400">


### Z_Stabilizers:
Their roles is to identify the $\textbf{Bit-Flip}$ errors, the kind of errors that occurs due to X-Pauli operation on a data qubit. To identify this error, Z-Syndromes are used. The steps to identify the bit-flip error on quantum circuits:

- Reset data/ancilla qubits' state to $\lvert 0 \rangle$
- Apply X (Bit-Flip) error to data qubit
- Apply CNOT gate by setting data as control and ancillas as target qubit: CNOT(D,A)
- Measure ancillas

<img src="z_stabilizer.png" width="400">

### Running X and Z Stabilizers:
To run the X and Z stabilizer:
```
python test/stabilizer.py A B C
```

Where `A` is the central qubit of the selected patch, `B` is the data qubit next to the Z-syndrome, and `C` is the data qubit adjacent to the X-syndrome. After running the above code the syndrome outcome will be shown as follow:

```
python test/stabilizer.py 23 4 16

X_stabilizer (Should show flips in bits next to X-Syndrome qubits):
Syndrome Outcomes: {'00010100': 1024}
Z_stabilizer (Should show flip in bit for Qubit 4 error):
Syndrome Outcomes: {'00000010': 1024}
```

### Simulating Crosstalk effect:
Crosstalk arises from microwave‑pulse interference between neighboring qubits during gate execution. In practice, it refers to how a two‑qubit gate applied to a target pair can unintentionally disturb the state or performance of a nearby spectator qubit.

To observe this effect cleanly, the target pair should be chosen from qubits with low intrinsic gate error, ensuring that crosstalk is not masked by unrelated noise sources.

In the simulation, the crosstalk effect is modeled as an unintended rotation $\theta$ induced on the spectator qubit. This rotation occurs because the spectator is coupled to the target pair through a tunable coupler with interaction frequency $\xi$. The induced rotation during a CZ‑gate of duration $t_{cz}$ is approximated by:

$$
\theta = 2 \times \pi \times \xi \times t_{cz} 
$$

IBM’s tunable couplers typically limit $\xi$ to around $1 kHz$, which makes the real crosstalk effect extremely small. To make the effect visible in simulation, the code scales $\xi$ by a factor of 100.

Because IBM’s open‑access accounts provide limited access to real hardware backends, the crosstalk analysis is performed entirely on a simulated backend. To run the simulation, use the following code:

``` python test/xtalk_error.py ```

The effect of crosstalk with $\xi = 100 kHz$ is shown below with Simultaneus being the Crosstalk included effect and the Isolated legend showing without the effect:

![My Figure](xt_effect.png)

### Simulating Leakage error:
Leakage occurs when a transmon qubit is excited to energy states beyond the computational subspace ($ket{0}$ and $\ket{1}$). While the anharmonicity—provided by the nonlinear inductance of the Josephson Junction—is designed to isolate the first two levels, high-speed microwave (MW) pulses can still inadvertently drive transitions to higher energy levels, such as the $\ket{2}$ state.

To analyze these effects, this simulation utilizes the following framework:

1 - Gate Context: Analysis is focused on leakage induced during single-qubit gate operations.

2 - Three-Level System: The transmon is modeled as a qutrit with states $\ket{0}$, $\ket{1}$, and $\ket{2}$.

3 - Target Operation: The study simulates a standard X-gate (bit-flip), requiring a precise $\pi$-pulse calibration.

4 - DRAG Compensation: Used Derivative Removal by Adiabatic Gate (DRAG) technique. This uses a parameter $\beta$ to scale the derivative of the Gaussian envelope, effectively canceling out-of-subspace (phase) transitions.

5 - Parameter Optimization: Both the MW-pulse duration ($T$) and the DRAG coefficient ($\beta$) are swept to identify the "Leakage Floor" of the system.

6 - Pulse Shaping: The drive utilizes a Gaussian envelope, where the complex component (phase shift) is determined by the $\beta$-scaled derivative.

To run the leakage estimation script, execute the following:

``` python test/leakage.py ```

The simulation generates a population plot (shown below) tracking the state evolution. The "Optimized" result minimizes the population in state $\ket{2}$ at the end of the gate duration.

 ![My Figure](leakage_population.png)


