#!/usr/bin/env python3

import pandas as pd
from qiskit_aer.noise import NoiseModel, QuantumError, thermal_relaxation_error, depolarizing_error
from qiskit.circuit import QuantumCircuit

# ------------------------------------------------------------
# 1. Load CSV
# ------------------------------------------------------------
df = pd.read_csv("datafiles/ibm_fez_calibrations_2026-01-06T18_10_35Z.csv")

# ------------------------------------------------------------
# 2. Parse CZ, RZZ, and gate times from your CSV
# ------------------------------------------------------------
cz_errors = {}
cz_gate_times = {}
rzz_errors = {}

def parse_pairwise(entry, source_qubit):
    out = {}
    if isinstance(entry, str) and entry.strip():
        for item in entry.split(";"):
            if ":" in item:
                tgt, val = item.split(":")
                pair = tuple(sorted((source_qubit, int(tgt))))
                out[pair] = float(val)
    return out

for _, row in df.iterrows():
    q = int(row["Qubit"])
    cz_errors.update(parse_pairwise(row["CZ error"], q))
    cz_gate_times.update(parse_pairwise(row["Gate length (ns)"], q))
    rzz_errors.update(parse_pairwise(row["RZZ error"], q))
print("done with parsing: CZ error, Gate length (ns), and RZZ error")
#print("   cz_errors: ",cz_errors," cz_gate_times: ",cz_gate_times, "rzz_errors: ", rzz_errors)
# ------------------------------------------------------------
# 3. Helper: T1/T2 extraction
# ------------------------------------------------------------
def get_T1_T2_eff(q):
    row_q = df.loc[df["Qubit"] == q].iloc[0]
    T1 = float(row_q["T1 (us)"]) * 1e-6 # to seconds
    T2 = float(row_q["T2 (us)"]) * 1e-6 # to seconds
    T2_eff = min(T2, 2*T1 - 1e-12)
    return T1, T2_eff

# ------------------------------------------------------------
# 4. Build identity noise channel
# ------------------------------------------------------------
qc_id = QuantumCircuit(1)
id_noise = QuantumError([(qc_id, 1.0)])
print("done with building identity noise channel")
# ------------------------------------------------------------
# 5. Build H-gate noise
# ------------------------------------------------------------
def build_h_noise(q):
    row_q = df.loc[df["Qubit"] == q].iloc[0]

    T1, T2_eff = get_T1_T2_eff(q)

    t1q_ns = float(row_q["Single-qubit gate length (ns)"])
    t1q = t1q_ns * 1e-9 # to seconds

    tr = thermal_relaxation_error(T1, T2_eff, t1q)

    p_sx = float(row_q["√x (sx) error"])
    p_rx = float(row_q["RX error"])
    p_x  = float(row_q["Pauli-X error"])

    depol_sx = depolarizing_error(p_sx, 1)
    depol_rx = depolarizing_error(p_rx, 1)
    depol_x  = depolarizing_error(p_x, 1)

    return tr.compose(depol_sx).compose(depol_rx).compose(depol_x)

# ------------------------------------------------------------
# 6. Build CNOT noise model
# ------------------------------------------------------------
noise_model = NoiseModel()

for idx, ((i, j), cz_err) in enumerate(cz_errors.items()):
    if idx >= 3:
        break

    t_cz = cz_gate_times.get((i, j), None)
    t_cz_s = t_cz * 1e-9 if t_cz else None

    T1_i, T2_i_eff = get_T1_T2_eff(i)
    T1_j, T2_j_eff = get_T1_T2_eff(j)

    if t_cz_s:
        tr_i = thermal_relaxation_error(T1_i, T2_i_eff, t_cz_s)
        tr_j = thermal_relaxation_error(T1_j, T2_j_eff, t_cz_s)
        tr_cz = tr_i.expand(tr_j)
    else:
        tr_cz = None

    depol_cz = depolarizing_error(cz_err, 2)
    cz_noise = tr_cz.compose(depol_cz) if tr_cz else depol_cz

    h_noise = build_h_noise(j)
    H_on_j = id_noise.tensor(h_noise)

    cnot_noise = H_on_j.compose(cz_noise).compose(H_on_j)

    noise_model.add_quantum_error(cnot_noise, 'cx', [i, j])
    print("done with building CNOT noise model")

# ------------------------------------------------------------
# 7. Print summary
# ------------------------------------------------------------
print("CNOT noise model built successfully.")
print(noise_model)
