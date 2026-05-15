from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit import QuantumCircuit
from qiskit_experiments.library import StateTomography

import numpy as np
import os

# connecting to ibm
service = QiskitRuntimeService()
backend = service.backend("ibm_fez")

def job_state(n1,n2,i,j):
    qc = QuantumCircuit(2)
    if i > 0: qc.x(0)
    if j > 0: qc.x(1)
    qc.measure_all()
    exp = StateTomography(qc)
    exp.set_transpile_options(translation_method='translator',optimization_level=1 )
    exp.set_run_options(meas_level=1, meas_return='single')
    job = exp.run(backend,qubits=[n1,n2],shots=15000)
    job.block_for_results()
    
    iq1, iq2 = job.data()[0]['memory'], job.data()[0]['memory']
    I1,Q1 = iq1[:, 0, 0], iq1[:, 0, 1]                # slot 0, I component
    I2,Q2 = iq2[:, 1, 0], iq2[:, 1, 1] 
    return I1,Q1, I2,Q2



I1_00, Q1_00, I2_00, Q2_00 = job_state(127, 128, 0, 0) # both qubits in state |0>
I1_11, Q1_11, I2_11, Q2_11 = job_state(127, 128, 1, 1) # both qubits in state |1>
I1_01, Q1_01, I2_01, Q2_01 = job_state(127, 128, 0, 1) # 127: |0>, 128: |1>
I1_10, Q1_10, I2_10, Q2_10 = job_state(127, 128, 1, 0) # 127: |1>, 128: |0>


X00 = np.column_stack([I1_00, Q1_00, I2_00, Q2_00])  # shape (shots_00, 4)
X01 = np.column_stack([I1_01, Q1_01, I2_01, Q2_01])  # shape (shots_01, 4)
X10 = np.column_stack([I1_10, Q1_10, I2_10, Q2_10])  # shape (shots_10, 4)
X11 = np.column_stack([I1_11, Q1_11, I2_11, Q2_11])  # shape (shots_11, 4)


np.savez(os.getcwd()+'/../files/datafiles/x.npz', X00=X00, X01=X01, X10=X10, X11=X11)