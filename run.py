# Runs the GENESIS simulation (Propagation/propagator.py) with the inputs in Inputs/ and writes the
# results to Outputs/.
#     python run.py
from Propagation.propagator import run

if __name__ == '__main__':
    run()
    print("Sim finished!")
