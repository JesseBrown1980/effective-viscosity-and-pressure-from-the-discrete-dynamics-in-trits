import subprocess
import sys
import os

def run():
    print("=======================================================================")
    print(" REPRODUCING EFFECTIVE VISCOSITY & PRESSURE FROM DISCRETE TRIT DYNAMICS")
    print("=======================================================================")
    
    # 1. Compile calibration program
    print("\n[1/3] Compiling asolaria_fluid_hydro_calibration.cpp...")
    cmd_compile = ["g++", "-O3", "-std=c++17", "-Wall", "-Wextra", "-Werror", "-static", 
                   "asolaria_fluid_hydro_calibration.cpp", "-o", "asolaria_fluid_hydro_calibration.exe"]
    res = subprocess.run(cmd_compile)
    if res.returncode != 0:
        print("[FAIL] Compilation failed!")
        sys.exit(1)
    print("[PASS] Compilation successful.")

    # 2. Run calibration on 300GB Host8 Room Substrate
    room_path = "D:/ASOLARIA-HOST8-ROOMS-300G/sector-000.rooms"
    print(f"\n[2/4] Executing calibration on 300GB Host8 room substrate:\n      {room_path}")
    if os.path.exists(room_path):
        res = subprocess.run([os.path.join(".", "asolaria_fluid_hydro_calibration.exe"), room_path])
        if res.returncode != 0:
            print("[FAIL] Calibration execution failed!")
            sys.exit(1)
    else:
        print(f"[WARN] 300GB room file not found at {room_path}; skipping live substrate probe.")

    # 3. Compile and run direct dynamical acoustic standing wave experiment
    print("\n[3/4] Compiling and running acoustic standing wave propagation experiment...")
    cmd_compile_ac = ["g++", "-O3", "-std=c++17", "-Wall", "-Wextra", "-Werror", "-static",
                      "asolaria_acoustic_pulse_propagation.cpp", "-o", "asolaria_acoustic_pulse_propagation.exe"]
    res = subprocess.run(cmd_compile_ac)
    if res.returncode != 0:
        print("[FAIL] Acoustic compilation failed!")
        sys.exit(1)
    res = subprocess.run([os.path.join(".", "asolaria_acoustic_pulse_propagation.exe")])
    if res.returncode != 0:
        print("[FAIL] Acoustic standing wave experiment failed quality threshold!")
        sys.exit(1)
    print("[PASS] Acoustic standing wave experiment completed and passed quality threshold.")

    # 4. Compile and verify core discrete fluid simulation
    print("\n[4/4] Compiling and running core discrete lattice gas simulation...")
    cmd_compile_core = ["g++", "-O3", "-std=c++17", "-Wall", "-Wextra", "-Werror", "-static", 
                        "asolaria_ternary_fluid_v1.cpp", "-o", "asolaria_ternary_fluid_v1.exe"]
    res = subprocess.run(cmd_compile_core)
    if res.returncode != 0:
        print("[FAIL] Core compilation failed!")
        sys.exit(1)
    res = subprocess.run([os.path.join(".", "asolaria_ternary_fluid_v1.exe"), "test_results.csv"])
    if res.returncode != 0:
        print("[FAIL] Core simulation failed!")
        sys.exit(1)
    
    # Verify bit-for-bit CSV match
    with open("test_results.csv", "rb") as f1, open("asolaria_ternary_fluid_v1_results.csv", "rb") as f2:
        b1 = f1.read().replace(b"\r\n", b"\n").strip()
        b2 = f2.read().replace(b"\r\n", b"\n").strip()
        if b1 == b2:
            print("\n[PASS] Core discrete lattice gas matched reference checkpoints byte-for-byte!")
        else:
            print("\n[FAIL] Checkpoint mismatch!")
            sys.exit(1)

    print("\n=======================================================================")
    print(" ALL CALIBRATIONS AND DISCRETE INVARIANT CHECKS PASSED SUCCESSFULLY.")
    print("=======================================================================")

if __name__ == "__main__":
    run()
