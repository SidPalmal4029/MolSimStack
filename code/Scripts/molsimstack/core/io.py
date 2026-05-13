from rdkit import Chem
import os
import glob
import pandas as pd
import numpy as np



def load_molecules(input_path, outdir=None):
    mols = []
    names = []
    failed_entries = []

    total = 0

    # Case 1: Directory input
    if os.path.isdir(input_path):
        sdf_files = sorted(glob.glob(os.path.join(input_path, "*.sdf")))

        if len(sdf_files) == 0:
            raise ValueError("No .sdf files found in directory")

        print(f"[IO] Detected directory input with {len(sdf_files)} SDF files")

        for file_idx, sdf_file in enumerate(sdf_files):
            supplier = Chem.SDMolSupplier(sdf_file)

            for mol_idx, mol in enumerate(supplier):
                total += 1

                if mol is None:
                    failed_entries.append(f"{sdf_file}::record_{mol_idx}")
                    continue

                mols.append(mol)

                # Name priority
                if mol.HasProp("_Name"):
                    names.append(mol.GetProp("_Name"))
                else:
                    names.append(os.path.basename(sdf_file))

    # Case 2: Single SDF file
    else:
        supplier = Chem.SDMolSupplier(input_path)

        print(f"[IO] Detected single SDF input")

        for i, mol in enumerate(supplier):
            total += 1

            if mol is None:
                failed_entries.append(f"{input_path}::record_{i}")
                continue

            mols.append(mol)

            if mol.HasProp("_Name"):
                names.append(mol.GetProp("_Name"))
            else:
                names.append(f"mol_{i}")

    # Summary
    valid = len(mols)
    failed = len(failed_entries)

    print(f"[IO] Total records     : {total}")
    print(f"[IO] Valid molecules  : {valid}")
    print(f"[IO] Failed molecules : {failed}")

    if failed > 0:
        print(f"[IO] Failed entries  : {failed_entries[:5]}{' ...' if failed > 5 else ''}")

    # Write failure log
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        log_path = os.path.join(outdir, "failed_molecules.txt")

        with open(log_path, "w") as f:
            for entry in failed_entries:
                f.write(f"{entry}\n")

        print(f"[IO] Failure log written to: {log_path}")

    # Critical failure
    if valid == 0:
        raise ValueError("No valid molecules found. Cannot proceed.")

    # Preview
    print(f"[IO] Proceeding with {valid} molecules")

    preview_n = min(5, valid)
    print("[IO] Example molecules:")
    for i in range(preview_n):
        print(f"   - {names[i]}")

    return mols, names


