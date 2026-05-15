import os
import csv
import pickle

from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    Crippen,
    Lipinski,
    rdMolDescriptors,
    AllChem,
    MACCSkeys
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.AtomPairs import Pairs, Torsions

# 1. DESCRIPTORS + SCAFFOLds

def compute_features(mol):
    """ Function for the computation of descriptors and other data
        for ech molecules that has bee loaded from the input directory"""
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "RotBonds": Lipinski.NumRotatableBonds(mol),
        "Charge": Chem.GetFormalCharge(mol),
        "Scaffold": Chem.MolToSmiles(
            MurckoScaffold.GetScaffoldForMol(mol)
        )
    }

def compute_fingerprints(mols):
  """Function for generation of the fingerprints of the molecules"""
    fps = {
        "morgan": [],
        "maccs": [],
        "atompair": [],
        "torsion": []
    }

    for mol in mols:
        # Morgan (ECFP4)
        fps["morgan"].append(
            AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        )

        # MACCS
        fps["maccs"].append(
            MACCSkeys.GenMACCSKeys(mol)
        )

        # Atom Pair
        fps["atompair"].append(
            Pairs.GetAtomPairFingerprintAsBitVect(mol)
        )

        # Topological Torsion
        fps["torsion"].append(
            Torsions.GetTopologicalTorsionFingerprintAsBitVect(mol)
        )

    return fps


def save_descriptors(features, names, outdir):
    outfile = os.path.join(outdir, "features.csv")

    headers = ["Name"] + list(features[0].keys())

    with open(outfile, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)

        for name, feat in zip(names, features):
            row = [name] + [feat[h] for h in headers[1:]]
            writer.writerow(row)


def save_fingerprints(fps_dict, names, outdir):
    for fp_type, fps in fps_dict.items():
        outfile = os.path.join(outdir, f"fingerprints_{fp_type}.pkl")

        with open(outfile, "wb") as f:
            pickle.dump(dict(zip(names, fps)), f)


def save_bitstrings(fps_dict, names, outdir):
    """
    Optional human-readable version (only for bit vectors)
    """
    for fp_type, fps in fps_dict.items():

        # Skip non-bitvect (safety, though current ones are bitvect)
        if not hasattr(fps[0], "ToBitString"):
            continue

        outfile = os.path.join(outdir, f"{fp_type}_bits.csv")

        with open(outfile, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Name", f"{fp_type}_bits"])

            for name, fp in zip(names, fps):
                writer.writerow([name, fp.ToBitString()])


def save_metadata(outdir):
    """
    Stores fingerprint parameters for reproducibility
    """
    meta_file = os.path.join(outdir, "fingerprint_metadata.csv")

    metadata = [
        ["type", "details"],
        ["morgan", "radius=2,nBits=2048"],
        ["maccs", "166-bit MACCS keys"],
        ["atompair", "RDKit default AtomPair bit vector"],
        ["torsion", "RDKit topological torsion bit vector"]
    ]

    with open(meta_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(metadata)

def compute_and_store_features(mols, names, outdir,
                               save_bits=True,
                               force=False):

    os.makedirs(outdir, exist_ok=True)

    feat_file = os.path.join(outdir, "features.csv")

    if os.path.exists(feat_file) and not force:
        print("[Features] Existing features found. Skipping computation.")
        return

    print("[Features] Computing descriptors")
    features = [compute_features(m) for m in mols]
    save_descriptors(features, names, outdir)


    print("[Features] Computing fingerprints")
    fps_dict = compute_fingerprints(mols)

    print("[Features] Saving fingerprints (PKL)")
    save_fingerprints(fps_dict, names, outdir)

    print("[Features] Saving bitstring versions")
    save_bitstrings(fps_dict, names, outdir)

    print("[Features] Writing metadata")
    save_metadata(outdir)

    print("[Features] Completed")
