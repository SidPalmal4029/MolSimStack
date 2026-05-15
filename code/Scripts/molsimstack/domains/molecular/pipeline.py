from core import io
import os
from core.comparison.compare import compare_methods, save_results
from .methods import (
    morgan, topo, maccs, shape3d,
    atompair, torsion, descriptors,
    pharmacophore, scaffold, mcs
)
from core.features import compute_and_store_features
from core.comparison.per_molecule import (
    per_molecule_comparison,
    save_per_molecule,
    summarize_per_molecule
)

DEFAULT_METHODS = [
    "morgan", "topo", "maccs", "shape3d"
]

ALL_METHODS = [
    "morgan", "topo", "maccs", "shape3d",
    "atompair", "torsion", "descriptors",
    "pharmacophore", "scaffold", "mcs"
]

METHOD_MAP = {
    "morgan": morgan,
    "topo": topo,
    "maccs": maccs,
    "shape3d": shape3d,
    "atompair": atompair,
    "torsion": torsion,
    "descriptors": descriptors,
    "pharmacophore": pharmacophore,
    "scaffold": scaffold,
    "mcs": mcs
}

def run(args):
    print("[Molecular] Loading molecules")

    mols, names = io.load_molecules(
        args.input,
        outdir=args.outdir
    )
    # NEW: Per-molecule feature generation
    print("[Molecular] Computing per-molecule features")

    feature_dir = os.path.join(args.outdir, "per_molecule_features")

    compute_and_store_features(
        mols,
        names,
        feature_dir,
        save_bits=True,     # default ON
        force=False         # change to True if recompute needed
    )
    # Resolve methods
    methods = args.methods if args.methods else DEFAULT_METHODS

    # Handle "all" keyword
    if "all" in methods:
        methods = ALL_METHODS

    # Remove duplicates while preserving order
    seen = set()
    methods = [m for m in methods if not (m in seen or seen.add(m))]

    # Validate methods
    for method in methods:
        if method not in METHOD_MAP:
            raise ValueError(f"Invalid method: {method}")

    print(f"[Molecular] Methods to run: {methods}")

    # Run methods
    for method in methods:
        print(f"[Molecular] Running method: {method}")

        method_outdir = os.path.join(args.outdir, method)
        METHOD_MAP[method].run(mols, names, method_outdir, args.threads)

    # Post-analysis (only if ≥2 methods)
    if len(methods) >= 2:
        print("[Molecular] Running cross-method comparison")

        results = compare_methods(args.outdir, methods)
        if results:
            save_results(results, args.outdir)

        print("[Molecular] Running per-molecule comparison")

        df = per_molecule_comparison(args.outdir, methods)
        if df is not None:
            save_per_molecule(df, args.outdir)

            # Optional summary (if implemented)
            try:
                summarize_per_molecule(df, args.outdir)
            except NameError:
                pass
