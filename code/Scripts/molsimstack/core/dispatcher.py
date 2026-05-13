from domains.molecular import pipeline as mol_pipeline

def main():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--domain", required=True)
    parser.add_argument("--methods", nargs="+")
    parser.add_argument("--threads", type=int, default=1)

    args = parser.parse_args()

    if args.domain == "mol":
        mol_pipeline.run(args)
    else:
        raise ValueError("Only molecular domain implemented for now")

if __name__ == "__main__":
    main()
