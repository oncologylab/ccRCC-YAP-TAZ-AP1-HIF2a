from pathlib import Path
import argparse
from figures.figure_6.process import verify_d as verify


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    result = verify()
    if a.output:
        if a.output.exists():
            raise FileExistsError(a.output)
        a.output.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(a.output, sep="\t", index=False)
    print("PASS: 45 values; 15 adjusted P values; three 13,833-gene Wald/BH families.")
