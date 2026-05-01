import argparse
import urllib.request
from pathlib import Path

KIRCHER_ELEMENTS_URL = (
    "https://raw.githubusercontent.com/kircherlab/"
    "MPRA_SaturationMutagenesis/master/data/elements.tsv.gz"
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="data/raw/kircher_elements.tsv.gz")
    parser.add_argument("--url", default=KIRCHER_ELEMENTS_URL)
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    print(f"Downloading real MPRA data from: {args.url}")
    print(f"Writing to: {output}")

    urllib.request.urlretrieve(args.url, output)

    if output.stat().st_size == 0:
        raise RuntimeError("Downloaded file is empty.")

    print(f"Downloaded {output.stat().st_size:,} bytes")

if __name__ == "__main__":
    main()
