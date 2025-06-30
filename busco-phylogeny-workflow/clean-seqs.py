import os
import shutil
from pathlib import Path
import argparse

def process_faa_files(base_dir):
    print(f"[INFO] Scanning subdirectories in: {base_dir}")
    for subdir in sorted(Path(base_dir).iterdir()):
        if subdir.is_dir():
            print(f"[INFO] Entering: {subdir}")
            faa_files = list(subdir.glob("*.faa"))
            if not faa_files:
                print(f"[WARNING] No .faa files found in {subdir}")
            for faa_file in faa_files:
                print(f"[INFO] Processing {faa_file.name}")
                faa2_file = faa_file.with_suffix(faa_file.suffix + "2")
                cut_file = faa2_file.with_name(faa2_file.name + ".cut")

                # Keep only first FASTA record
                with faa_file.open() as infile, faa2_file.open("w") as outfile:
                    n = 0
                    for line in infile:
                        if line.startswith(">"):
                            n += 1
                            if n > 1:
                                break
                        outfile.write(line)

                # Modify headers: add subdir name, remove slashes
                with faa2_file.open("r") as infile:
                    content = infile.read()
                content = content.replace(">", f">{subdir.name} ").replace("/", "")
                with faa2_file.open("w") as outfile:
                    outfile.write(content)

                # Keep only first word, remove '*'
                with faa2_file.open("r") as infile, cut_file.open("w") as outfile:
                    for line in infile:
                        if line.startswith(">"):
                            line = line.split()[0].replace("*", "") + "\n"
                        outfile.write(line)

                # Replace original file
                shutil.move(str(cut_file), faa_file)
                faa2_file.unlink()
                print(f"[INFO] Finished: {faa_file.name}")

def concatenate_files(base_dir, script_dir):
    busco_file = script_dir / "busco-seqs.txt"
    if not busco_file.exists():
        print(f"[ERROR] Missing busco-seqs.txt at {busco_file}")
        return

    with busco_file.open() as f:
        for line in f:
            gene_file = line.strip()
            if not gene_file:
                continue
            print(f"[INFO] Concatenating: {gene_file}")
            with open(base_dir / gene_file, "w") as outfile:
                for subdir in sorted(base_dir.iterdir()):
                    if subdir.is_dir():
                        file_path = subdir / gene_file
                        if file_path.exists():
                            with file_path.open() as infile:
                                shutil.copyfileobj(infile, outfile)
                        else:
                            print(f"[WARNING] {file_path} not found")

def move_final_faa_files(base_dir):
    outdir = Path(base_dir) / "prot-seqs"
    outdir.mkdir(exist_ok=True)
    for faa in Path(base_dir).glob("*.faa"):
        print(f"[INFO] Moving {faa.name} to prot-seqs/")
        shutil.move(str(faa), outdir / faa.name)

def main():
    parser = argparse.ArgumentParser(description="Process and concatenate BUSCO .faa files.")
    parser.add_argument("-d", "--directory", required=True, help="Base directory containing BUSCO output folders")
    args = parser.parse_args()

    base_dir = Path(args.directory).resolve()
    script_dir = Path(__file__).resolve().parent

    if not base_dir.exists():
        raise FileNotFoundError(f"Directory not found: {base_dir}")

    process_faa_files(base_dir)
    concatenate_files(base_dir, script_dir)
    move_final_faa_files(base_dir)
    print(f"[DONE] Processing complete.")

if __name__ == "__main__":
    main()
