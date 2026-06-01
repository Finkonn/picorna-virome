#!/usr/bin/env python3

import argparse
import pandas as pd
import glob
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert DIAMOND results into readable summary table"
    )

    parser.add_argument(
        "-d",
        "--diamond",
        required=True,
        nargs="+",  # Accept one or more files
        help="DIAMOND best hits table(s). Can specify multiple files or use wildcards (requires --pattern flag)"
    )

    parser.add_argument(
        "-p",
        "--pattern",
        action="store_true",
        help="Treat --diamond as a file pattern (e.g., 'results_*.tsv')"
    )

    parser.add_argument(
        "-m",
        "--metadata",
        required=True,
        help="SRA metadata table"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="summary_table.tsv",
        help="Output table"
    )

    parser.add_argument(
        "--merge-strategy",
        choices=["union", "intersection", "per_file"],
        default="union",
        help="How to combine multiple DIAMOND files: 'union' (all hits, default), 'intersection' (hits in all files), 'per_file' (add filename column)"
    )

    return parser.parse_args()


def load_diamond_files(diamond_inputs, use_pattern=False):
    """Load one or more DIAMOND files"""
    all_dfs = []
    file_sources = []
    
    if use_pattern:
        # Expand wildcard pattern
        if len(diamond_inputs) > 1:
            print("Warning: With --pattern, only the first argument is used as pattern")
        pattern = diamond_inputs[0]
        files = glob.glob(pattern)
        if not files:
            raise ValueError(f"No files found matching pattern: {pattern}")
        print(f"Found {len(files)} files matching pattern: {pattern}")
        file_sources = files
    else:
        file_sources = diamond_inputs
    
    # Load each file
    for filepath in file_sources:
        try:
            df = pd.read_csv(filepath, sep="\t")
            df["_source_file"] = os.path.basename(filepath)
            all_dfs.append(df)
            print(f"Loaded: {filepath} ({len(df)} rows)")
        except Exception as e:
            print(f"Error loading {filepath}: {e}")
    
    if not all_dfs:
        raise ValueError("No DIAMOND files could be loaded")
    
    return all_dfs


def merge_diamond_data(all_dfs, merge_strategy):
    """Combine multiple DIAMOND dataframes based on strategy"""
    if merge_strategy == "union":
        # Combine all hits
        combined = pd.concat(all_dfs, ignore_index=True)
        print(f"Union: {len(combined)} total rows from all files")
        
    elif merge_strategy == "intersection":
        # Keep only qseqid that appear in ALL files
        qseqid_sets = [set(df["qseqid"]) for df in all_dfs]
        common_qseqids = set.intersection(*qseqid_sets)
        
        # Combine dataframes but keep only common qseqids
        combined = pd.concat(all_dfs, ignore_index=True)
        combined = combined[combined["qseqid"].isin(common_qseqids)]
        print(f"Intersection: {len(common_qseqids)} common qseqids across all files")
        print(f"  {len(combined)} total rows after intersection")
        
    elif merge_strategy == "per_file":
        # Keep all data but add source file info (already added during loading)
        combined = pd.concat(all_dfs, ignore_index=True)
        print(f"Per file: {len(combined)} total rows (source info preserved)")
    
    # Remove duplicate hits (keep best based on bitscore if duplicates exist)
    if merge_strategy != "per_file":
        combined = combined.sort_values(["qseqid", "bitscore"], ascending=[True, False])
        combined = combined.drop_duplicates(subset=["qseqid"], keep="first")
        print(f"After deduplication: {len(combined)} rows")
    
    return combined


def first_nonempty(row, columns):
    for col in columns:
        if col not in row:
            continue
        value = row[col]
        if pd.notna(value) and str(value).strip() != "":
            return value
    return None


def extract_info(row):
    info_parts = []
    for col in [
        "BioProject",
        "SRA Study",
        "Experiment",
        "Submitter_Id",
        "project_name",
        "manuscript_id",
        "_source_file",  # Include source file if present
    ]:
        if col in row and pd.notna(row[col]) and str(row[col]).strip() != "":
            info_parts.append(f"{col}: {row[col]}")
    return " | ".join(info_parts)


def main():
    args = parse_args()
    
    # Load all DIAMOND files
    all_diamond_dfs = load_diamond_files(args.diamond, args.pattern)
    
    # Merge based on strategy
    diamond = merge_diamond_data(all_diamond_dfs, args.merge_strategy)
    
    # Load metadata
    try:
        if args.metadata.endswith(('.xlsx', '.xls')):
            meta = pd.read_excel(args.metadata)
        else:
            meta = pd.read_csv(args.metadata, sep="\t")
        print(f"Loaded metadata: {args.metadata} ({len(meta)} rows)")
    except Exception as e:
        print(f"Error loading metadata: {e}")
        return
    
    # Process DIAMOND results
    diamond["Sample"] = diamond["qseqid"].apply(lambda x: x.split("|")[0])
    diamond["Profile"] = diamond["qseqid"].apply(
        lambda x: x.split("|")[1].split("/")[1] if len(x.split("|")) > 1 and len(x.split("|")[1].split("/")) > 1 else ""
    )
    diamond["Contig_length"] = diamond["qseqid"].apply(
        lambda x: x.split("|")[1].split("/")[2] if len(x.split("|")) > 1 and len(x.split("|")[1].split("/")) > 2 else ""
    )
    
    diamond["Blast_hit"] = diamond["stitle"].apply(
        lambda x: (x.split(" ", 1)[1].replace(";", " ")) if " " in str(x) else x
    )
    
    # Merge with metadata
    merged = diamond.merge(
        meta,
        left_on="Sample",
        right_on="Run",
        how="left"
    )
    
    print(f"Merged data: {len(merged)} rows")
    
    # Extract additional fields
    merged["Country_final"] = merged.apply(
        lambda row: first_nonempty(
            row,
            [
                "geo_loc_name_country",
                "geo_loc_name",
                "geographic_location_(country_and/or_sea)",
                "geo_loc_name_country_continent",
            ]
        ),
        axis=1
    )
    
    merged["Host_final"] = merged.apply(
        lambda row: first_nonempty(
            row,
            [
                "Host",
                "host_scientific_name",
                "Scientific_Name",
                "Organism",
            ]
        ),
        axis=1
    )
    
    merged["Year_raw"] = merged.apply(
        lambda row: first_nonempty(
            row,
            [
                "Collection_Date",
                "ReleaseDate",
                "create_date",
            ]
        ),
        axis=1
    )
    merged["Year"] = merged["Year_raw"].astype(str).str.extract(r"(\d{4})")[0]
    
    merged["Info"] = merged.apply(extract_info, axis=1)
    
    # Prepare final dataframe
    final_columns = {
        "Образец": merged["Sample"],
        "Страна": merged["Country_final"],
        "Год": merged["Year"],
        "Хозяин": merged["Host_final"],
        "Длина контига": merged["Contig_length"],
        "Профиль": merged["Profile"],
        "Схожесть (%)": merged["pident"],
        "Нач. контига": merged["qstart"],
        "Кон. контига": merged["qend"],
        "Нач. реф.": merged["sstart"],
        "Кон. реф.": merged["send"],
        "Дл. выравн.": merged["align_len"],
        "Покр. контига (%)": merged["qcovhsp"],
        "Покр. реф. (%)": merged["scovhsp"],
        "Бласт хит": merged["Blast_hit"],
        "Инфо": merged["Info"]
    }
    
    # Add source file column if per_file strategy was used
    if args.merge_strategy == "per_file" and "_source_file" in merged.columns:
        final_columns["Источник"] = merged["_source_file"]
    
    final = pd.DataFrame(final_columns)
    
    # Save output
    final.to_csv(
        args.output,
        sep="\t",
        index=False,
        encoding="utf-8-sig"
    )
    
    print(f"\nSaved: {args.output}")
    print(f"Total rows in output: {len(final)}")
    print(f"Unique samples: {final['Образец'].nunique()}")


if __name__ == "__main__":
    main()