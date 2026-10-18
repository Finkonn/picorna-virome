#!/usr/bin/env python3

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert DIAMOND results into readable summary table"
    )

    parser.add_argument(
        "-d",
        "--diamond",
        required=True,
        help="DIAMOND best hits table"
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

    return parser.parse_args()


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
    ]:

        if col in row and pd.notna(row[col]) and str(row[col]).strip() != "":
            info_parts.append(f"{col}: {row[col]}")

    return " | ".join(info_parts)


def main():

    args = parse_args()

    diamond = pd.read_csv(args.diamond, sep="\t")

    meta = pd.read_csv(args.metadata, sep=",", low_memory=False)

    diamond["Sample"] = diamond["qseqid"].apply(
        lambda x: x.split("|")[0]
    )

    diamond["Profile"] = diamond["qseqid"].apply(
        lambda x: x.split("|")[1].split("/")[1]
    )

    diamond["Contig_length"] = diamond["qseqid"].apply(
        lambda x: x.split("|")[1].split("/")[2]
    )

    diamond["Blast_hit"] = diamond["stitle"].apply(
        lambda x: (
            x.split(" ", 1)[1]
            .replace(";", " ")
        )
        if " " in str(x)
        else x
    )

    merged = diamond.merge(
        meta,
        left_on="Sample",
        right_on="Run",
        how="left"
    )

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

    merged["Year"] = (
        merged["Year_raw"]
        .astype(str)
        .str.extract(r"(\d{4})")[0]
    )

    merged["Info"] = merged.apply(
        extract_info,
        axis=1
    )


    final = pd.DataFrame({
        "Образец": merged["Sample"],
        "Страна": merged["Country_final"],
        "Год извлечения": merged["Year"],
        "Хозяин": merged["Host_final"],
        "Длина контига": merged["Contig_length"],
        "Профиль": merged["Profile"],
        "На что бластанулось": merged["Blast_hit"],
        "Bitscore": merged["bitscore"],
        "Информация": merged["Info"]
    })

    final.to_csv(
        args.output,
        sep="\t",
        index=False,
        encoding="utf-8-sig"
    )

    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()