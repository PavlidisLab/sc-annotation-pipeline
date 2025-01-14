#!/bin/python

import argparse
import json
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Extract metadata from a JSON file and save to a TSV file.")
    parser.add_argument("json_file", type=str, help="Path to the JSON file.")
   # parser.add_argument("output_file", type=str, help="Path to save the output metadata TSV file.")
    args = parser.parse_args()
    
    # Process the JSON file and save metadata
    process_json_file(args.json_file)

def process_sample(sample):
    sample_dict = {}
    
    # Extract basic details using .get()
    sample_id = sample.get("id")
    accession = sample.get("accession", {}).get("accession")
    organism = sample.get("arrayDesign", {}).get("taxon", {}).get("scientificName", "").replace(" ", "_").lower()
    
    # Add basic details to the dictionary
    sample_dict["sample_id"] = sample_id
    sample_dict["accession"] = accession
    sample_dict["organism"] = organism
    
    # Extract characteristics
    for characteristic in sample.get("sample", {}).get("characteristics", []):
        category = characteristic.get('category', '')
        value = characteristic.get('value', '')
        if category:
            sample_dict[category] = value
    
    # Extract additional details if available
    bio_source = sample.get("sample", {}).get("BioSource")
    gene = sample.get("Gene")
    
    if bio_source:
        sample_dict["BioSource"] = bio_source
    if gene:
        sample_dict["Gene"] = gene

    for item in sample.get("sample", {}).get("factorValues", []):
        characteristics = item.get("characteristics", [])
        for characteristic in characteristics:
            if characteristic.get("category") == "developmental stage":
                sample_dict["developmental_stage"] = characteristic.get("value")
        
    # Convert to DataFrame
    return pd.DataFrame([sample_dict])

def process_json_file(json_file, output_file):
    sample_meta = pd.DataFrame()
    
    # Read and process the JSON file
    with open(json_file) as f:
        print(f"Processing file: {json_file}")
        data = json.load(f)
        for sample in data["data"]:
            # Process each sample and collect its metadata
            sample_data = process_sample(sample)
            # Append to the DataFrame
            sample_meta = pd.concat([sample_meta, pd.DataFrame(sample_data)], ignore_index=True)
            
    output_file = json_file.replace(".json", ".tsv")
    # Save the metadata to a TSV file
    sample_meta.to_csv(output_file, sep="\t", index=False)
    print(f"Saved metadata to {output_file}")


if __name__ == "__main__":
    main()
