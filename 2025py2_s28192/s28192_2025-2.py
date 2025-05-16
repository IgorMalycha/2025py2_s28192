#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever
Podstawowy skrypt do łączenia się z NCBI i pobierania rekordów sekwencji genetycznych dla danego identyfikatora taksonomicznego.
"""

from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid):
        print(f"Searching for records with taxID: {taxid}")
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            search_term = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print(f"No records found for {organism_name}")
                return None

            print(f"Found {count} records")
            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count

            return count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=10):
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            batch_size = min(max_records, 500)
            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )
            return list(SeqIO.parse(handle, "genbank"))
        except Exception as e:
            print(f"Error fetching records: {e}")
            return []

def save_to_csv(records, filename, min_len, max_len):
    data = []
    for rec in records:
        length = len(rec.seq)
        if min_len <= length <= max_len:
            data.append((rec.id, length, rec.description))
    df = pd.DataFrame(data, columns=["Accession", "Length", "Description"])
    df.to_csv(filename, index=False)
    return df

def plot_lengths(df, output_file):
    df = df.sort_values(by="Length", ascending=False)
    plt.figure(figsize=(12, 6))
    plt.plot(df["Accession"], df["Length"], marker="o")
    plt.xticks(rotation=90, fontsize=6)
    plt.xlabel("GenBank Accession")
    plt.ylabel("Sequence Length")
    plt.title("Filtered Sequence Lengths")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def main():
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")
    min_len = int(input("Enter minimum sequence length: "))
    max_len = int(input("Enter maximum sequence length: "))

    retriever = NCBIRetriever(email, api_key)
    count = retriever.search_taxid(taxid)

    if not count:
        print("No records found. Exiting.")
        return

    print("\nFetching records...")
    records = retriever.fetch_records(start=0, max_records=50)
    csv_file = f"taxid_{taxid}_filtered.csv"
    df = save_to_csv(records, csv_file, min_len, max_len)
    print(f"Saved filtered records to {csv_file}")

    plot_file = f"taxid_{taxid}_plot.png"
    plot_lengths(df, plot_file)
    print(f"Saved plot to {plot_file}")

if __name__ == "__main__":
    main()
