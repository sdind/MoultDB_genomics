import csv

# Read input file path from Snakemake
input_file = snakemake.input[0]

# Read output file path from Snakemake
output_file = snakemake.output[0]

with open(input_file, newline='') as csvfile, open(output_file, 'w') as outfile:
    reader = csv.DictReader(csvfile, delimiter='\t')  # Assuming the file is tab-delimited
    for row in reader:
        ncbi_id = row['ncbi_id']
        taxid = row['taxid']
        combined_string = f"{taxid}_{ncbi_id}"
        outfile.write(combined_string + '\n')



