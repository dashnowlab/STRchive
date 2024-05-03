import json

def generate_gene_pages(json_file, template_file, output_dir):
    # Read the JSON file
    genes = set()
    with open(json_file, 'r') as f:
        data = json.load(f)
        for x in data:
            genes.add(x['gene'])

    # Read the template file
    with open(template_file, 'r') as f:
        template = f.read()

    # Generate HTML page for each gene
    for gene in genes:
        gene_html = template.replace('<<gene>>', gene)

        # Write the HTML page to the output directory
        output_file = f"{output_dir}/{gene}.html"
        with open(output_file, 'w') as f:
            f.write(gene_html)