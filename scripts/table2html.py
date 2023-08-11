import csv
from jinja2 import Template

def read_csv(csvfile):
    with open(csvfile) as fh:
        data = []
        reader = csv.DictReader(fh, )
        # List of dicts, where each dict is a row from the csv
        for row in reader:
            if row['id'] != '':
              row['pos_hg38'] = row['chrom'] + ':' + row['start_hg38'] + '-' + row['stop_hg38']
              data.append(dict(row))
    #print(data[0])
    return data

def html_table(data):

  titles = ['Position (hg38)', 'Gene', 'Disease', 'Inheritance', 'Motif']
  cols = ['pos_hg38', 'gene', 'disease', 'Inheritance', 'repeatunit_pathogenic_geneorientation']

  table_html = '''
  <table id="table" class="ui sortable celled table" >
    <thead>
      <tr>
      {% for title in titles %}
        <th>{{title}}</th>
      {% endfor %}
      </tr>
      
    </thead>
    <tbody>
      {% for row in data %}
      <tr>
        {% for col in cols %}
        <td>{{row[col]}}</td>
        {% endfor %}
      </tr>
      {% endfor %}
    </tbody>
  </table>
  '''

  final_html = f'''
  <div class="ui middle aligned container">
  <h1 class="ui massive header">STRchive</h1>
  <p>STRchive (Short Tandem Repeat Archive) is a database of STRs associated with disease in humans.</p>
  <p>The full dataset can be viewed on 
    <a class="ui icon button label" href="https://github.com/hdashnow/STRchive">
      <i class="github icon"> GitHub</i>
    </a>
  </p>
  <br>
  
  {table_html}
  </div>
    '''
  my_templ = Template(final_html)
  return my_templ.render(titles=titles, data=data, cols=cols)

def html_surrounds(html_text, template_file):
  with open(template_file) as fh:
    template_html = fh.read()

    return template_html.replace('insert_html_here', html_text)


def main():
    table = html_table(read_csv('../STR-disease-loci.csv'))
    out_html = html_surrounds(table, 'minimal-template.html')
    with open('../table.html', 'w') as fh:
        fh.write(out_html)
    

if __name__ == '__main__':
    main()