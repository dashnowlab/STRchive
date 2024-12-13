import plotly.graph_objects as go
import json
import jsbeautifier
import argparse
import doctest

def parse_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Generate plotly data (JSON) for pathogenic size ranges")
    parser.add_argument("input", help="JSON file of STRchive loci database")
    parser.add_argument("output", help="Output JSON file for plotly library")
    return parser.parse_args()


def get_minimum_value(info):
    """
    Get pathogenic min for sorting
    """
    values = []
    try: values.append(int(info['pathogenic_min']))
    except TypeError: pass
    try: values.append(int(info['pathogenic_max']))
    except TypeError: pass

    #values = list(filter(lambda x: x>0, values))
    values = [x * info['motif_len'] for x in values if x > 0]
    if len(values) == 0: return 0

    return min(values)


def valid_range(range_min, range_max, scale):
    """
    Ensure the range is valid, multiply by scale, and return a list of two integers

    >>> valid_range(None, None, 1)
    [None, None]
    >>> valid_range(None, 10, 1)
    [10, 10]
    >>> valid_range(10, None, 1)
    [10, 10]
    >>> valid_range(10, 20, 2)
    [20, 40]
    >>> valid_range(0, 10, 1)
    [1, 10]
    """
    new_range = [range_min, range_max]
    if range_min is None and range_max is None: return [0, 0]
    if range_min is None: new_range = [int(range_max), int(range_max)]
    if range_max is None: new_range =  [int(range_min), int(range_min)]
    for i in range((len(new_range))):
        new_range[i] = int(new_range[i]) * scale
        if new_range[i] < 1: new_range[i] = 1
    return new_range


def build_JSON(args):
    json_txt = ""
    with open(args.input) as fh:
        for line in fh: json_txt += line.strip()
    strchive_info = json.loads(json_txt)

    strchive_info = list(filter(lambda x: "conflicting_evidence" not in x["locus_tags"], strchive_info))
    strchive_info = sorted(strchive_info, key=get_minimum_value, reverse=True)

    diseases = []
    bar_data = {'Benign': {'base': [], 'x': []}, 'Intermediate': {'base': [], 'x': []}, 'Pathogenic': {'base': [], 'x': []}}
    scatter_data = {'Benign': {'y': [], 'x': []}, 'Intermediate': {'y': [], 'x': []}, 'Pathogenic': {'y': [], 'x': []}}
    line_data = {'y': [], 'xmin': [], 'xmax': []}
    colors = {'Benign': '#00afbb', 'Intermediate': '#e7b800', 'Pathogenic': '#fc4e07'}
    for info in strchive_info:
        diseases.append(info['disease_id'])
        motif_len = info['motif_len']
        benign_min, benign_max = valid_range(info['benign_min'], info['benign_max'], motif_len)
        pathogenic_min, pathogenic_max = valid_range(info['pathogenic_min'], info['pathogenic_max'], motif_len)
        intermediate_min, intermediate_max = valid_range(info['intermediate_min'], info['intermediate_max'], motif_len)
        
        values = [[benign_min, benign_max], [intermediate_min, intermediate_max], [pathogenic_min, pathogenic_max]]
        values = sorted(values, key=lambda x: x[0])
        if values[0][1] != 0 and values[1][0] > values[0][1]:
                line_data['y'].append(info['disease_id'])
                line_data['xmin'].append(values[0][1])
                line_data['xmax'].append(values[1][0])
        if values[1][1] != 0 and values[2][0] > values[1][1]:
                line_data['y'].append(info['disease_id'])
                line_data['xmin'].append(values[1][1])
                line_data['xmax'].append(values[2][0])
        
        if benign_min==benign_max and benign_min>0: scatter_data['Benign']['x'].append(benign_min); scatter_data['Benign']['y'].append(info['disease_id'])
        if intermediate_min==intermediate_max and intermediate_min>0: scatter_data['Intermediate']['x'].append(intermediate_min); scatter_data['Intermediate']['y'].append(info['disease_id'])
        if pathogenic_min==pathogenic_max and pathogenic_min>0: scatter_data['Pathogenic']['x'].append(pathogenic_min); scatter_data['Pathogenic']['y'].append(info['disease_id'])
        
        # print(info['disease_id'], benign_min, benign_max, intermediate_min, intermediate_max, pathogenic_min, pathogenic_max)
        bar_data['Benign']['base'].append(benign_min); bar_data['Benign']['x'].append(benign_max-benign_min)
        bar_data['Intermediate']['base'].append(intermediate_min); bar_data['Intermediate']['x'].append(intermediate_max-intermediate_min)
        bar_data['Pathogenic']['base'].append(pathogenic_min); bar_data['Pathogenic']['x'].append(pathogenic_max-pathogenic_min)

    fig = go.Figure(layout=go.Layout())
    fig.update_layout( yaxis=dict( tickmode="array", tickvals=diseases))
    fig = go.Figure(layout=go.Layout()) #width = 600, height=len(yvalues_line)*10)

    fig.update_layout(yaxis=dict( tickmode="array", tickvals=diseases))

    for key in bar_data:
        fig.add_trace(go.Bar(
            y=diseases,
            x=bar_data[key]['x'],
            base=bar_data[key]['base'],
            name= key,
            legendgroup=key,
            orientation = 'h',
            marker=dict(color=colors[key]),
            width=0.6,
            offset=-0.3,
            hovertemplate="Disease: %{y} <br> Range: %{base} - %{x} bp"
        ))

    for i,key in enumerate(line_data['y']):
        fig.add_trace(go.Scatter(
            y=[key, key],
            x=[line_data['xmin'][i], line_data['xmax'][i]],
            mode='lines',
            line=dict(dash='dot', color='#aeaeae'),
            showlegend=False,
            hoverinfo='skip'
        ))

    for key in scatter_data:
        fig.add_trace(go.Scatter(
            y=scatter_data[key]['y'],
            x=scatter_data[key]['x'],
            mode='markers',
            name=key,
            legendgroup=key,
            marker=dict(color=colors[key], size=8),
            hovertemplate="Disease: %{y} <br> Value: %{x} bp",
            showlegend=False
        ))

    fig.update_layout(autosize=False, width=800, height=len(diseases)*20,
                    margin=dict(l=10, r=10, t=10, b=10))
    fig.update_layout(template="plotly_white", legend_title='Allele size')
    fig.update_xaxes(title_text="Allele size in base pairs", type="log")
    fig.update_yaxes(title_text="Disease")
    #fig.show()
    fig_json = fig.to_json()

    with open(args.output, "w") as file:
        options = jsbeautifier.default_options()
        options.indent_size = 2
        options.brace_style="expand"
        file.write(jsbeautifier.beautify(fig_json, options))

if __name__=="__main__":
    args = parse_args()
    build_JSON(args)