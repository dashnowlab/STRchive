import plotly.graph_objects as go
import json
import jsbeautifier
import argparse

def parse_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Generate plotly data (JSON) for pathogenic size ranges")
    parser.add_argument("input", help="JSON file of STRchive loci database")
    parser.add_argument("output", help="Output JSON file for plotly library")
    return parser.parse_args()



def valid_range(amin, amax):
    if amin is None and amax is None: return [0,0]
    if amin is None: return [amax,amax+1]
    if amax is None: return [amin, amin+1]
    return [amin, amax]

def get_min_age(info, i = 0):
    """ Get the minimum age of onset from the info dictionary.
    If i is provided, it will return the i-th smallest value from the list of ages.
    For example, if i = 0, it will return the minimum age of onset.
    If i = 1, it will return the second smallest age of onset, and so on.
    If no ages are available, it will return 0.
    """
    values = []
    try: values.append(int(info['age_onset_min']))
    except TypeError: pass
    try: values.append(int(info['age_onset_max']))
    except TypeError: pass
    try: values.append(int(info['typ_age_onset_min']))
    except TypeError: pass
    try: values.append(int(info['typ_age_onset_max']))
    except TypeError: pass

    # values = list(filter(lambda x: x>0, values))
    if len(values) == 0: return 0

    # ensure that the values are unique and sorted
    values = sorted(set(values))

    try:
        return values[i]
    except IndexError:
        return 0


def build_JSON(args):
    json_txt = ""
    with open(args.input) as fh:
        for line in fh: json_txt += line.strip()
    strchive_info = json.loads(json_txt)

    strchive_info = list(filter(lambda x: "conflicting_evidence" not in x["locus_tags"], strchive_info))
    # sort by minimum age of onset, then by typical then by maximum
    # sort by multiple outputs of get_min_age
    strchive_info = sorted(strchive_info, key=lambda x: (get_min_age(x, 0), get_min_age(x, 1), get_min_age(x, 2)), reverse=True)

    plot_data = {}
    diseases = []
    colors = {'AD': '#33a02c', 'AR': '#b2df8a', 'XD': '#1f78b4', 'XR': '#a6cee3'}
    for info in strchive_info:
        disease = info['disease_id']
        if len(info['inheritance']) == 0: continue
        
        onset_min, onset_max = valid_range(info['age_onset_min'], info['age_onset_max'])
        typonset_min, typonset_max = valid_range(info['typ_age_onset_min'], info['typ_age_onset_max'])
        if max([onset_min, onset_max, typonset_min, typonset_max]) == 0: continue
        
        diseases.append(disease)
        inheritance = info['inheritance'][0]
        for x in colors:
            if x == inheritance:
                if inheritance in plot_data:
                    plot_data[inheritance]["Onset"]["y"].append(disease)
                    plot_data[inheritance]["Onset"]["base"].append(onset_min)
                    plot_data[inheritance]["Onset"]["x"].append(onset_max-onset_min)
                    plot_data[inheritance]["TypOnset"]["y"].append(disease)
                    plot_data[inheritance]["TypOnset"]["base"].append(typonset_min)
                    plot_data[inheritance]["TypOnset"]["x"].append(typonset_max-typonset_min)
                else:
                    plot_data[inheritance] = {"Onset": {'base': [onset_min], 'x': [onset_max-onset_min], 'y': [disease], 'color': [inheritance]}, 
                                            "TypOnset": {'base': [typonset_min], 'x': [typonset_max-typonset_min], 'y': [disease], 'color': inheritance}}
            else:
                if x in plot_data:
                    plot_data[x]["Onset"]["y"].append(disease)
                    plot_data[x]["Onset"]["base"].append(0)
                    plot_data[x]["Onset"]["x"].append(0)
                    plot_data[x]["TypOnset"]["y"].append(disease)
                    plot_data[x]["TypOnset"]["base"].append(0)
                    plot_data[x]["TypOnset"]["x"].append(0)
                else:
                    plot_data[x] = {"Onset": {'base': [0], 'x': [0], 'y': [disease]}, 
                                    "TypOnset": {'base': [0], 'x': [0], 'y': [disease]}}

    fig = go.Figure(layout=go.Layout()) #width = 600, height=len(yvalues_line)*10)

    fig.update_layout( yaxis=dict( tickmode="array", tickvals=diseases))

    for key in plot_data:
        fig.add_trace(go.Bar(
            y=plot_data[key]['Onset']['y'],
            x=plot_data[key]['Onset']['x'],
            base=plot_data[key]['Onset']['base'],
            name= key,
            legendgroup=key,
            orientation = 'h',
            marker=dict(color=colors[key]),
            width=0.2,
            offset=-0.1,
            hovertemplate="Disease: %{y} <br> Range onset: %{base} - %{x} years"
        ))
    for key in plot_data:
        fig.add_trace(go.Bar(
            y=plot_data[key]['TypOnset']['y'],
            x=plot_data[key]['TypOnset']['x'],
            base=plot_data[key]['TypOnset']['base'],
            name= key,
            legendgroup=key,
            orientation = 'h',
            marker=dict(color=colors[key]),
            width=0.6,
            offset=-0.3,
            hovertemplate="Disease: %{y} <br> Typical onset: %{base} - %{x} years",
            showlegend=False
        ))

    fig.add_trace(go.Scatter(
        y=[diseases[0], diseases[-1]],
        x=[18, 18],
        mode="lines",
        line=dict(dash='dot', color='#aeaeae'),
        showlegend=False,
        hoverinfo='text',
        hovertemplate="Age 18",
        name="Age 18"
    ))

    fig.update_layout(height=max(len(diseases)*20, 1600),
                    template="plotly_white",
                    legend_title="Inheritance",
                    margin=dict(l=10, r=10, t=10, b=10),
                    ) 
    fig.update_xaxes(title_text="Age of onset (years)")
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