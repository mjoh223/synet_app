import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_cytoscape as cyto
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import csv
import re
import random
import json
import plotly.express as px
import socket
cyto.load_extra_layouts()
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server
file_name = '/app/net_files/040920_1.tsv'
network_file = '/app/net_files/040920_1_network.tsv'
node_attr_file = '/app/net_files/040920_1_network_node_attr.tsv'
locus_width = 1500
h = .5
s = 1
#Purpose: Parse the locus attributes file to be drawn on plotly
#input: filename
#output: returns dictionary with assembly/organism/marker as keys and list of orfs attributes as values
def file_parse(fname):
    loci_dict = {}
    with open(fname) as tsvfile:
        tsvreader = csv.reader(tsvfile, delimiter="\t")
        for line in tsvreader:
            if line[0] != '\n':
                anchor = False
                if line[3] == line[4]:
                    anchor = True
                key = '|'.join(line[0:4])
                if key in loci_dict.keys():
                    loci_dict[key].append( line[3:]+[anchor] )
                else:
                    loci_dict[key] =  [ line[3:]+[anchor] ,]
    tsvfile.close()
    return loci_dict
#Purpose: align locus around orf of interest (anchor) and draw in plotly
#input: filtered dictionary, height for each orf to be drawn, v is a scale for
#the spacing between loci, total locus width, and anchor. Accepts in dash core componenet checklist value as query.
#output: returns a figure
def organize_map(filt_dict, h, locus_width, anchor, neighbor_nodes):
    fig = go.Figure()
    filt_arch_set = dict()
    for i, [assembly, locus] in enumerate(filt_dict.items()):
        #Get anchor start
        for orf in locus:
            clan = orf[8]
            if clan == anchor:
                anchor_begin = int(orf[2])
                strand_of_anchor = orf[4]
                shift_width = int(orf[3])-int(orf[2]) #only needed for negative anchor strand
                updated_locus = {}
                orfs_in_locus = []
                key = ''
                for orf in locus:
                    key+=orf[8]
                    #Fix all orfs based on anchor start
                    updated_coord = []
                    if strand_of_anchor == '+':
                        for orf in locus:
                            start, stop = int(orf[2])-anchor_begin, int(orf[3])-anchor_begin
                            updated_coord.append([start]+[stop]+orf[4:])
                    else:
                        for orf in locus:
                            start, stop = int(orf[2])-anchor_begin, int(orf[3])-anchor_begin
                            updated_coord.append([start]+[stop]+orf[4:])
                        updated_coord_flip = []
                        for x in updated_coord:
                            if x[2] == "+":
                                strand = '-'
                            else:
                                strand = '+'
                            updated_coord_flip.append([-x[1]+shift_width]+[-x[0]+shift_width]+[strand]+x[3:])
                            updated_coord = updated_coord_flip
                orfs_in_locus.append(updated_coord)
                if updated_locus.get(assembly):
                    updated_locus[assembly].append(orfs_in_locus)
                else:
                    updated_locus[assembly] = [orfs_in_locus,]
                if filt_arch_set.get(key):
                    filt_arch_set[key].append(updated_locus)
                else:
                    filt_arch_set[key] = [updated_locus,]
    #testme = filt_arch_set['WP_031668916.1WP_039194677.1WP_060869021.1YP_009622344.1WP_009928181.1WP_060869024.1WP_031667948.1WP_047933338.1']
    #print(filt_arch_set)
    plot_dict = {}
    for i, k in enumerate(sorted(filt_arch_set, key=lambda k: len(filt_arch_set[k]), reverse=True)):
        group_locus = filt_arch_set[k]
        group_size = len(group_locus)
        orf_dict={}
        for assembly in group_locus: #assembly is a dict
            for assembly_name, locus in assembly.items():
                locus = [item for sublist in locus for item in sublist][0]
                for j, orf in enumerate(locus):
                    product = orf[5].replace('[','').replace(']','').replace("'", '')
                    clan = orf[6]
                    element = orf[-2]
                    if not orf_dict.get(j):
                        orf_dict[j] = [set(),set(),set(),set()]
                    orf_dict[j][0].add(clan)
                    orf_dict[j][1].add(product)
                    orf_dict[j][2].add(element)
                    orf_dict[j][3].add(assembly_name)
        plot_dict[k] = [orf_dict, group_size, locus]
    n = len(plot_dict)
    fig.update_layout(
        #width = locus_width/10,
        #height = 10+n*h,
        yaxis = dict(
            showgrid = False,
            zeroline = False,
            tickmode = 'array',
            tickvals = list(range(0,n)),
            ticktext = list(range(0,n))
        ),
        xaxis = dict(
            showgrid = False,
            zeroline = False)
        )
    sorted_plot_dict = {k: v for k, v in sorted(plot_dict.items(), key=lambda item: str(item[1][0][0][3]).split('|')[1])}
    print(sorted_plot_dict)
    for i, [k,v] in enumerate(sorted_plot_dict.items()):
        group_key = k
        group_data = v[0]
        group_size = v[1]
        #assembly = group_data[0][3]
        locus = v[2]
        for j, orf in enumerate(locus):
            clan = orf[6]
            start, stop = orf[0], orf[1]
            color = clan_dict[clan]
            line_color = color
            neighbor_ids = [n['data']['id'] for n in neighbor_nodes]
            neighbor_id_name_ = [[n['data']['id'], n['data']['name_']] for n in neighbor_nodes]
            opacity = 1
            product = ''.join(group_data[j][1]).replace("'", "").replace("[","").replace("]","")
            if clan in neighbor_ids:
                line_color = 'black'
                opacity = 1
                product = [x[1] for x in neighbor_id_name_ if clan in x[0]]
                product = ''.join(product).replace("'", "").replace("[","").replace("]","")
            if orf[2] == '+':
                fig.add_trace(go.Scatter(
                                x=(start, start+50, start, stop-50, stop, stop-50, start),
                                y=(i, i+h/2, i+h, i+h, i+h/2, i, i),
                                fill='toself',
                                fillcolor=color,
                                line_color=line_color,
                                hovertemplate='{}<br>'.format(product)+
                                              '{}<br>'.format(str(group_data[j][2]))+
                                              '{}<br>'.format(str(group_size))+
                                              '{}<br>'.format('<br>'.join(list(group_data[j][3]))),
                                opacity=opacity,
                                marker=dict(size=1)))
            else:
                fig.add_trace(go.Scatter(
                                x=(start+50, start, start+50, stop, stop-50, stop, start+50),
                                y=(i, i+h/2, i+h, i+h, i+h/2, i, i),
                                fill='toself',
                                fillcolor=color,
                                line_color=line_color,
                                hovertemplate='{}<br>'.format(product)+
                                              '{}<br>'.format(str(group_data[j][2]))+
                                              '{}<br>'.format(str(group_size))+
                                              '{}<br>'.format('<br>'.join(list(group_data[j][3]))),
                                opacity=opacity,
                                marker=dict(size=1)))
            fig.layout.plot_bgcolor = 'white'
            fig.layout.paper_bgcolor = 'white'
    fig.update_layout(showlegend=False)
    print('<br>'.join(list(group_data[j][3])))
    return fig, n

#Style for network
#purpose:
#input:
#output:
def csv_parse(fname, header):
    return pd.read_csv(fname, sep="\t", header=header)
#purpose:
#input:
#output:
def attr_dict_(df_attr):
    attr_dict = {}
    for i, row in df_attr.iterrows():
        attr_dict[row['cluster_rep']] = row
    return attr_dict
#purpose:
#input:
#output:
def filter_net(df_net, filter):
    df = df_net.query('{} <= weight <= {}'.format(filter[0],filter[1]))
    return df
#purpose:
#input:
#output:
#def plot_network(filtered_net):

df_attr = csv_parse(node_attr_file,header=0)
attr_dict = attr_dict_(df_attr)
df_net = csv_parse(network_file,header=None)
df_net.columns = ['source','target','weight']
max_weight = max(df_net.weight)
cy_edges = []
cy_nodes = []
nodes = set()
neighbor_node_di = {}
neighbor_edge_di = {}
for i, edges in df_net.iterrows():
    source, target, weight = edges
    if attr_dict[source]['name_'] != '[]' and attr_dict[target]['name_'] != '[]':
        s_name_ = ''.join(attr_dict[source]['name_']).replace("'","").replace("]","").replace("[","")
        cy_source = {'data': {'id': source,
                                   'name_': s_name_,
                                   'member_ids': attr_dict[source]['member_ids'],
                                   'locus': attr_dict[source]['locus'],
                                   'assembly': attr_dict[source]['assembly'],
                                   'size': attr_dict[source]['n_members'],
                                   'plasmid': attr_dict[source]['n_plasmids'],
                                   'virus': attr_dict[source]['virus_percent'],
                                   'product': attr_dict[source]['product'],
                                   'iso': attr_dict[source]['avg_iso'],
                                   'phylum': attr_dict[source]['phylum'],
                                   'genus': attr_dict[source]['species'],
                                   'mean_length': attr_dict[source]['mean_length']},
                    }
        t_name_ = ''.join(attr_dict[target]['name_']).replace("'","").replace("]","").replace("[","")
        cy_target = {'data': {'id': target,
                                   'name_': t_name_,
                                   'member_ids': attr_dict[target]['member_ids'],
                                   'locus': attr_dict[target]['locus'],
                                   'assembly': attr_dict[target]['assembly'],
                                   'size': attr_dict[target]['n_members'],
                                   'plasmid': attr_dict[target]['n_plasmids'],
                                   'virus': attr_dict[target]['virus_percent'],
                                   'product': attr_dict[target]['product'],
                                   'iso': attr_dict[target]['avg_iso'],
                                   'phylum': attr_dict[target]['phylum'],
                                   'genus': attr_dict[target]['species'],
                                   'mean_length': attr_dict[target]['mean_length']},
                    }
        cy_edge = {"data": {'id': source+target, "source": source, "target": target, 'weight': weight}, }

        if source not in nodes:
            nodes.add(source)
            cy_nodes.append(cy_source)
        if target not in nodes:
            nodes.add(target)
            cy_nodes.append(cy_target)
       # Process dictionary of following
        if not neighbor_node_di.get(source):
            neighbor_node_di[source] = []
        if not neighbor_node_di.get(target):
            neighbor_node_di[target] = []
        if not neighbor_edge_di.get(source):
            neighbor_edge_di[source] = []
        if not neighbor_edge_di.get(target):
            neighbor_edge_di[target] = []

        neighbor_node_di[source].append(cy_target)
        neighbor_node_di[target].append(cy_source)
        neighbor_edge_di[source].append(cy_edge)
        neighbor_edge_di[target].append(cy_edge)

        cy_edges.append(cy_edge)
default_elements = cy_edges+cy_nodes

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    },
    'width': '50%'
}
default_stylesheet = [
    # Group selectors
    {
        "selector": 'node',
        'style': {
            "opacity": 1,
            'z-index': 9999,
            "label": "data(name_)"
        }
    },
    {
        "selector": 'edge',
        'style': {
            "curve-style": "bezier",
            "opacity": 1,
            'z-index': 5000
            }
    },
       {
        'selector': ':selected',
        "style": {
            "border-width": 2,
            "border-color": "black",
            "border-opacity": 1,
            "background-color": "white",
            "opacity": 1,
            'z-index': 1
        }
    }
]

app.layout = html.Div([
    html.H1('SYNET server'),
    html.Div([
        html.Div([
            html.H4('Gene co-occurence threshold'),
            dcc.RangeSlider(
                id='cytoscape_range',
                min=1,
                max=max_weight,
                step=1,
                value=[1, max_weight],
                marks={
                    1: {'label': '1'},
                    200: {'label': '200'},
                    300: {'label': '300'},
                    400: {'label': '400'},
                    600: {'label': '600'},
                    }),
            html.H4('Selected gene info'),
            html.Pre(id='cytoscape-tapNodeData-json', style=styles['pre'])
        ],className = 'four columns'),

        html.Div([
            html.H4('Association network'),
            html.P('Query genes are colored orange'),
            cyto.Cytoscape(
                    id='cytoscape',
                    elements=default_elements,
                    style={'width': '100%', 'height': '800px'},
                    layout={'name':'dagre'},)
        ],className = 'eight columns')
    ],className="row"),

    html.Div([
        html.H3('orf map'),
        dcc.Graph(id='orf_map',
                  figure=go.Figure())
    ])
])

loci_dict = file_parse(file_name)
clan_keys = set(sum([ [orf[8] for orf in locus] for locus in loci_dict.values()], []))
rgb_values = ['rgb{}'.format(tuple(np.random.choice(range(256), size=3))) for i in range(len(clan_keys))]
clan_dict = dict(zip(clan_keys,rgb_values))

@app.callback(Output('cytoscape', 'elements'),
              [Input('cytoscape_range', 'value'),
               Input('cytoscape', 'tapNodeData')])
def plot_network(range, nodeData):
    if not nodeData:
        return default_elements
    neighbor_nodes = neighbor_node_di.get(nodeData['id'])
    filtered_nodes = set()
    updated_edges = []
    updated_net_stylesheet = []
    for edge in cy_edges:
        if range[0] <= edge['data']['weight'] <= range[1]: #for speed to try to find a filter option
            source, target = edge['data']['source'], edge['data']['target']
            filtered_nodes.add(source)
            filtered_nodes.add(target)
            updated_edges.append(edge)
    updated_nodes = []
    updated_net_stylesheet = []
    for node in cy_nodes:
        if node['data']['id'] in filtered_nodes:
            node['classes'] = ''
            updated_nodes.append(node)
            if node['data']['id'] in [n['data']['id'] for n in neighbor_nodes]:
                node['classes'] = 'neighborNode'
                updated_nodes.append(node)
                default_stylesheet.append({
                        'selector': '[id = "{}"]'.format(node['data']['id']),
                        'style': {
                            'background-color': '{}'.format(clan_dict[node['data']['id']])
                        }})
    return updated_nodes+updated_edges

@app.callback(Output('cytoscape', 'stylesheet'),
              [Input('cytoscape', 'elements')])
def plot_network(elements):
    return default_stylesheet

@app.callback([Output('cytoscape-tapNodeData-json', 'children'),
              Output('orf_map', 'figure'),
              Output('orf_map', 'style')],
              [Input('cytoscape', 'tapNodeData')])
def plot_orfmap(nodeData):
    if not nodeData:
        return json.dumps('select node', indent=2), go.Figure(), {}
    anchor = nodeData['id']
    anchor_name = nodeData['name_'].replace('[', '').replace(']', '').replace("'", '')
    neighbor_nodes = neighbor_node_di.get(nodeData['id'])
    filt_loci_dict = dict(filter(lambda elem: anchor in [x[8] for x in elem[1]], loci_dict.items()))
    fig, n = organize_map(filt_loci_dict, h, locus_width, anchor, neighbor_nodes)
    style = {
        'height': '{}vh'.format(n*5),
        'width': '100%',
    }
    return json.dumps(nodeData, indent=2), fig, style

if __name__ == '__main__':
    app.run_server(debug=True)
