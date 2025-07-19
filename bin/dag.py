#!/bin/python

import graphviz

# Define the DAG
dag = graphviz.Digraph(name="Nextflow_Pipeline", format='png')
dag.attr(rankdir='LR', size='12,8')
dag.attr('node', fontname='Helvetica', fontsize='10')

# Color groups
color_map = {
    'start': '#8dd3c7',
    'download': '#80b1d3',
    'setup': '#fdb462',
    'process': '#b3de69',
    'classify': '#fb8072',
    'visualize': '#bebada',
    'finalize': '#d9d9d9'
}

# Define nodes (merged logic for both per-sample and combined paths)
fancy_nodes = {
    'save_params_to_file':       ('start', 'ellipse'),
    'runSetup':                  ('setup', 'box'),
    'DOWNLOAD_STUDIES_SUBWF':    ('download', 'box'),
    'processQuery':              ('process', 'box'),  # abstracted name
    'getCensusAdata':            ('process', 'box'),
    'getMeta':                   ('process', 'box'),
    'rfClassify':                ('classify', 'box'),
    'loadCTA':                   ('finalize', 'box'),
    'processQC':                 ('visualize', 'component'),
    'loadCLC':                   ('finalize', 'box'),
    'runMultiQC':                ('visualize', 'component'),
    'publishMultiQC':            ('finalize', 'ellipse')
}

# Add styled nodes
for node, (group, shape) in fancy_nodes.items():
    dag.node(node, style='filled,setlinewidth(1)', fillcolor=color_map[group], shape=shape)

# Define edges (omit combineCTA, combineCLC, combineQC)
edges = [
    ('save_params_to_file', 'runSetup'),
    ('DOWNLOAD_STUDIES_SUBWF', 'processQuery'),
    ('runSetup', 'processQuery'),
    ('processQuery', 'rfClassify'),
    ('getCensusAdata', 'rfClassify'),
    ('rfClassify', 'loadCTA'),
    ('getMeta', 'processQC'),
    ('rfClassify', 'processQC'),
    ('processQuery', 'processQC'),
    ('processQC', 'loadCLC'),
    ('processQC', 'runMultiQC'),
    ('runMultiQC', 'publishMultiQC')
]

# Add edges to the DAG
for src, dst in edges:
    dag.edge(src, dst)

# Render the DAG
dag.render('nextflow_pipeline_dag', cleanup=False)
