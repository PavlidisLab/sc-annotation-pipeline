#!/bin/python

import graphviz

# Define the DAG
dag = graphviz.Digraph(name="Nextflow_Pipeline", format='png')
dag.attr(rankdir='LR', size='10,6')
dag.attr('node', fontname='Helvetica', fontsize='10')

# Optional: Define color palette for steps
color_map = {
    'start': '#8dd3c7',
    'download': '#80b1d3',
    'setup': '#fdb462',
    'process': '#b3de69',
    'classify': '#fb8072',
    'visualize': '#bebada',
    'finalize': '#d9d9d9'
}

# Define nodes with categories and custom styles
fancy_nodes = {
    'save_params_to_file':        ('start', 'ellipse'),
    'runSetup':                   ('setup', 'box'),
    'DOWNLOAD_STUDIES_SUBWF':     ('download', 'box'),
    'processQuery':              ('process', 'box'),
    'getCensusAdata':            ('process', 'box'),
    'rfClassify':                ('classify', 'box'),
    'getMeta':                   ('process', 'box'),
    'plotQC':                    ('visualize', 'component'),
    'runMultiQC':                ('visualize', 'component'),
    'publishMultiQC':           ('finalize', 'ellipse'),
    'loadResults':              ('finalize', 'box')
}

# Add styled nodes
for node, (group, shape) in fancy_nodes.items():
    dag.node(node, style='filled,setlinewidth(1)', fillcolor=color_map[group], shape=shape)

# Define edges
edges = [
    ('DOWNLOAD_STUDIES_SUBWF', 'processQuery'),
    ('runSetup', 'processQuery'),
    ('processQuery', 'rfClassify'),
    ('getCensusAdata', 'rfClassify'),
    ('rfClassify', 'loadResults'),
    ('rfClassify', 'plotQC'),
    ('processQuery', 'plotQC'),
    ('getMeta', 'plotQC'),
    ('plotQC', 'runMultiQC'),
    ('runMultiQC', 'publishMultiQC'),
    ('save_params_to_file', 'runSetup')  # Optional start node
]

# Add edges
for src, dst in edges:
    dag.edge(src, dst)

# Render the DAG
dag.render('nextflow_pipeline_dag', cleanup=False)
