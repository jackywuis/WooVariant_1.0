import argparse
import plotly
import plotly.graph_objs as go
from plotly.graph_objs import *


parser = argparse.ArgumentParser(description='WooGraph - A tool to draw heatmap for nt/aa variant vs position. '
                                             'For research team of Dr. Liu in iSynBio, SIAT only. '
                                             'Usage: python WooGraph.py -i file/name.csv '
                                             'Output: file/name.html '
                                             'Copyleft: Jacky Woo from ZHK Research team, iSynBio, SIAT.')
parser.add_argument('-i', '--input_file', type=str, help='Input csv file')
parser.add_argument('-m', '--max', type=float, default=1.0, required=False,
                    help='Maximum threshold of mutation abundance (Default: 1.00 [100%])')
parser.add_argument('-o', '--output_file', type=str, required=False,
                    help='Output html file prefix (Default: input name and file)')

args = parser.parse_args()
if args.output_file is not None:
    outfile = args.output_file + '.html'
else:
    outfile = args.input_file.rsplit('.csv')[0] + '.html'
reader = open(args.input_file)
x = []
y = []
z = []
for line in reader:
    if 'Pos' == line.split('\t')[0]:
        x = line.split('\n')[0].split('\t')[2:]
    else:
        y += [line.split('\t')[0]]
        out = []
        for item in line.split('\n')[0].split('\t')[2:]:
            if item is not '':
                if float(item) > args.max:
                    out += [args.max]
                else:
                    out += [float(item)]
            else:
                out += [item]
        z += [out]
trace = go.Heatmap(x=x, y=y, z=z, colorscale='Jet')
title = outfile.split('/')[1].split('.csv')[0] + ' abundance heatmap'
if '*' in x:
    xtitle = 'Amino acids'
else:
    xtitle = 'Nucleotides'
fig=go.Figure(data=[trace], layout=go.Layout(title=title, xaxis=dict(title=xtitle), yaxis=dict(title='Position')))
plotly.offline.plot(fig, filename=outfile)
