from bokeh.plotting import figure, output_notebook, show, ColumnDataSource, output_file
from bokeh.models.tools import HoverTool, BoxSelectTool, BoxZoomTool, PanTool, WheelZoomTool, SaveTool, ResetTool
from bokeh.models.tickers import FixedTicker
from bokeh.models.widgets import CheckboxGroup, RangeSlider, Tabs, TextInput, Button
from bokeh.layouts import column, row, WidgetBox
from bokeh.models import Panel
from bokeh.io import show, curdoc
import os
import pandas as pd
import numpy as np
from math import pi
from copy import deepcopy

os.chdir('/Users/Kimberly/Documents/projects/ecoli_promoters/endo/bokeh_viz')


global frag_plus_pileup_norm
global frag_minus_pileup_norm
global colors
global genes

def norm_pileup(pileup):
    # exclude data below median expression. However, for graphing we need a continuous line
    # for all positions. We will set all data below median to the median value, and then scale
    # data so minimum value is 1. On log-scale, log(1) = 0 so any data <= median will be set to
    # zero on plot
    median = np.median(pileup.expression)
    norm_pileup = deepcopy(pileup)
    norm_pileup.expression = np.where(pileup.expression <= median,
                                     median, pileup.expression)
    scalar = 1.0/median
    norm_pileup.expression = norm_pileup.expression * scalar
    # make sure everything is equal to 1, due to decimal multiplication some may be 0.9999
    norm_pileup.expression = np.where(norm_pileup.expression < 1,
                                     1, norm_pileup.expression)
    return norm_pileup


def make_pileup_dataset(conditions, start, end, colors, log_transform=True):
    
    # subset relevant condition, conditions must be list
    plus_pileup_subset = frag_plus_pileup_norm[(frag_plus_pileup_norm.condition.isin(conditions)) &
                                         (frag_plus_pileup_norm.position.isin(range(start,end)))]
    
    minus_pileup_subset = frag_minus_pileup_norm[(frag_minus_pileup_norm.condition.isin(conditions)) &
                                           (frag_minus_pileup_norm.position.isin(range(start,end)))]
    
    plus_pileup_subset = plus_pileup_subset.merge(colors[colors.strand == '+'], on='condition', how='left')
    minus_pileup_subset = minus_pileup_subset.merge(colors[colors.strand == '-'], on='condition', how='left')
    
    # patch accepts list of lists, make one list for each combination of strand and condition
    # x-axis position will be same for every condition and strand, define coordinates for edges at the top (pileup line) and the 
    position_unique = plus_pileup_subset.sort_values('position').position.unique()
    position = np.hstack((position_unique, position_unique[::-1]))
    # bottom edge along zero (positions in reverse to create closed shape)
    bottom_edge = np.zeros(len(position))
    # take negative of reverse coverage to display on graph
    counts = []
    position_all = []
    condition_colors = []
    for condition in conditions:
        # plus strand
        position_all.append(position)
        if log_transform:
            counts.append(np.hstack((np.log(plus_pileup_subset.expression[plus_pileup_subset.condition == condition]),
                               bottom_edge)))
        else:
            counts.append(np.hstack((plus_pileup_subset.expression[plus_pileup_subset.condition == condition],
                                   bottom_edge)))
        condition_colors.append(colors.color[(colors.condition == condition) & (colors.strand == '+')].to_string(index=False))
        
        # minus strand
        position_all.append(position)
        if log_transform:
            counts.append(np.hstack((-np.log(minus_pileup_subset.expression[minus_pileup_subset.condition == condition]),
                               bottom_edge)))
        else:
            counts.append(np.hstack((-minus_pileup_subset.expression[minus_pileup_subset.condition == condition],
                                   bottom_edge)))
        condition_colors.append(colors.color[(colors.condition == condition) & (colors.strand == '-')].to_string(index=False))
    
    patches_src = pd.DataFrame({'position' : position_all, 
                               'count' : counts,
                               'color' : condition_colors})    
    
    return ColumnDataSource(patches_src)


def make_region_genes(genes, start, end):
    # grab genes within coordinates, with buffer of 100bp
    start_buffer = start - 100
    end_buffer = end + 100
    region_genes = genes[(genes.start >= start_buffer) & (genes.end <= end_buffer)]

    # center is midpoint of gene
    gene_center = region_genes.start + abs(region_genes.end - region_genes.start)/2.0
    gene_center = gene_center.tolist()
    gene_color = '#e2b306'
    gene_width = region_genes.end - region_genes.start
    gene_width = gene_width.tolist()

    src_gene = ColumnDataSource(data=dict(
    gene_center = gene_center,
    # y-center is on x-axis, i.e y = 0
    gene_center_y = np.zeros(len(region_genes)),
    gene_color = [gene_color] * len(region_genes),
    gene_width = gene_width,
    # set center of triangle to start or end depending on strand
    tri_x = [region_genes.end.iloc[i] if region_genes.strand.iloc[i] == '+' else region_genes.start.iloc[i] \
           for i in range(len(region_genes))],
    angle = [-90 if strand == '+' else 90 for strand in region_genes.strand],
    gene_name = region_genes.name.tolist()))
    
    return src_gene



def make_region_plot(src, src_gene):
    '''
    Construct pileup plot based on src
    '''

    # output_file(html_output)

    # flatten list of lists in count column of src, find max value of absolute expression
    count_range = max(map(abs, [x for count in src.data['count'] for x in count]))
    
    # draw blank figure of correct size with tools
    p = figure(y_range=(-count_range, count_range), plot_width=900, plot_height=700, 
               tools=[BoxSelectTool(), BoxZoomTool(), PanTool(), WheelZoomTool(), 
                      SaveTool(), ResetTool()],
               toolbar_location='above')
    
    # format axis and colors
    p.xaxis.axis_label = 'position'
    p.xaxis.major_label_orientation = pi/4
    p.xaxis[0].formatter.use_scientific = False
    # p.xaxis[0].ticker=FixedTicker(ticks=range(start, end, 100))
    p.yaxis.axis_label = 'log normalized expression (RNA/DNA)'
        
    p.patches(source=src, xs='position', ys='count', fill_color='color', line_color=None, alpha=0.50)

    # plot genes
    p.rect(x='gene_center', y='gene_center_y', width='gene_width',
                 color='gene_color', height=10, height_units='screen', 
                  alpha=0.75, source=src_gene)
    p.triangle(x='tri_x', y=0, size=20, angle='angle', angle_units='deg',
                     fill_color='gene_color', line_color=None, alpha=0.75, source=src_gene)
    p.text(x='gene_center', y='gene_center_y', text='gene_name', text_color='black',
          text_align='center', text_baseline='middle', text_font_size='10pt', source=src_gene)
        
    return p
        

def update_region_plot(attr, old, new):
    '''
    Update plot based on interactive widgets
    '''
    # get list of conditions for graph
    conditions_to_plot = [condition_selection.labels[i] for i in condition_selection.active]
    
    # # get selected position, value for range slider is tuple (start, end)
    # position_start = position_selection.value[0]
    # position_end = position_selection.value[1]

    # get selected position from input text box
    # position = position_selection.value
    # # split comma separated string, remove any leading or trialing whitespaces,
    # # convert to int
    # position_start, position_end = map(int, map(unicode.strip, position.split(',')))

    position_start = int(start.value)
    position_end = int(end.value)
    
    # make new subset based on selected conditions
    new_src = make_pileup_dataset(conditions_to_plot, position_start, position_end, colors)
    new_src_gene = make_region_genes(genes, position_start, position_end)
    
    # update source used in patch glyphs
    src.data.update(new_src.data)
    # update source used in gene rectangle and triangle glyphs
    src_gene.data.update(new_src_gene.data)


def update_axis():
    count_range = max(map(abs, [x for count in src.data['count'] for x in count]))
    p.y_range.start = -count_range
    p.y_range.end = count_range
    
############################ Read in data ######################################

# endo TSS expression
endo_tss_lb = pd.read_table('../processed_data/endo_tss/lb/rLP5_Endo2_expression.txt',
                           sep=' ')

# LB genomic shearing fragment pileup
frag_lb_plus = pd.read_table('../processed_data/frag_peak_calling/lb/plus_frag_pileup.wig',
                            sep='\t', skiprows=1, names=['position', 'expression'])

frag_lb_minus = pd.read_table('../processed_data/frag_peak_calling/lb/minus_frag_pileup.wig',
                            sep='\t', skiprows=1, names=['position', 'expression'])

# M9 minimal genomic shearing fragment pileup
frag_m9_plus = pd.read_table('../processed_data/frag_peak_calling/m9/plus_frag_pileup_M9.wig',
                            sep='\t', skiprows=1, names=['position', 'expression'])

frag_m9_minus = pd.read_table('../processed_data/frag_peak_calling/m9/minus_frag_pileup_M9.wig',
                            sep='\t', skiprows=1, names=['position', 'expression'])


# format frag pileup data
frag_plus_pileup = frag_lb_plus.merge(frag_m9_plus, on='position', how='outer', suffixes=['_lb', '_m9'])
frag_plus_pileup.columns = ['position', 'LB', 'M9']
frag_minus_pileup = frag_lb_minus.merge(frag_m9_minus, on='position', how='outer', suffixes=['_lb', '_m9'])
frag_minus_pileup.columns = ['position', 'LB', 'M9']

frag_plus_pileup = pd.melt(frag_plus_pileup, id_vars=['position'], value_vars=['LB', 'M9'],
             var_name='condition', value_name='expression')
frag_minus_pileup = pd.melt(frag_minus_pileup, id_vars=['position'], value_vars=['LB', 'M9'],
             var_name='condition', value_name='expression')

# normalize and scale
frag_plus_pileup_norm = norm_pileup(frag_plus_pileup)
frag_minus_pileup_norm = norm_pileup(frag_minus_pileup)


# read in gene annotation
genes = pd.read_table('U00096.2_genes_clean.bed', sep = '\t', header = None,
                      names=['chrom', 'start', 'end', 'name', 'score', 'strand', 
                            'thick_start', 'thick_end', 'item_rgb', 'block_count',
                            'block_sizes', 'block_start'])
# drop unnecessary columns and simplify name
genes = genes[['start', 'end', 'name', 'strand']]
genes.name = [x.split(':')[1] for x in genes.name.tolist()]

# create colors for each condition and strand
colors = pd.DataFrame({'condition' : ['LB', 'LB', 'M9', 'M9'],
                      'strand' : ['+', '-', '+', '-'],
                      'color' : ['#8daec9', '#edac80', '#528ecb', '#ef8137']})


############################ Define widgets ####################################

condition_selection = CheckboxGroup(labels=['LB', 'M9'], active = [0])
# link change in selected buttons to update function
condition_selection.on_change('active', update_region_plot)

# # RangeSlider to change positions of region
# position_selection = RangeSlider(start = 2238500 - 300, end = 2238500 + 300,
#                           value = (2238500 - 300, 2238500 + 300),
#                           step = 100, title = 'genomic position')

# TextInput to define start and end position
# position_selection = TextInput(value='362455, 365529', title='genome U00096.2 start position , end position (Comma-separated)')
start = TextInput(value='362455', title='genome U00096.2 start position')
end = TextInput(value='365529', title='genome end position')

# update plot when value is changed
# position_selection.on_change('value', update_region_plot)
start.on_change('value', update_region_plot)
end.on_change('value', update_region_plot)

# update axis range on click
axis_button = Button(label='Update Axis')
axis_button.on_click(update_axis)


############################## Initialize ######################################

# find initial conditions and position
initial_conditions = [condition_selection.labels[i] for i in condition_selection.active]
# split comma separated string, remove any leading or trialing whitespaces, convert to int
# initial_start, initial_end = map(int, map(unicode.strip, position_selection.value.split(',')))
initial_start = int(start.value)
initial_end = int(end.value)

src = make_pileup_dataset(initial_conditions,
    initial_start,
    initial_end,
    colors)

src_gene = make_region_genes(genes, initial_start, initial_end)

p = make_region_plot(src, src_gene)



################################# Layout ######################################

# combine all elements onto one page by defining layout

# put controls in single element
# controls = WidgetBox(condition_selection, position_selection)
controls = WidgetBox(condition_selection, start, end, axis_button)


# create row layout
layout = row(controls, p)

# make tab with layout
tab = Panel(child=layout, title = 'Genomic fragment shearing pileup')
tabs = Tabs(tabs=[tab])

# add to current document (displays plot)
curdoc().add_root(tabs)

