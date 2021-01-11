### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python selfie_heatmap.py
					--exp <EXPRESSION_FILE>
					--genes <GENE_INFO_FILE>
					--samples <SAMPLE_INFO_FILE>
					--out <OUTPUT_FILE>
					"""

import matplotlib.pyplot as plt
import os, sys, re, math
import numpy as np
import seaborn as sns
from pandas import DataFrame

# --- end of imports --- #

def load_data( data_file ):
	
	data = {}
	with open( data_file, "r" ) as f:
		samples = f.readline().strip().split('\t')[1:]
		vals_per_sample = {}
		for sample in samples:
			vals_per_sample.update( { sample: [] } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			exp_per_gene = {}
			for idx, val in enumerate( parts[1:] ):
				exp_per_gene.update( { samples[ idx ]: float( val ) } )
				vals_per_sample[ samples[ idx ] ].append( float( val ) )
			data.update( { parts[0]: exp_per_gene } )
			line = f.readline()
	return data, samples, vals_per_sample


def calculate_z_scores( values ):
	"""! @brief calculate z-scores per gene """
	
	if sum( values ) > 0:
		avg = np.mean( values )
		std = np.std( values )
		zscores = []
		for val in values:
			if val < avg:
				zscores.append( -1 * math.log( ( abs( val - avg ) / std ), 2 ) )
			else:
				zscores.append( math.log( ( ( val - avg ) / std ), 2 ) )
		return zscores
	else:
		return values


def construct_data_output_file( data, candidate_genes, candidate_samples, outputfile, gene_name_mapping_table, sample_order ):
	"""! @brief write expression values of all candidate genes into output file """
	
	datamatrix = []
	genes = []
	with open( outputfile, "w" ) as out:
		new_line = [ "gene" ] + sample_order	#sorted( candidate_samples.keys() )
		tissues = new_line[1:]
		out.write( "\t".join( new_line ) + '\n' )
		for gene in candidate_genes:
			try:
				new_line = [ gene_name_mapping_table[ gene ] ]
			except KeyError:
				new_line = [ gene ]
				
			for tissue in tissues:
				tmp_value = []
				for sample in candidate_samples[ tissue ]:
					try:
						tmp_value.append( data[ gene ][ sample ] )
					except KeyError:
						pass
				if len( tmp_value ) > 0:
					new_line.append( sum( tmp_value ) / len( tmp_value ) )
				else:
					new_line.append( 0 )
			zscores = calculate_z_scores( new_line[1:] )
			out.write( "\t".join( [ new_line[0] ] + map( str, zscores ) ) + '\n' )
			datamatrix.append( zscores )
			try:
				genes.append( gene_name_mapping_table[ gene ] )
			except KeyError:
				genes.append( gene )
	return genes, tissues, datamatrix


def construct_heatmap( datamatrix, genes, tissues, heatmap_file ):
	"""! @brief construct heatmap from given data matrix """
	
	print "number of genes for heatmap construction: " + str( len( genes ) )
	print "number of samples for heatmap construction: " + str( len( tissues ) )
	
	df = DataFrame( datamatrix, index=genes[::-1], columns=tissues)
	
	fig, ax = plt.subplots( figsize=(7,7) )
	
	x = sns.heatmap( 	df, ax=ax, linewidths=0.3, annot=True, annot_kws={'fontsize':3}, cbar=True, cmap='seismic',  center=0,	#cool
									cbar_kws= { 'label': 'log2 of gene expression z-scores', 'shrink': .2 } 
								)
	#cmap='YlGnBu'  #binary	
	#set min and max values: vmin=0, vmax= my_vmax
	# change annotated values to integer: fmt="g",
	#square=True,
	x.figure.axes[-1].yaxis.label.set_size(5)
	x.figure.axes[-1].tick_params(axis='both', which='major', labelsize=5)
	x.figure.axes[-1].tick_params(axis='both', which='minor', labelsize=5)
	
	for idx, gene in enumerate( genes ):
		ax.text( -2.75, idx+0.6, gene, fontsize=5 )
	
	for idx, tissue in enumerate( tissues ):
		ax.text( idx+0.4, len( genes )+1, tissue, fontsize=5, rotation=90 )	#, rotation=90
	
	ax.set_yticklabels( [], rotation=0, fontsize=2 )
	ax.set_xticklabels( [] , rotation=90, fontsize=3  )
	
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('none')
	
	plt.yticks( rotation=0 )
	plt.subplots_adjust( left=0.15, right=0.99, top=0.99, bottom=0.12, wspace=0.2 )
	
	plt.savefig( heatmap_file, dpi=300  )


def load_sample_infos( sample_info_file ):
	"""! @brief load sample info from given file """
	
	candidate_samples = {}
	sample_order = []
	
	with open( sample_info_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			sample_order.append( parts[0] )
			candidate_samples.update( { parts[0]: parts[1].split(',') } )
			line = f.readline()
	
	return candidate_samples, sample_order


def load_candidate_genes( gene_file ):
	"""! @brief load candidate gene names """
	
	candidate_gene_order = []
	gene_name_mapping_table = {}
	
	with open( gene_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			candidate_gene_order.append( parts[0] )
			gene_name_mapping_table.update( { parts[0]: parts[1] } )
			line = f.readline()
	
	return candidate_gene_order, gene_name_mapping_table


def main( arguments ):
	
	data_file = arguments[ arguments.index( '--exp' )+1 ]
	gene_file = arguments[ arguments.index( '--genes' )+1 ]
	heatmap_file = arguments[ arguments.index( '--out' )+1 ]
	sample_info_file = arguments[ arguments.index( '--samples' )+1 ]

	data, samples, vals_per_sample = load_data( data_file )
	
	outputfile = heatmap_file + "plotted_values.txt"
	candidate_gene_order, gene_name_mapping_table = load_candidate_genes( gene_file )
	candidate_samples, sample_order = load_sample_infos( sample_info_file )
	
	genes, tissues, datamatrix = construct_data_output_file( data, candidate_gene_order, candidate_samples, outputfile, gene_name_mapping_table, sample_order )
	
	construct_heatmap( datamatrix, genes, tissues, heatmap_file )


if "--exp" in sys.argv and '--genes' in sys.argv and '--samples' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
