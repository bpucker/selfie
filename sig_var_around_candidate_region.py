### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python sig_var_around_candidate_region.py
					--sig <VCF_WITH_SIGNIFICANT_VARIANTS>
					--all <VCF_FILE_WITH_ALL_VARIANTS>
					--seq <NAME_OF_SEQ_OF_INTEREST>
					--fig <FIGURE_FILENAME>
					
					optional:
					--window <SIZE_OF_SLIDING_WINDOW>
					--step <STEP_SIZE>
					--x1 <INTERVALL_STAR_POSITION>
					--x2 <INTERVALL_END_POSITION>
					--score <ACTIVATES_P_VALUE_BASED_SCORE>
					
					bug reports and feature requests:
					bpucker@cebitec.uni-bielefeld.de
					"""


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os, sys, math

# --- end of imports --- #

def load_var_pos( vcf, seq_of_interest ):
	"""! @brief load all variants on sequence of interest """
	
	var_pos = []
	adj_p_values = []
	with open( vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[0] == seq_of_interest:
					var_pos.append( int( parts[1] ) )
					try:
						adj_p_values.append( float( parts[7] ) )	#this only works for the sig variant file
					except:
						pass
			line = f.readline()
	return var_pos, adj_p_values


def generate_figure( fig_file, general_var_pos, sig_var_pos, seq_of_interest, window, step, xstart, xend, adj_p_values, score_status ):
	"""! @brief visualize distribution of all variants and the significant variants """
	
	# --- calculation of values --- #
	if xend:
		max_var_pos = 0 + xend	#use defined region
	else:
		max_var_pos = max( general_var_pos )	#use last variant position
	
	y1 = []		#general variants
	y2 = []		#significant variants
	y3 = []		#ratio
	x = []
	y4 = []	#adjusted p-values
	
	if xstart:
		start = 0 + xstart
	else:
		start = 0
		xstart = 0
	end = start + window
	while start < max_var_pos:
		general = 0
		sig = 0
		p = 0
		for pos in general_var_pos:
			if start < pos < end:
				general += 1
		for k, pos in enumerate( sig_var_pos ):
			if start < pos < end:
				sig += 1
				p += 1 / adj_p_values[ k ]
		
		y1.append( general )
		y2.append( sig )
		if general > 0:
			y3.append( sig / float( general ) )
		else:
			y3.append( 0 )
		if p > 0:
			y4.append( math.log( p / float( general ) ) )
		else:
			y4.append( 0 )
		x.append( (start+end) / 2000000.0 )
		
		start += step
		end += step
	
	# --- generating figure --- #
	fig, ax = plt.subplots()
	
	ax.plot( x, y1, marker="o", linewidth=0, markersize=1, color="grey" )	#
	
	ax2 = ax.twinx()
	ax2.plot( x, y2, marker="o", linewidth=0, markersize=1, color="blue")
	
	ax3 = ax.twinx()
	ax3.plot( x, y3, marker="o", linewidth=0, markersize=1, color="magenta" )
	
	
	ax4 = ax.twinx()
	if score_status:
		ax4.plot( x, y4, marker="o", linewidth=0, markersize=1, color="black" )
	
	ax.set_xlabel( "position on sequence [Mbp]" )
	#ax.set_title( seq_of_interest )
	
	ax.set_ylabel( "all variants" )
	ax2.set_ylabel( "significant variants" )
	
	
	ax.set_xlim( xstart/1000000.0, max_var_pos/1000000.0 )
	ax.set_ylim( 0, max(  y1 )*1.01 )
	ax2.set_ylim( 0, max( y2 )*1.01 )
	ax3.set_ylim( 0, max( y3 )*1.01 )
	ax4.set_ylim( 0, max( y4 )*1.01 )
	
	my_legend = [ 	mpatches.Patch(color='grey', label='all variants'),
								mpatches.Patch(color='blue', label='significant variants'),
								mpatches.Patch(color='magenta', label='normalized sig. variants'),
								mpatches.Patch(color='black', label='region score')	#region score is sum of inverse adj. p-values
							]
	
	ax.legend( handles=my_legend, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.15) )
	
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
	ax2.spines['top'].set_visible(False)
	ax2.spines['left'].set_visible(False)
	
	ax3.spines['top'].set_visible(False)
	ax3.spines['left'].set_visible(False)
	ax3.spines['right'].set_visible(False)
	
	ax3.get_yaxis().set_ticks([])
	
	ax4.spines['top'].set_visible(False)
	ax4.spines['left'].set_visible(False)
	ax4.spines['right'].set_visible(False)
	
	ax4.get_yaxis().set_ticks([])
	
	fig.savefig( fig_file )


def main( arguments ):
	"""! @brief run everything """

	sig_var_vcf = arguments[ arguments.index( '--sig' )+1 ]
	general_vcf = arguments[ arguments.index( '--all' )+1 ]
	seq_of_interest = arguments[ arguments.index( '--seq' )+1 ]
	fig_file = arguments[ arguments.index( '--fig' )+1 ]
	
	if '--window' in arguments:
		window = int( arguments[ arguments.index( '--window' )+1 ] )
	else:
		window = 500000
	
	if '--step' in arguments:
		step = int( arguments[ arguments.index( '--step' )+1 ] )
	else:
		step = 100000
	
	if '--x1' in arguments:
		xstart = int( arguments[ arguments.index( '--x1' )+1 ] )
	else:
		xstart = False

	if '--x' in arguments:
		xend = int( arguments[ arguments.index( '--x2' )+1 ] )
	else:
		xend = False
	
	if '--score' in arguments:
		score_status = True
	else:
		score_status = False

	general_var_pos, waste = load_var_pos( general_vcf, seq_of_interest )
	sig_var_pos, adj_p_values = load_var_pos( sig_var_vcf, seq_of_interest )
	generate_figure( fig_file, general_var_pos, sig_var_pos, seq_of_interest, window, step, xstart, xend, adj_p_values, score_status )


if '--all' in sys.argv and '--sig' in sys.argv and '--seq' in sys.argv and '--fig' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
