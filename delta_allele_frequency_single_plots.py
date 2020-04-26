### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.3 ###

__usage__ = """
					python delta_allele_frequency.py\n
					--input_vcf <FILENAME>
					--reference_file <FILENAME>
					--output_dir <DIRECTORY_NAME>[will be generated if required]
					--pool1 <sample name in VCF; multiple samples names can be provided comma-seperated>
					--pool2 <sample name in VCF; multiple samples names can be provided comma-seperated>
					
					optional:
					--minP1cov <FLOAT, lower coverage cutoff for pool1>
					--maxP1cov <FLOAT, upper coverage cutoff for pool1>
					--minP2cov <FLOAT, lower coverage cutoff for pool2>
					--maxP2cov <FLOAT, upper coverage cutoff for pool2>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import matplotlib.pyplot as plt
import sys, os
import numpy as np

# --- end of imports --- #

def get_coverage( input_vcf, high_name, low_name ):
	"""! @brief get average coverage for both samples of interest """
	
	coverage_high = []
	coverage_low = []	
	
	with open( input_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[6] == "PASS":
					sample_data = parts[9:]
					
					tmp_high = []
					for index in high_sample_idx:
						hd = sample_data[ index ]	#variant data for high sample
						if hd[:3] != "./.":
							hd_parts = hd.split(':')
							tmp_high.append( int( hd_parts[2] ) )
					if tmp_high > 0:
						coverage_high.append( sum( tmp_high ) )
					
					tmp_low = []
					for index in low_sample_idx:
						ld = sample_data[ index ]	#variant data for low sample
						if ld[:3] != "./.":
							ld_parts = ld.split(':')
							tmp_low.append( int( ld_parts[2] ) )
					if tmp_low > 0:
						coverage_low.append( sum( tmp_low ) )
			else:
				try:
					samples = line.strip().split('\t')[9:]
					high_sample_idx = []
					for each in high_name:
						try:
							high_sample_idx.append( samples.index( each ) )
						except ValueError:
							print "ERROR: sample not detected - " + each
					low_sample_idx = []
					for each in low_name:
						try:
							low_sample_idx.append( samples.index( each ) )
						except ValueError:
							print "ERROR: sample not detected - " + each
				except IndexError:
					pass
			line = f.readline()
	
	#plt.hist( coverage_high, bins = 20000 )
	#plt.title( "coverage high" )
	#plt.xlim( [ 0, 200 ] )
	#plt.show()
	
	#plt.hist( coverage_low, bins = 20000 )
	#plt.title( "coverage low" )
	#plt.xlim( [ 0, 200 ] )
	#plt.show()
	
	return sorted( coverage_high )[ len( coverage_high ) / 2 ], sorted( coverage_low )[ len( coverage_low ) / 2 ]


def get_delta_allel_frequencies( input_vcf, output_file, unique_variant_vcf, high_name, low_name, avg_cov_high, avg_cov_low, min_high_cov, max_high_cov, min_low_cov, max_low_cov ):
	"""! @brief calculate all possible delta allele frequencies
	@note reference allele frequency is devided by alternative allele frequency
	@note triallelic variants (and higher) are skipped
	"""
	
	#coverage cofidence intervall borders for pool1 (high)
	if not min_high_cov:
		min_high_cov = 0.75*avg_cov_high	#0.75
	if not max_high_cov:
		max_high_cov = 1.5*avg_cov_high	#1.5
	
	#coverage cofidence intervall borders for pool2 (low)
	if not min_low_cov:
		min_low_cov = 0.75 * avg_cov_low	#0.75
	if not max_low_cov:
		max_low_cov = 1.5*avg_cov_low	#1.5
	
	with open( unique_variant_vcf, "w" ) as unique_out:
		with open( output_file, "w" ) as out:
			with open( input_vcf, "r" ) as f:
				line = f.readline()
				while line:
					if line[0] != '#':
						parts = line.strip().split('\t')
						if parts[6] == "PASS" and not "," in parts[4]:	#avoid triallelic variants:
							sample_data = parts[9:]
							
							# --- combine information of all pool1 samples --- #
							h_cov = 0
							h_allele1 = 0
							h_allele2 = 0
							for each in high_sample_idx:
								hd = sample_data[ each ]	#variant data for high sample
								hd_parts = hd.split(':')
								try:
									try:
										h_cov += int( hd_parts[2] )
									except IndexError:
										pass
								except ValueError:
									pass
								if hd[:3] != "./.":
									h_allele1 += int( hd_parts[1].split(',')[0] )
									h_allele2 += int( hd_parts[1].split(',')[1] )
							
							# --- combine information of all pool2 samples --- #
							l_cov = 0
							l_allele1 = 0
							l_allele2 = 0
							for each in low_sample_idx:
								ld = sample_data[ each ]	#variant data for low sample
								ld_parts = ld.split(':')
								try:
									try:
										l_cov += int( ld_parts[2] )
									except IndexError:
										pass
								except ValueError:
									pass
								if ld[:3] != "./.":	
									l_allele1 += int( ld_parts[1].split(',')[0] )
									l_allele2 += int( ld_parts[1].split(',')[1] )
							
							# --- check if coverage of locus is within expected range ---- #
							if h_cov >= min_high_cov and h_cov <= max_high_cov and l_cov >= min_low_cov and l_cov <= max_low_cov:
								af_high_status = False
								try:
									af_high =  ( float( h_allele2 ) / float( h_allele1+h_allele2 ) )
								except ZeroDivisionError:
									af_high = 1
									af_high_status = True
								
								af_low_status = False
								try:
									af_low =  ( float( l_allele2 ) / float( l_allele1+l_allele2 ) )
								except ZeroDivisionError:
									af_low = 1
									af_low_status = True
								
								if af_high_status + af_low_status < 2:
									af_delta = af_high - af_low
									out.write( "\t".join( parts[:7] +  map( str, [ ".", ".", h_cov, h_allele1, h_allele2, l_cov, l_allele1, l_allele2, af_delta ] ) ) + '\n' )
							elif h_cov >= min_high_cov and h_cov <= max_high_cov and l_cov == 0:
								unique_out.write( "\t".join( parts[:7] + map( str, [ ".", ".", h_cov, h_allele1, h_allele2, l_cov, l_allele1, l_allele2, 1 ] ) ) + '\n' )
							elif l_cov >= min_low_cov and l_cov <= max_low_cov and h_cov == 0:
								unique_out.write( "\t".join( parts[:7] + map( str, [ ".", ".", h_cov, h_allele1, h_allele2, l_cov, l_allele1, l_allele2, -1 ] ) ) + '\n' )
					else:
						try:
							try:
								samples = line.strip().split('\t')[9:]
								high_sample_idx = []
								for each in high_name:
									high_sample_idx.append( samples.index( each ) )
								low_sample_idx = []
								for each in low_name:
									low_sample_idx.append( samples.index( each ) )
								out.write( "\t".join( line.strip().split('\t')[:9] + [ "Pool1Coverage\tPool1RefCov\tPool1AltCov\tPool2Coverage\tPool2RefCov\tPool2AltCov\tdelta_AF" ] ) + '\n' )
								unique_out.write( "\t".join( line.strip().split('\t')[:9] + [ "Pool1Coverage\tPool1RefCov\tPool1AltCov\tPool2Coverage\tPool2RefCov\tPool2AltCov\tdelta_AF" ] ) + '\n' )
							except ValueError:
								pass
						except IndexError:
							pass
					line = f.readline()


def plot_genome_wide_delta_allele_frequencies( af_frequency_vcf, unique_variant_vcf, chr_lengths, window_size, step_size ):
	"""! @brief show genome wide distribution of AF frequencies """
	
	# --- load information about chromosomes --- #
	chr_names = sorted( chr_lengths.keys() )
	raw_data_x = [ ]
	raw_data_y = [ ]
	af_data_y_pool1 = []	#same X values can be used for all
	af_data_y_pool2 = []
	for  k in chr_names:
		raw_data_x.append( [] )
		raw_data_y.append( [] )
		af_data_y_pool1.append( [] )
		af_data_y_pool2.append( [] )
	
	# --- load normal variant information (variants detected in both pools) --- #
	with open( af_frequency_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				raw_data_x[ chr_names.index( parts[0] ) ].append( int( parts[1] ) )
				raw_data_y[ chr_names.index( parts[0] ) ].append( abs( float( parts[-1] ) ) )
				af_data_y_pool1[ chr_names.index( parts[0] ) ].append( float( parts[10] ) / float( parts[9] ) )
				af_data_y_pool2[ chr_names.index( parts[0] ) ].append( float( parts[13] ) / float( parts[12] ) )
			line = f.readline()
	
	# # --- load information about unique variants --- #
	# additional_raw_data_x = []
	# additional_raw_data_y = []
	# for  k in chr_names:
		# additional_raw_data_x.append( [] )
		# additional_raw_data_y.append( [] )
	# with open( unique_variant_vcf, "r" ) as f:
		# line = f.readline()
		# while line:
			# if line[0] != '#':
				# parts = line.strip().split('\t')
				# additional_raw_data_x[ chr_names.index( parts[0] ) ].append( int( parts[1] ) )
				# additional_raw_data_y[ chr_names.index( parts[0] ) ].append( abs( float( parts[-1] ) ) )
			# line = f.readline()
	
	# --- prepare normal data (variants detected in both pools) for construction of plot --- #
	data_to_plot_x = [ ]
	data_to_plot_y = [ ]
	af_pool1 = []
	af_pool2 = []
	intervalls = []
	for k in chr_names:
		data_to_plot_x.append( [] )
		data_to_plot_y.append( [] )
		intervalls.append( [] )
		af_pool1.append( [] )
		af_pool2.append( [] )
	
	for idx, chr_data in enumerate( raw_data_x ):
		start = 0
		end = 0 + window_size
		while end < len( chr_data ):	
			intervalls[ idx ].append( ( chr_data[ start ], chr_data[ end ] ) )	
			if start > len( chr_data )-window_size:	#end analysis of this chromosome
				try:
					tmp_x = chr_data[ start: ]
					tmp_y = raw_data_y[ idx ][ start: ] 
					tmp_y_pool1 = af_data_y_pool1[ idx ][ start: ] 
					tmp_y_pool2 = af_data_y_pool2[ idx ][ start: ] 
					data_to_plot_x[ idx ].append( sum( tmp_x ) / ( len( tmp_x ) * 1000000.0 ) )
					data_to_plot_y[ idx ].append( np.median( tmp_y ) )		#sum( tmp_y ) / float( len( tmp_y ) )
					af_pool1[ idx ].append( np.median( tmp_y_pool1 ) )
					af_pool2[ idx ].append( np.median( tmp_y_pool2 ) )
				except ZeroDivisionError:
					data_to_plot_x[ idx ].append( 0 )
					data_to_plot_y[ idx ].append( 0 )
				break
				
			else:
				tmp_x = chr_data[ start:end ]
				tmp_y = raw_data_y[ idx ][ start:end ] 
				tmp_y_pool1 = af_data_y_pool1[ idx ][ start:end ] 
				tmp_y_pool2 = af_data_y_pool2[ idx ][ start:end ] 
				data_to_plot_x[ idx ].append( sum( tmp_x ) / ( len( tmp_x ) * 1000000.0 ) )
				data_to_plot_y[ idx ].append( np.median( tmp_y ) )		#sum( tmp_y ) / float( len( tmp_y ) )
				af_pool1[ idx ].append( np.median( tmp_y_pool1 ) )
				af_pool2[ idx ].append( np.median( tmp_y_pool2 ) )
				start += step_size
				end += step_size
	
	# # --- prepare unique variant data for plot construction --- #
	# additional_data_to_plot_x = [ ]
	# additional_data_to_plot_y = [ ]
	# for k in chr_names:
		# additional_data_to_plot_x.append( [] )
		# additional_data_to_plot_y.append( [] )
	
	# for idx, chr_data in enumerate( additional_raw_data_x ):
		# start = 0
		# end = 0 + window_size
		
		# while end < len( chr_data ):			
			# if start > len( chr_data )-window_size:	#end analysis of this chromosome
				# try:
					# tmp_x = chr_data[ start: ]
					# tmp_y = raw_data_y[ idx ][ start: ] 
					# additional_data_to_plot_x[ idx ].append( sum( tmp_x ) / ( len( tmp_x ) * 1000000.0 ) )
					# additional_data_to_plot_y[ idx ].append( len( tmp_y ) / float( tmp_x[-1]-temp_x[0] ) )	#density of unique variant in current interval
				# except ZeroDivisionError:
					# additional_data_to_plot_x[ idx ].append( 0 )
					# additional_data_to_plot_y[ idx ].append( 0 )
				# break
				
			# else:
				# tmp_x = chr_data[ start:end ]
				# tmp_y = raw_data_y[ idx ][ start:end ] 
				# additional_data_to_plot_x[ idx ].append( sum( tmp_x ) / ( len( tmp_x ) * 1000000.0 ) )
				# additional_data_to_plot_y[ idx ].append( len( tmp_y ) / float( tmp_x[-1]-tmp_x[0] ) )
				# start += step_size
				# end += step_size
	
	
	# --- construct single plots --- #
	for idx, each in enumerate( chr_names ):
		fig_output_file = af_frequency_vcf + "." + str( window_size ) + "." + each + ".genome_wide.png"
	
		fig, ax = plt.subplots(  figsize=(20, 3) )		
		# --- adding chromosomes --- #
		ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0, 0 ] , color="black", linewidth=.5 )
		
		# --- add helper lines --- #
		ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0.25, 0.25 ] , color="black", linewidth=.1 )
		ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0.5, 0.5 ] , color="black", linewidth=.1 )
		ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0.75, 0.75 ] , color="black", linewidth=.1 )
		
		# --- adding variant information --- #
		ax.scatter( data_to_plot_x[ idx ], data_to_plot_y[ idx ], c=map( abs, data_to_plot_y[ idx ] ), s=1, cmap="cool" )
		
		#ax.scatter( data_to_plot_x[ idx ], af_pool1[ idx ], s=1, color="red" )
		#ax.scatter( data_to_plot_x[ idx ], af_pool2[ idx ], s=1, color="lime" )
		
		
		#hot = high values invisible
		#cool = Hanna's choice
		#winter = better than hot
		
		## --- adding unique variant information --- #
		#ax.scatter( additional_data_to_plot_x[ idx ], additional_data_to_plot_y[ idx ], c=map( abs, additional_data_to_plot_y[ idx ] ), s=1, cmap="winter", marker="+" )		#hot
		
		ax.set_xlabel( "chromosome position [Mbp]" )
		ax.set_ylabel( "delta Allele Frequency" )
		ax.set_title( each )
		
		ax.spines["top"].set_visible(False)
		#ax.spines["left"].set_visible(False)
		ax.spines["right"].set_visible(False)
		#ax.set_frame_on(False)
		#ax.axes.get_yaxis().set_visible(False)
		
		ax.set_ylim( 0, 1 )
		ax.set_xlim(  0, chr_lengths[ each ]/1000000.0 )
		
		start, end = ax.get_xlim()
		ax.xaxis.set_ticks( np.arange( start, end, 10 ) )
		
		plt.subplots_adjust( left=0.03, right=0.98, top=0.98, bottom=0.2 )
		fig.savefig( fig_output_file, dpi=300 )
		plt.close('all')
	return data_to_plot_x, data_to_plot_y, intervalls	#additional_data_to_plot_x, additional_data_to_plot_y,


def plot_genome_wide_single_pos_dAF( af_frequency_vcf, unique_variant_vcf, chr_lengths ):
	"""! @brief show genome wide distribution of AF frequencies """
	
	# --- load information about chromosomes --- #
	chr_names = sorted( chr_lengths.keys() )
	raw_data_x = [ ]
	raw_data_y = [ ]
	additional_raw_data_x = []
	additional_raw_data_y = []
	for  k in chr_names:
		raw_data_x.append( [] )
		raw_data_y.append( [] )
		additional_raw_data_x.append( [] )
		additional_raw_data_y.append( [] )
	
	# --- load normal variant information (variants detected in both pools) --- #
	with open( af_frequency_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				raw_data_x[ chr_names.index( parts[0] ) ].append( int( parts[1] )/1000000.0 )
				raw_data_y[ chr_names.index( parts[0] ) ].append( abs( float( parts[-1] ) ) )
			line = f.readline()
	
	# --- load information about unique variants --- #
	with open( unique_variant_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				additional_raw_data_x[ chr_names.index( parts[0] ) ].append( int( parts[1] )/1000000.0 )
				additional_raw_data_y[ chr_names.index( parts[0] ) ].append( float( parts[-1] ) )
			line = f.readline()
	
	# --- construct single plots --- #
	for idx, each in enumerate( chr_names ):
		fig_output_file = af_frequency_vcf + "." + each + ".genome_wide.single_variants.png"
	
		fig, ax = plt.subplots(  figsize=(20, 3) )		
		# --- adding chromosomes --- #
		ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0, 0 ] , color="black", linewidth=.5 )
		for i in range(10):
			ax.plot( [ 0,  chr_lengths[ each ]/1000000.0 ], [ 0.1*i, 0.1*i ] , color="grey", linewidth=.1 )
		
		
		# --- adding variant information --- #
		ax.scatter( raw_data_x[ idx ], raw_data_y[ idx ], c=map( abs, raw_data_y[ idx ] ), s=1, cmap="cool" )	#hot, binary
		
		# --- adding unique variant information --- #
		ax.scatter( additional_raw_data_x[ idx ], additional_raw_data_y[ idx ], c=map( abs, additional_raw_data_y[ idx ] ), s=1, cmap="cool", marker="+" )
		
		ax.set_xlabel( "chromosome position [Mbp]" )
		ax.set_ylabel( "delta Allele Frequency" )
		ax.set_title( each )
		
		ax.spines["top"].set_visible(False)
		#ax.spines["left"].set_visible(False)
		ax.spines["right"].set_visible(False)
		#ax.set_frame_on(False)
		#ax.axes.get_yaxis().set_visible(False)
		
		ax.set_ylim( 0, 1 )
		ax.set_xlim(  0, chr_lengths[ each ]/1000000.0 )
		
		start, end = ax.get_xlim()
		ax.xaxis.set_ticks( np.arange( start, end, 10 ) )
		
		plt.subplots_adjust( left=0.03, right=0.98, top=0.98, bottom=0.2 )
		fig.savefig( fig_output_file, dpi=300 )
		plt.close('all')


def load_seq_lengths( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	seq_lens = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					seq_lens.update( { header: len( seq ) } )
					header = line.strip()[1:].split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		seq_lens.update( { header: len( seq ) } )
	return seq_lens


def construct_delta_AF_frequency_hist( af_frequency_vcf ):
	"""! @brief construct histogram of delta AF distribution """
	
	delta_AFs = []
	
	with open( af_frequency_vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				delta_AFs.append( float( line.strip().split('\t')[-1] ) )
			line = f.readline()
	
	fig_output_file = af_frequency_vcf + ".AF_hist.png"
	
	fig, ax = plt.subplots()
	ax.hist( delta_AFs, bins=10000 )
	ax.set_ylim( 0, 1000 )
	ax.set_ylabel( "number of variants" )
	ax.set_xlabel( "delta allele frequency (pool1 - pool2)" )
	
	fig.savefig( fig_output_file, dpi=300 )
	print "total number of calculated delta allele frequencies: " + str( len( delta_AFs ) )


def normalized_euclidean_distance( x, y ):
	"""! @brief calculate normalized euclidean distance """
	
	x = map( float, x )
	y = map( float, y )
	
	differences = []
	for i in range( len( x ) ):
		differences.append( x[i]-y[i] )
	
	distance = 0.5 * np.var( differences ) / ( np.var( x ) + np.var( y ) )
	return distance
	

def main( arguments ):
	"""! @brief run all parts of this script """
	
	# ---- collecting all inputs --- #
	input_vcf = arguments[ arguments.index( '--input_vcf' ) + 1 ]
	fasta_ref_file = arguments[ arguments.index( '--reference_file' ) + 1 ]
	output_dir = arguments[ arguments.index( '--output_dir' ) + 1 ]
	
	if not output_dir[-1] == "/":
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	raw_high_name = arguments[ arguments.index( '--pool1' ) + 1 ]
	high_name = []
	for each in raw_high_name.split(','):
		if len( each ) > 0:
			high_name.append( each )
	
	raw_low_name = arguments[ arguments.index( '--pool2' ) + 1 ]
	low_name = []
	for each in raw_low_name.split(','):
		if len( each ) > 0:
			low_name.append( each )
	
	if '--minP1cov' in arguments:
		min_high_cov = float( arguments[ arguments.index( '--minP1cov' ) + 1 ] )
	else:
		min_high_cov = False
	if '--maxP1cov' in arguments:
		max_high_cov = float( arguments[ arguments.index( '--maxP1cov' ) + 1 ] )
	else:
		max_high_cov = False
	
	if '--minP2cov' in arguments:
		min_low_cov = float( arguments[ arguments.index( '--minP2cov' ) + 1 ] )
	else:
		min_low_cov = False
	if '--maxP2cov' in arguments:
		max_low_cov = float( arguments[ arguments.index( '--maxP2cov' ) + 1 ] )
	else:
		max_low_cov = False
	
	
	# --- calling all functions and running detection of loci --- #
	print "pool1: " + str( high_name )
	print "pool1 size: " + str( len( high_name ) )
	print "pool2: " + str( low_name )
	print "pool2 size: "+ str( len( low_name ) )
	
	af_frequency_vcf = output_dir + "allele_frequencies.vcf"
	unique_variant_vcf = output_dir + "unique_variant.vcf"
	window_sizes = [ 100 ]	#25, 50, 100, 500 
	step_size = 5
	
	cov_high, cov_low = get_coverage( input_vcf, high_name, low_name )
	print "high cov: " + str( cov_high )
	print "low cov: " + str( cov_low )
	
	get_delta_allel_frequencies( input_vcf, af_frequency_vcf, unique_variant_vcf, high_name, low_name, cov_high, cov_low, min_high_cov, max_high_cov, min_low_cov, max_low_cov )
	
	construct_delta_AF_frequency_hist( af_frequency_vcf )
	
	
	chr_lengths = load_seq_lengths( fasta_ref_file )
	#print chr_lengths
	
	value_output_file = output_dir + "value_output_file.txt"
	with open( value_output_file, "w" ) as out:
		plot_genome_wide_single_pos_dAF( af_frequency_vcf, unique_variant_vcf, chr_lengths )
		for window_size in window_sizes:
			data_to_plot_x, data_to_plot_y, intervalls = plot_genome_wide_delta_allele_frequencies( af_frequency_vcf, unique_variant_vcf, chr_lengths, window_size, step_size )	#additional_data_to_plot_x, additional_data_to_plot_y, 
			chr_names = sorted( chr_lengths.keys() )
			for idx, k in enumerate( chr_names ):
				intervalls_of_interest = intervalls[ idx ]
				y_values_of_interest = data_to_plot_y[ idx ]
				for i, intervall in enumerate( intervalls_of_interest ):
					out.write( k + '\t' + str( y_values_of_interest[ i ] ) + '\t' + str( intervall[0] ) + '\t' + str( intervall[1] ) + '\n' )

if __name__ == '__main__':
	
	main( sys.argv )
	
	print "all done!"
