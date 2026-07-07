#!/usr/bin/env ruby

# get_pi_sites.rb from mPCRselect version 0.3.2
# Michael G. Campana, 2022-2024
# Smithsonian Institution

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford
# University have waived all copyright and related or neighboring rights to mPCRselect;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see
# <http://creativecommons.org/publicdomain/zero/1.0/>.

# We politely request that this work be cited as:
# Armstrong EE, Li C, Campana MG, Ferrari T, Kelley JL, Petrov DA, Solari KA, Mooney JA.
# In prep. Recommendations for population and individual diagnostic SNP selection in non-
# model species.

# Script to get highest pi sites from VCFtools site-pi file
# Sites are ranked by cluster density (number of other high-pi sites within 400 bp on
# the same chromosome), then by pi as a tiebreaker, so that co-amplifiable SNP groups
# are prioritised for multiplex PCR panel design.
# Usage: ruby get_pi_sites.rb <file> <minimum-Pi> <hard-maximum-number-of-sites> <soft-percent-sites-to-retain>

WINDOW = 400  # bp window for cluster scoring

@header = true # Switch to ignore header line
@pos_count = 0.0 # Total number of sites
@sites = [] # Array of [pi, line] for candidate sites above minimum pi
File.open(ARGV[0]) do |f1|
	while line = f1.gets
		if @header
			@header = false
		elsif line.strip != "" # Discount blank lines
			@pos_count += 1
			pi = line.strip.split[2].to_f
			@sites.push([pi, line]) if pi >= ARGV[1].to_f
			@pos_count += 1.0
		end
	end
end

# Determine actual number of sites to retain
soft_limit = (@pos_count * ARGV[3].to_f).to_i
@limit = 0
ARGV[2].to_i < soft_limit ? @limit = ARGV[2].to_i : @limit = soft_limit
@limit = @sites.size if @limit > @sites.size

if @limit > 0
	# Build per-chromosome sorted [pos_int, key] arrays from all candidates above minPi.
	# Using all candidates (not just the top @limit) ensures cluster scoring reflects
	# the full local SNP density available for multiplex PCR amplicon design.
	chrom_data = {}  # chr -> array of [pos_int, "chr\tpos_str"]
	@sites.each do |pi, line|
		parts = line.strip.split
		chr     = parts[0]
		pos_str = parts[1]
		chrom_data[chr] ||= []
		chrom_data[chr] << [pos_str.to_i, "#{chr}\t#{pos_str}"]
	end
	chrom_data.each_value { |arr| arr.sort_by! { |x| x[0] } }

	# Compute cluster score for each site using an O(N) sliding window per chromosome.
	# cluster_score["chr\tpos"] = number of OTHER high-pi sites within WINDOW bp.
	cluster_score = {}
	chrom_data.each do |chr, pos_arr|
		n     = pos_arr.size
		left  = 0
		right = 0
		pos_arr.each_with_index do |(pos, key), i|
			right = i if right < i
			right += 1 while right + 1 < n && pos_arr[right + 1][0] <= pos + WINDOW
			left  += 1 while pos_arr[left][0]  <  pos - WINDOW
			cluster_score[key] = right - left  # window count minus self
		end
	end

	# Sort by cluster score descending, then by pi descending as tiebreaker
	@sorted_sites = @sites.sort do |a, b|
		a_parts = a[1].strip.split
		b_parts = b[1].strip.split
		a_key   = "#{a_parts[0]}\t#{a_parts[1]}"
		b_key   = "#{b_parts[0]}\t#{b_parts[1]}"
		cs_cmp  = (cluster_score[b_key] || 0) <=> (cluster_score[a_key] || 0)
		cs_cmp != 0 ? cs_cmp : b[0] <=> a[0]
	end

	for i in 0 ... @limit
		puts @sorted_sites[i][1]
	end
end
