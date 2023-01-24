#!/usr/bin/ruby

# hoelzer.martin@gmail.com

# Input: a tab-separated output file from covsonar (e.g. from --lineage XBB.1.5)
# Task: calculate the "consensus mutation profile" where all mutations are above a certain frequency cutoff (e.g. 75%)
# Output: Consensus mutation profile 

require 'bio'

freq_cutoff = 0.75

covsonar_tsv = File.open(ARGV[0],'r')

hits = 0
lineage = ''
nt_mutations = {}
aa_changes = {}
covsonar_tsv.each do |line|
    next if line.start_with?('accession')
    hits += 1
    s = line.split("\t")
    # get nt mutations
    s[19].split(' ').each do |nt_mutation|
        if nt_mutations[nt_mutation]
            nt_mutations[nt_mutation] += 1
        else
            nt_mutations[nt_mutation] = 1
        end
    end
    # get aa changes
    s[20].split(' ').each do |aa_change|
        if aa_changes[aa_change]
            aa_changes[aa_change] += 1
        else
            aa_changes[aa_change] = 1
        end
    end
    # get lineage
    lineage = s[17].chomp.strip if lineage == ''
end
covsonar_tsv.close

# filter based on frequency
nt_mutations_filtered = {}
aa_changes_filtered = {}

nt_mutations.each do |nt_mutation, count|
    freq = count.to_f / hits
    nt_mutations_filtered[nt_mutation] = freq.round(3) if freq >= freq_cutoff
end

aa_changes.each do |aa_change, count|
    freq = count.to_f / hits
    aa_changes_filtered[aa_change] = freq.round(3) if freq >= freq_cutoff
end

# write out again a covsonar-styled file
output = File.open(ARGV[0].sub('.tsv','.consensus.tsv'),'w')
output << %w(accession	description	lab	source	collection	technology	platform	chemistry	material	ct	software	software_version	gisaid	ena	zip	date	submission_date	lineage	seqhash	dna_profile	aa_profile	fs_profile).join("\t")
output << "\n"
dna_profile = nt_mutations_filtered.keys.join(' ')
aa_profile = aa_changes_filtered.keys.join(' ')
output << "CONSENSUS-PROFILE\tCONSENSUS-PROFILE\t\tmerge-profiles.rb\t\t\t\t\t\t\t\t\t\t\t\t\t#{lineage}\t\t\t#{dna_profile}\t#{aa_profile}\t\n"
output.close


# read in reference genome and replace INDELs and substitutions - experimental feature bc/ covsonar will maybe better do that...
consensus_seq = ''
Bio::FastaFormat.open('NC_045512.2.fasta').each do |entry|
    ref_seq = entry.seq
    ref_seq_a = ref_seq.split('')
    nt_mutations_filtered.each do |nt, position|
        if nt.start_with?('del:')
            # deletion, e.g. "del:29734:26"
            ref_pos = nt.split(':')[1].to_i
            del_length = nt.split(':')[2].to_i
            # iterate over the positions that should be deleted and mask them
            # delete them later to have no index problems
            ref_pos = ref_pos - 1 
            del_length.times do |i|
                ref_seq_a[ref_pos+i] = '?'
            end
        else
            # substitution or insertion
            ref_pos = nt.scan(/\d+/).first
            ref_nt = nt.split(ref_pos)[0]
            alternative_nt = nt.split(ref_pos)[1]
            #puts "#{ref_nt}\t#{ref_pos}\t#{alternative_nt}"
            ref_pos = ref_pos.to_i
            if ref_seq_a[ref_pos-1] == ref_nt
                ref_seq_a[ref_pos-1] = alternative_nt
            else
                abort 'STOP! This should not happen. We check if the correct reference nucleotide is replaced at the extracted position. Please debug...'
            end
        end
    end
    # now finally create the consensus and also remove the DELETIONS
    consensus_seq = ref_seq_a.join('').gsub('?','')
end

output = File.open(ARGV[0].sub('.tsv','.consensus.fasta'),'w')
output << ">#{lineage}-Consensus\n#{consensus_seq.chomp}\n"
output.close

puts aa_changes_filtered