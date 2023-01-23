#!/usr/bin/ruby

# hoelzer.martin@gmail.com

# Input: a tab-separated output file from covsonar (e.g. from --lineage XBB.1.5)
# Task: calculate the "consensus mutation profile" where all mutations are above a certain frequency cutoff (e.g. 75%)
# Output: Consensus mutation profile 

freq_cutoff = 0.50

covsonar_tsv = File.open(ARGV[0],'r')

hits = 0
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
output = File.open(ARGV[0].sub('.tsv','.consensus-profile.tsv'),'w')
output << %w(accession	description	lab	source	collection	technology	platform	chemistry	material	ct	software	software_version	gisaid	ena	zip	date	submission_date	lineage	seqhash	dna_profile	aa_profile	fs_profile).join("\t")
output << "\n"
lineage = ''
dna_profile = nt_mutations_filtered.keys.join(' ')
aa_profile = aa_changes_filtered.keys.join(' ')
output << "CONSENSUS-PROFILE\tCONSENSUS-PROFILE\t\tmerge-profiles.rb\t\t\t\t\t\t\t\t\t\t\t\t\t#{lineage}\t\t\t#{dna_profile}\t#{aa_profile}\t\n"
output.close