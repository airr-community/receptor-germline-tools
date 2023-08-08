add_germline_annotations rearrangement.tsv rearrangement_output.tsv germline_set.json
add_germline_annotations rearrangement_output.tsv rearrangement_undummified.tsv germline_set.json -u -c v_call_dummy,d_call_dummy,j_call_dummy
convert_fasta_labels labelled_fasta.fasta dummified_fasta.fasta -g fasta_germline_set.json
convert_fasta_labels dummified_fasta.fasta undummified_fasta.fasta -g fasta_germline_set.json -u
create_fasta Human_IGH_VDJ_rev_3_ex.json Human_IGH_VDJ_rev_3_ex.fasta -n -f gs_release_version,release_version,aliases,paralogs
create_fasta "Mouse_129S1_SvlmJ IGLV_rev_2.json" "Mouse_129S1_SvlmJ IGLV_rev_2.fasta" -n -f gs_release_version,release_version,aliases,paralogs
