#functions to build multiple template-based model using Modeller

from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb


def build_model(seq, tempid, ch):
    # Step 1: Create a query .ali file
    heading = ">P1;query" + "\n" + "sequence:query" + ":::::::0.00: 0.00"
    query_seq = heading + "\n" + seq + "*"
    file2 = open('query.ali', 'w')
    file2.write(query_seq)
    file2.close()
    
    # Step 2: Align query sequence against the template
    env = environ()
    aln = alignment(env)
    mdl = model(env, file=tempid, model_segment=('FIRST:A', 'LAST:A'))
    aln.append_model(mdl, atom_files=tempid, align_codes=tempid + ch)
    aln.append(file='query.ali', align_codes='query')
    aln.align2d()
    aln.write(file='query_' + tempid + '.ali', alignment_format='PIR')
    
    # Step 3: Model the structure
    env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)
    env.io.atom_files_directory = ['.']
    a = automodel(env, alnfile='query_' + tempid + '.ali', knowns=tempid + ch,
                  sequence='query', assess_methods=(assess.DOPE, assess.GA341))
    a.starting_model = a.ending_model = 1
    a.library_schedule = autosched.slow
    a.max_var_iterations = 300
    a.md_level = refine.slow
    a.repeat_optimization = 2
    a.max_molpdf = 1e6
    a.make()
    
    # Additional assessments
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    mdl = complete_pdb(env, 'query.B99990001.pdb')
    s = selection(mdl)
    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='query.profile',
                  normalize_profile=True, smoothing_window=15)

def multi_model(seq, tmp, new_al, strc_det):
    # Create a query .ali file for multiple models
    heading = ">P1;query" + "\n" + "sequence:query" + ":::::::0.00: 0.00"
    query_seq = heading + "\n" + seq + '*'
    file2 = open('PDB-dataset/MULTI/query.ali', 'w')
    file2.write(query_seq)
    file2.close()

    log.verbose()
    env = environ()
    env.io.atom_files_directory = ['PDB-dataset/MULTI']
    knwns = [str(t) + 'A' for t in tmp]

    # Align multiple templates
    aln = alignment(env)
    for code in tmp:
        chain = 'A'
        mdl = model(env, file=code, model_segment=('FIRST:A', 'LAST:A'))
        aln.append_model(mdl, atom_files=code, align_codes=code + chain)

    # Perform structure alignment
    for (weights, write_fit, whole) in [((1., 0., 0., 0., 1., 0.), False, True),
                                        ((1., 0.5, 1., 1., 1., 0.), False, True),
                                        ((1., 1., 1., 1., 1., 0.), True, False)]:
        aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                   rr_file='$(LIB)/as1.sim.mat', overhang=30,
                   gap_penalties_1d=(-450, -50),
                   gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
                   dendrogram_file='fm00495.tree',
                   alignment_type='tree',
                   feature_weights=weights,
                   improve_alignment=True, fit=True, write_fit=write_fit,
                   write_whole_pdb=whole, output='ALIGNMENT QUALITY')

    aln.write(file='PDB-dataset/MULTI/temp_alignment.ali', alignment_format='PIR')

    # Perform additional alignments
    aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
               gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
               alignment_type='progressive', feature_weights=[0]*6,
               improve_alignment=False, fit=False, write_fit=True,
               write_whole_pdb=False, output='QUALITY')

    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    log.verbose()
    env = environ()
    aln = alignment(env)

    # Read aligned structures and sequences
    aln.append(file='PDB-dataset/MULTI/temp_alignment.ali', align_codes='all')
    aln_block = len(aln)
    aln.append(file='PDB-dataset/MULTI/query.ali', align_codes='query')

    # Perform structure-sensitive variable gap penalty sequence-sequence alignment
    aln.salign(output='', max_gap_length=20,
               gap_function=True,
               alignment_type='PAIRWISE', align_block=aln_block,
               feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
               gap_penalties_1d=(-450, 0),
               gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
               similarity_flag=True)

    aln.write(file='PDB-dataset/MULTI/final_alignment.ali', alignment_format='PIR')

    # Edit alignment file
    old_aln = open('PDB-dataset/MULTI/final_alignment.ali', 'r')
    new_aln = open('PDB-dataset/MULTI/new_final_alignment.ali', 'a+')
    old_aln_lines = old_aln.readlines()
    headers = list()
    strc = list()

    # Extract headers and structures from old alignment
    for lin in old_aln_lines:
        if lin.startswith('>'):
            headers.append(lin)
    for x in range(len(new_al)):
        new_aln.write(str(headers[x]))
        new_aln.write(str(strc_det[x]))
        new_aln.write(str(new_al[x])+'*\n')

    # Add query sequence to new alignment
    new_aln.write('>P1;query\nsequence:query:     : :     : ::: 0.00: 0.00\n' + str(seq) + '*')

    new_aln.close()
    old_aln.close()

    env.schedule_scale = physical.values(default=1.0, soft_sphere=0.
