# -*- coding: utf-8 -*-
from __future__ import division
import sys
import urllib, urllib2
import os
import math
import re
import string
import time
from Bio.Blast import NCBIXML, NCBIWWW
import traceback
import shutil
import logging



aa_dict={'ALA' :  'A' , 'ARG' : 'R' , 'ASN' : 'N' , 'ASP' : 'D' , 'CYS' : 'C', 'GLU' : 'E' , 'GLN' : 'Q' , 'GLY' : 'G',
                 'HIS' : 'H', 'ILE' : 'I' , 'LEU' : 'L' , 'LYS' : 'K' , 'MET' : 'M', 'PHE' : 'F' , 'PRO' : 'P', 'SER' : 'S', 'THR': 'T',
                 'TRP' : 'W', 'TYR' : 'Y', 'VAL' : 'V'}

def modelit(seq):

    out_dir=''
    #read config file
    with open('config.txt','r') as config:
        for path in config:
            out_dir=str(path)
    os.chdir(out_dir)
    pdb_dataset="PDB_dataset" #set path of pdb_dataset generated in 1st step
    #get logger
    logger=fragMod_logger()
   

            
    #identifier=raw_input("Enter unique identifier : ")
    #in_seq=raw_input("Enter Primary sequence of your query protein: \n")
    in_seq=seq
    #print in_seq
    
    in_seq=re.sub('\n','',in_seq)
    in_length=len(in_seq)
    ############################################### BLAST against swissprot #############################################################
    domain_seq=[' ' for i in xrange(len(in_seq))]
    logger.info("running BLASTp")
    result_handle = NCBIWWW.qblast("blastp", "swissprot", str(in_seq))
    #result_handle=open('swissprot.xml','r')
    blast_records = NCBIXML.parse(result_handle)
    ranges=dict()
    acc=list()
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
                    
            for hsp in alignment.hsps:
                              
                if hsp.expect < 1:
                    if int(hsp.identities) > int(in_length/2):
                        print('\n\n****Alignment****')
                                     
                        accessin=alignment.accession
                        acc.append(str(accessin))            
                        ranges.update({str(accessin) : str(hsp.sbjct_start)+'-'+str(hsp.sbjct_end)})
                                     
                        '''print ('Accession:',alignment.accession)
                        print hsp.sbjct_start
                        print hsp.sbjct_end
                        print hsp.query[0:]+"\n"+hsp.match[0:]+"\n"+hsp.sbjct[0:]'''
                        
                     
    logger.info("Fetching Domains from uniprot")

    query_domains_seq=list()
    domains_name=list()
    domain_ranges=list()
    for a in acc:
        try:
            rang=str(ranges[str(a)]).split('-')
            
            swispro=urllib.urlopen('http://www.uniprot.org/uniprot/'+str(a)+'.txt')
            data=swispro.readlines()
            
            for r in data:
                k=r.split()
               
                try:
                    if str(k[1]).startswith('DOMAIN'):     
                            start=int(k[2])
                            end=int(k[3])
                            domain=' '.join(k[4:])
                            domain=re.sub(' ','_',domain)
                            domain=domain.upper()
                            #domain=re.sub(".","",domain)
                            
                            #check if domain exists in the aligned region
                            if start > int(rang[0])-1:
                                if end < int(rang[1])+1:
                                    
                                    try:
                                            
                                            seq_index = [i for i, item in enumerate(query_domains_seq) if re.search(str(in_seq[start-int(rang[0]):end-int(rang[0])]), item)]
                                            
                                            domain_name_index=[i for i, item in enumerate(domains_name) if re.search(str(domain), item)]
                                            #print domain_name_index
                                            #print seq_index
                                            
                                            if len(seq_index) == 0 and len(domain_name_index) == 0:
                                                
                                                
                                                #print domain
                                                query_domains_seq.append(str(in_seq[start-int(rang[0]):end-int(rang[0])]))
                                                domain_ranges.append([int(start-int(rang[0])),int(end-int(rang[0]))])
                                                #print in_seq[start-int(rang[0]):end-int(rang[0])]
                                                logger.info(domain+" found")
                                                domains_name.append(str(domain)) #save name of domain relative to the sequence
                                                for i in xrange(start-int(rang[0]),end-int(rang[0])):
                                                    domain_seq[i] = 'D'
                                            else:
                                                traceback.print_exc()
                                    except ValueError:
                                        traceback.print_exc()
                except IndexError:
                    traceback.print_exc()
                                            
                                        
                                    
                                    
                                
                        
        except:
            logger.warning("No domain exists in "+str(a))
            pass

    #remove overlapping domain boundaries

    logger.info("removing overlapping domain boundaries")
    new_ranges=list()
    new_dom_names=list()

    for i in range(len(domain_ranges)):
        k=domain_ranges[i]
        print(k)
        bad=0
        for d in domain_ranges:
            if int(d[0]) <= int(k[0]) and int(k[1]) <= int(d[1]):
                #print 'caught',k
                bad=bad+1
                
            
            

        if bad < 2:
            #print 'GOOD: ',k
            new_dom_names.append(str(domains_name[i]))
            new_ranges.append([int(k[0]),int(k[1])])
            logger.info("registering domain sequence in query")
            query_domains_seq.append(str(in_seq[int(k[0]):int(k[1])])) #register domain sequences in query 
            #print str(in_seq[int(k[0]):int(k[1])])
    #for i in range(len(new_ranges)):
        #logger.info("Domain name: "+new_dom_names[i]+"\t Range:"+new_ranges[i])
    coil_ranges=list()

    #save coil ranges
    for j in range(len(new_ranges)-1):
        first=new_ranges[j]
        second=new_ranges[j+1]
        diff=int(second[0])-int(first[1])
        if diff > 15:
            coil_ranges.append([int(first[1]+1),int(second[0]-1)])

    if len(coil_ranges) > 0:
        logger.info("Listing coil range(s)...")
    for c in coil_ranges:        
        logger.info('Coil: ',c)


    
    #model coils
    
   
    for i in coil_ranges:
        seq=in_seq[int(i[0]):int(i[1])]
        name=str(i[0])+'-'+str(i[1])
        os.mkdir(out_dir+'/'+name+'_files', 0777)
        #print name
        logger.info("runnig BLASTp on coil sequence: "+seq)
        ok=1
        #definitions
        
        coil_query=dict()
        coil_pdbid=list()
        coil_chain=list()
        coil_sb_strt=dict()
        coil_sbjct_en=dict()
        coil_query=dict()
        coil_match=dict()
        coil_sbjct=dict()
        coil_q_start=dict()
        coil_q_end=dict()
        
        while ok > 0:
            try:
                result_handle = NCBIWWW.qblast("blastp", "pdb", str(seq))
                ok=0
                logger.warning("trying to reconnect..")
            except:
                ok=1
        logger.info("connected. Fetching results..")
        #result_handle=open('coil.xml','r')
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            
            for alignment in blast_record.alignments:
                        
                for hsp in alignment.hsps:
                                  
                    
                    #if int(hsp.identities) > int(in_length/2):
                        
                    #print '\n\n****Alignment****'          
                            
                    coil_acc=str(alignment.accession).split('_')
                        
                    coil_pdbid.append(str(coil_acc[0]))   #save pdbid of best aligned results
                        
                    #print str(coil_acc[0])
                        
                    #alignment_file.write(str(coil_acc[0])+'\n\n')
                        
                    coil_chain.append(str(coil_acc[1])) #save chain
                        
                    coil_q_start.update({str(coil_acc[0]) : hsp.query_start})
                        
                    coil_q_end.update({str(coil_acc[0]) : hsp.query_end})
                        
                    coil_sb_strt.update({str(coil_acc[0]) : hsp.sbjct_start}) #save start range of alignment
                            
                    coil_sbjct_en.update({str(coil_acc[0]) : hsp.sbjct_end}) #save end range of alignment
                            
                    #print 'length:', alignment.length
                        
                    #print str(coil_acc[0])
                        
                    #print 'Aligned Length:', hsp.align_length
                         
                    #print 'e_value:', hsp.expect
                        
                    #print 'score:', hsp.score
                        
                    #print 'Identities:',hsp.identities
                        
                    coil_query.update({str(coil_acc[0]) : hsp.query}) #save query
                        
                    coil_match.update({str(coil_acc[0]) : hsp.match}) #save match
                        
                    coil_sbjct.update({str(coil_acc[0]) : hsp.sbjct}) #save template
                        
                    #alignment_file.write(str(hsp.query[0:])+"\n"+str(hsp.match[0:])+"\n"+str(hsp.sbjct[0:])+'\n\n########################################\n\n')
                        
                    #print hsp.query[0:]+"\n"+hsp.match[0:]+"\n"+hsp.sbjct[0:]
        for k in range(len(coil_pdbid)):
            pro=str(coil_pdbid[k])
            #print pro
            #ch=str(coil_chain[k])
            TmpStrt=coil_sb_strt[pro]
            TmpEnd=coil_sbjct_en[pro]
            coil_query_seq=str(coil_query[pro])
            coil_match_seq=list(coil_match[pro])
            coil_sbjct_seq=str(coil_sbjct[pro])
            coil_query_start=int(coil_q_start[pro])
            coil_query_end=int(coil_q_end[pro])
            
            req = urllib.urlopen('http://rcsb.org/pdb/files/%s.pdb' %pro)
            
            
            the_page = req.readlines()
            
            logger.info("Extracting Secondary Structure of template ...")
            helix_initial=list()
            helix_terminal=list()
            sheet_initial=list()
            sheet_terminal=list()
            for line in the_page:
                            
                if line.startswith('HELIX'):
                                    
                                
                    if line[19:20].strip() == str(ch):
                        ini=int(line[22:25].strip())
                        ter=int(line[34:37].strip())
                        helix_initial.append(ini)
                        helix_terminal.append(ter)  

                if line.startswith('SHEET'):
                    if line[21:22].strip() == str(ch):
                        ini=int(line[23:26].strip())
                        ter=int(line[34:37].strip())
                        sheet_initial.append(ini)
                        sheet_terminal.append(ter)
         
                        
            #check how many secondary structure elements are present in aligned region
            NoOfHel=0
            NoOfSheet=0
            logger.info("checking how many secondary structure elements are present in aligned region")
            for k in xrange(len(helix_initial)):
                    
                    
                if (TmpStrt <= int(helix_initial[k]) <= TmpEnd) or (TmpStrt <= int(helix_terminal[k]) <= TmpEnd):
                        
                    NoOfHel=NoOfHel+1
                        
                else:
                    pass
            for k in xrange(len(sheet_initial)):
                    
                if (TmpStrt <= int(sheet_initial[k]) <= TmpEnd) or (TmpStrt <= int(sheet_terminal[k]) <= TmpEnd):                    
                    NoOfSheet=NoOfSheet+1
                        
                else:
                    pass
                
            logger.info("Helices:"+NoOfHel+"\t Sheets:"+NoOfSheet)
            
            coil_missings=list()
            if NoOfHel < 2 and NoOfSheet < 2:
                    
                #finding gaps:
                logger.info("finding gap regions")
                coil_blanks=list()
                
                for i in xrange(len(coil_match_seq)):
                    m=str(coil_match_seq[i])
                    q=str(coil_query_seq[i])
                    if q != '-':
                        
                        if m != q:
                            
                            if m == ' ' or m == '+':
                                
                                
                                en=0
                                st=0
                                #extend range forward $ backward
                                
                                if i > 3 and i < (len(seq)-3):
                                    st=i-1
                                    en=i+1
                                else:
                                    pass
                                
                                
                                if en != 0 and st !=0:
                                    #print int(st+coil_query_start),'-',int(en+coil_query_start+1)
                                    coil_blanks.append([st+coil_query_start,en+(coil_query_start+1)]) #save ranges of gaps
                                else:
                                    pass
                if coil_blanks:
                    
                    for ran in coil_blanks:
                        
                        
                        try:
                            logger.info(coil_query_seq[ran[0]:ran[1]]+' is a LOOP')
                                #search in LOOPS Dataset
                            for root, dirs, files in os.walk(pdb_dataset+"\LOOPS\\"):
                                for f in files:
                                        

                                    coil_loop_seq=coil_query_seq[int(ran[0]):int(ran[1])]
                                    if re.search("("+coil_loop_seq+".*)",str(f)):
                                        #print f
                                            
                                            
                                        for three, one in aa_dict.iteritems():   #convert the substring into three letter codes of AA
                                            if one == coil_loop_seq[0]:
                                                first=three
                                            if one == coil_loop_seq[1]:
                                                second=three
                                            if one == coil_loop_seq[2]:
                                                third=three

                                            

                                                    
                                        coil_loop_temp=open(pdb_dataset+"\LOOPS\\"+str(f),'r')
                                        coil_loop_temp_lines=coil_loop_temp.readlines()
                                        a_resid=0
                                        b_resid=0
                                        c_resid=0
                                            
                                        for a in coil_loop_temp_lines:        #check first residue and save its residue ID
                                                
                                                
                                            a_cols=a.split()
                                                
                                                
                                            if str(a_cols[3]) == first:
                                                a_resid=float(a[23:26])
                                                break
                                                    
                                                
                                        for c in coil_loop_temp_lines:    #check third residue and save its residue ID

                                                    
                                            c_cols=c.split()
                                                            
                                            if str(c_cols[3]) == third and int(c[23:26]) == a_resid+2:
                                                c_resid=float(c[23:26])
                                                break
                                                                    
                                        if a_resid+2 == c_resid:    #if first and last residue are same then check middle (main) residue
                                                
                                            for b in coil_loop_temp_lines:
                                                                               
                                                b_cols=b.split()
                                                        
                                                if str(b_cols[3]) == second and int(b[23:26]) == a_resid+1:
                                                                
                                                    b_resid=float(b[23:26])
                                                        
                                        if b_resid: #if substring matches then fetch the coordinates 
                                            coil_gap_name=int(ran[0])+(TmpStrt-coil_query_start)+1
                                            coil_missings.append(coil_gap_name)
                                            coil_gap_fil=open(out_dir+'/'+name+'_files/'+str(coil_gap_name)+'.pdb','a+')
                                            CA_x=0
                                            CA_y=0
                                            CA_z=0
                                                
                                            for lin in coil_loop_temp_lines:
                                                cols=lin.split()
                                                if int(lin[23:26]) == b_resid:
                                                    if str(cols[2]).strip() == 'CA':
                                                        CA_x=float(cols[6])   #getting coordinates of carbon alpha
                                                        CA_y=float(cols[7])
                                                        CA_z=float(cols[8])

                                            for lin in coil_loop_temp_lines:
                                                cols=lin.split()
                                                if int(lin[23:26]) == b_resid:
                                                    
                                                    if str(cols[2]).strip() == 'CA':
                                                        continue
                                                    
                                                    elif str(cols[2]).strip() == 'C':
                                                        continue
                                                    elif str(cols[2]).strip() == 'N':
                                                        continue
                                                    elif str(cols[2]).strip() == 'O':
                                                        continue
                                                    else:
                                                            
                                                        x=CA_x-float(cols[6])     #getting differenc between CA and others
                                                        y=CA_y-float(cols[7])
                                                        z=CA_z-float(cols[8])
                                                        #print lin           
                                                        coil_gap_fil.write(str(cols[2])+'\t'+str(cols[3])+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\n')

                                            coil_missings.append(coil_gap_name)
                                            coil_gap_fil.close()
                                            coil_loop_temp.close()
                                            del a_resid
                                            del b_resid
                                            del c_resid
                                        break
                        except:
                            pass
                
                logger.info("scanning template...")
                ini_res=0
                for l in the_page:
                    if l.startswith('ATOM'):
                        if l[21:22].strip() =='A':
                            ini_res=int(l[23:26])
                            break
                logger.info("writing PDB file..")
                atom_no=1
                coil_pdb_out=open(out_dir+'/'+name+'_files/'+name+'_'+pro+'.pdb','a+')
                if int(len(coil_missings)) > 0:
                    
                    for i in xrange(TmpStrt,TmpEnd):
                        
                        try:
                                        
                            if coil_missings.index(int(i)):
                                #print "############################"
                                #print i
                                #print "############################"
                                
                                ref_x=0
                                ref_y=0
                                ref_z=0
                                f=''
                                for line in the_page:
                                    
                                    if line.startswith('ATOM'):
                                        if str(line[21:22]) == str(ch):
                                            if int(line[23:26].strip()) == int(i+ini_res):
                                                cols=line.split()
                                                if str(cols[2]).strip() == 'N':
                                                    atom_name=str(cols[2])[0]
                                                    
                                                    f=f+str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+'AAA'.ljust(3)+' '+str(cols[4]).rjust(1)+str(cols[5]).rjust(4)+'    '+str(cols[6]).rjust(8)+str(cols[7]).rjust(8)+str(cols[8]).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n'
                                                    atom_no=atom_no+1
                                                elif str(cols[2]).strip() == 'CA':
                                                    ref_x=float(cols[6])   #getting coordinates of carbon alpha
                                                    ref_y=float(cols[7])
                                                    ref_z=float(cols[8])
                                                    atom_name=str(cols[2])[0]
                                                    
                                                    f=f+str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+'AAA'.ljust(3)+' '+str(cols[4]).rjust(1)+str(cols[5]).rjust(4)+'    '+str(cols[6]).rjust(8)+str(cols[7]).rjust(8)+str(cols[8]).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n'
                                                    atom_no=atom_no+1
                                                elif str(cols[2]).strip() == 'C':
                                                    atom_name=str(cols[2])[0]
                                                    
                                                    f=f+str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+'AAA'.ljust(3)+' '+str(cols[4]).rjust(1)+str(cols[5]).rjust(4)+'    '+str(cols[6]).rjust(8)+str(cols[7]).rjust(8)+str(cols[8]).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n'
                                                    atom_no=atom_no+1
                                                elif str(cols[2]).strip() == 'O':
                                                    atom_name=str(cols[2])[0]
                                                    
                                                    f=f+str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+'AAA'.ljust(3)+' '+str(cols[4]).rjust(1)+str(cols[5]).rjust(4)+'    '+str(cols[6]).rjust(8)+str(cols[7]).rjust(8)+str(cols[8]).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n'
                                                    atom_no=atom_no+1
                                                else:
                                                    pass
                                            
                                                
                                coil_mising_file=open(out_dir+'/'+name+'_files/'+str(i)+'.pdb','r')
                                coil_mising_lines=coil_mising_file.readlines()
                                for m in coil_mising_lines:
                                    m_cols=m.split()
                                    aa=str(m_cols[1])
                                    f=re.sub('AAA',str(aa),f)
                                    xx=float(m_cols[2])
                                    yy=float(m_cols[3])
                                    zz=float(m_cols[4])
                                    new_x=xx+ref_x
                                    new_y=yy+ref_y
                                    new_z=zz+ref_z
                                    atom_name=str(m_cols[0])[0]
                                    #f=f+str('ATOM').ljust(6)+str(atom_no).rjust(5)+'  '+str(m_cols[0]).ljust(4)+' '+str(m_cols[1]).ljust(3)+' '+str(ch).rjust(1)+str(i).rjust(4)+'    '+str(new_x).ljust(8)+str(new_y).ljust(8)+str(new_z).ljust(8)+str('1.00').ljust(6)+str('0.00').ljust(6)+str(atom_name).ljust(12)+'\n'
                                    
                                    f=f+str('ATOM').ljust(6)+str(atom_no).rjust(5)+'  '+str(m_cols[0]).ljust(4)+str(m_cols[1]).ljust(3)+' '+str(ch).rjust(1)+str(i+ini_res).rjust(4)+'    '+str(new_x).rjust(8)+str(new_y).rjust(8)+str(new_z).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n'
                                    atom_no=atom_no+1
                                print(f)
                                coil_pdb_out.write(str(f))
                                
                        except ValueError:
                            logger.error("Value Error at index",i)
                            for line in the_page:
                                
                                if line.startswith('ATOM'):
                                    if str(line[21:22]) == str(ch):
                                        
                                            
                                            if int(line[23:26].strip()) == int(i+ini_res):
                                                cols=line.split()
                                                atom_name=str(cols[2])[0]
                                                    
                                                
                                                #pdb_out.write(str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+' '+str(cols[3]).ljust(3)+' '+str(cols[4]).rjust(1)+str(cols[5]).rjust(4)+'    '+str(cols[6]).ljust(8)+str(cols[7]).ljust(8)+str(cols[8]).ljust(8)+str('1.00').ljust(6)+str('0.00').ljust(6)+str(atom_name).ljust(12)+'\n')
                                                coil_pdb_out.write(str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+str(cols[3]).ljust(3)+' '+str(cols[4]).rjust(1)+str(i+ini_res).rjust(4)+'    '+str(cols[6]).rjust(8)+str(cols[7]).rjust(8)+str(cols[8]).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n')

                                                atom_no=atom_no+1
                else:
                    for i in xrange(TmpStrt,TmpEnd):
                        for line in the_page:
                                    
                            if line.startswith('ATOM'):
                                if str(line[21:22]) == str(ch):
                                            
                                                
                                    if int(line[23:26].strip()) == int(i+ini_res):
                                        coil_pdb_out.write(line)

                coil_pdb_out.close()
               
                break
            else:
                pass
    
    #model domains one by one    
    for i in xrange(len(query_domains_seq)):
        
      
        seq=str(query_domains_seq[i])    
        dom=str(domains_name[i])    
        domain_directory_name=str(dom[0:4])
        if(os.path.isdir(domain_directory_name) == False):
            os.mkdir(domain_directory_name, 0777)
        logger.info("running BLASTp on PDB\n seq:"+seq+"\ndomain: "+str(dom))
        
        ok=1
        
        while ok > 0:
            try:
                result_handle = NCBIWWW.qblast("blastp", "pdb", str(seq))
                ok=0
            except:
                logger.exception("Exception occured while running BLASTp. Trying again...")
                ok=1
        
        blast_results=result_handle.read()
        #Next, we save the output
        save_file = open(domain_directory_name+'/'+domain_directory_name+'_BLAST.xml', "w")
        save_file.write(blast_results)
        save_file.close()
        
        logger.info("Fetching Results ")
        
        result_handle=open(domain_directory_name+'/'+domain_directory_name+'_BLAST.xml','r')
        blast_records = NCBIXML.parse(result_handle)
        
        #definitions
        
        query=dict()
        pdbid=list()
        chain=list()
        sb_strt=dict()
        sbjct_en=dict()
        query=dict()
        match=dict()
        sbjct=dict()
        q_start=dict()
        q_end=dict()
        identities=list()
        logger.info("creating alignment file...")
        
        alignment_file=open(str(dom)+'_alignment.txt','a+')
        
        for blast_record in blast_records:
            
            for alignment in blast_record.alignments:
                        
                for hsp in alignment.hsps:
                                  
                    
                    if hsp.identities:
                        
                        #print '\n\n****Alignment****'
                        
                        acc=str(alignment.accession).split('_')
                        
                        pdbid.append(str(acc[0]))   #save pdbid of best aligned results
                        
                        #print str(acc[0])
                        
                        alignment_file.write(str(acc[0])+'\n\n')
                        
                        chain.append(str(acc[1])) #save chain
                        
                        q_start.update({str(acc[0]) : hsp.query_start})
                        
                        q_end.update({str(acc[0]) : hsp.query_end})
                        
                        sb_strt.update({str(acc[0]) : hsp.sbjct_start}) #save start range of alignment
                            
                        sbjct_en.update({str(acc[0]) : hsp.sbjct_end}) #save end range of alignment
                            
                        #print 'length:', alignment.length
                        
                        #print str(acc[0])
                        
                        #print 'Aligned Length:', hsp.align_length
                         
                        #print 'e_value:', hsp.expect
                        
                        #print 'score:', hsp.score
                        
                        #print 'Identities:',hsp.identities
                        identities.append(int(hsp.identities))
                        
                        query.update({str(acc[0]) : hsp.query}) #save query
                        
                        match.update({str(acc[0]) : hsp.match}) #save match
                        
                        sbjct.update({str(acc[0]) : hsp.sbjct}) #save template
                        
                        alignment_file.write(str(hsp.query[0:])+"\n"+str(hsp.match[0:])+"\n"+str(hsp.sbjct[0:])+'\n\n########################################\n\n')
                        
                        #print hsp.query[0:]+"\n"+hsp.match[0:]+"\n"+hsp.sbjct[0:]

        max_id=max(identities) #get index of maximum identity temp
        max_index=identities.index(max_id)
        print(max_index)
        #for i in xrange(len(pdbid)):
        try:
            pro=str(pdbid[max_index])
           # pro='2ETM'
            ch=str(chain[max_index])
            #ch='A'
            query_seq=str(query[pro])
            match_seq=list(match[pro])
            sbjct_seq=str(sbjct[pro])
            temp_start=int(sb_strt[pro])
            #print temp_start
            
            temp_end=int(sbjct_en[pro])
            ##print temp_end
            query_start=int(q_start[pro])
           
            
            query_end=int(q_end[pro])
            #print query_seq
            #print match_seq
            
            #print sbjct_seq
            
            
            req = urllib.urlopen('http://rcsb.org/pdb/files/%s.pdb' %pro)
            #req=open('GAPS/P'+str(dom)+"_new.pdb",'r')
            the_page = req.readlines()
            #new_seq=''
            #for new_lines in the_page:
              #  if new_lines[13:15] == 'CA':
              #      res=str(new_lines[17:20])
              #      try:
              #          new_seq=new_seq+aa_dict[res]
              #      except:
              #          new_seq=new_seq+'Z'
            #print "######\n"+new_seq
            
            logger.info("finding gaps")
            
            blanks=list()
            
            for i in xrange(len(match_seq)):
                m=str(match_seq[i])
                q=str(query_seq[i])
                if q != '-':
                    
                    if m != q:
                        
                        if m == ' ' or m == '+':
                            
                            
                            en=0
                            st=0
                            #extend range forward $ backward
                            
                            if i > 3 and i < (len(seq)-3):
                                st=i-1
                                en=i+1
                            else:
                                pass
                            
                            
                            if en != 0 and st !=0:
                                print(int(st+query_start),'-',int(en+query_start+1))
                                blanks.append([st+query_start,en+(query_start+1)]) #save ranges of gaps
                            else:
                                pass
                      


            logger.info("Extracting Secondary Structure of template ...")
            helix_initial=list()
            helix_terminal=list()
            sheet_initial=list()
            sheet_terminal=list()
            for line in the_page:
                            
                if line.startswith('HELIX'):
                                    
                                
                    if line[19:20].strip() == str(ch):
                        ini=int(line[22:25].strip())
                        ter=int(line[34:37].strip())
                        helix_initial.append(ini)
                        helix_terminal.append(ter)  

                if line.startswith('SHEET'):
                    if line[21:22].strip() == str(ch):
                        ini=int(line[23:26].strip())
                        ter=int(line[34:37].strip())
                        sheet_initial.append(ini)
                        sheet_terminal.append(ter)
         
                        
            logger.info("searching gaps coordinates in dataset...")
            missings=list()
            for ran in blanks:
                
                
                ss=0
                
                #check if gap sequence is a part of HELIX or SHEET
                for k in xrange(len(helix_initial)):
                    
                    
                    if (helix_initial[k] <= (int(ran[0])+(temp_start-query_start)) <= helix_terminal[k]) or (helix_initial[k] <= (int(ran[1])+(temp_start-query_start)) <= helix_terminal[k]):
                        
                        ss=1
                        
                    else:
                        pass
                for k in xrange(len(sheet_initial)):
                    
                    if (sheet_initial[k] <= (int(ran[0])+(temp_start-query_start)) <= sheet_terminal[k]) or (sheet_initial[k] <= (int(ran[1])+(temp_start-query_start)) <= sheet_terminal[k]):
                        
                        ss=2
                        
                    else:
                        pass

                  
                if ss == 0:
                    logger.info(query_seq[ran[0]:ran[1]]+' is a LOOP')
                    #search in LOOPS Dataset
                    for root, dirs, files in os.walk(pdb_dataset+"\LOOPS\\"):
                        for f in files:
                            

                            loop_seq=query_seq[int(ran[0]):int(ran[1])]
                            if re.search("("+loop_seq+".*)",str(f)):
                                print f
                                
                                
                                for three, one in aa_dict.iteritems():   #convert the substring into three letter codes of AA
                                    if one == loop_seq[0]:
                                        first=three
                                    if one == loop_seq[1]:
                                        second=three
                                    if one == loop_seq[2]:
                                        third=three

                                

                                        
                                loop_temp=open(pdb_dataset+"\LOOPS\\"+str(f),'r')
                                loop_temp_lines=loop_temp.readlines()
                                a_resid=0
                                b_resid=0
                                c_resid=0
                                
                                for a in loop_temp_lines:        #check first residue and save its residue ID
                                    
                                    
                                    a_cols=a.split()
                                    
                                    
                                    if str(a_cols[3]) == first:
                                        a_resid=float(a[23:26])
                                        break
                                        
                                    
                                for c in loop_temp_lines:    #check third residue and save its residue ID

                                        
                                    c_cols=c.split()
                                                
                                    if str(c_cols[3]) == third and int(c[23:26]) == a_resid+2:
                                        c_resid=float(c[23:26])
                                        break
                                                        
                                if a_resid+2 == c_resid:    #if first and last residue are same then check middle (main) residue
                                    
                                    for b in loop_temp_lines:
                                                                   
                                        b_cols=b.split()
                                            
                                        if str(b_cols[3]) == second and int(b[23:26]) == a_resid+1:
                                                    
                                            b_resid=float(b[23:26])
                                            
                                if b_resid: #if substring matches then fetch the coordinates 
                                    gap_name=int(ran[0])+(temp_start-query_start)+1
                                    missings.append(gap_name)
                                    gap_fil=open(str(dom)+'_files/'+str(loop_seq[1])+str(gap_name)+'.pdb','a+')
                                    CA_x=0
                                    CA_y=0
                                    CA_z=0
                                    
                                    for lin in loop_temp_lines:
                                        cols=lin.split()
                                        if int(lin[23:26]) == b_resid:
                                            if str(cols[2]).strip() == 'CA':
                                                CA_x=float(cols[6])   #getting coordinates of carbon alpha
                                                CA_y=float(cols[7])
                                                CA_z=float(cols[8])

                                    for lin in loop_temp_lines:
                                        cols=lin.split()
                                        if int(lin[23:26]) == b_resid:
                                            
                                            if str(cols[2]).strip() == 'CA':
                                                continue
                                            
                                            elif str(cols[2]).strip() == 'C':
                                                continue
                                            elif str(cols[2]).strip() == 'N':
                                                continue
                                            elif str(cols[2]).strip() == 'O':
                                                continue
                                            else:
                                                
                                                x=CA_x-float(cols[6])     #getting differenc between CA and others
                                                y=CA_y-float(cols[7])
                                                z=CA_z-float(cols[8])
                                                #print lin           
                                                gap_fil.write(str(cols[2])+'\t'+str(cols[3])+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\n')
                                    gap_fil.close()
                                    loop_temp.close()
                                    del a_resid
                                    del b_resid
                                    del c_resid
                                break            
                    
                elif ss == 1:
                    logger.info(query_seq[ran[0]:ran[1]]+' is a HELIX')
                    #search in HELIX Dataset
                    for root, dirs, files in os.walk(pdb_dataset+"\HELIX\\"):
                        for f in files:
                            

                            hel_seq=query_seq[int(ran[0]):int(ran[1])]
                            if re.search("("+hel_seq+".*)",str(f)):
                                #print f
                                
                                
                                for three, one in aa_dict.iteritems():   #convert the substring into three letter codes of AA
                                    if one == hel_seq[0]:
                                        first=three
                                    if one == hel_seq[1]:
                                        second=three
                                    if one == hel_seq[2]:
                                        third=three

                                

                                        
                                hel_temp=open(pdb_dataset+"\HELIX\\"+str(f),'r')
                                hel_temp_lines=hel_temp.readlines()
                                a_resid=0
                                b_resid=0
                                c_resid=0
                                
                                for a in hel_temp_lines:        #check first residue and save its residue ID
                                    
                                    
                                    a_cols=a.split()
                                    
                                    
                                    if str(a_cols[3]) == first:
                                        a_resid=float(a[23:26])
                                        break
                                        
                                    
                                for c in hel_temp_lines:    #check third residue and save its residue ID

                                        
                                    c_cols=c.split()
                                                
                                    if str(c_cols[3]) == third and int(c[23:26]) == a_resid+2:
                                        c_resid=float(c[23:26])
                                        break
                                                        
                                if a_resid+2 == c_resid:    #if first and last residue are same then check middle (main) residue
                                    
                                    for b in hel_temp_lines:
                                                                   
                                        b_cols=b.split()
                                            
                                        if str(b_cols[3]) == second and int(b[23:26]) == a_resid+1:
                                                    
                                            b_resid=float(b[23:26])
                                            
                                if b_resid: #if substring matches then fetch the coordinates 
                                    gap_name=int(ran[0])+(temp_start-query_start)+1
                                    missings.append(gap_name)
                                    gap_fil=open(str(dom)+'_files/'+str(hel_seq[1])+str(gap_name)+'.pdb','a+')
                                    CA_x=0
                                    CA_y=0
                                    CA_z=0
                                    
                                    for lin in hel_temp_lines:
                                        cols=lin.split()
                                        if int(lin[23:26]) == b_resid:
                                            if str(cols[2]).strip() == 'CA':
                                                CA_x=float(cols[6])   #getting coordinates of carbon alpha
                                                CA_y=float(cols[7])
                                                CA_z=float(cols[8])

                                    for lin in hel_temp_lines:
                                        cols=lin.split()
                                        if int(lin[23:26]) == b_resid:
                                            
                                            if str(cols[2]).strip() == 'CA':
                                                continue
                                            
                                            elif str(cols[2]).strip() == 'C':
                                                continue
                                            elif str(cols[2]).strip() == 'N':
                                                continue
                                            elif str(cols[2]).strip() == 'O':
                                                continue
                                            else:
                                                
                                                x=CA_x-float(cols[6])     #getting differenc between CA and others
                                                y=CA_y-float(cols[7])
                                                z=CA_z-float(cols[8])
                                                #print lin           
                                                gap_fil.write(str(cols[2])+'\t'+str(cols[3])+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\n')
                                    gap_fil.close()
                                    hel_temp.close()
                                break              
                                                
                                        
                elif ss == 2:
                    logger.info(query_seq[ran[0]:ran[1]]+' is a SHEET')
                    #search in SHEET dataset
                    for root, dirs, files in os.walk(pdb_dataset+"\SHEET\\"):
                        for f in files:

                            sheet_seq=query_seq[int(ran[0]):int(ran[1])]    
                            if re.search("("+sheet_seq+".*)",str(f)):
                                print f
                                
                                                
                                for three, one in aa_dict.iteritems():   #convert the substring into three letter codes of AA
                                    if one == sheet_seq[0]:
                                        first=three
                                    if one == sheet_seq[1]:
                                        second=three
                                    if one == sheet_seq[2]:
                                        third=three


                                                        
                                sheet_temp=open(pdb_dataset+"\SHEET\\"+str(f),'r')
                                sheet_temp_lines=sheet_temp.readlines()
                                a_resid=0
                                b_resid=0
                                c_resid=0
                                                
                                for a in sheet_temp_lines:        #check first residue and save its residue ID
                                                    
                                                    
                                    a_cols=a.split()
                                                    
                                                    
                                    if str(a_cols[3]) == first:
                                        a_resid=float(a[23:26])
                                        break
                                                        
                                                    
                                for c in sheet_temp_lines:    #check third residue and save its residue ID

                                                        
                                    c_cols=c.split()
                                                                
                                    if str(c_cols[3]) == third and int(c[23:26]) == a_resid+2:
                                        c_resid=float(c[23:26])
                                        break
                                                                        
                                if a_resid+2 == c_resid:    #if first and last residue are same then check middle (main) residue
                                                    
                                    for b in sheet_temp_lines:
                                                                                   
                                        b_cols=b.split()
                                                            
                                        if str(b_cols[3]) == second and int(b[23:26]) == a_resid+1:
                                                                    
                                            b_resid=float(b[23:26])
                                                            
                                if b_resid: #if substring matches then fetch the coordinates 
                                    gap_name=int(ran[0])+(temp_start-query_start)+1
                                    missings.append(gap_name)
                                    gap_fil=open(str(dom)+'_files/'+str(sheet_seq[1])+str(gap_name)+'.pdb','a+')
                                    CA_x=0
                                    CA_y=0
                                    CA_z=0
                                    
                                    for lin in sheet_temp_lines:
                                        cols=lin.split()
                                        if int(lin[23:26]) == b_resid:
                                            if str(cols[2]).strip() == 'CA':
                                                CA_x=float(cols[6])   #getting coordinates of carbon alpha
                                                CA_y=float(cols[7])
                                                CA_z=float(cols[8])

                                    for lin in sheet_temp_lines:
                                        cols=lin.split()
                                        if int(lin[23:26]) == b_resid:
                                            
                                            if str(cols[2]).strip() == 'CA':
                                                continue
                                            
                                            elif str(cols[2]).strip() == 'C':
                                                continue
                                            elif str(cols[2]).strip() == 'N':
                                                continue
                                            elif str(cols[2]).strip() == 'O':
                                                continue
                                            else:
                                                
                                                x=CA_x-float(cols[6])     #getting differenc between CA and others
                                                y=CA_y-float(cols[7])
                                                z=CA_z-float(cols[8])
                                                            
                                                gap_fil.write(str(cols[2])+'\t'+str(cols[3])+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\n')
                                            
                                    gap_fil.close()
                                    sheet_temp.close()
                                break
                 
            logger.info("scanning template...")
            ini_res=0
            for l in the_page:
                if l.startswith('ATOM'):
                    if l[21:22].strip() =='A':
                        ini_res=int(l[23:26])
                        break
            logger.info("writing PDB file..")
            atom_no=1
            pdb_out=open(domain_directory_name+'/'+domain_directory_name+'.pdb','a+')
            if int(len(missings)) > 0:
                
                for i in xrange(temp_start,temp_end):
                    
                    try:
                                    
                        if missings.index(int(i)):
                            #print "############################"
                            #print i
                            #print "############################"
                            
                            ref_x=0
                            ref_y=0
                            ref_z=0
                            f=''
                            for line in the_page:
                                
                                if line.startswith('ATOM'):
                                    if str(line[21:22]) == str(ch):
                                        if int(line[23:26].strip()) == int(i+ini_res):
                                            cols=line.split()
                                            if str(cols[2]).strip() == 'N':
                                                atom_name=str(cols[2])[0]
                                                
                                                f=f+str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+'AAA'.ljust(3)+' '+str(cols[4]).rjust(1)+str(cols[5]).rjust(4)+'    '+str(cols[6]).rjust(8)+str(cols[7]).rjust(8)+str(cols[8]).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n'
                                                atom_no=atom_no+1
                                            elif str(cols[2]).strip() == 'CA':
                                                ref_x=float(cols[6])   #getting coordinates of carbon alpha
                                                ref_y=float(cols[7])
                                                ref_z=float(cols[8])
                                                atom_name=str(cols[2])[0]
                                                
                                                f=f+str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+'AAA'.ljust(3)+' '+str(cols[4]).rjust(1)+str(cols[5]).rjust(4)+'    '+str(cols[6]).rjust(8)+str(cols[7]).rjust(8)+str(cols[8]).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n'
                                                atom_no=atom_no+1
                                            elif str(cols[2]).strip() == 'C':
                                                atom_name=str(cols[2])[0]
                                                
                                                f=f+str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+'AAA'.ljust(3)+' '+str(cols[4]).rjust(1)+str(cols[5]).rjust(4)+'    '+str(cols[6]).rjust(8)+str(cols[7]).rjust(8)+str(cols[8]).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n'
                                                atom_no=atom_no+1
                                            elif str(cols[2]).strip() == 'O':
                                                atom_name=str(cols[2])[0]
                                                
                                                f=f+str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+'AAA'.ljust(3)+' '+str(cols[4]).rjust(1)+str(cols[5]).rjust(4)+'    '+str(cols[6]).rjust(8)+str(cols[7]).rjust(8)+str(cols[8]).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n'
                                                atom_no=atom_no+1
                                            else:
                                                pass
                                        
                                            
                            mising_file=open(domain_directory_name+'/'+str(i)+'.pdb','r')
                            mising_lines=mising_file.readlines()
                            for m in mising_lines:
                                m_cols=m.split()
                                aa=str(m_cols[1])
                                f=re.sub('AAA',str(aa),f)
                                xx=float(m_cols[2])
                                yy=float(m_cols[3])
                                zz=float(m_cols[4])
                                new_x=xx+ref_x
                                new_y=yy+ref_y
                                new_z=zz+ref_z
                                atom_name=str(m_cols[0])[0]
                                #f=f+str('ATOM').ljust(6)+str(atom_no).rjust(5)+'  '+str(m_cols[0]).ljust(4)+' '+str(m_cols[1]).ljust(3)+' '+str(ch).rjust(1)+str(i).rjust(4)+'    '+str(new_x).ljust(8)+str(new_y).ljust(8)+str(new_z).ljust(8)+str('1.00').ljust(6)+str('0.00').ljust(6)+str(atom_name).ljust(12)+'\n'
                                
                                f=f+str('ATOM').ljust(6)+str(atom_no).rjust(5)+'  '+str(m_cols[0]).ljust(4)+str(m_cols[1]).ljust(3)+' '+str(ch).rjust(1)+str(i+ini_res).rjust(4)+'    '+str(new_x).rjust(8)+str(new_y).rjust(8)+str(new_z).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n'
                                atom_no=atom_no+1
                            #print f
                            pdb_out.write(str(f))
                            
                    except ValueError:
                        logger.error("ValueError at ",i)
                        for line in the_page:
                            
                            if line.startswith('ATOM'):
                                if str(line[21:22]) == str(ch):
                                    
                                        
                                        if int(line[23:26].strip()) == int(i+ini_res):
                                            cols=line.split()
                                            atom_name=str(cols[2])[0]
                                                
                                            
                                            #pdb_out.write(str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+' '+str(cols[3]).ljust(3)+' '+str(cols[4]).rjust(1)+str(cols[5]).rjust(4)+'    '+str(cols[6]).ljust(8)+str(cols[7]).ljust(8)+str(cols[8]).ljust(8)+str('1.00').ljust(6)+str('0.00').ljust(6)+str(atom_name).ljust(12)+'\n')
                                            pdb_out.write(str(cols[0]).ljust(6)+str(atom_no).rjust(5)+'  '+str(cols[2]).ljust(4)+str(cols[3]).ljust(3)+' '+str(cols[4]).rjust(1)+str(i+ini_res).rjust(4)+'    '+str(cols[6]).rjust(8)+str(cols[7]).rjust(8)+str(cols[8]).rjust(8)+str('1.00').rjust(6)+str('0.00').rjust(6)+str(atom_name).rjust(12)+'\n')

                                            atom_no=atom_no+1
            else:
                for i in xrange(temp_start,temp_end):
                    for line in the_page:
                                
                        if line.startswith('ATOM'):
                            if str(line[21:22]) == str(ch):
                                        
                                            
                                if int(line[23:26].strip()) == int(i+ini_res):
                                    pdb_out.write(line)

            pdb_out.close()
            logger.info("Job Completed")
           
        except:
            
            traceback.print_exc()

    
#modelit(in_seq)
def fragMod_logger():
    logging.basicConfig(filename="fragmod.log",
                            filemode='a',
                            format='%(asctime)s, %(name)s \t %(levelname)s: %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG)
    return logging.getLogger('FragMode')
