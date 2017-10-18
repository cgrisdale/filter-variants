#!/usr/bin/python

import sys
import re
from collections import Counter

#Take in BT- and PT-*transcriptome.somatic.annotated.synon-filtered.txt files, combine events with same eventID from same file into
#dict, then combine dict's into a single one by merging events with same ID; increase count (v[0]) for every sample that same event
#is found in, this makes up part of total count, unique count is number of samples that a particular event is found in. 

#Global variables
FileNames = {'PT-AB0029':"BT-106",'PT-VO7089':"BT-100",'PT-GC1519':"BT-108",'PT-AB6372':"BT-119",'PT-HN6692':"BT-126",'PT-BM772':"BT-127",'PT-BK0248':"BT-134",'PT-RL7940':"BT-140",'PT-EV3071':"BT-142",'PT-MB9777':"BT-143",'PT-DS9789':"BT-147",'PT-MD9088':"BT-156",'PT-GE7528':"BT-157",'PT-SK0976':"BT-161",'PT-JW6420':"BT-166",'PT-LR9369':"BT-167",'PT-RD1291':"BT-169",'PT-KM5291':"BT-172",'PT-JP2405':"BT-181",'PT-CM1209':"BT-189",'PT-SJ5453':"BT-191",'PT-AK7565':"BT-194",'PT-LS4891':"BT-198",'PT-PV2594':"BT-206",'PT-AH1410':"BT-208",'PT-SS3647':"BT-220",'PT-GJ3716':"BT-238",'PT-JE6375':"BT-241",'PT-CA2271':"BT-245",'PT-GR2309':"BT-246",'PT-CM1209':"BT-248",'PT-ML9537':"BT-275",'PT-TM5196':"BT-280",'PT-RV2286':"BT-282",'PT-WP9124':"BT-284",'PT-GB9483':"BT-287",'PT-AR3050':"BT-301",'PT-PD6881':"BT-41",'PT-GI2070':"BT-42",'PT-LC6372':"BT-47",'PT-CM3220':"BT-48",'PT-PR5617':"BT-50",'PT-DF5919':"BT-53",'PT-JB1730':"BT-63",'PT-FB6711':"BT-67",'PT-HO0394':"BT-69",'PT-RW9277':"BT-73",'PT-RL5404':"BT-75",'PT-WT4796':"BT-76",'PT-AL4257':"BT-84",'PT-GB9186':"BT-85",'PT-SO0258':"BT-89",'PT-WR7927':"BT-90",'PT-AR5365':"BT-92",'PT-LC3356':"BT-94"}

def EventNames(astring):
  '''Event names in different SV files sometimes differ, this will map to common event names'''
  m={'del':'deletion', 'ins':'insertion', 'dup':'duplication', 'inv':'inversion', 'ITD':'ITD'} #ITD=internal tandem duplication

  try:
    return m[astring]

  except:
    print astring
    raise ValueError('Three letter abbreviation not present')

def Read_in_files(File):
  '''Read in file, create features, add to dictionary with duplicates (based on event positions) being combined as single dict entry'''
  fnameBT,fnamePT,fname,newlist='','','',[]

  with open(File, 'r') as f:

    if str(File)[0:2]=="BT":
      fnameBT=re.search('BT-[0-9]*',str(File)).group()
      fname=fnameBT

    elif str(File)[0:2]=="PT":
      fnamePT=re.search('PT-[A-Z]{0,4}[0-9]*',str(File)).group()
      if fnamePT in FileNames:
        fname=FileNames[fnamePT]
      else:
        sys.exit("PT- file name does not exist in dictionary")

    else:
      print str(File)
      sys.exit("Not BT or PT file")

    svd={}
    header_line=f.readline() #save header

    for line in f:
      myfeatures=line.strip().split('\t') #myfeatures[1:6]= BT-xx, gene, indel, evtype, positions, GBM_xx;GBM_yy
      event=myfeatures[0].split('_')
      readsup=event[-1]
      eventid=myfeatures[0].split('_(')[0]

      gene=myfeatures[6]

      #If it's a mitochondrial gene, skip the rest of the while loop so it doesn't get added
      if gene[0:3]=="MT-":
        continue
      
      #te=[fname]+[gene]+myfeatures[2:3]+myfeatures[4:6]+[readsup]
      #print te
      #if event ID already in dict add additional features for current event to same event ID
      if eventid in svd:

        #myfeatures[4] is breakpoint positions
        if str(myfeatures[4])==str(svd[eventid][3]) and str(myfeatures[2])==str(svd[eventid][2]): #same breakpoint and event type

          if gene==svd[eventid][1]: #same gene name in same file, do not add any info

            print "Same eventID, positions, event type, and gene name. This should not happen.",myfeatures,svd[eventid]
            pass

          else: #additional gene name, add to values unless it's already contained therein

            if str(gene) in str(svd[eventid][1]) or str(svd[eventid][1]) in str(gene): #if there are combined gene names, don't duplicate names
              pass

            else:
              newgene=str(svd[eventid][1]+"_"+gene)
              newlist=[svd[eventid][0]]+[newgene]+svd[eventid][2:]
              svd[eventid]=newlist #only add gene name

        else:
          print "Position identifiers don't match",str(File),myfeatures[4],svd[eventid][4],gene,str(svd[eventid][1])
          sys.exit(0)

      else:
        svd[eventid]=[fname]+[gene]+myfeatures[2:3]+myfeatures[4:6]+[readsup]
        #print svd[eventid]
        #sys.exit(0)

  #print len(svd),svd
  return svd

def Combine_d(dl):
  '''Combine list of dict's into one'''
  newd,templ={},[]

  for x in dl: #list of dict's

    for k,v in x.items(): #one dict per SV run

      if k in newd:

        #if identical event from same BT-/PT- then combine! Multiple files per sample/patient in some cases where many libraries sequenced
        if newd[k][2]==v[0]: #same BT- sampleID
          #print "Same key and BT-: ",newd[k],'\n',v,'\n'
          newd[k]=newd[k] #[int(newd[k][0])+1]+newd[k][1:] #Add 1 to "total" count ie newd[k][0]

        else:
          #print "Same key different BT-: ",newd[k],'\n',v,'\n'
          templ=newd[k][2:]+v[0:3]+v[3:6] # v[0] is  is file/ID name (eg BT-248), only want it in v once
          newd[k]=[int(newd[k][0])+1]+[int(newd[k][1])+1]+templ #+1 to Total and Unique counts b/c event has same eventID but not same sampleID/gene

      else: #put 1 as first part of dict value, this will be used to get counts of events

        if len(v) > 7: #too many fields
          #tot=len(v)-7 #each entry past 7 is another gene name for the same event
          #newd[k]=[tot]+[1]+v #[0:8]
          print "Dict values >7, exiting", k,v
          sys.exit(0)
        else:
          newd[k]=[1]+[1]+v #

  #print newd
  #sys.exit(0)
  return newd

def Load_dict(ad,name1,gname1,evtype1):
  #name=sampleID (eg. BT-241), gname=gene name, evtype=event type (eg. deletion)

  if gname1 in ad:
    temp=[]

    for item in ad[gname1]:
      temp.append(item)
    temp.append([name1,evtype1])
    ad[gname1]=temp

  else:
    ad[gname1]=[[name1,evtype1]]

  return ad

def Get_stats(nd):
  '''Get counts of gene names, event types, and sample names'''
  clist,gened,eventd,sampled=[],{},{},{} #list of numbers to count occurences of SV events

  for k,v in nd.items(): #eventID [1, 1, BT-xy, fusion, inversion, ...]
    #print len(v),k,v

    for x in range(0,(len(v)-2)/6): #first 3 are nums and filename, then repeats of 6 items per event
      name,gname,evtype,totc,unqc=v[(x*6)+2],v[(x*6)+3],v[(x*6)+4],v[0],v[1] #file/patient ID, gene name, event type

      #unique event names for gened
      try:
        eventpos=v[(x*6)+5].split(':')[1] #eg. MT:7882-7882 => 7882-7882
        eventuq=evtype+"_"+eventpos #eg. del_7882-7882 these need to be unique to get correct counts for genes
      except IndexError:
        print "IndexError:",k,v
        sys.exit(0)
      #print name,gname,evtype,eventuq,totc,unqc

      if len(evtype)==3: #some files have event names del for delete etc.
        evtype=EventNames(evtype)

      if len(gname)>30: #I think this statement is to avoid the wrong type of data in the gname field
        print gname,name,'\n',v,'\n'
        sys.exit("Gene name should not be more than 15 characters")

      gened=Load_dict(gened,name,gname,evtype)
      eventd=Load_dict(eventd,name,evtype,gname) #use event type as key instead of gene name
      sampled=Load_dict(sampled,evtype,name,gname) #use sampeID as key instead of gene name

  myout=[gened,eventd,sampled]
  #print '\n'
  #print myout
  return myout #list of three dict's

def Output(outd,outn):
  outF=open(outn,'w')

  for l,m in outd.items():
    cnames,cevents=Counter(),Counter()
    names,events=[y[0] for y in m],[y[1] for y in m]

    try:
      cnames,cevents=Counter(names),Counter(events)

    except TypeError:
      print l,m,'\n',names,'\n',events

    ntotal,etotal=sum(cnames.values()),sum(cevents.values()) #total times event occurs, including multiple times in same sample
    nsamples,esamples=len(cnames.values()),len(cevents.values()) #in how many different samples does it occur
    #print l,ntotal,nsamples,cnames,etotal,esamples,cevents#' '.join(['{0}:{1}'.format(k,v) for k,v in cnames.items()])
    outF.write(l+"\t"+str(ntotal)+"\t"+str(nsamples)+"\t"+str(esamples)+'\n')
  outF.close()
    

if __name__ == "__main__":

  if len(sys.argv)>2:
    Files=sys.argv[1:] #

  else:
    sys.exit("Script works for multiple files only; exiting")

  SV=[]

  for i in range(len(Files)):
    sv=Read_in_files(Files[i])
    SV.append(sv)
  #SV is list of dict's with K=eventID V=event info
  c=Combine_d(SV)
  outf=open('SV.outfile.tsv', 'w')

  for k,v in c.items():
    #print k,"\t",v[0],"\t",v[1:]
    outf.write(k+"\t"+str(v[0])+"\t"+str(v[1])+"\t"+"\t".join(v[2:])+'\n')
  outf.close()

  outlist=Get_stats(c)
  z=0
  for x in outlist:
    z+=1
    outname="file."+str(z)+".txt"
    Output(x,outname)




  #for l,m in sorted(cc.items(), key=lambda x:x[1], reverse=True):
  #  if m > 5:
  #    print l,m


#      if len(event)<5: #some indels have only one gene name
      #can either use both gene names with '_' between, or use single gene names from col6
        #gene1,gene2=event[2].split(':')[0],"NA"
#      else:
        #gene1,gene2=event[2].split(':')[0],event[3]
