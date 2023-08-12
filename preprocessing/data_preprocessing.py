# Validate correlation between genes of the predicted ppis 
# Filtering experiments of NCBI GEO DataSets database using (cancer) AND "Homo sapiens"[porgn:__txid9606] and selecting the "Expression profiling by high throughput sequencing" as study type
# The filter applied returned 5291 experiments

from selenium import webdriver
from selenium.webdriver.common.by import By
options = webdriver.ChromeOptions()
options.add_argument("--incognito")
options.add_argument("--headless")

import os, time
from Bio import Entrez


class Process_counts_file:
    def _init__(self, folder, output_folder):
        self.folder = folder
        self.fout = output_folder
        
    # filter files with "count" in the name
    # build dictionaries of mapping for hgnc and ensembl
    # verify the column that corresponds to one of the mappings above to find gene identifier
    # verify normalization
    # test if at least three lines splitted by , or ; or \t or space in txt files have the same number of columns to discover separator in txt
    def build_dictionaries_mapping(self):
        hgnc={}
        f=open("mapping_hgnc_uniprot.txt","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            hgnc[l[0]] = l[1]
        f.close()

        ensembl={}
        f=open("mapping_ensembl-uniprot.txt","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(l[1]!=""):
                ensembl[l[0]] = l[1]
        f.close()

        return [hgnc, ensembl]

    def get_separator(self, first_4_lines):
        separator=""
        sep_candidates=[",",";","\t"]
        gt=0
        for s in sep_candidates:
            count=len(first_4_lines[0].split(s))
            if(gt<count):
                separator = s
                gt=count
        return separator

    def get_position_geneName(self, first_4_lines, separator):
        mapping=self.build_dictionaries_mapping()
        hgnc=mapping[0]
        ensembl=mapping[1]

        #indexs=[]
        #final_dict=[]
        #for i in range(1,len(first_4_lines)):
        l=first_4_lines[1].split(separator)
        finalIndex=0
        dict_=""
        index=0
        for w in l:
            if(w.find("|")!=-1):
                w=w.split("|")[0]
            if(w.find("ENSG")!=-1):
                w=w.split(".")[0]
            
            if(w in hgnc.keys()):
                finalIndex=index
                dict_='hgnc'
            if(w in ensembl.keys()):
                finalIndex=index
                dict_='ensembl'

            index+=1

        return [finalIndex, dict_]

    def get_positions_counts(self, first_4_lines, separator):
        header=first_4_lines[0]
        h=header.split(separator)
        i=0
        positions_array=[]
        for i in range(1,len(first_4_lines)):
            positions=[]
            sentence=first_4_lines[i].split(separator)
            index=0
            for s in sentence:
                try:
                    a=int(s)
                    if(h[index].lower().find("strand")==-1 and h[index].lower().find("fpkm")==-1 and h[index].lower().find("start")==-1 and h[index].lower().find("end")==-1 and h[index].lower().find("transcript")==-1):
                        positions.append(index)
                except:
                    pass
                index+=1
            positions_array.append(positions)
        
        equals=0
        positions=[]
        final_positions=positions_array[0]
        for pa in positions_array:
            if(final_positions==pa):
                equals+=1
        if(equals==len(positions_array)):
            positions=final_positions
        
        return positions

    def rewrite_matrix_count(self, file, acc, filename):
        if( not os.path.isdir( f"{self.fout}data_matrices_count") ):
            os.system( f"mkdir {self.fout}data_matrices_count")
        
        first_4_lines=[]
        c=0
        print(file)
        try:
            f=open(file,"r", encoding="utf-8")
            for line in f:
                l=line.replace('"',"").replace("\n","")
                if(c<4):
                    first_4_lines.append(l)
                else:
                    break
                c+=1
            f.close()
        except:
            f=open(file,"r", encoding="cp1252")
            for line in f:
                l=line.replace('"',"").replace("\n","")
                if(c<4):
                    first_4_lines.append(l)
                else:
                    break
                c+=1
            f.close()

        if (file.split(".")[-1].lower()=="tsv"):
            separator="\t"
        elif (file.split(".")[-1].lower()=="csv"):
            separator=","
        else:
            separator="\t"
            #separator=self.get_separator(first_4_lines)
            
        positions=[]
        mapping=self.build_dictionaries_mapping()

        mapinfo=self.get_position_geneName(first_4_lines,separator)
        positions.append(mapinfo[0])
        # treat mapping identifier later
        #if(mapinfo[1]=='hgnc'):
        #    mapping=mapping[0]
        #else:
        #    mapping=mapping[1]

        count_positions=self.get_positions_counts(first_4_lines, separator)
        positions+=count_positions

        matrix_name= f"{self.fout}data_matrices_count/"+acc+"/newmcount-"+filename.replace("."+filename.split(".")[-1],"")
        
        c=0
        try:
            f=open(file,"r", encoding="utf-8")
            for line in f:
                l=line.replace('"',"").replace("\n","").split(separator)
                if(c>0):
                    text=[]
                    i=0
                    ok=True
                    for info in l:
                        if(i in positions):
                            if(i==0):
                                info=info.split(".")[0]
                                cm=0
                                for m in mapping:
                                    info=info.split(".")[0]
                                    if(info in m.keys()):
                                        info=m[info]
                                        cm+=1
                                        break
                                if(cm==0):
                                    ok=False
                            text.append(info)
                        i+=1
                    #print(l[0],ok)
                    if(ok):
                        if(len(text)>2):
                            #print(text)
                            if( not os.path.isdir( f"{self.fout}data_matrices_count/"+acc) ):
                                os.system( f"mkdir {self.fout}data_matrices_count/"+acc)
                            with open(matrix_name+".tsv","a") as g:
                                g.write(("\t".join(text))+"\n")
                #else:
                #    print(l)
                c+=1
            f.close()
            print("ok")
        except:
            f=open(file,"r", encoding="cp1252")
            for line in f:
                l=line.replace('"',"").replace("\n","").split(separator)
                if(c>0):
                    text=[]
                    i=0
                    ok=True
                    for info in l:
                        if(i in positions):
                            if(i==0):
                                info=info.split(".")[0]
                                cm=0
                                for m in mapping:
                                    if(info in m.keys()):
                                        info=m[info]
                                        cm+=1
                                        break
                                if(cm==0):
                                    ok=False
                            text.append(info)
                        i+=1
                    #print(l[0],ok)
                    if(ok):
                        if(len(text)>2):
                            if( not os.path.isdir( f"{self.fout}data_matrices_count/"+acc) ):
                                os.system( f"mkdir {self.fout}data_matrices_count/"+acc)
                            #print(text)
                            with open(matrix_name+".tsv","a") as g:
                                g.write(("\t".join(text))+"\n")
                #else:
                #    print(l)
                c+=1
            f.close()
            
class Retrieve_rnaseq_experiments:
    def _init__(self, folder, output_folder, filter_file_folder):
        self.folder = folder
        self.fout = output_folder
        self.filter_file_folder = filter_file_folder

    def get_info_gds_result(self):
        filter_file_folder = self.filter_file_folder
        
        if( not os.path.isdir( f"{self.fout}geo_data") ):
            os.system(f"mkdir {self.fout}geo_data")
        
        for fi in os.listdir(self.folder+'/'+filter_file_folder):
            filter_file = self.folder+'/'+filter_file_folder+'/'+fi
            
            appendix=filter_file.replace(filter_file.split(".")[-1],"").replace(".","")
            f=open( f"{self.fout}data_experiments-"+appendix+".tsv","w")
            f.close()

            info=[]
            c=0
            f=open(filter_file,"r")
            for line in f:
                if(line=="\n"):
                    if(len(info)>=5):
                        text=""
                        for i in info:
                            text+=i+"\t"
                        
                        links = self.get_link_supplementary(info[4])

                        if(len(links)>0):
                            text+=" ".join(links)
                            with open( f"{self.fout}data_experiments-"+appendix+".tsv","a") as g:
                                g.write(text+"\n")
                    info=[]
                    c=0
                else:
                    l=line.replace("\n","")
                    if c==0 :
                        print(l)
                        info.append(l.split(". ")[0])
                        info.append(l.split(". ")[1])
                    
                    if c==4 or c==5:
                        if(l.find("Platform")!=-1):
                            platforms=[]
                            things=l.split(" ")
                            for t in things:
                                if(t.find("GPL")!=-1):
                                    platforms.append(t)
                            if(len(platforms)>0):
                                info.append(" ".join(platforms))

                            i=0
                            for t in things:
                                try:
                                    samples=int(t)
                                    if(i!=0):
                                        info.append(str(samples))
                                except:
                                    pass
                                i+=1
                    
                    if c==5 or c==6 or c==7 :
                        if(l.find("Series\t\tAccession:")!=-1):
                            info.append(l.split("\t")[2].replace("Accession: ",""))
                            info.append(l.split("\t")[3].replace("ID: ",""))
                    
                    c+=1

            f.close()

    def get_link_supplementary(self, accession):
        links=[]
        try:
            b = webdriver.Chrome(executable_path='./chromedriver',chrome_options=options)
            b.get("https://www.ncbi.nlm.nih.gov/geo/download/?acc="+accession)
            for f in b.find_elements_by_css_selector("div.section a"):
                a=f.get_attribute("href")
                if(a.find(accession+"/suppl/"+accession+"_")!=-1 and a.lower().find("count")!=-1 and (a.lower().find(".csv")!=-1 or a.lower().find(".txt")!=-1 or a.lower().find(".tsv")!=-1)):
                    if( not os.path.isdir( f"{self.fout}geo_data/"+accession) ):
                        os.system( "mkdir {self.fout}geo_data/"+accession)
                    if( not os.path.isfile( f"{self.fout}geo_data/"+accession+"/"+a.split("/")[-1]) ):
                        os.system("wget "+a+" -P {self.fout}geo_data/"+accession+"/")
                    links.append(a)
        except:
            pass
        return links

    def get_summary(self):
        
        ids=[]
        for fi in os.listdir(self.fout):
            if( fi.startswith('data_experiments') ):
                f=open(f"{self.fout}{fi}.tsv","r")
                for line in f:
                    l=line.replace("\n","").split("\t")
                    ids.append(l[5])
                f.close()

        Entrez.email = 'ycfrenchgirl2@gmail.com'
        Entrez.api_key = "4543094c8a41e6aecf9a1431bff42cfac209"
        f=open( f"{self.fout}entries_summary.tsv","w")
        f.close()
        c=0
        while(c<len(ids)):
            if((c+100) > len(ids)):
                middle=ids[c:]
            else:
                middle=ids[c:(c+100)]
            query=",".join(middle)
            handle = Entrez.esummary(db='gds', id=query, rettype="summary", retmode="text" )
            results = Entrez.read(handle)
            c+=100

            for r in results:
                with open( f"{self.fout}entries_summary.tsv","a") as g:
                    g.write(r['Id']+"\t"+r['summary']+"\n")  
            time.sleep(5)  

    def read_process_matrix_count(self):
        if( not os.path.isdir( f'{self.fout}geo_data/') ):
            os.system(f'mkdir {self.fout}geo_data/')
            
        os.system( f"rm -rf {self.fout}data_matrices_count/*")

        for fi in os.listdir(self.fout):
            if( fi.startswith('data_experiments') ):
                ids=[]
                f=open(f"{self.fout}{fi}.tsv","r")
                for line in f:
                    l=line.replace("\n","").split("\t")
                    accession = l[4]
                    ids.append(accession)
                    
                    if( not os.path.isdir( f"{self.fout}geo_data/"+accession) ):
                        os.system( "mkdir {self.fout}geo_data/"+accession)
                        for count_file in l[-1].split(" "):
                            if( count_file != ''):
                                os.system("wget "+a+" -P {self.fout}geo_data/"+accession+"/")
                f.close()

                p = Process_counts_file(self.folder, self.fout)
                for id in ids:
                    for root, folders, files in os.walk("geo_data/"+id):
                        for fi in files:
                            filename = f"{self.fout}geo_data/"+id+"/"+fi.replace(".gz","")
                            if(fi.split(".")[-1]=="gz" and fi.lower().find("count")!=-1):
                                os.system( f"gunzip {self.fout}geo_data/"+id+"/"+fi)
                            if(os.path.isfile(filename) and fi.lower().find("count")!=-1):
                                p.rewrite_matrix_count(filename,id,fi.replace(".gz",""))
        
        
class Run_preprocessing:
    def run(self, folder, folder_filter_gds):
        if( os.path.isdir(folder) ):
            if(folder[-1]=='/'):
                folder=folder[:-1]
            if( os.path.isdir(folder+'/'+folder_filter_gds) ):
                fout = folder+'/resultsAggExp/'
                if( not os.path.isdir(fout) )
                    a = Retrieve_rnaseq_experiments(folder, fout, filter_gds)
                    a.get_info_gds_result()
                    
                a = a.read_process_matrix_count()
            else:
                print("Error - Folder with the lists of filtered ncbi results for geo datasets not found") 
        else:
            print("Error - working directory does not exist")     

# uncompress the files and see those that have a csv or tsv pattern as a table and rewrite normalizing
# except those that already are normalized

import sys
working_directory = sys.argv[1]
name_folder_filter_gds = sys.argv[2]

a = Run()
a.Run_preprocessing(working_directory, name_folder_filter_gds)

