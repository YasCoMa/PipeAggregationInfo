import pandas as pd
from sklearn import preprocessing
import scipy
import os
from Bio import Entrez
import numpy as np

class Expression_pattern_evaluation:

    def check_existence_pair(self, pair, file):
        p=pair.split(",")
        foundp1=False
        foundp2=False
        
        os.system("grep -ic "+p[0]+" "+file+" > count.txt")
        f=open("count.txt","r")
        for line in f:
            l=line.replace("\n","")
            if(int(l)>0):
                foundp1= True
        f.close()

        os.system("grep -ic "+p[1]+" "+file+" > count.txt")
        f=open("count.txt","r")
        for line in f:
            l=line.replace("\n","")
            if(int(l)>0):
                foundp2= True
        f.close()

        return (foundp1 and foundp2)

    def build_normalize_correlation(self, file):
        # Loading matrix
        m=pd.DataFrame.from_csv(file, sep='\t', header=None)
        
        # Getting list of proteins
        proteins=list(m.index)

        # Scaling
        x = m.values #returns a numpy array
        min_max_scaler = preprocessing.MinMaxScaler()
        x_scaled = min_max_scaler.fit_transform(x)
        df = pd.DataFrame(x_scaled)

        # Transposing
        df=df.T

        # Calculating correlation
        mc=df.corr(method='pearson')

        return mc, proteins

    def getindex_protein(self, protein, listp):
        i=0
        for p in listp:
            if(p==protein):
                return i
            i+=1

        return -1

    def process_pairs(self, folder, interactome_file):
        f=open(folder+"expPattern_evaluation_pairs.tsv","w")
        f.close()

        pairs=[]
        f=open(folder+interactome_file, "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(not [l[0],l[1]] in pairs):
                pairs.append([l[0],l[1]])
        f.close()

        for exps in os.listdir('./data_matrices_count'):
            for matfiles in os.listdir('./data_matrices_count/'+exps):
                if(matfiles.find(".tsv")!=-1):
                    #print('./data_matrices_count/'+exps+'/'+matfiles)
                    if(self.check_existence_pair(l[0]+","+l[1], './data_matrices_count/'+exps+'/'+matfiles)):
                        try:
                            mc, proteins=self.build_normalize_correlation('./data_matrices_count/'+exps+'/'+matfiles)
                            for p in pairs:
                                indexp1=self.getindex_protein(p[0],proteins)
                                indexp2=self.getindex_protein(p[1],proteins)
                                correlation=mc.iloc[indexp1, indexp2]
                                #print(correlation)

                                body=[exps+'/'+matfiles,p[0],p[1],str(correlation)]
                                with open(folder+"expPattern_evaluation_pairs.tsv","a") as gf:
                                    gf.write(("\t".join(body))+"\n")
                        except:
                            pass

    def process_correlations(self, folder):
        d={}
        f=open(folder+"expPattern_evaluation_pairs.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            
            key=l[1]+","+l[2]
            if( not key in d.keys() ):
                d[key]={}
            
            value=0
            if(l[3]!="nan"):
                value=float(l[3])
            d[key][l[0]]=value
        f.close()

        if(os.path.isdir(folder+"ranking_expression_profile")):
            os.system("rm "+folder+"ranking_expression_profile/*")
        else:
            os.system("mkdir "+folder+"ranking_expression_profile")

        for p in d.keys():
            a=d[p]
            sorted_ = sorted(a, key=a.get, reverse=True)
            f=open(folder+"ranking_expression_profile/"+p.replace(",","-")+"_top50.tsv","w")
            for v in sorted_:
                f.write("%s\t%.2f\n" %(v, a[v]) )
            f.close()

    def select_search_subset_pairs(self, folder, selected_pairs):
        if(os.path.isdir(folder+"filtered_pairs")):
            os.system("rm "+folder+"filtered_pairs/*")
        else:
            os.system("mkdir "+folder+"filtered_pairs")

        f=open(folder+selected_pairs, "r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(os.path.isfile(folder+"ranking_expression_profile/"+l[0]+"-"+l[1]+"_top50.tsv")):
                os.system("cp "+folder+"ranking_expression_profile/"+l[0]+"-"+l[1]+"_top50.tsv "+folder+"filtered_pairs/"+l[0]+"-"+l[1]+"_top50.tsv")
            
            if(os.path.isfile(folder+"ranking_expression_profile/"+l[1]+"-"+l[0]+"_top50.tsv")):
                os.system("cp "+folder+"ranking_expression_profile/"+l[1]+"-"+l[0]+"_top50.tsv "+folder+"filtered_pairs/"+l[1]+"-"+l[0]+"_top50.tsv")
        f.close()

        print("Results for the selected pairs can be found in the folder filtered_pairs/")

class Running_config:
    def run(self, args):
        if(args.running_type in [1,2]):
            a=Expression_pattern_evaluation()

            if(args.running_type==1):
                if(args.folder!="" and args.interactome_file!=""):
                    if(os.path.isdir(args.folder) and os.path.isfile(args.folder+args.interactome_file) ):
                        a.process_pairs(args.folder, args.interactome_file)
                        a.process_correlations(args.folder)
                    else:
                        print("Error: There are invalid folder or files, some of them were not found")
                else:
                    print("Error: There are missing parameters")

            if(args.running_type==2):
                if(args.folder!="" and args.selected_pairs_file!=""):
                    if(os.path.isdir(args.folder) and os.path.isfile(args.folder+args.selected_pairs_file) ):
                        a.select_search_subset_pairs(args.folder, args.selected_pairs_file)
                    else:
                        print("Error: There are invalid folder or files, some of them were not found")
                else:
                    print("Error: There are missing parameters")
        else:
            print("Error: This option is invalid")

if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description='Pipeline to find gene expression experiments where PPIs are highly correlated', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-rt", "--running_type", action="store", help="Required - 1 - Make the process of finding the experiments and ranking them by correlation\n2 - Select pairs that were already processed and ranked making a separated folder of interest", type=int)
    parser.add_argument("-fo", "--folder", action="store", help="Folder to store the files (use the folder where the other files required can be found, ex.: /home/user/experiment/ )")
    parser.add_argument("-if", "--interactome_file", action="store", help="File with all PPIs (two columns with uniprot identifiers in tsv format) ")
    parser.add_argument("-spf", "--selected_pairs_file", action="store", help="File with PPIs of interest (two columns with uniprot identifiers in tsv format) ")
    args = parser.parse_args()
    r=Running_config()
    r.run(args)
