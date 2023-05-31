from matplotlib.pyplot import axis
# from consumers.datasets.scraping.chembl_scarping import chembl
# from consumers.datasets.scraping.BindingDB_scarping import bindingDB
# from consumers.datasets.scraping.pubchem_scraping import pubchem
# from consumers.datasets.utils.ID_conversion import Conversion
import pandas as pd, requests, time, traceback
from tqdm import tqdm
import pandas as pd,numpy as np
from consumers.datasets.processing.salt_remover import clean_smiles
import structure_curation as cur
from rdkit import Chem
import process_smiles as ps
from rdkit import Chem
##### key words used to distingush the assay as cell based ######
Assay_classification_keywords=['cell growth','growth inhibit','proliferation','viability','MTT','MTS','cytotoxicity','cells']

def automate_ligand_scrape(Uniport_ID,activity,assay_type):
    import urllib.parse
    import urllib.request
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
    'from': 'ACC+ID',
    'to': 'GENENAME',
    'format': 'tab',
    'query': Uniport_ID
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    print(req)
    with urllib.request.urlopen(req) as f:
       response = f.read()

    decoded_data= response.decode('utf-8')
    decode_list=decoded_data[8:-1].split('\t')
    #print(decode_list)
    protein_name=decode_list[1]

    Chembl_data,Chembl_info=chembl(protein_name,Uniport_ID)
    BindingDB_data,bindingDB_info=bindingDB(Uniport_ID)
    pubchem_data,pubchem_info=pubchem(Uniport_ID)

    info={"Chembl":Chembl_info,"BindingDB":bindingDB_info,"Pubchem":pubchem_info}
    print(info)
    
    combine_data=BindingDB_data.merge(pubchem_data,how='outer',on=['PubChem CID','std_type','std_value'])
    combine_data['assay']=combine_data['assay'].astype(str)


    for sim in range(0,5):
        for i,r in tqdm(combine_data.iterrows()):
            if not pd.isna(r['SMILES']) and not pd.isna(r['assay']):
                continue
            elif pd.isna(r['SMILES']) and not r['assay']=='nan':
                try:
                    cid_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+str(r['PubChem CID'])+"/property/CanonicalSMILES/JSON"
                    re=requests.get(cid_url)
                    output_cid=re.json()
                    s=output_cid['PropertyTable']['Properties'][0]['CanonicalSMILES']
                    combine_data.loc[i,'SMILES']=s
                except:
                    time.sleep(15)
                    traceback.print_exc()
                    continue
            elif not pd.isna(r['SMILES']) and r['assay']=='nan':
                try: 
                    Pubchem_url='https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22cid%22:%22'+str(r['PubChem CID'])+'%22}]},%22order%22:[%22acvalue,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22CID_'+str(r['PubChem CID'])+'_bioactivity%22,%22nullatbottom%22:1}'
                    pubchem_data=pd.read_csv(Pubchem_url)
                    combine_data.loc[i,'assay']=list(pubchem_data[(pubchem_data["repacxn"]==Uniport_ID)&(pubchem_data['acname']==r['std_type'])]['aidname'])[0]
                except:
                    continue
        if combine_data['SMILES'].count()!=len(combine_data):
            break
    for sim in range(0,5):
        for i,r in tqdm(combine_data.iterrows()):
            if not pd.isna(r['SMILES']) and not pd.isna(r['assay']):
                continue
            elif pd.isna(r['SMILES']) and not r['assay']=='nan':
                    continue
            elif not pd.isna(r['SMILES']) and r['assay']=='nan':
                try: 
                    Pubchem_url='https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22cid%22:%22'+str(r['PubChem CID'])+'%22}]},%22order%22:[%22acvalue,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22CID_'+str(r['PubChem CID'])+'_bioactivity%22,%22nullatbottom%22:1}'
                    pubchem_data=pd.read_csv(Pubchem_url)
                    combine_data.loc[i,'assay']=list(pubchem_data[(pubchem_data["repacxn"]==Uniport_ID)&(pubchem_data['acname']==r['std_type'])]['aidname'])[0]
                except:
                    time.sleep(10)
                    continue
        if combine_data['assay'].count()!=len(combine_data):
            break


    combine_data=combine_data.drop('PubChem CID',axis=1)
    combine_data=combine_data.dropna(axis=0)

    All_three_DB_data=pd.concat([Chembl_data,combine_data],axis=0)


    All_three_DB_data['assay']=[np.nan if i=='nan' else i for i in All_three_DB_data['assay'] ]

    All_three_DB_data=All_three_DB_data.dropna(axis=0)


    ##### Classifying the data as Mutant/Wild based on the assay info ######
    All_three_DB_data['Mutant Data']=["Wild" if i.find('mutant')  else "Mutant" for i in All_three_DB_data['assay'] ]


    ##### Classifying the data as Biochemical or cell based based on the assay info ######
    All_three_DB_data['Type of assay']=["Cell assay" if any(map(i.__contains__, Assay_classification_keywords)) else "Biochemical Assay" for i in All_three_DB_data['assay'] ]

    
    ##### Canonical smiles ######
    All_three_DB_data['SMILES']=[Chem.CanonSmiles(i) for i in All_three_DB_data['SMILES']]


    All_three_DB_data=All_three_DB_data[All_three_DB_data['std_type']==activity]

    All_three_DB_data_cell_based=All_three_DB_data[All_three_DB_data['Type of assay']=='Cell assay']
    All_three_DB_data_bio_based=All_three_DB_data[All_three_DB_data['Type of assay']=='Biochemical Assay']
    print("data after spliting  cell based :",All_three_DB_data_cell_based.shape)
    print("data after spliting  bio based :",All_three_DB_data_bio_based.shape)


    smiles=[]
    activity_data=[]
    All_three_DB_data_cell_based=All_three_DB_data_cell_based.drop_duplicates(['SMILES','std_value'])
    All_three_DB_data_cell_based['std_value']=All_three_DB_data_cell_based['std_value'].astype(float)
    print("data after spliting  cell based after drop_duplicates:",All_three_DB_data_cell_based.shape)
    for i in All_three_DB_data_cell_based.groupby('SMILES').groups.items():
        sum1=[]
        for inx in list(i[1]):
            
            sum1.append(All_three_DB_data_cell_based.loc[inx,'std_value'])
            
    #     mean_1=mean(sum1)
        mean_1=sum(sum1) / len(sum1)
        smiles.append(i[0])
        activity_data.append(mean_1)
    All_three_DB_data_cell_based=pd.DataFrame(zip(smiles,activity_data),columns=['SMILES','ACTIVITY'])
    print("data after spliting  cell based after new:",All_three_DB_data_bio_based.shape)


    smiles,activity_data=[],[]
    All_three_DB_data_bio_based=All_three_DB_data_bio_based.drop_duplicates(['SMILES','std_value'])
    All_three_DB_data_bio_based['std_value']=All_three_DB_data_bio_based['std_value'].astype(float)

    print("data after spliting  bio based after drop_duplicates:",All_three_DB_data_bio_based.shape)

    for i in All_three_DB_data_bio_based.groupby('SMILES').groups.items():
        sum1=[]
        for inx in list(i[1]):
            
            sum1.append(All_three_DB_data_bio_based.loc[inx,'std_value'])
            
    #     mean_1=mean(sum1)
        mean_1=sum(sum1) / len(sum1)
        smiles.append(i[0])
        activity_data.append(mean_1)
    All_three_DB_data_bio_based=pd.DataFrame(zip(smiles,activity_data),columns=['SMILES','ACTIVITY'])
    

    All_three_DB_data_bio_based=curate_smiles(All_three_DB_data_bio_based)
    All_three_DB_data_bio_based['Assay_type']='Biochemical'
    All_three_DB_data_cell_based=curate_smiles(All_three_DB_data_cell_based)
    All_three_DB_data_cell_based['Assay_type']='Cell Based'
    info[f"assay_info for_{activity}"]={"num of biochemical assay":len(All_three_DB_data_bio_based),"num of cell_based assay":len(All_three_DB_data_cell_based)}
    if assay_type.lower()=='biochemical':
        return All_three_DB_data_bio_based,info
    if assay_type.lower()=='cell_based':
        return All_three_DB_data_cell_based,info
    else:
        return pd.concat([All_three_DB_data_bio_based,All_three_DB_data_cell_based],axis=0),info

def curate_smiles(data):
    data['SMILES']=clean_smiles(data['SMILES'].tolist())


    data_cur = cur.Curator()

    for i, row in data.iterrows():
        smi = row['SMILES']
        data_cur.get_rdkit_mol(smi)
        substance_type, sanitized_smiles = data_cur.filter_smiles()   
        data.loc[data.index == i,'structure_curated'] = sanitized_smiles
        data.loc[data.index == i,'substance_type_name'] = substance_type
    
    return data

def input_data_type(input_data,data,assay):
    
    try:
      #conversion from ensemble id to uniport id  
      if (input_data[0:3]=="ENS"):
          i = Conversion('gkiju@gmail.com')
          entrez_id=i.convert_ensembl_to_entrez(input_data)
          UniProtId=i.convert_entrez_to_uniprot(entrez_id)
        #print(UniProtId)
          data_return =automate_ligand_scrape(UniProtId,data,assay)
          
       #conversion from gene id to uniport id
      elif input_data.isdigit():
          i = Conversion('gkiju@gmail.com')
        
          UniProtId = i.convert_entrez_to_uniprot(input_data)
        #print(UniProtId)
          data_return =automate_ligand_scrape(UniProtId,data,assay)
          
        
      else:
          data_return = automate_ligand_scrape(input_data,data,assay)
        #print(UniProtId)
      return data_return
    except Exception as e:
        traceback.print_exc()
        print("Enter valid input-ensembl,entrez and uniport ID")
        error="Enter valid input-ensembl,entrez and uniport ID"
        return (error)

if __name__=="__main__":
    c=automate_ligand_scrape('Q13546','IC50','biochemical')
    print(c)
    