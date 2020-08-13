import pandas as pd
import requests
import os.path
from os import path

#Data loader functions belong here. This is where
#  information about the data files is found.

def load_max_quant():
    #Takes a file and returns a dataframe.
    #    file: the file path to read from
    #    The rest of the paramters are used to select the columns.
    #    By default, it will look for ones starting with 'Reporter intensity'
    #        that do not contain 'count' or 'corrected' and use the 'Protein IDs'
    #        column as the indecies. These will be the raw intensity values.
    file = download_file(download_to_path="data/proteinGroups-Sub1-ungrouped-not_LFQ_norm.txt", url_file_path="data/proteinGroups_url.txt")
        
    prefix="Intensity"
    contains=["_"]

    with open(file, 'r') as _file:
        line = _file.readline().strip()
        headings = line.split('\t')
    headings = [i.strip('"') for i in headings]
    if prefix:#filter by columns beginning in prefix
        headings = [i for i in headings if i.startswith(prefix)]
    for req in contains:
        headings = [i for i in headings if req in i]
    
    final_col = ["Protein IDs"]
    for i in headings: final_col.append(i)
    df = pd.read_csv(file, sep='\t', header=0, index_col=0)
    df = df.drop(df[df['Potential contaminant'] == '+'].index)
    df = df.drop(df[df.Reverse == '+'].index)
    df = df.drop(df[df['Only identified by site'] == '+'].index)
    df = df[headings]

    return df

def load_FragPipe(month='June', contains=['Subject1']):
    #Takes a file and returns a dataframe.
    #    file: the file path to read from
    #    The rest of the paramters are used to select the columns.
    #    By default, it will look for ones ending with 'Total intensity'
    #        that do not contain 'count' or 'corrected' and use the 'Protein IDs'
    #        column as the indecies. These will be the raw intensity values.
    file_name="data/combined_protein_{0}_FP.tsv".format(month)
    url_file_path="data/combined_protein_{0}_FP_url.txt".format(month)
    file = download_file(download_to_path=file_name, url_file_path=url_file_path)
    if file==1:
        print("Error with file download.")
        return False
        
    suffix="Total Intensity"

    with open(file, 'r') as _file:
        line = _file.readline().strip()
        headings = line.split('\t')
    headings = [i.strip('"') for i in headings]
    if suffix:#filter by columns beginning in prefix
        headings = [i for i in headings if i.endswith(suffix)]
    for req in contains:
        headings = [i for i in headings if req in i]
    
    final_col = ["Protein IDs"]
    for i in headings: final_col.append(i)
    df = pd.read_csv(file, sep='\t', header=0, index_col=3)
    #df = df.drop(df[df['Potential contaminant'] == '+'].index)
    #df = df.drop(df[df.Reverse == '+'].index)
    #df = df.drop(df[df['Only identified by site'] == '+'].index)
    df = df[headings]

    return df

def download_file(download_to_path="data/datafile.txt", url_file_path="data/url.txt", 
                  password_file_path="data/password.txt", redownload=False):
    """Download a file from a given url to the specified location.
    Parameters:
    path (str): The path to the file to save the file to on the local machine.
    Returns:
    str: The path the file was downloaded to.
    """
        
    if redownload or path.exists(download_to_path) == False: #If the file has been downloaded, or the user wants to update, download the file
        if path.exists(url_file_path):
            url_file = open(url_file_path, 'r')
            url = url_file.read().strip()
            url_file.close()
        else: 
            print("MISSING URL FILE")
            return 1
        
        if path.exists(password_file_path):
            password_file = open(password_file_path, 'r')
            password = password_file.read().strip()
            password_file.close()
        else:
            print("MISSING PASSWORD FILE")
            return 1
        
        for i in range(2):

            with requests.Session() as session: # Use a session object to save cookies
                # Construct the urls for our GET and POST requests
                get_url = url
                post_url = get_url.replace("https://byu.box.com/shared", "https://byu.app.box.com/public")

                # Send initial GET request and parse the request token out of the response
                get_response = session.get(get_url)
                soup = bs4.BeautifulSoup(get_response.text, "html.parser")
                token_tag = soup.find(id="request_token")
                token = token_tag.get("value")

                # Send a POST request, with the password and token, to get the data
                payload = {
                    'password': password,
                    'request_token': token}
                response = session.post(post_url, data=payload)

                with open(download_to_path, 'wb') as dest:
                    dest.write(response.content)
    return download_to_path

def load_fasta():
    file="data/uniprot-filtered-proteome_3AUP000005640_reviewed_human.fasta"
    #file is formated:
    #>sp|Q96IY4|CBPB2_HUMAN Carboxypeptidase B2 OS=Homo sapiens OX=9606 GN=CPB2 PE=1 SV=2
    #MKLCS...
    headings = {}
    with open(file) as f:
        for line in f:
            if line.startswith('>'):#header line
                ID = line.split('|')[1]
                name=line.split('|')[2].split('=')[0].strip('OS')
                headings[ID]=name
    headings = pd.Series(list(headings.values()), index=headings.keys())
    
    return headings

        
def names_max_quant():
    file = download_file(download_to_path="data/proteinGroups.txt", url_file_path="data/proteinGroups_url.txt")
    df = pd.read_csv(file, sep='\t', header=0, index_col=0, usecols=['Protein IDs','Gene names','Fasta headers'])
    
    """by_group={}
    for ID in df.index:
        if ';' in ID:
            headers = df['Fasta headers'][ID].split(';')
            names=''
            for i in headers:
                if i !='':
                    n = i.split('|')[2].split('=')[0].strip('OS')
                    names=names+n+';'
            names=names[:-1]
            by_group[ID]=names
        else: 
            i = df['Fasta headers'][ID]
            if str(i) != 'nan':
                i.split('|')[2].split('=')[0].strip('OS')
                by_group[ID]=i
    by_group = pd.Series(list(by_group.values()), index=by_group.keys())

    df['Names']=by_group"""

    return df
    
    
                