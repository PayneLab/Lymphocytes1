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
    file = download_file(download_to_path="data/proteinGroups.txt", url_file_path="data/proteinGroups_url.txt")
        
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
    try: df = pd.read_csv(file, sep='\t', header=0, index_col=0, usecols=final_col)
    except: df = pd.read_excel(file, sep='\t', header=0, index_col=0, usecols=final_col)

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
        url_file = open(url_file_path, 'r')
        url = url_file.read().strip()
        url_file.close()
        
        if path.exists(password_file_path):
            password_file = open(password_file_path, 'r')
            password = password_file.read().strip()
            password_file.close()
        else:
            print("MISSING PASSWORD FILE")
        
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