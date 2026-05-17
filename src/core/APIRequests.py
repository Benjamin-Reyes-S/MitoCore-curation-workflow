import pandas as pd
import time
from urllib import request
import requests
    

def get_uniprot_from_hgnc(df):


    session = requests.Session()

    ids = {'model_id': [], 'hgnc': [], 'uniprot': []}

    for _, row in df.iterrows():

        model_id = row['model_id']
        hgnc_id = row['hgnc'][0]
        uniprot_acc = None

        url = f"https://rest.genenames.org/fetch/hgnc_id/{hgnc_id}"

        try:
            resp = session.get(
                url,
                headers={"Accept": "application/json"},
                timeout=10
            )

            if resp.ok:
                docs = resp.json()["response"]["docs"]
                if docs and "uniprot_ids" in docs[0]:
                    uniprot_acc = docs[0]["uniprot_ids"][0]

        except requests.exceptions.RequestException:
            print(f"Request failed for {hgnc_id}")

        ids['model_id'].append(model_id)
        ids['hgnc'].append(hgnc_id)
        ids['uniprot'].append(uniprot_acc)

    return pd.DataFrame(ids)


def get_hsa_from_uniprot(df):
    # ensure no whitespace or version suffixes
    df['uniprot'] = df['uniprot'].astype(str).str.strip()

    output_dict = {'uniprot': [], 'hsa': []}
    batch_size = 10
    df = list(df['uniprot'])

    for i in range(0, len(df), batch_size):
        batch = df[i:i + batch_size]
        url = "https://rest.kegg.jp/conv/hsa/" + "+".join([f"uniprot:{uid}" for uid in batch]) #/conv/genes/

        try:
            with request.urlopen(url) as f:
                response = f.read()

            lines = response.decode('utf-8').splitlines()
            print(f"Requested {len(batch)} IDs; got {len(lines)} mappings")

            if lines:
                for line in lines:
                    if '\t' in line:
                        uniprot_entry, hsa_id = line.strip().split('\t')
                        uniprot = uniprot_entry.replace("up:", "").strip()
                        hsa = hsa_id.replace("hsa:", "").strip()

                        output_dict['uniprot'].append(uniprot)
                        output_dict['hsa'].append(hsa)

        except Exception as e:
            print(f"Failed on batch starting with {batch[0]}: {e}")
            for uid in batch:
                output_dict['uniprot'].append(uid)
                output_dict['hsa'].append(None)

        time.sleep(3)
    return  pd.DataFrame.from_dict(output_dict) 


def get_kegg_from_hsa():
    # Download the full KEGG gene-to-KO mapping for human
    url = "http://rest.kegg.jp/link/ko/hsa"
    response = request.urlopen(url)
    response_content = response.read().decode('utf-8')  # Decode the response content

    # Parse the mapping
    mapping = []
    for line in response_content.strip().split("\n"):
        gene_id, ko_id = line.split("\t")
        mapping.append((gene_id, ko_id.split(":")[1]))  # remove 'ko:' prefix

    # Convert to DataFrame
    return pd.DataFrame(mapping, columns=["hsa", "kegg.genes"])


