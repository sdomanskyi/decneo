import hashlib
import json
import os
import requests
from requests.exceptions import HTTPError

BASE_URL = 'https://api.figshare.com/v2/{endpoint}'
PRE_REQ = 'account/articles'

with open('/mnt/home/domansk6/figshare_upload_token.txt', 'r') as tempFile:
    TOKEN = tempFile.readline()

CHUNK_SIZE = 1048576

def raw_issue_request(method, url, data=None, binary=False):
    headers = {'Authorization': 'token ' + TOKEN}
    if data is not None and not binary:
        data = json.dumps(data)
    response = requests.request(method, url, headers=headers, data=data)
    try:
        response.raise_for_status()
        try:
            data = json.loads(response.content)
        except ValueError:
            data = response.content
    except HTTPError as error:
        print('Caught an HTTPError: {}'.format(error))
        print('Body:\n', response.content.decode())
        raise

    return data

def issue_request(method, endpoint, *args, **kwargs):

    return raw_issue_request(method, BASE_URL.format(endpoint=endpoint), *args, **kwargs)

def list_articles():
    result = issue_request('GET', PRE_REQ)
    print('Listing current articles:')
    if result:
        for item in result:
            print(u'  {url} - {title}'.format(**item))
    else:
        print('  No articles.')
    print()

def create_article(title):
    data = {
        'title': title,  # You may add any other information about the article here as you wish.
    }
    result = issue_request('POST', PRE_REQ, data=data)
    print('Created article:', result['location'], '\n')

    result = raw_issue_request('GET', result['location'])

    return result['id']

def list_files_of_article(article_id):
    result = issue_request('GET', PRE_REQ + '/{}/files'.format(article_id))
    print('Listing files for article {}:'.format(article_id))
    if result:
        for item in result:
            print('  {id} - {name}'.format(**item))
    else:
        print('  No files.')

    print()

def get_file_check_data(file_name):
    with open(file_name, 'rb') as fin:
        md5 = hashlib.md5()
        size = 0
        data = fin.read(CHUNK_SIZE)
        while data:
            size += len(data)
            md5.update(data)
            data = fin.read(CHUNK_SIZE)
        return md5.hexdigest(), size

def initiate_new_upload(article_id, file_name):

    md5, size = get_file_check_data(file_name)
    data = {'name': os.path.basename(file_name),
            'md5': md5,
            'size': size}

    result = issue_request('POST', PRE_REQ + '/{}/files'.format(article_id), data=data)
    print('Initiated file upload:', result['location'], '\n')

    result = raw_issue_request('GET', result['location'])

    return result

def complete_upload(article_id, file_id):

    issue_request('POST', PRE_REQ + '/{}/files/{}'.format(article_id, file_id))

def upload_parts(file_info, FILE_PATH):
    url = '{upload_url}'.format(**file_info)
    result = raw_issue_request('GET', url)

    print('Uploading parts:')
    with open(FILE_PATH, 'rb') as fin:
        for part in result['parts']:
            upload_part(file_info, fin, part)
    print()

def upload_part(file_info, stream, part):
    udata = file_info.copy()
    udata.update(part)
    url = '{upload_url}/{partNo}'.format(**udata)

    stream.seek(part['startOffset'])
    data = stream.read(part['endOffset'] - part['startOffset'] + 1)

    raw_issue_request('PUT', url, data=data, binary=True)
    print('  Uploaded part {partNo} from {startOffset} to {endOffset}'.format(**part))

def upload(filePath, fileName, article_id = None):

    if article_id is None:
        article_id = create_article(fileName)

    list_articles()
    list_files_of_article(article_id)

    # Upload the file
    file_info = initiate_new_upload(article_id, filePath)

    # Use the figshare upload service API
    upload_parts(file_info, filePath)

    # Complete the file upload process
    complete_upload(article_id, file_info['id'])

    list_files_of_article(article_id)

    return

if __name__ == '__main__':

    fileName = 'VoightChoroid4567RemappedData.h5'

    upload('/mnt/home/domansk6/Projects/Endothelial/scripts/demo/' + fileName, fileName, '13179818')