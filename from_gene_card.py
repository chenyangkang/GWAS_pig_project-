# 导入所需的第三方库
import requests
from lxml import etree
import openpyxl
import xlsxwriter
import time
import pandas as pd

# 设置base_url是为了得到url的cookies，网站设置了反爬机制
base_url = "https://www.genecards.org"
url ="https://www.genecards.org/cgi-bin/carddisp.pl?gene=A2BP1"

headers = {
    'upgrade-insecure-requests': '1',
    'cookie': 'rvcn=rq6mqJwBYsnKc12jFOkjzDrkwS1f_r0rymXE7GstdDFLf_7VLy_P7gtjDD5hP4PcrWNpiQYjdhYoJ8heXepq9X_Mch41; visid_incap_146342=xLquDkzqSwC8QGN5SXqmqZsuo2EAAAAAQUIPAAAAAAArFuKZxo0xC0Y6fIhku3Sh; ASP.NET_SessionId=ebs5luykhdmgit0w5xjjrnua; ARRAffinity=5b37b8ccdf277074fbf566cec0b5c40d45b5fae0727ca76b5409c1b6da6b189c; ARRAffinitySameSite=5b37b8ccdf277074fbf566cec0b5c40d45b5fae0727ca76b5409c1b6da6b189c; nlbi_146342=VW4kVx/4GFBlMkmBRXA9UQAAAADEF6rKa63Xq2x1ND1kYbYO; incap_ses_637_146342=parRIdmWrCallpB9jRTXCO4jpmEAAAAAPiggnhtk5TOPj5i+KTnB3A==; _ga=GA1.2.937150690.1638278134; _gid=GA1.2.690157045.1638278134',
    'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/96.0.4664.45 Safari/537.36',
}


# 获取cookies
def get_cookies(url):
    try:
        requests.session()
        sessions = requests.get(url, headers=headers)
        cookie = sessions.cookies
    except:
        cookie = ''
    return cookie


# 获取所需详情页信息
def get_search_response(url, cookies):
    try:
        response = requests.get(url, headers=headers, cookies=cookies)
        content = response.content.decode()
        data = etree.HTML(content)
        entrez_gene_summary= data.xpath('//section[@id="summaries"]/div/ul/li/p/text()')[0]
    except:
        entrez_gene_summary = 'NA'
        
    return entrez_gene_summary


# 读取保存基因名称的Excel文件
def read_excel():
    file = openpyxl.open('C:/Users/86151/Documents/国科大/考试/统计基因组/T4modify_result_df.xlsx')
    table = file.get_sheet_by_name('Gene')
    ncols = table.columns
    gene_list = []
    for col in ncols:
        line = [row.value for row in col]
        line = list(filter(None, line))
        gene_list.append(line)
    return gene_list



whole_gene_list = read_excel()
res ={}
res_list=[]
for index, li in enumerate(whole_gene_list):
        name_lis = []
        content_lis = []
        for item in li[1:]:
            time.sleep(1)
            print(item)
            new_url = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + item
            cookies = get_cookies(base_url)
            print(get_search_response(new_url, cookies))
            content = get_search_response(new_url, cookies)
            res_list.append([item, ''.join(content[0::])])
            res[item] = content

data_df = pd.DataFrame(res_list)
data_df.to_csv("C:/Users/86151/Documents/国科大/考试/统计基因组/T4gene_summary.csv",sep='\t', header=True)
