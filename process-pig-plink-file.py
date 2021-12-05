#coding = utf-8
import pandas as pd
import pymysql
import time
import os
import sys
import shutil

dir = "C:/Users/86151/Documents/国科大/考试/统计基因组/"
ori_file_path="C:/Users/86151/Documents/国科大/考试/统计基因组/T4.assoc.qassoc.adjusted"
backup_dir = "C:/Users/86151/Documents/国科大/考试/统计基因组/process-file/"
def file_rename(old_name_path):
  try:
    fpath,fname=os.path.split(old_name_path)
    # shutil.copyfile(old_name_path,backup_dir+fname)
    # os.rename(old_name_path,old_name_path.replace("adjusted","csv"))
    return old_name_path.replace("adjusted","csv")
  except Exception as e:
      print(e)
      print('rename file fail\r\n')

# GC值小于1.0e-6
def filter_P_value(p_value):
  return True if float(p_value) < float('1.0e-6') else False


def process(file_path):
  info = []
  for line in open(ori_file_path):
    info_list = [i for i in line.split("\t") if i and i != '\n']
    if len(info_list):
      info.append(info_list)
  info_all = pd.DataFrame(info[1::], columns=info[0])
  info_all = info_all.loc[info_all['GC'].apply(filter_P_value),:]
  snp_list = info_all['SNP'].values.tolist()
  df = pd.DataFrame()
  conn = pymysql.connect(host='lilaboc.dynv6.net',user = "xiehb",passwd = "weallfly", port=7016, db="speciation")
  cursor=conn.cursor()
  for snp in snp_list:
    (chr, snp_site)  = snp.split(".")[1:3]
    select_sql = f"select genename,start, end, (cast(end as signed)-cast(start as signed)) as length from speciation.genespan where chr={chr} and start< {snp_site} and end >{snp_site} order by length desc;"
    cursor.execute(select_sql)
    result = cursor.fetchall()
    result_df = pd.DataFrame(list(result),columns=['genename','start','end','length'])
    result_df['site'] = snp_site
    result_df['chr'] = chr
    # print(result_df)
    if not df.shape[0] >0:
      df = result_df
      continue
    df = df.append(result_df)
  df= (df.reset_index())[['genename','start','end','length','site','chr']]
  modify_result_df = df.loc[df['length']< 1000000,:]
  modify_result_df = (modify_result_df.reset_index())[['genename','start','end','length','site','chr']]
  df.to_csv(dir+"all_df.csv",sep='\t', header=True)
  modify_result_df.to_csv(dir+"T4modify_result_df.csv",sep='\t', header=True)
def main():
  file = file_rename(ori_file_path)
  process(file)
if __name__ == "__main__":
    main()








