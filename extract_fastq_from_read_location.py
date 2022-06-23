#! /usr/bin/python
import re
f=open("unclassified_trim800_1000.fastq","r")
# for id_search in id_text:
#     print(id_search)
#     id_search=id_search.strip()
#     id_list=id_search.split("_")
#     id_search=id_list[0]
#     start=id_list[1]
#     end=id_list[2]
#     print(id_search)
#     if id_search == ">bb9d24a1-5947-4b1b-b111-e0040f90f46b":
#         print("OK")
# id_text.close()
out=open("test.fastq","+w")
for line in f:
    line=line.strip()
    if re.match("@",line):
        info = line.split(" ")
        id = str(info[0].replace("@",">"))
        seq=f.readline().strip()
        plus=f.readline().strip()
        quality=f.readline().strip()
        id_text=open("id_location.txt","r")
        for id_line in id_text:
            id_line=id_line.strip()
            id_list=id_line.split("_")
            id_search=id_list[0]
            start=int(id_list[1])
            end=int(id_list[2])
            if id == id_search:
                seq_true=seq[start:end]
                print(seq_true)
        id_text.close()
out.close()
f.close()
