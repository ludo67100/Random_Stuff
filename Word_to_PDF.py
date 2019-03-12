# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 16:57:50 2019

J'ai la flemme de convertir toutes les partoches au format PDF pour impression 
du coup c'est ce bon vieux serpent qui va le faire pour moi 

@author: ludov
"""


def convert_pdf():
    
    """
    Converts list of .docx documents in .pdf
    Specifiy directories
    """
    
    import comtypes.client, os, re

    basedir = "D:\\Mygales/14_03/"
    workdir = "D:\\Mygales/14_03/00_PDF/"
    
    partoches = os.listdir(basedir)
    
    wdFormatPDF = 17 #Obviously the pdf format 
    
    for partoch in partoches:
        if partoch.endswith('.docx'):
            input_file = basedir + partoch
            
            new_partoch = re.sub('.docx$', '', partoch)
            
            output_file = workdir + new_partoch
            print(partoch, 'will be converted to pdf')
            
            word = comtypes.client.CreateObject('Word.Application')
            doc = word.Documents.Open(input_file)
            doc.SaveAs(output_file, FileFormat=wdFormatPDF)
            doc.Close()
            word.Quit()
            
            
        else:
            continue
        
def merge_pdf():
    
    """
    Merges pdfs into a single pdf file
    Specigy workdir
    """
    
    from PyPDF2 import PdfFileMerger
    import os 
    
    basedir = "D:\\Mygales/14_03/00_PDF/"
    
    pdf_list = os.listdir(basedir)
    
    pdfs = []
    
    for pdf in pdf_list:

        pdfs.append(basedir+pdf)
            
    merger = PdfFileMerger()
    
    for pdf in pdfs:
        merger.append(pdf)
    
    merger.write('{}/00_ALL_PDFS.pdf'.format(basedir))
    
     
        
if __name__ == '__main__':
    
    convert_pdf()
    merge_pdf()
    




