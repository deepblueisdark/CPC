import locale
locale.setlocale(locale.LC_ALL,'pt_BR.UTF-8')
from datetime import datetime
import pandas as pd
from statistics import mean 
import openpyxl
from openpyxl.workbook import Workbook




filename = "CHUVA_CPC_COMPLETA_SIN.xlsx"
wb=openpyxl.load_workbook(filename)

sheet=wb['testechuva']


#ANO MES DIA grande iguacu jacui mosaico_SIN panema paraibadosul paranaiba parana saofrancisco sistema_sudeste tiete tocantins uruguai
#2022 11 20  3.50  0.06  0.00  5.04  0.01  0.64  6.92  1.50  6.47  4.29  2.92 11.15  0.08
#grande iguacu jacui panema paraibadosul parana paranaiba saofrancisco sistema_sudeste tiete tocantins uruguai



header=['MES/DIA',
        'DATA',
		'ANO',
        'MES',
        'DIA',
        'HORA',
        'GRANDE',
        'IGUACU',
        'JACUI',
        'PARANAPANEMA',
        'PARAIBASUL',
        'PARANÁ',
        'PARANAIBA',
        'SAOFRANCISCO',
        'SIN',
        'TIETE',
        'TOCANTINS',
        'URUGUAI']

# ANO MES DIA SIN SIN SIN SIN SIN SIN SIN SIN SIN SIN SIN SIN
# 2022 11 28  1.98  2.35  0.03  4.13  6.09  3.48  0.38 10.16  5.97  4.19  5.89  1.46
# 2022 11 29  8.16  3.05  0.00  6.15 12.16  8.64  0.85 18.65  8.64  9.40 13.71  0.78
# 2022 11 30 11.47  3.48  0.32  3.29 10.88  8.32  3.05 11.81  9.40 15.46 13.45  2.57
# 2022 12  1  8.39  2.59  4.63  6.72 12.18  0.51  2.84  3.81  4.86  8.75  5.06  3.89
# 2022 12  2  3.84 13.81  0.00  5.90  6.41  3.70  3.09  8.47  5.52  1.53  7.58  2.77
# 2022 12  3  6.93  3.04 15.42  4.73  6.99  7.59  8.75  6.04  6.60  8.55 10.22 10.67
# 2022 12  4 15.61 18.66  5.44  7.67 10.20 11.12 11.78  4.19 11.99 14.18 20.75  5.24

col_names=['ano',
           'mes',
           'dia',
           'URUGUAI', ##eira 
           'IGUACU',  ## cacho
           'PARANA',  ##Madre
           'PARANAPANEMA',  ## Beni
           'TIETE',   ##mamor
           'GRANDE',   ##inc ate JRB
           'PARANAIBA',  ##guapore
           'SAOFRANCISCO',   ##JRB jirau
           'TOCANTINS',   ##Jirau-SAE
           'MADEIRA',    ##JRB-SAE
           'SUDESTE',   ##montJP
           'PARAIBADOSUL', 
           'JACUI']    ## Inc_JPVILA

 
coluna=("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","Y","Z")   
 
#df = pd.read_csv("testechuva.csv", delim_whitespace=True,names=col_names,header=None,skiprows=1)
df = pd.read_excel("chuva_SIN.xlsx", sheet_name="chuva_media", names=col_names, header=0)
num_itens=len(df)
print(df)

#print(df.ano)
#df['SAE01']=df['SAE01'].astype(float)
#print(df.dtypes) 


num_cols=len(header)
print(num_cols)
for n in range(0,num_cols):
    celula=coluna[n]+"1"
    sheet[celula]=header[n]









for dia in range(0,num_itens):
    start0=datetime(1979,1,1,0,0)
    atual=datetime(df.ano[dia], df.mes[dia], df.dia[dia], 0, 0)
    num_dias_time=atual-start0
    #print(num_dias_time)
    #print(type(num_dias_time))
    #print(type(atual))
    num_dias=str(num_dias_time)
    row=2+int(num_dias[0:6],10) 
	

    celula=coluna[2]+str(row)
    sheet[celula]=df.ano[dia]   

    celula=coluna[3]+str(row)
    sheet[celula]=df.mes[dia]   

    celula=coluna[4]+str(row)
    sheet[celula]=df.dia[dia]

    celula=coluna[5]+str(row)
    sheet[celula]=0.0

    celula=coluna[6]+str(row)
    sheet[celula]=df.GRANDE[dia] 
	
    celula=coluna[7]+str(row)
    sheet[celula]=df.IGUACU[dia]    ## mamore
    
    celula=coluna[8]+str(row)
    sheet[celula]=df.JACUI[dia]     ## mont_jpvilla
    
    celula=coluna[9]+str(row)
    sheet[celula]=df.PARANAPANEMA[dia]     ### JRB-JIRAU
    
    celula=coluna[10]+str(row)
    sheet[celula]=df.PARAIBADOSUL[dia]    ### jirau-SAE

    celula=coluna[11]+str(row)
    sheet[celula]=df.PARANA[dia] 
   
	
    celula=coluna[12]+str(row)
    sheet[celula]=df.PARANAIBA[dia] 
    ##Inc_JRB	###Beni
    celula=coluna[13]+str(row)
    sheet[celula]=df.SAOFRANCISCO[dia]   ###guapore

    celula=coluna[14]+str(row)
    sheet[celula]=df.SUDESTE[dia]    ###jrb_sae

    celula=coluna[15]+str(row)
    sheet[celula]=df.TIETE[dia]   ###JRB

    celula=coluna[16]+str(row)
    sheet[celula]=df.TOCANTINS[dia]    ###madredios

    celula=coluna[17]+str(row)
    sheet[celula]=df.URUGUAI[dia]      ### inc_JVILA     
 
	

wb.save(filename)