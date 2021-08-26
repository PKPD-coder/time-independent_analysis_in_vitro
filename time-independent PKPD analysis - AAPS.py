#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday Jun 18 16:11:17 2019

@author: ovacika

"""



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from sklearn.metrics import auc
from scipy.optimize import curve_fit
import xlsxwriter
mpl.rc('figure', max_open_warning = 0) ## this is to turn off an error for plotting  



def plotsObsOverTime(DrugName,ReadoutName,CellLineName,Data):
    # this functions plots the observatsions over time at ecach concentration
    # Figure name is Drug name + Readout + Cell line name, each concentraiton is plotted on the same Figure
    fig= plt.figure()
    Get_ConcandTime=Data.loc[DrugName,ReadoutName,CellLineName]
    Conc_Id=Get_ConcandTime.index.codes[0]
    my_List=[]
    for x in Conc_Id:
        my_List.append(x)
    my_List = list(set(my_List))         
    colors=['purple','blue','teal','green','yellow','orange','red','darkred','maroon']
    color_index=0
    for y in my_List:  
        Get_Conc=Get_ConcandTime.index.levels[0][y] 
        Final_Data=Get_ConcandTime.loc[Get_Conc]     
        label_name=str(Get_Conc)+' pM'
        Final_Data.plot( x='Time',y = 'Observation', label=label_name, kind='line',color='black',  marker='o',markerfacecolor=colors[color_index])
        color_index=color_index+1
    plt.xlabel('Time(hr)')
    plt.ylabel('Observation')
    fig_name=str(drug) + ' '+  cell + '_'+ Read 
    fig.suptitle(fig_name)
    plt.legend()
    plt.savefig('./Figures/'+Read+ ' Over Time ' + str(drug)+' '+cell) 
    return Final_Data


def Sigmoidal_model(Concentration,Top,Bottom,EC50,hill):
    # this function simulates Emax-EC50 curve for a given concentration range
    return Bottom + pow(Concentration,hill)*(Top-Bottom)/(pow(EC50,hill) + pow(Concentration,hill))

def Hockey_StickModel(Concentration,Bottom,Threshold_conc,S):
    # this function simulates Hockey_Stick curves  for a given concentration range
    P=len(Concentration)*[1]
    for i in range(len(Concentration)):
        if Concentration[i] <= Threshold_conc:
            P[i]=0   
    return Bottom + S*(Concentration-Threshold_conc)*P

def Calculate_AUCs(Drug_id,Readout_id,Cellline_id,SeriesData):
   # this function calculates AUCs for a given experiment
   # the  inputs are the drug name , Readout name and cell line  name, dataset and the timepoints
   # the outputs are the AUC valuues at each concentration 
    AUC=[]
    Get_ConcandTime=SeriesData.loc[Drug_id,Readout_id,Cellline_id]
    Conc_Id=Get_ConcandTime.index.codes[0] 
    Time_Id=Get_ConcandTime.index.codes[1]            
    Time_List=[]
    for a in Time_Id:
        Time_List.append(a)
    Time_List=list(set(Time_List))  
    Time_newList=[]            
    for b in Time_List:
        Get_Time=Get_ConcandTime.index.levels[1][b] 
        Time_newList.append(Get_Time)     
    my_List=[]
    for x in Conc_Id:
         my_List.append(x)
    Conc_List=[]
    my_List = list(set(my_List))                   
    for y in my_List:  
        Get_Conc=Get_ConcandTime.index.levels[0][y] 
        Conc_List.append(Get_Conc)
        Final_Data=Get_ConcandTime.loc[Get_Conc]  
        Time_data=[]
        for a in range(len(Final_Data.index)):   
            Time_data.append(Final_Data.index[a])
        AUC.append(auc(Time_data,Final_Data.to_numpy()))
    return AUC,Conc_List

def Calculate_Rmax(Drug_id,Readout_id,Cellline_id,SeriesData):
   # this function calculates AUCs for a given experiment
   # the  inputs are the drug name , Readout name and cell line  name, dataset and the timepoints
   # the outputs are the AUC valuues at each concentration 
    Rmax=[]
    Get_ConcandTime=SeriesData.loc[Drug_id,Readout_id,Cellline_id]
    Conc_Id=Get_ConcandTime.index.codes[0] 
    Time_Id=Get_ConcandTime.index.codes[1]            
    Time_List=[]
    for a in Time_Id:
        Time_List.append(a)
    Time_List=list(set(Time_List))  
    Time_newList=[]            
    for b in Time_List:
        Get_Time=Get_ConcandTime.index.levels[1][b] 
        Time_newList.append(Get_Time)     
    my_List=[]
    for x in Conc_Id:
         my_List.append(x)
    Conc_List=[]
    my_List = list(set(my_List))                   
    for y in my_List:  
        Get_Conc=Get_ConcandTime.index.levels[0][y] 
        Conc_List.append(Get_Conc)
        Final_Data=Get_ConcandTime.loc[Get_Conc]  
        Time_data=[]
        for a in range(len(Final_Data.index)):   
            Time_data.append(Final_Data.index[a])        
        Rmax.append(max(Final_Data.to_numpy()))
    return Rmax,Conc_List


def AUCdata_toFit(AUC_data,Conc_List):
   # this function fit the concentration AUC data to Emax and EC50 model. 
   Dict={}
   List=[]
   init_guess =[max(AUC_data),min(AUC_data),(max(Conc_List)+min(Conc_List))/2,1] 
   Fitted_params,covariates = curve_fit(Sigmoidal_model,Conc_List ,AUC_data, p0=init_guess ,bounds=((-np.inf,-np.inf,min(Conc_List),0.01), (np.inf,np.inf,max(Conc_List),10)),maxfev = 500000)
   SE_perc = np.sqrt(np.diag(covariates))/Fitted_params
   SE_Top,SE_Bottom,SE_EC50,SE_hill=100*SE_perc      
   fit_Top,fit_Bottom,fit_EC50,fit_hill = Fitted_params 
   dict_name= str(drug) + ' '+  Read + ' '+ cell
   Dict.update({dict_name :fit_EC50})
   List.append([drug,Read,cell,fit_EC50])
   return fit_Top,fit_Bottom,fit_EC50,fit_hill,SE_Top,SE_Bottom,SE_EC50,SE_hill

def AUCdata_toFitHockey_stick(AUC_data,Conc_List):
   # this function fit the concentration AUC data to Hockey stick model
   Dict={}
   List=[]
   init_guess =[min(AUC_data),1,1] 
   Fitted_params,covariates = curve_fit(Hockey_StickModel,Conc_List ,AUC_data, p0=init_guess ,bounds=((-np.inf,-np.inf,-np.inf),(np.inf,np.inf,np.inf)),maxfev = 500000)
   SE_perc = np.sqrt(np.diag(covariates))/Fitted_params
   SE_Bottom,SE_Threshold,SE_S=100*SE_perc      
   fit_Bottom,fit_Threshold,fit_S,= Fitted_params 
   dict_name= str(drug) + ' '+  Read + ' '+ cell
   Dict.update({dict_name :fit_Threshold})
   List.append([drug,Read,cell,fit_Threshold])
   return fit_Bottom,fit_Threshold,fit_S,SE_Bottom,SE_Threshold,SE_S


def Simulation_SigmoidalModel(Conc_List,T,B,EC,h):
   Simulation_Conc=np.linspace(min(Conc_List),max(Conc_List),100000)
   Stim_output=[]
   for x in Simulation_Conc:
        Stim_output.append(Sigmoidal_model(x,T,B,EC,h))    
   return  Simulation_Conc,Stim_output   


def Simulation_Hockey_StickModel(Conc_List,Bo,Thres,S):
   Simulation_Conc=np.linspace(min(Conc_List),max(Conc_List),100000)
   Stim_output=[]
   Stim_output=Hockey_StickModel(Simulation_Conc,Bo,Thres,S)   
   return  Simulation_Conc,Stim_output 

def plotConcvsAUC(Sim_x,Sim_y,Data_x,Data_y,DrugName,Readout_name,CellLine_Name):
   fig = plt.figure()
   plt.semilogx(Sim_x,Sim_y,label="Sigmoidal Model Fit ",color='black')                       
   plt.semilogx(Data_x,Data_y,'o',label=str(drug) + ' '+  Read + ' '+ cell, color='black',  marker='o',markerfacecolor='red')
   plt.xlabel('Concentration (pM)')
   plt.ylabel('AUCE')
   fig_name= 'Model Fit ' + str(DrugName) + ' '+ CellLine_Name + '_'+  Readout_name  # Figure name is ID + Readout + Treatment, each concentraiton is plotted on the same Figure
   fig.suptitle(fig_name)
   plt.legend()
   plt.savefig('./Figures/'+ fig_name +'.jpeg')    
   plt.close()
   return str(DrugName) + ' '+  CellLine_Name + '_'+  Readout_name 

def plotConcvsRmax(Sim_x,Sim_y,Data_x,Data_y,DrugName,Readout_name,CellLine_Name):
   fig = plt.figure()
   plt.semilogx(Sim_x,Sim_y,label="Sigmoidal Model Fit ",color='black')                       
   plt.semilogx(Data_x,Data_y,'o',label=str(drug) + ' '+ cell + '_'+ Read, color='black',  marker='o',markerfacecolor='red')
   plt.xlabel('Concentration (pM)')
   plt.ylabel('Rmax')
   fig_name= 'Model Fit ' + str(DrugName) + ' '+ CellLine_Name + '_'+  Readout_name  # Figure name is ID + Readout + Treatment, each concentraiton is plotted on the same Figure
   fig.suptitle(fig_name)
   plt.legend()
   plt.savefig('./Figures/'+ fig_name +'_Rmax.jpeg')    
   plt.close()
   return str(DrugName) + ' '+  CellLine_Name + '_'+  Readout_name 
 
# def plotConcvsAUCHockeyStick(Sim_x,Sim_y,Data_x,Data_y,DrugName,Readout_name,CellLine_Name):
#    fig = plt.figure()
#    plt.semilogx(Sim_x,Sim_y,label="Sigmoidal Model Fit ",color='black')                       
#    plt.semilogx(Data_x,Data_y,'o',label=str(drug) + ' '+  Read + ' '+ cell, color='black',  marker='o',markerfacecolor='red')
#    plt.xlabel('Concentration (pM)')
#    plt.ylabel('AUCE')
#    fig_name= 'Model Fit ' + str(DrugName) + ' '+ CellLine_Name + '_'+  Readout_name+ '_HOCKEYSTICK'  # Figure name is ID + Readout + Treatment, each concentraiton is plotted on the same Figure
#    fig.suptitle(fig_name)
#    plt.legend()
#    plt.savefig('./Figures/'+ fig_name +'.jpeg')    
#    plt.close()
#    return str(DrugName) + ' '+  CellLine_Name + '_'+  Readout_name 

      
def File_write(Filename,FigureName,T,B,EC,h,err_T,err_B,err_EC50,err_h): 
      Filename.write(FigureName+"\t"+ str(round(T))+"\t" + str(round(err_T,2))+"\t" + str(round(B)) +"\t"+str(round(err_B,3))+"\t" + str(round(EC,3))+ "\t"+str(round(err_EC50,3))+"\t" + str(round(h,3))+ "\t"+str(round(err_h,3))+"\t" + '\n' )   
      
def File_write2(Filename,FigureName,Bott,Threshold,S,err_Bott,err_Thr,err_Sv): 
      Filename.write(FigureName+"\t"+ str(round(Bott))+"\t" + str(round(err_Bott,2))+"\t" + str(round(Threshold,3)) +"\t"+str(round(err_Thr,3))+"\t" + str(round(S))+ "\t"+str(round(err_Sv,3))+"\t" + '\n' )



# In[Initialization] :
## the code starts here
plt.ioff()
 
      
File= open('Sigmoidal Model Parameters.txt','w') ## This line creates a text file and write the below line on top of the txt line
File.write('Figure name '+"\t"+'Top'+ "\t" +'Top SE %' +"\t" + 'Bottom' +"\t"+ 'Bottom SE %' +"\t" 'EC50' +"\t" +'EC50 SE %' + "\t" 'Hill' +"\t" +'Hill SE %' +'\n' )
File2= open('Hockey Stick Model Parameters.txt','w') ## This line creates a text file and write the below line on top of the txt line
File2.write('Figure name '+"\t"+'Bottom'+ "\t" +'Bottom SE %' +"\t" + 'Threshold value' +"\t"+ 'Threshold value SE %' +"\t" 'S' +"\t" +'S SE %'  +'\n' )
File3= open('Sigmoidal Model Parameters Rmax.txt','w') ## This line creates a text file and write the below line on top of the txt line
File3.write('Figure name '+"\t"+'Top'+ "\t" +'Top SE %' +"\t" + 'Bottom' +"\t"+ 'Bottom SE %' +"\t" 'EC50' +"\t" +'EC50 SE %' + "\t" 'Hill' +"\t" +'Hill SE %' +'\n' )
File4= open('Hockey Stick Model Parameters Rmax.txt','w') ## This line creates a text file and write the below line on top of the txt line
File4.write('Figure name '+"\t"+'Bottom'+ "\t" +'Bottom SE %' +"\t" + 'Threshold value' +"\t"+ 'Threshold value SE %' +"\t" 'S' +"\t" +'S SE %'  +'\n' )


Data_read = pd.read_excel('CEA-TCB dataset for dynamic analysis - AAPS.xlsx', sheet_name='Sheet1')  ##This line uploads excelDataFile to DataFrame
df_Drug=Data_read['Drug'] 
Drug_Ids=df_Drug.unique(); # unique Drug Ids
df_CellLine=Data_read['CellLine']
CellLine=df_CellLine.unique(); # unique Cell Lone descriptions
df_PD_Readout=Data_read['Readout']
Readouts=df_PD_Readout.unique(); # unique Readouts - cell killing or cytokines
Groupd_ser_data=Data_read.groupby(['Drug','Readout','CellLine','Conc','Time']).Obs.mean() # this transforms dataframe to a series


# In[Initialization] :

#This loop creates 
# 1-) Oberservation over time at each concentration figures
# 2-) a txt file with EMax-EC50 model parameters and associated SE % 
# 3-) figures for each model fit
# 4-) a txt file with simulations for each model fit

Dict_f={}
List_f=[]
Flag_f=[]

workbook=xlsxwriter.Workbook('Sigmoidal Model Simulations.xlsx')
workbook2=xlsxwriter.Workbook('Calculated AUCE.xlsx')
for drug in Drug_Ids: #for each unique ID
    for cell in CellLine: 
      for Read in Readouts:  #for each unique REadout
            Mero=plotsObsOverTime(drug,Read,cell,Groupd_ser_data)   
            AUC_f=[]          
            AUC_f,Conc_List_f=Calculate_AUCs(drug,Read,cell,Groupd_ser_data)           
            Top,Bottom,EC50,hill,err_Top,err_Bottom,err_EC50,err_hill= AUCdata_toFit(AUC_f,Conc_List_f)
            flag=0
            if abs(err_EC50) > 1.5 or EC50<0 :
                flag=1
            List_f.append([drug,Read,cell,EC50])      
            Flag_f.append([drug,Read,cell,flag])
            Sim_C_f,Simulation_f=Simulation_SigmoidalModel(Conc_List_f,Top,Bottom,EC50,hill)         
            Fig_name=plotConcvsAUC(Sim_C_f,Simulation_f,Conc_List_f, AUC_f,drug,Read,cell)           
            File_write(File,Fig_name,Top,Bottom,EC50,hill,err_Top,err_Bottom,err_EC50,err_hill)  
            Report_Sheet= workbook.add_worksheet(Fig_name)            
            for row_ind, row_value in enumerate(zip(Sim_C_f,Simulation_f)):
                for col_ind, col_value in enumerate(row_value):
                    Report_Sheet.write(row_ind + 1, col_ind, col_value)
            Report_Sheet.write(0, 0, 'Concentration pM')
            Report_Sheet.write(0, 1, 'AUCE')
            Report_Sheet= workbook2.add_worksheet(Fig_name)
            for row_ind, row_value in enumerate(zip(Conc_List_f,AUC_f)):
                for col_ind, col_value in enumerate(row_value):
                    Report_Sheet.write(row_ind + 1, col_ind, col_value)
            Report_Sheet.write(0, 0, 'Concentration pM')
            Report_Sheet.write(0, 1, 'AUCE')
File.close()           
workbook.close()
workbook2.close()


workbook3=xlsxwriter.Workbook('Hockey Stick Model Simulations.xlsx')
for drug in Drug_Ids: #for each unique ID
    for cell in CellLine: 
      for Read in Readouts:  #for each unique REadout
            AUC_f=[]          
            AUC_f,Conc_List_f=Calculate_AUCs(drug,Read,cell,Groupd_ser_data)           
            #Top,Bottom,EC50,hill,err_Top,err_Bottom,err_EC50,err_hill= AUCdata_toFit(AUC_f,Conc_List_f)
            if abs(err_Top) > .2 or abs(err_Bottom) > .2 or abs(err_EC50) > .2:
                Bottom_HS,Threshold,S,err_Bottom_HS,err_Threshold,err_S= AUCdata_toFitHockey_stick(AUC_f,Conc_List_f)
                Sim_Conc_HS,Simulation_HS=Simulation_Hockey_StickModel(Conc_List_f,Bottom_HS,Threshold,S)
                #if abs(err_Threshold) < 5:
                Fig_name=plotConcvsAUC(Sim_Conc_HS,Simulation_HS,Conc_List_f, AUC_f,drug,Read,cell+'_HOCKEYSTICK')     
                Report_Sheet= workbook3.add_worksheet(str(drug)+' '+ cell + '_' + Read)            
                for row_ind, row_value in enumerate(zip(Sim_Conc_HS,Simulation_HS)):
                        for col_ind, col_value in enumerate(row_value):
                            Report_Sheet.write(row_ind + 1, col_ind, col_value)
                            Report_Sheet.write(0, 0, 'Concentration pM')
                            Report_Sheet.write(0, 1, 'AUCE')
                    
                File_write2(File2,str(drug)+' '+ cell + '_' + Read,Bottom_HS,Threshold,S,err_Bottom_HS,err_Threshold,err_S)  
                    
File2.close()           
workbook3.close()

# In[Initialization] : 

workbook4=xlsxwriter.Workbook('Sigmoidal Model Simulations for Rmax.xlsx')
workbook5=xlsxwriter.Workbook('Observed Rmax.xlsx')

for drug in Drug_Ids: #for each unique ID
    for cell in CellLine: 
      for Read in Readouts:  #for each unique REadout
            Rmax_f=[]          
            Rmax_f,Conc_List_f=Calculate_Rmax(drug,Read,cell,Groupd_ser_data)          
            Top,Bottom,EC50,hill,err_Top,err_Bottom,err_EC50,err_hill= AUCdata_toFit(Rmax_f,Conc_List_f)       
            flag=0
            if abs(err_EC50) > 0.5 or EC50<0 :
                flag=1
            List_f.append([drug,Read,cell,EC50])      
            Flag_f.append([drug,Read,cell,flag])
            Sim_C_f,Simulation_f=Simulation_SigmoidalModel(Conc_List_f,Top,Bottom,EC50,hill)         
            Fig_name=plotConcvsRmax(Sim_C_f,Simulation_f,Conc_List_f, Rmax_f,drug,Read,cell)           
            File_write(File3,Fig_name,Top,Bottom,EC50,hill,err_Top,err_Bottom,err_EC50,err_hill)  
            Report_Sheet= workbook4.add_worksheet(Fig_name)            
            for row_ind, row_value in enumerate(zip(Sim_C_f,Simulation_f)):
                for col_ind, col_value in enumerate(row_value):
                    Report_Sheet.write(row_ind + 1, col_ind, col_value)
            Report_Sheet.write(0, 0, 'Concentration pM')
            Report_Sheet.write(0, 1, 'Rmax')
            Report_Sheet= workbook5.add_worksheet(Fig_name)
            for row_ind, row_value in enumerate(zip(Conc_List_f,AUC_f)):
                for col_ind, col_value in enumerate(row_value):
                    Report_Sheet.write(row_ind + 1, col_ind, col_value)
            Report_Sheet.write(0, 0, 'Concentration pM')
            Report_Sheet.write(0, 1, 'Rmax')
File3.close()           
workbook4.close()
workbook5.close()


workbook6=xlsxwriter.Workbook('Hockey Stick Model Simulations Rmax.xlsx')
for drug in Drug_Ids: #for each unique ID
    for cell in CellLine:
      for Read in Readouts:  #for each unique REadout
            Rmax_f=[]          
            Rmax_f,Conc_List_f=Calculate_Rmax(drug,Read,cell,Groupd_ser_data)           
            Top,Bottom,EC50,hill,err_Top,err_Bottom,err_EC50,err_hill= AUCdata_toFit(Rmax_f,Conc_List_f)
#            if abs(err_Top) > .5:
            Bottom_HS,Threshold,S,err_Bottom_HS,err_Threshold,err_S= AUCdata_toFitHockey_stick(Rmax_f,Conc_List_f)
            Sim_Conc_HS,Simulation_HS=Simulation_Hockey_StickModel(Conc_List_f,Bottom_HS,Threshold,S)         
#                if abs(err_Threshold) < 5:
            Fig_name=plotConcvsRmax(Sim_Conc_HS,Simulation_HS,Conc_List_f, Rmax_f,drug,Read,cell+'_HOCKEYSTICK')           
            Report_Sheet= workbook6.add_worksheet(str(drug)+' '+ cell + '_' + Read)            
            for row_ind, row_value in enumerate(zip(Sim_Conc_HS,Simulation_HS)):
                  for col_ind, col_value in enumerate(row_value):
                            Report_Sheet.write(row_ind + 1, col_ind, col_value)
                            Report_Sheet.write(0, 0, 'Concentration pM')
                            Report_Sheet.write(0, 1, 'Rmax')
            File_write2(File4,str(drug)+' '+ cell + '_' + Read,Bottom_HS,Threshold,S,err_Bottom_HS,err_Threshold,err_S)  
                    
File4.close()           
workbook6.close()


# In[Initialization] :


EC50_data = pd.DataFrame(List_f,columns = ['Drug' , 'Y_type', 'CL','EC50']) 
Flag_data=pd.DataFrame(Flag_f,columns = ['Drug' , 'Y_type', 'CL','Flag']) 
New_EC50_data=EC50_data.groupby(['Drug','Y_type','CL']).EC50.mean()
New_Flag_data=Flag_data.groupby(['Drug','Y_type','CL']).Flag.mean()


#This loop creates bar charts  of EC50 for each readout across drugs and expression levels
for Read in Readouts:
    Fig2=plt.figure()
    dat=[]
    fla=[]
    str_dat=[]
    for drug in Drug_Ids: 
        for cell in CellLine:
            dat.append(New_EC50_data.loc[drug][Read][cell])
            fla.append(New_Flag_data.loc[drug][Read][cell])
            str_dat.append(str(drug)+' '+cell)
    for i in range(len(dat)):
        if fla[i]==1:
            dat[i]=0   
    plt.bar(str_dat,dat,color=['limegreen','darkgreen','blue','darkblue']) 
    plt.gcf().subplots_adjust(bottom=0.35)
    plt.xticks(range(len(str_dat)), str_dat, rotation=45)
    plt.xlabel(' Drug and Cell line')
    plt.ylabel('EC50 (pM)')
    fig2_name=Read # Figure name is Readout +, each concentraiton is plotted on the same Figure
    Fig2.suptitle(fig2_name + '\n')
    plt.savefig('./Figures/'+ 'Cumulative Potency ' + fig2_name +'.jpeg')  



